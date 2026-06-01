#ifndef MOLRESPONSE_V3_CALC_CALC_EXECUTOR_HPP
#define MOLRESPONSE_V3_CALC_CALC_EXECUTOR_HPP

// ===========================================================================
// CalcManager executor — the WORLD-BOUND half of the calc manager (doc 15,
// 15a). The scheduling core (calc_manager.hpp) decides what runs and in what
// order; this header actually drives the solves.
//
//   FdResponseExecutor : ICalcExecutor
//       run_protocol(WorkItem) -> one protocol step of an FD/NuclearFD
//       solve. Dispatches Static/Full × Closed/Open at runtime, builds the
//       perturbation source, honors the reconcile action for seed loading,
//       runs solvers::iterate_protocol over the single protocol step threshold (with
//       a reproject-prepare and a save post_step), and reports convergence.
//
//   CalcManager
//       build(n_atoms) -> dag from the plan; run(world, exec) -> the
//       one-wave-per-pass reschedule loop. Each pass reloads the aggregate
//       metadata and re-schedules, so every protocol step's action is computed against
//       up-to-date disk state (the protocol-step-level barrier the design decided on).
//
// Scope (15a): FD (dipole + nuclear) closed-shell, dipole open-shell. ES
// bundles and DerivedFD expansion are recognized but not yet solved here —
// run_protocol returns not-converged with a clear message and run() stops cleanly
// once it can make no further progress. The ES executor + derived expansion
// are the next increment; the expansion hook below is already wired so only
// the ES solve has to land.
//
// This header reuses, unchanged, exactly the machinery the FD skeleton test
// exercises: GroundState::prepare, build_response_ground_state_*,
// dipole/nuclear_perturbation, FDSolver, iterate_protocol, try_load_fd_state,
// save_fd_state.
// ===========================================================================

#include "calc_manager.hpp"

#include "../GroundState.hpp"
#include "../Perturbations.hpp"
#include "../ResponseProtocol.hpp"
#include "../kernels/full.hpp"
#include "../kernels/static.hpp"
#include "../kernels/tags.hpp"
#include "../solvers/build_response_ground_state.hpp"
#include "../solvers/convergence_policy.hpp"
#include "../solvers/fd_problem.hpp"
#include "../solvers/fd_save_load.hpp"
#include "../solvers/fd_solver.hpp"
#include "../solvers/iterate_protocol.hpp"
#include "../solvers/response_metadata.hpp"
#include "../solvers/response_state.hpp"

#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace molresponse_v3 {

/// Everything an FD protocol step solve needs beyond the node itself. References to
/// `world` / `gs` are borrowed — the executor lives within the solve scope.
struct ExecutorContext {
  madness::World   &world;
  GroundState      &gs;            // mutable: re-prepared at each protocol
  double            L = 0.0;       // cubic cell half-edge (from archive header)
  std::string       fock_json;     // moldft.fock.json path ("" -> none)
  ConvergencePolicy policy;
  PrintLevel        print_level = PrintLevel::Normal;
  std::string       calc_dir;      // holds response_metadata.json + archives
  int               max_iters = 25;
};

// ---------------------------------------------------------------------------
namespace detail_exec {

inline bool is_static_freq(double freq) { return std::abs(freq) < 1e-12; }

/// Alpha-spin perturbation source vector for one (pert) at the active
/// protocol. Throws for sources 15a does not build.
inline vector_real_function_3d
pert_source_alpha(madness::World &world, GroundState &gs,
                  const Perturbation &pert) {
  switch (pert.kind) {
    case Perturbation::Kind::Dipole:
      return dipole_perturbation(world, gs, pert.axis);
    case Perturbation::Kind::NuclearDisplacement:
      return nuclear_perturbation(world, gs, pert.atom, pert.axis);
    case Perturbation::Kind::Magnetic:
      throw std::runtime_error(
          "calc executor: magnetic perturbation source not implemented (15a)");
  }
  throw std::runtime_error("calc executor: unknown perturbation kind");
}

/// Beta-spin source (open shell). Only dipole has a beta builder today.
inline vector_real_function_3d
pert_source_beta(madness::World &world, GroundState &gs,
                 const Perturbation &pert) {
  if (pert.kind == Perturbation::Kind::Dipole)
    return dipole_perturbation_beta(world, gs, pert.axis);
  throw std::runtime_error(
      "calc executor: open-shell source for this perturbation not "
      "implemented (15a) — only dipole has a beta builder");
}

/// Fill the FD source carrier (== the solver Storage shape) for (Type, Shell).
template <typename Type, typename Shell, typename Carrier>
void fill_source(madness::World &world, GroundState &gs,
                 const Perturbation &pert, Carrier &src) {
  const auto pa = pert_source_alpha(world, gs, pert);
  src.x_alpha = pa;
  if constexpr (std::is_same_v<Type, Full>) src.y_alpha = pa;
  if constexpr (std::is_same_v<Shell, OpenShell>) {
    const auto pb = pert_source_beta(world, gs, pert);
    src.x_beta = pb;
    if constexpr (std::is_same_v<Type, Full>) src.y_beta = pb;
  }
}

/// Build the (Type, Shell)-correct GS-side kernel inputs.
template <typename Shell>
ResponseGroundState build_gs(madness::World &world, GroundState &gs) {
  if constexpr (std::is_same_v<Shell, ClosedShell>)
    return build_response_ground_state_closed_shell(
        world, gs, gs.hf_exchange_coefficient(), gs.params().lo());
  else
    return build_response_ground_state_open_shell(
        world, gs, gs.hf_exchange_coefficient(), gs.params().lo());
}

/// Re-project every in-flight component to a new wavelet basis. Generic over
/// (Type, Shell) via Storage::blocks().
template <typename State>
void reproject_state(State &st, int k, double thresh) {
  for (auto &ch : st.responses)
    for (auto *blk : ch.blocks())
      for (auto &fn : *blk) fn = madness::project(fn, k, thresh);
  for (auto &rho : st.rho_alpha_prev)
    rho = madness::project(rho, k, thresh);
}

} // namespace detail_exec

// ---------------------------------------------------------------------------
// One protocol step of one FD node.
// ---------------------------------------------------------------------------

/// Solve (pert, freq) at a single protocol `thresh`. `action` is the reconcile
/// verdict: Fresh starts from the perturbation guess (never loads, so a
/// diverged archive can't poison the solve); Restart/Resume load the nearest
/// converged-or-partial seed via try_load_fd_state and re-project it in the
/// first prepare(). Saves the result (+ metrics) through save_fd_state.
template <typename Type, typename Shell>
NodeResult solve_fd_protocol(ExecutorContext &ctx, const Perturbation &pert,
                         double freq, double thresh, NodeAction action) {
  using namespace madness;
  using Solver = FDSolver<Type, Shell>;
  World &world = ctx.world;
  GroundState &gs = ctx.gs;

  // Bring FunctionDefaults<3> + the ground state to this protocol so the
  // source builders and try_load (which keys on protocol_key()) are correct.
  set_response_protocol(world, ctx.L, thresh);
  const double t0 = FunctionDefaults<3>::get_thresh();
  {
    auto coulop = poperatorT(
        CoulombOperatorPtr(world, gs.params().lo(), 0.001 * t0));
    gs.prepare(world, 0.001 * t0, coulop, ctx.fock_json);
  }

  // Target at this protocol (rebuilt inside prepare for the ramp's single step).
  FDProblem<Type, Shell> tgt;
  tgt.gs = detail_exec::build_gs<Shell>(world, gs);
  tgt.responses.resize(1);
  tgt.responses[0].omega = freq;
  detail_exec::fill_source<Type, Shell>(world, gs, pert,
                                        tgt.responses[0].source);

  // Initial guess: the source carrier IS a valid Storage, so copy it.
  typename Solver::State s0;
  s0.responses.resize(1);
  s0.responses[0] = tgt.responses[0].source;

  // Seed from disk unless the reconcile said Fresh.
  if (action != NodeAction::Fresh) {
    auto loaded =
        try_load_fd_state<Type, Shell>(world, ctx.calc_dir, pert, freq);
    if (loaded) s0 = std::move(loaded->state);
  }

  Solver solver(world, tgt, ctx.policy, ctx.print_level);

  // The protocol + ground state + target are already set up above for `thresh`.
  // iterate_protocol calls prepare() before the (single) step; re-doing the
  // expensive gs.prepare()/build_gs() there would DOUBLE the per-step prep cost
  // (ground-state projection dominates — ~half of wall time per solve). So only
  // re-prepare if the protocol actually changed; otherwise just re-project the
  // (possibly coarser-loaded) state into the active basis.
  double prepared_t = t0;
  auto prepare = [&](double th, Solver &solv, typename Solver::State &st) {
    const double cur_t = FunctionDefaults<3>::get_thresh();
    const bool changed =
        std::abs(th - prepared_t) > 1e-15 * std::max(th, prepared_t) ||
        std::abs(cur_t - th)      > 1e-15 * std::max(cur_t, th);
    if (changed) {
      set_response_protocol(world, ctx.L, th);
      const double new_t = FunctionDefaults<3>::get_thresh();
      auto coulop = poperatorT(
          CoulombOperatorPtr(world, gs.params().lo(), 0.001 * new_t));
      gs.prepare(world, 0.001 * new_t, coulop, ctx.fock_json);

      FDProblem<Type, Shell> nt;
      nt.gs = detail_exec::build_gs<Shell>(world, gs);
      nt.responses.resize(1);
      nt.responses[0].omega = freq;
      detail_exec::fill_source<Type, Shell>(world, gs, pert,
                                            nt.responses[0].source);
      solv.set_target(std::move(nt));
      prepared_t = new_t;
    }
    detail_exec::reproject_state(st, FunctionDefaults<3>::get_k(),
                                 FunctionDefaults<3>::get_thresh());
  };

  // Real convergence = not diverged AND residuals within this protocol's
  // targets. A max-iters STALL must read as unconverged so it is neither
  // skipped on restart nor fed to property assembly as if good.
  auto converged_now = [](const typename Solver::State &st, const Solver &sv) {
    if (st.diverged) return false;
    double mb = 0.0, md = 0.0;
    for (double r : st.last_bsh_residual)     mb = std::max(mb, r);
    for (double r : st.last_density_residual) md = std::max(md, r);
    const auto &t = sv.targets();
    return mb <= t.bsh_residual && md <= t.density_residual;
  };

  auto post_step = [&](double, Solver &solv, typename Solver::State &st) {
    save_fd_state<Type, Shell>(world, st, ctx.calc_dir, pert, freq,
                               /*converged=*/converged_now(st, solv));
  };

  solvers::IterateProtocolPolicy pp;
  pp.max_iters_per_step = ctx.max_iters;
  const std::vector<double> one_protocol = {thresh};
  auto sf = solvers::iterate_protocol(solver, std::move(s0), one_protocol,
                                      prepare, post_step, pp);

  NodeResult r;
  r.converged = converged_now(sf, solver);
  r.reached_protocol_key = protocol_key();  // active defaults reflect this protocol step
  return r;
}

// ---------------------------------------------------------------------------
// FdResponseExecutor — dispatches the runtime (Type × Shell) for FD nodes.
// ---------------------------------------------------------------------------

class FdResponseExecutor : public ICalcExecutor {
public:
  explicit FdResponseExecutor(ExecutorContext ctx) : ctx_(ctx) {}

  NodeResult run_protocol(const WorkItem &item) override {
    const CalcNode &node = *item.node;

    if (node.kind == CalcKind::ES || node.kind == CalcKind::DerivedFD ||
        node.kind == CalcKind::VBC) {
      if (ctx_.world.rank() == 0) {
        madness::print("[CALC] run_protocol: node", node.id, "kind not handled by "
                       "the 15a FD executor (ES/DerivedFD/VBC) — skipping");
      }
      return NodeResult{/*converged=*/false, /*reached_protocol_key=*/"", {}};
    }

    const bool restricted = ctx_.gs.is_spin_restricted();
    const bool is_static   = detail_exec::is_static_freq(node.freq);

    if (ctx_.world.rank() == 0) {
      madness::print("[CALC] run_protocol: node", node.id,
                     " pert=", node.pert.description(),
                     " freq=", node.freq, " thresh=", item.thresh,
                     " type=", (is_static ? "static" : "full"),
                     " shell=", (restricted ? "closed" : "open"),
                     " action=", action_name(item.action));
    }

    if (restricted && is_static)
      return solve_fd_protocol<Static, ClosedShell>(ctx_, node.pert, node.freq,
                                                 item.thresh, item.action);
    if (restricted && !is_static)
      return solve_fd_protocol<Full, ClosedShell>(ctx_, node.pert, node.freq,
                                               item.thresh, item.action);
    if (!restricted && is_static)
      return solve_fd_protocol<Static, OpenShell>(ctx_, node.pert, node.freq,
                                               item.thresh, item.action);
    return solve_fd_protocol<Full, OpenShell>(ctx_, node.pert, node.freq,
                                          item.thresh, item.action);
  }

private:
  static const char *action_name(NodeAction a) {
    switch (a) {
      case NodeAction::Skip:    return "skip";
      case NodeAction::Restart: return "restart";
      case NodeAction::Resume:  return "resume";
      case NodeAction::Fresh:   return "fresh";
    }
    return "?";
  }
  ExecutorContext ctx_;
};

// ---------------------------------------------------------------------------
// CalcManager — top-level driver.
// ---------------------------------------------------------------------------

/// Seed-selection strategy (nearest converged; the knob is the metric, not a
/// fixed topology — see doc 15). Hoisted to namespace scope so it is a
/// complete type when used as a `= {}` default argument inside CalcManager.
enum class SeedStrategy { NearestConverged };
struct CalcManagerPolicy {
  SeedStrategy seed = SeedStrategy::NearestConverged;
  int          max_iters_per_step = 25;
};

// ===========================================================================
// Execution model (15a — single group). WHY all ranks must agree on the plan.
//
// A MADNESS Function is distributed across ALL ranks of the communicator, so
// solving one response state is a COLLECTIVE operation: every rank holds a
// piece and they collaborate on every step (BSH apply, fences, ...). 15a runs
// single-group — the whole communicator solves ONE state at a time; a "wave"
// of independent states is simply looped over (all ranks per state), NOT fanned
// out across ranks. Distributing a wave's states to rank SUBGROUPS sized to fit
// memory is the 15c design (STATE_PARALLEL / subworlds); the wave structure
// exists precisely so that later layer can partition it.
//
// Consequence — DETERMINISM. Every rank runs schedule() independently, then all
// ranks collectively solve waves.front()[i] together. If two ranks computed
// different waves they would issue mismatched collective calls and DEADLOCK.
// So the schedule must come out byte-identical on every rank: same
// response_metadata.json in (all ranks read the same file after a fence),
// deterministic logic out (sorted ramp + std::map perturbation order + pure
// reconcile). The std::map (not unordered_map) is a correctness requirement,
// not a style choice.
//
// run() in one sentence: forever — re-read what's done from disk, compute the
// next batch of independent work, solve it (all ranks together, one state at a
// time), save each; stop when nothing is left, or when a batch repeats
// unchanged (no progress possible).
// ===========================================================================
class CalcManager {
public:
  using Policy = CalcManagerPolicy;

  CalcManager(ResponsePlan plan, std::string calc_dir, Policy policy = {})
      : plan_(std::move(plan)), calc_dir_(std::move(calc_dir)),
        policy_(policy) {}

  /// Build the node DAG. n_atoms expands the nuclear_all sentinel.
  void build(int n_atoms) { dag_ = build_dag(plan_, n_atoms); }

  const std::vector<CalcNode> &dag() const { return dag_; }
  const std::string &calc_dir() const { return calc_dir_; }

  /// Drive every node to its target protocol against `exec`. One wave per
  /// pass: each pass reloads the aggregate metadata, re-schedules (so actions
  /// reflect current disk), runs the first non-empty wave, fences, then
  /// expands any newly-converged ES bundle into concrete DerivedFD nodes.
  /// Terminates when the schedule is empty (all Skip) or the same wave repeats
  /// with no progress (a node that cannot converge does not spin the loop).
  void run(madness::World &world, ICalcExecutor &exec) {
    const std::vector<double> ramp = global_ramp();
    const std::string meta_path = calc_dir_ + "/response_metadata.json";

    std::string last_sig;
    for (;;) {
      world.gop.fence();
      auto meta = ResponseMetadata::load_or_create(meta_path);
      auto waves = schedule(dag_, ramp, meta.json());
      if (waves.empty()) {
        if (world.rank() == 0)
          madness::print("[CALC] run: nothing left to schedule — done");
        break;
      }

      const std::string sig = wave_signature(waves.front());
      if (sig == last_sig) {
        if (world.rank() == 0)
          madness::print("[CALC] run: no progress on wave {", sig,
                         "} — stopping (unconverged or unhandled nodes)");
        break;
      }
      last_sig = sig;

      if (world.rank() == 0)
        madness::print("[CALC] run: wave of", (int)waves.front().size(),
                       "protocol step(s): {", sig, "}");
      for (const auto &item : waves.front()) {
        NodeResult res = exec.run_protocol(item);
        maybe_record_es_roots(item, res);
      }
      world.gop.fence();

      expand_converged_es(world);   // dormant until the ES executor lands
    }
  }

private:
  /// Coarse->fine union of every node's protocol ladder.
  std::vector<double> global_ramp() const {
    std::vector<double> all;
    for (const auto &n : dag_)
      for (double t : n.protocols) all.push_back(t);
    return detail_calc::sorted_ramp(std::move(all));
  }

  /// Stable, rank-deterministic signature of a wave (sorted id@thresh list).
  static std::string wave_signature(const std::vector<WorkItem> &wave) {
    std::vector<std::string> parts;
    parts.reserve(wave.size());
    for (const auto &it : wave)
      parts.push_back(it.node->id + "@" + ResponseMetadata::freq_key(it.thresh));
    std::sort(parts.begin(), parts.end());
    std::string s;
    for (const auto &p : parts) { s += p; s += ';'; }
    return s;
  }

  /// Stash ES root frequencies from a converged ES protocol step for expansion. (ES
  /// solves aren't wired in 15a, so this stays empty; kept so the ES executor
  /// only has to populate NodeResult::es_root_freqs.)
  void maybe_record_es_roots(const WorkItem &item, const NodeResult &res) {
    if (item.node->kind == CalcKind::ES && res.converged &&
        !res.es_root_freqs.empty()) {
      es_root_freqs_[item.node->id] = res.es_root_freqs;
    }
  }

  /// Turn each symbolic DerivedFD "*" whose ES bundle has converged into one
  /// concrete FD node per root (freq = root energy). Appends to dag_. No-op
  /// when no ES roots are known yet (15a). Idempotent via expanded_.
  void expand_converged_es(madness::World &world) {
    std::vector<CalcNode> additions;
    for (const auto &sym : dag_) {
      if (!sym.is_symbolic()) continue;
      if (expanded_.count(sym.id)) continue;
      // The symbolic node's single hard_dep is its ES bundle id.
      if (sym.prerequisites.empty()) continue;
      const std::string &es_id = sym.prerequisites.front();
      auto it = es_root_freqs_.find(es_id);
      if (it == es_root_freqs_.end()) continue;  // ES not converged yet
      int root = 0;
      for (double w : it->second) {
        CalcNode c;
        // Promote: a derived FD at ω = root energy is an ordinary FD point.
        // Solving it as CalcKind::FD lets the FD executor run it and the normal
        // skip / restart / nearest-frequency-seed logic apply uniformly.
        c.kind      = CalcKind::FD;
        c.pert      = sym.pert;
        c.freq      = w;
        c.protocols = sym.protocols;
        c.es_root_id = make_es_root_label(root++);  // provenance: source ES root
        c.id        = fd_node_id(c.pert, w);  // concrete -> fd_states subtree
        additions.push_back(std::move(c));
      }
      expanded_.insert(sym.id);
    }
    if (!additions.empty()) {
      if (world.rank() == 0)
        madness::print("[CALC] expanded", (int)additions.size(),
                       "derived FD node(s) from converged ES roots");
      for (auto &c : additions) dag_.push_back(std::move(c));
    }
  }

  static std::string make_es_root_label(int i) {
    char buf[24];
    std::snprintf(buf, sizeof buf, "es_root_%04d", i);
    return buf;
  }

  ResponsePlan          plan_;
  std::string           calc_dir_;
  Policy                policy_;
  std::vector<CalcNode> dag_;
  std::map<std::string, std::vector<double>> es_root_freqs_;
  std::set<std::string> expanded_;
};

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_CALC_CALC_EXECUTOR_HPP
