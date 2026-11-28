#ifndef MOLRESPONSE_V3_CALC_CALC_EXECUTOR_HPP
#define MOLRESPONSE_V3_CALC_CALC_EXECUTOR_HPP

// ===========================================================================
// CalcManager executor — the WORLD-BOUND half of the calc manager (doc 15,
// 15a). The scheduling core (calc_manager.hpp) decides what runs and in what
// order; this header actually drives the solves.
//
//   FdResponseExecutor : ICalcExecutor
//       run_protocol(WorkItem) -> one protocol step of a driven state, an
//       excited-state bundle or a VBC source, each through its node in the
//       response session (solvers/shared_response_session.hpp). Dispatches
//       Static/Full x Closed/Open at runtime for driven states.
//
//   CalcManager
//       build(n_atoms) -> dag from the plan; run(world, exec) -> the
//       one-wave-per-pass reschedule loop. Each pass reads what is saved from
//       the node folders (node_snapshot.hpp) and re-schedules from the
//       session's plan for every state.
//
// The assemble_* functions contract the qualified states of the session into
// the requested properties.
// ===========================================================================

#include "../solvers/shared_response_session.hpp"
#include "../solvers/response_property_metadata.hpp"
#include "calc_manager.hpp"
#include "node_snapshot.hpp"

#include "../GroundState.hpp"
#include "../Perturbations.hpp"
#include "../ResponseProtocol.hpp"
#include "../kernels/full.hpp"
#include "../kernels/static.hpp"
#include "../kernels/tags.hpp"
#include "../kernels/tpa.hpp"
#include "../solvers/build_response_ground_state.hpp"
#include "../solvers/convergence_policy.hpp"
#include "../solvers/fd_problem.hpp"
#include "../solvers/dalton_import.hpp"   // dalton_fd_guess (seed.freq_tol)
#include "../kernels/beta.hpp"
#include "../kernels/tpa_source_spec.hpp"
#include "../solvers/es_analysis.hpp"
#include "../solvers/es_seed_guard.hpp"
#include "../solvers/response_seed.hpp"
#include "../solvers/es_solver.hpp"
#include "../solvers/fd_solver.hpp"
#include "../solvers/iterate_protocol.hpp"
#include "../solvers/node_subworlds.hpp"
#include "../solvers/response_metadata.hpp"
#include "../solvers/response_state.hpp"
#include "../solvers/saved_states.hpp"

#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace molresponse_v3 {

/// Executor configuration carried by ResponseWorkflowInput. Every response
/// run has a session (shared_response); it owns the nodes of all states.
struct ExecutorSettings {
  std::shared_ptr<response_shared::Session> shared_response;
  // The session's plan for every state (D-016).
  StatePlans        state_plans;
  ConvergencePolicy policy;
  PrintLevel        print_level = PrintLevel::Normal;
  // Guards the shared print functions (IterationOutput.h); from the same print_level.
  madness::iteration::Verbosity verbosity = madness::iteration::Verbosity::detail;
  std::string       calc_dir;      // holds response_metadata.json (output only)
  // DALTON nearest-frequency seed for FD legs absent from the RSPVEC (the
  // derived two-photon legs at omega_f/2). Set from the deck (`dalton.dir`,
  // `seed.freq_tol`) or --dalton-dir/--seed-freq-tol. 0 = off (2026-09-09).
  std::string       dalton_dir;
  double            seed_freq_tol         = 0.0;
  // Deck `diagnostics.properties`: cheap property values in each history row
  // (polarizability rows of dipole states against these axes, excitation
  // energies). Information only; never a result, never a second-order step.
  bool              diagnostic_properties = false;
  std::vector<int>  diagnostic_dipole_axes;
  // Excited-state (Full / TDA-warmup) solve settings — defaults for the Full
  // closed-shell path (random guess, 10 warmup iters, oversampled warmup, KAIN).
  ESGuessMode       es_guess              = ESGuessMode::SolidHarmonics;  // sweep-validated default
  // AO basis projected for the VirtualAO guess (CLI --es-guess-basis / deck
  // response.excited.guess_basis). Radial rank per l-sector = #shells(l) −
  // #occupied(l), so a larger basis reaches further up the radial ladder
  // (Rydberg series want d-aug-cc-pvtz or similar). Ignored by other modes.
  std::string       es_guess_basis        = "aug-cc-pvdz";
  int               es_tda_warmup_iters   = 10;
  // Warmup oversampling: keep the lowest n_roots of ceil(factor*n_roots)
  // partially-converged warmup trials. 3.0 (was 2.0) because on the PRODUCTION
  // ladder (moldft climbing to 1e-6/k8) the extra cut-line margin lets the
  // downselect keep h2o's true 4th state 0.4097 au and match Dalton
  // d-aug-cc-pVQZ to 0.03% (ladder job 2084168), vs a 2.4% missing-state error
  // at 2.0 (job 2083837, which reported the 5th state 0.420 in its place).
  // KNOWN LIMITATION (follow-up): the underlying cause is the TDA warmup
  // MISORDERING the near-degenerate 0.40/0.41 cluster at the coarse rung — a
  // 1e-4/k6-only run still misorders it (the true 4th sits just above the
  // 4-root cut) and does not fully tighten that root in 25 iters. Oversampling
  // only buys margin, it can't fix an ordering inversion at the cut; the real
  // fix is an energy-ordered guess (ESGuessMode::VirtualAO, the NWChem
  // CIS-diagonal path) or more warmup iters. Cost of 3.0 is warmup-only
  // (coarse rung, es_tda_warmup_iters iters).
  double            es_warmup_oversample  = 3.0;
  int               es_kain_maxsub        = 8;
  double            es_maxrotn            = 0.5;
  // KAIN in the MAIN ES solve (deck `kain`, when set). The warmup is KAIN-free
  // regardless: its whole length is the tda_warmup_iters no-KAIN window.
  bool              es_kain               = true;
  // Delay KAIN onset in the MAIN ES solve by this many iters (pure BSH +
  // step-restriction first, so roots stabilize before KAIN starts recording
  // history). 0 = KAIN from iter 1 (previous behaviour). Distinct from
  // es_tda_warmup_iters, which is the separate oversampled-warmup pre-pass.
  int               es_main_kain_delay    = 0;
  // Full-deflation locking of converged ES roots (workstream B). ON by default
  // (sweep-validated: required to converge near-degenerate h2o/c2h4 clusters).
  bool              es_lock_converged     = true;
  // D-122: the deck wrote kain.min_residual (policy.kain_min_residual is then
  // the written value, which wins over the cold-start rule below).
  bool              es_kain_min_residual_written = false;
  // D-133: the bundle's plateau detector (deck stall.window / stall.ratio).
  // Off unless written; the driven legs' policy keeps it off.
  int               es_stall_window       = 0;
  double            es_stall_ratio        = 0.10;
  // Inc-3c tensor-layer ES γ (--es-tensor; TDA/ClosedShell only): the bundle's
  // exchange runs off the per-protocol cached g0 tensor + per-root contractions
  // (ESSolver::set_gamma_tensor / tda_batch::compute_gamma_flat) instead of one
  // Exchange operator per root per iter. A/B-to-floor vs the reference (~1e-3
  // rel in θ → ~1e-6 in a converged ω), NOT bitwise. Off by default.
  bool              es_gamma_tensor   = false;
  // F2 (doc 32 §5): fan independent states out across subworlds. This is the
  // number of subworlds PER PHYSICAL NODE (F2f): 0 = single-World reference
  // (byte-identical); 1 = one subworld per node (node-aligned); >1 = sub-node /
  // NUMA packing (φ replicated P× per node, more state-parallelism — small-system
  // regime). Total subworlds G = n_nodes × this. G≤1 short-circuits to the G=0
  // path. Fans FD / NuclearFD / VBC (F2g); ES stays single-World.
  int               fd_subworlds = 0;
  // Large-system regime (2026-09-10): ranks PER SUBWORLD. > 0 replaces the
  // per-node packing above with a universe-level contiguous split into
  // G = min(items, universe/fd_subworld_ranks) subworlds that may span nodes
  // (valinomycin, 300 MOs: one SCF wants 6-8 nodes x 8 ranks, so each state
  // needs a 48-64-rank subworld). 0 = per-node packing (fd_subworlds).
  int               fd_subworld_ranks = 0;
  // F2d (doc 32 §5.6): per-subworld log attribution. Empty in single-World
  // mode ⇒ output byte-identical to the G=0 reference. Set by the fan-out
  // closure. A non-empty log_prefix is ALSO the "I am in a
  // subworld" signal that suppresses the redundant GS/protocol banners (F2d-ii).
  std::string       log_prefix;     // "[g2/3 xm005] " prepended to per-state human lines
};

/// Everything an FD protocol step solve needs beyond the node itself. The
/// World-bound part (`world` / `gs` borrowed — the executor lives within the
/// solve scope) plus the plain `ExecutorSettings` (inherited). The orchestrator
/// builds this from a loaded GroundState + the Input's settings.
struct ExecutorContext : ExecutorSettings {
  madness::World   &world;
  GroundState      &gs;            // mutable: re-prepared at each protocol
  double            L = 0.0;       // cubic cell half-edge (from archive header)
  std::string       fock_json;     // the ground-state node's Fock file ("" -> none)
  bool final_protocol = true; // role of the current work item in its full ladder

  ExecutorContext(madness::World &w, GroundState &g, double L_,
                  std::string fock, ExecutorSettings s = {})
      : ExecutorSettings(std::move(s)), world(w), gs(g), L(L_),
        fock_json(std::move(fock)) {}
};

// ---------------------------------------------------------------------------
namespace detail_exec {

/// Collective: rank 0 builds the scheduler's snapshot with `read` (run on
/// rank 0 only) and broadcasts it, so every rank schedules from the same
/// bytes rather than from whatever its own filesystem view shows after the
/// fence. A rank-0 failure is broadcast too, so all ranks throw together
/// instead of the others waiting in the broadcast.
inline nlohmann::json snapshot_rank0(madness::World &world,
                                     const std::function<nlohmann::json()> &read) {
  std::string text, err;
  if (world.rank() == 0) {
    try {
      text = read().dump();
    } catch (const std::exception &e) {
      err = e.what();
      if (err.empty()) err = "unknown error";
    }
  }
  world.gop.broadcast_serializable(err, 0);
  if (!err.empty()) throw std::runtime_error(err);
  world.gop.broadcast_serializable(text, 0);
  return nlohmann::json::parse(text);
}

/// The ES warmup and main-solve policies, from the caller's settings. Shared by
/// the TDA and Full executors. The warmup runs es_tda_warmup_iters iterations
/// inside a no-KAIN window of the same length, so it is KAIN-free whatever
/// `kain` says; the main solve delays KAIN by es_main_kain_delay.
/// D-122: the KAIN hold-off (kain_min_residual) is for a main solve that
/// starts from a seed or a saved bundle, whose first residual can make KAIN
/// extrapolate garbage. A main solve that starts from this solver's own TDA
/// warm-up (`own_warmup`: a fresh warm-up or the cached one) gets none: its
/// plain BSH step need not contract (C2H6 cold, 2026-11-04: max_res grew from
/// 0.11 every iteration, KAIN never engaged, the roots left for the
/// continuum). A written kain.min_residual wins in both cases.
struct EsIterationPolicies {
  ConvergencePolicy warm;
  ConvergencePolicy main;
};
inline EsIterationPolicies es_iteration_policies(const ExecutorSettings &s,
                                                 bool own_warmup = false) {
  EsIterationPolicies p;
  p.warm                          = s.policy;
  p.warm.kain                     = true;
  p.warm.kain_maxsub              = s.es_kain_maxsub;
  p.warm.maxrotn                  = s.es_maxrotn;
  p.warm.tda_warmup_iters         = s.es_tda_warmup_iters;
  p.warm.warmup_oversample_factor = s.es_warmup_oversample;
  p.warm.stall_window             = s.es_stall_window;
  p.warm.stall_ratio              = s.es_stall_ratio;
  p.main                  = p.warm;
  p.main.kain             = s.es_kain;
  p.main.tda_warmup_iters = s.es_main_kain_delay;
  p.main.lock_converged   = s.es_lock_converged;  // warmup never locks
  if (own_warmup && !s.es_kain_min_residual_written)
    p.main.kain_min_residual = 0.0;
  return p;
}

/// D-122: the main solve's KAIN hold-off and why, at Verbose.
inline void print_es_kain_hold_off(madness::World &world,
                                   const ConvergencePolicy &main,
                                   const ExecutorSettings &s, bool own_warmup) {
  if (world.rank() == 0 && s.print_level >= PrintLevel::Verbose)
    madness::print("[KAIN] main solve hold-off kain_min_residual =",
                   main.kain_min_residual,
                   s.es_kain_min_residual_written ? "(written in the deck)"
                   : own_warmup ? "(none: starts from its own TDA warm-up)"
                                : "(starts from a seed or a saved bundle)");
}


inline bool is_static_freq(double freq) { return std::abs(freq) < 1e-12; }

/// Rank 0 only. The scheduler's snapshot of the session's run folder: every
/// completely saved node under states/response/, read without changing it,
/// each judged against the current request (Session::node_meets_request).
inline nlohmann::json session_snapshot(const response_shared::Session &session) {
  const auto saved = saved_states::completed_response_nodes(session.root());
  std::map<std::string, madness::iteration::NodeDescription> by_name;
  for (const auto &node : saved)
    by_name.emplace(node.folder.filename().string(), node.description);
  std::vector<SnapshotNode> nodes;
  nodes.reserve(saved.size());
  for (const auto &[name, d] : by_name)
    nodes.push_back({name, d, session.node_meets_request(name, d, by_name)});
  return node_snapshot(nodes);
}

/// The saved state another node holds, as the starting guess of a shared
/// driven state with no node of its own. A static state (X only) becomes
/// (X, Y = X) for a dynamic target; a dynamic state gives a static target
/// its X. The other node's archive is only read.
template <typename Shell, typename Storage, typename Binding>
Storage load_start(madness::World &world, Binding &binding) {
  const auto layout = binding.source_settings().at("layout").template get<std::string>();
  const bool source_static = layout.ends_with("_static");
  constexpr bool target_static = std::is_same_v<Storage, ResponseStateX<Shell>>;
  if (source_static == target_static) return binding.template load<Storage>();
  Storage result;
  if constexpr (target_static) {
    auto xy = binding.template load<ResponseStateXY<Shell>>();
    result.x_alpha = std::move(xy.x_alpha);
    if constexpr (std::is_same_v<Shell, OpenShell>) result.x_beta = std::move(xy.x_beta);
  } else {
    auto x = binding.template load<ResponseStateX<Shell>>();
    result.x_alpha = madness::copy(world, x.x_alpha);
    result.y_alpha = std::move(x.x_alpha);
    if constexpr (std::is_same_v<Shell, OpenShell>) {
      result.x_beta = madness::copy(world, x.x_beta);
      result.y_beta = std::move(x.x_beta);
    }
  }
  return result;
}

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

/// Solve (pert, freq) at a single protocol `thresh` through the session's
/// node for the state: it starts from its own saved state, from the nearest
/// saved frequency, from an imported guess, or from zero, and the session
/// decides when it stops and saves.
template <typename Type, typename Shell>
NodeResult solve_fd_protocol(ExecutorContext &ctx, const Perturbation &pert,
                             double freq, double thresh) {
  using namespace madness;
  using Solver = FDSolver<Type, Shell>;
  World &world = ctx.world;
  GroundState &gs = ctx.gs;
  {
    auto& timing=ctx.shared_response->timing();
    timing.begin("protocol",{{"protocol",checkpoint::records::protocol_record({thresh,default_k_for_thresh(thresh)})}});
    timing.begin("task",{{"task_id",fd_node_id(pert,freq)}});
    timing.begin("prepare");
  }

  // Bring FunctionDefaults<3> + the ground state to this protocol so the
  // source builders and try_load (which keys on protocol_key()) are correct.
  const bool gs_verbose = ctx.log_prefix.empty();  // F2d-ii: quiet in subworlds
  set_response_protocol(world, ctx.L, thresh, /*override_k=*/-1, gs_verbose);
  const double t0 = FunctionDefaults<3>::get_thresh();
  {
    auto coulop = poperatorT(
        CoulombOperatorPtr(world, gs.params().lo(), 0.001 * t0));
    gs.prepare(world, 0.001 * t0, coulop, ctx.fock_json, gs_verbose);
  }

  // Target at this protocol (rebuilt inside prepare for the ramp's single step).
  FDProblem<Type, Shell> tgt;
  tgt.gs = detail_exec::build_gs<Shell>(world, gs);
  tgt.responses.resize(1);
  tgt.responses[0].omega = freq;
  detail_exec::fill_source<Type, Shell>(world, gs, pert,
                                        tgt.responses[0].source);

  // Initial guess: zero, replaced below by a saved or imported state.
  typename Solver::State s0;
  s0.responses.resize(1);
  // Fresh guess x0 = 0 (matches molresponse_v2, which zero-allocates the response):
  // the first BSH step then gives G(-2*V_p*phi), the bounded uncoupled first-order
  // response. Seeding x0 = the source directly blows up iter-1 for a large/peaked
  // source (nuclear displacement) — V0*x0 and the BSH shift*x0 are O(||source||),
  // which for nuclear exceeds the explosion guard. Same converged fixed point
  // either way (linear response is unique), so alpha/beta values are unchanged;
  // this only fixes the transient / divergence.
  // DEEP-copy the source's shape then zero it. A shallow copy (operator=) shares
  // the Function handles, so scaling it would also zero the perturbation source
  // itself -> no driving term -> x stays 0 (the "converged in 1 iter, res=0,
  // value=0" bug). madness::copy gives an independent buffer.
  s0.responses[0] = tgt.responses[0].source;    // shape (shallow)
  {
    auto zflat = madness::copy(world, s0.responses[0].flatten());
    madness::scale(world, zflat, 0.0);
    s0.responses[0].from_flat(zflat);           // independent zeros
  }
  namespace it=madness::iteration;
  // With no saved state of its own, the state starts from the nearest
  // frequency saved now (section 2 of step 4d).
  auto& binding=ctx.shared_response->start_nearby(fd_node_id(pert,freq));
  ctx.shared_response->timing().end();
  binding.activate({t0,FunctionDefaults<3>::get_k()});
  // Suppress the legacy stopping table; report the shared controller's checks.
  Solver solver(world,std::move(tgt),ctx.policy,PrintLevel::Silent,ctx.log_prefix);
  solver.set_print_level(ctx.print_level);
  if (binding.has_saved()) {
    it::ElapsedTimer<> load_timer;
    const auto started_from=binding.started_from();
    if (started_from) s0.responses.front()=detail_exec::load_start<Shell,typename Solver::Storage>(world,binding);
    else s0.responses.front()=binding.template load<typename Solver::Storage>();
    const double load_seconds=load_timer.seconds();
    ctx.shared_response->timing().append({{"name","load"},{"seconds",load_seconds}});
    if (started_from) {
      // Only a starting guess: no measured value or update count is carried.
      if (world.rank()==0 && it::prints(ctx.verbosity,it::needed::start_source))
        it::write_start_source(std::cout,pert.description(),freq,
                               binding.source_settings().at("frequency").template get<double>());
    } else if (world.rank()==0 && it::prints(ctx.verbosity,it::needed::restart_summary)) {
      const auto saved=binding.saved_metrics();
      it::write_restart_summary(std::cout,"response",saved,binding.controls().protocol,
          it::assess(binding.controls(),binding.expected(),saved),binding.reusable(),response_shared::metric_labels());
    }
    if (world.rank()==0 && it::prints(ctx.verbosity,it::needed::checkpoint_timing))
      it::write_checkpoint_timing(std::cout,"response",binding.info().task_id,load_seconds);
  }
  // With no saved state to start from, a seed of this node (seeds/), else the
  // nearest Dalton record computed now (seed.freq_tol), else the zero guess.
  std::optional<response_seed::Seed> seed;
  if(!binding.has_saved())
    seed=response_seed::find_seed(world,ctx.shared_response->root(),binding.name(),
        response_seed::seed_settings(binding.scientific_settings()),ctx.verbosity);
  if(!binding.has_saved() && (seed || !ctx.dalton_dir.empty())) {
    auto& timing=ctx.shared_response->timing();timing.begin("imported_guess");
    bool imported=false;
    if(seed) {
      s0.responses.front()=response_seed::load_driven_seed<typename Solver::Storage>(world,*seed);
      s0.iter=0;imported=true;
      binding.note_guess(seed->guess());
      if(world.rank()==0 && it::prints(ctx.verbosity,it::needed::start_source))
        it::write_seed_start(std::cout,binding.name(),seed->shown());
    }
    if constexpr(std::is_same_v<Type,Full> && std::is_same_v<Shell,ClosedShell>) {
      if(!imported && !ctx.dalton_dir.empty() && ctx.seed_freq_tol>0.0 && pert.kind==Perturbation::Kind::Dipole)
        if(auto guess=dalton_fd_guess(world,gs,ctx.dalton_dir,ctx.calc_dir,pert.axis,freq,ctx.seed_freq_tol)) {
          s0.responses.front()=std::move(*guess);imported=true;
          binding.note_guess({{"kind","dalton"},{"directory",ctx.dalton_dir}});
        }
    }
    if(imported && world.rank()==0) print("response: imported driven guess; fresh convergence checks");
    timing.end();
  }
  std::size_t evaluations=0;
  if (!binding.reusable()) {
    // Diagnostic polarizability row: the assembly's own inner products, at
    // this protocol, from the state just updated. One per requested axis.
    std::map<int,vector_real_function_3d> dipoles;
    if constexpr (std::is_same_v<Shell,ClosedShell>)
      if (ctx.diagnostic_properties && pert.kind==Perturbation::Kind::Dipole)
        for (int column:ctx.diagnostic_dipole_axes) dipoles[column]=dipole_perturbation(world,gs,column);
    const auto observe=[&](const typename Solver::State& state,const it::IterationMetrics& metrics,
                           const it::Decision& decision,const it::Progress& progress) {
      ++evaluations;
      std::vector<it::Diagnostic> diagnostics;
      if constexpr (std::is_same_v<Shell,ClosedShell>)
        for (const auto& [column,v]:dipoles) {
          const auto& x=state.responses.front();
          double value=0.0;
          if constexpr (requires {x.y_alpha;}) value=-2.0*(inner(x.x_alpha,v)+inner(x.y_alpha,v));
          else value=-2.0*(inner(x.x_alpha,v)+inner(x.x_alpha,v));
          diagnostics.push_back({"polarizability",metrics.states.front().id,
                                 std::string{char('x'+pert.axis),char('x'+column)},value});
        }
      if (world.rank()==0 && it::prints(ctx.verbosity,it::needed::iteration_metrics))
        it::write_iteration_metrics(std::cout,"response",metrics,decision,
            response_shared::metric_labels(),binding.final_level());
      binding.record_evaluation(metrics,decision,diagnostics);
      binding.publish(metrics,decision,progress,[&](const std::filesystem::path& directory) {
        MADNESS_CHECK_THROW(state.responses.size()==1,"driven task must contain one response state");
        state.responses.front().save(world,madness::native_archive_name((directory/"response").string()),IoBackend::Native);
      });
      if (decision.checkpoint && world.rank()==0 && it::prints(ctx.verbosity,it::needed::checkpoint_saved))
        it::write_checkpoint_saved(std::cout,decision);
      if (world.rank()==0 && decision.action==it::Action::stop && decision.reason=="iteration_limit" &&
          it::prints(ctx.verbosity,it::needed::unconverged_diagnostics))
        it::write_unconverged_diagnostics(std::cout,metrics,decision,diagnostics,response_shared::metric_labels());
    };
    const auto decision=solvers::iterate_shared(solver,s0,binding.controls(),binding.final_level(),binding.expected(),
        binding.info(),binding.previous_updates(),observe,
        [&](const it::IterationInfo& info,const it::IterationTiming& timing) {
          ctx.shared_response->record_iteration(info,timing);
          if (world.rank()==0 && it::prints(ctx.verbosity,it::needed::iteration_timing))
            it::write_iteration_timing(std::cout,"response",info,timing);
        });
    MADNESS_CHECK_THROW(decision.checkpoint,"driven response stopped without a valid current checkpoint");
  }
  ctx.shared_response->record_summary(binding,evaluations);
  ctx.shared_response->timing().end();
  ctx.shared_response->timing().end();
  ctx.shared_response->publish_timings();
  return {binding.complete(),protocol_key(),{}};
}

// ---------------------------------------------------------------------------
// Excited-state bundle solve (closed shell).
// ---------------------------------------------------------------------------

/// Solve the closed-shell excited-state bundle of `n_roots` at one protocol
/// step through the session's node for the set.
// Shared ES binding; the existing task DAG still schedules the whole bundle.
template<class Type>
NodeResult solve_shared_excited(ExecutorContext& ctx,int n_roots,double thresh) {
  using namespace madness;
  namespace it=madness::iteration;
  using Solver=ESSolver<Type,ClosedShell>;
  auto& world=ctx.world;auto& gs=ctx.gs;auto& session=*ctx.shared_response;
  auto& binding=session.excited_at(es_node_id(std::is_same_v<Type,TDA>,n_roots));
  auto& timing=session.timing();
  timing.begin("protocol",{{"protocol",checkpoint::records::protocol_record({thresh,default_k_for_thresh(thresh)})}});
  timing.begin("task",{{"task_id",binding.info().task_id}});timing.begin("prepare");
  set_response_protocol(world,ctx.L,thresh);
  binding.activate({thresh,FunctionDefaults<3>::get_k()});
  auto coulop=poperatorT(CoulombOperatorPtr(world,gs.params().lo(),0.001*thresh));
  gs.prepare(world,0.001*thresh,coulop,ctx.fock_json);
  ESProblem<Type,ClosedShell> problem{detail_exec::build_gs<ClosedShell>(world,gs),n_roots};
  auto policy=detail_exec::es_iteration_policies(ctx,!binding.has_saved()).main;
  policy.lock_converged=false;policy.tda_warmup_iters=0;
  Solver solver(world,std::move(problem),policy,PrintLevel::Silent);
  solver.set_print_level(ctx.print_level);solver.set_gamma_tensor(ctx.es_gamma_tensor);
  typename Solver::State state;
  timing.end();
  if(binding.has_saved()) {
    timing.begin("load");
    state=binding.load([&](const auto& directory,const auto& description) {
      return response_shared::load_excited_functions<Type,ClosedShell>(world,directory,description,binding.controls().protocol);
    });
    timing.end();
    if(world.rank()==0 && it::prints(ctx.verbosity,it::needed::restart_summary)) {
      const auto saved=binding.saved_metrics();
      it::write_restart_summary(std::cout,"response",saved,binding.controls().protocol,
          it::assess(binding.controls(),binding.expected(),saved),binding.reusable(),response_shared::excited_metric_labels());
    }
  }
  // With no saved set, the seed of this set (seeds/) when it holds as many roots.
  std::optional<response_seed::Seed> imported_seed;
  if(!binding.has_saved()) {
    const auto name=std::filesystem::path(binding.folder()).filename().string();
    auto seed=response_seed::find_seed(world,session.root(),name,
        response_seed::seed_settings(binding.scientific_settings()),ctx.verbosity);
    if(seed && seed->description.details.at("n_roots").template get<int>()!=n_roots) {
      if(world.rank()==0 && it::prints(ctx.verbosity,it::needed::seed_refused))
        print("Seed ignored:",seed->shown(),"holds",seed->description.details.at("n_roots").template get<int>(),
              "roots; this set has",n_roots);
      seed.reset();
    }
    if(seed || !ctx.dalton_dir.empty()) {
      timing.begin("imported_guess");
      if(seed) {
        state=response_seed::load_excited_seed<Type,ClosedShell>(world,*seed);
        binding.note_guess(seed->guess());
        if(world.rank()==0) print("response: imported excited bundle; fresh convergence checks");
        if(world.rank()==0 && it::prints(ctx.verbosity,it::needed::start_source))
          it::write_seed_start(std::cout,name,seed->shown());
        imported_seed=std::move(seed);
      }
      timing.end();
    }
  }
  if(state.roots.size()<static_cast<std::size_t>(n_roots)) {
    timing.begin("guess");
    // Preserve the legacy candidate space, but select its lowest Ritz roots
    // without the retired warmup updates. A single raw solid-harmonic trial
    // can otherwise pin a higher symmetry sector (e.g. transverse H2).
    const int trials=std::max(n_roots,static_cast<int>(std::ceil(ctx.es_warmup_oversample*n_roots)));
    auto tda=build_initial_guess_tda_closed_shell(world,gs,trials,ctx.es_guess,ctx.es_guess_basis);
    ESSolver<TDA,ClosedShell> guess_solver(world,
        build_es_problem_tda<ClosedShell>(world,gs,trials,gs.hf_exchange_coefficient(),gs.params().lo()),
        policy,PrintLevel::Silent);
    guess_solver.set_gamma_tensor(ctx.es_gamma_tensor);
    tda=guess_solver.prepare_initial_guess(std::move(tda));
    if(state.roots.empty()) tda=slice_state_lowest<ClosedShell>(tda,n_roots);
    else {
      // Added roots: take the new guesses from outside the saved roots' span.
      std::vector<ResponseStateX<ClosedShell>> saved_x(state.roots.size());
      for(std::size_t i=0;i<state.roots.size();++i) saved_x[i].x_alpha=state.roots[i].x_alpha;
      tda=select_roots_outside_saved(world,tda,saved_x,n_roots);
    }
    if(world.rank()==0 && ctx.print_level>=PrintLevel::Normal)
      print("response: selected",n_roots-static_cast<int>(state.roots.size()),"initial roots from",trials,"trials; no warmup updates");
    typename Solver::State guess;
    if constexpr(std::is_same_v<Type,TDA>) guess=std::move(tda);
    else guess=promote_tda_to_full_closed_shell(world,tda);
    // Retain every saved root; append only the newly requested guess functions.
    for(std::size_t i=0;i<state.roots.size();++i) guess.roots.at(i)=std::move(state.roots[i]);
    guess.stable_index=state.stable_index;
    for(std::size_t i=guess.stable_index.size();i<guess.roots.size();++i)
      guess.stable_index.push_back(static_cast<int>(i));
    state=std::move(guess);timing.end();
  }
  std::size_t evaluations=0;
  std::optional<double> guard_seconds;
  if(!binding.reusable()) {
    const auto observer=[&](const auto& current,const it::IterationMetrics& metrics,const it::Decision& decision,const it::Progress& progress) {
      ++evaluations;
      if(world.rank()==0 && it::prints(ctx.verbosity,it::needed::iteration_metrics))
        it::write_iteration_metrics(std::cout,"response",metrics,decision,response_shared::excited_metric_labels(),binding.final_level());
      // Diagnostic excitation energies: already measured, nothing computed.
      std::vector<it::Diagnostic> diagnostics;
      if(ctx.diagnostic_properties)
        for(const auto& root:metrics.states)
          if(const auto found=root.metrics.find(it::Metric::energy_change);
             found!=root.metrics.end() && found->second.evidence.current)
            diagnostics.push_back({"excitation_energy",root.id,"",*found->second.evidence.current});
      // One row in the node's history per decision, before any save can fail.
      binding.record_evaluation(metrics,decision,diagnostics);
      // A set that started from a seed: at the level's last decision, how
      // closely each root still matches its seed root. It is saved with the
      // roots it describes.
      nlohmann::json facts=nlohmann::json::object();
      if(imported_seed && decision.checkpoint && decision.action!=it::Action::update) {
        it::ElapsedTimer<> guard_timer;
        const auto guard=evaluate_es_seed_guard<Type,ClosedShell>(world,current,decision.action==it::Action::complete,
            [&] {return response_seed::load_excited_seed<Type,ClosedShell>(world,*imported_seed);},
            imported_seed->shown(),protocol_key(imported_seed->description.protocol.thresh,imported_seed->description.protocol.k));
        guard_seconds=guard_timer.seconds();
        print_es_seed_guard(world,guard);
        facts["seed_guard"]=guard.to_json();
        timing.append({{"name","seed_comparison"},{"seconds",*guard_seconds}});
      }
      binding.publish(metrics,decision,progress,current.stable_index,[&](const auto& directory) {
        response_shared::save_excited_functions(world,current,directory);
      },facts);
      if(decision.checkpoint && world.rank()==0 && it::prints(ctx.verbosity,it::needed::checkpoint_saved))
        it::write_checkpoint_saved(std::cout,decision);
      if(world.rank()==0 && decision.action==it::Action::stop && decision.reason=="iteration_limit" &&
         it::prints(ctx.verbosity,it::needed::unconverged_diagnostics))
        it::write_unconverged_diagnostics(std::cout,metrics,decision,diagnostics,response_shared::excited_metric_labels());
    };
    const auto decision=solvers::iterate_shared(solver,state,binding.controls(),binding.final_level(),binding.expected(),binding.info(),binding.previous_updates(),observer,
      [&](const it::IterationInfo& info,const it::IterationTiming& measured) {
        session.record_iteration(info,measured);
        if(world.rank()==0 && it::prints(ctx.verbosity,it::needed::iteration_timing)) it::write_iteration_timing(std::cout,"response",info,measured);
      });
    MADNESS_CHECK_THROW(decision.checkpoint,"excited bundle stopped without a valid current checkpoint");
  }
  session.record_summary(binding,evaluations,guard_seconds);
  timing.end();timing.end();session.publish_timings();
  return {binding.complete(),protocol_key(),{}};
}

template<class Type>
void assemble_shared_excited(ExecutorContext& ctx,int n_roots,double thresh) {
  using namespace madness;
  using json=nlohmann::json;
  auto& session=*ctx.shared_response;
  auto& binding=session.excited_at(es_node_id(std::is_same_v<Type,TDA>,n_roots));
  // One item per root holds every excited-state value of that root. They all
  // come from the root's own state, so the item has one convergence: the
  // root's metrics. A value that was not computed is absent and named in
  // availability; status is decided here, from what was computed.
  json excited={{"value_type","real"},
                {"units",{{"excitation_energy","hartree"},{"oscillator_strength","dimensionless"},
                          {"transition_dipole","atomic_units"},{"transition_quadrupole","atomic_units"}}},
                {"definitions",{{"transition_quadrupole","nontraceless second moment"}}},
                {"values",json::array()}};
  const auto set_id=es_node_id(std::is_same_v<Type,TDA>,n_roots);
  const auto set_name=session.node_name(set_id);
  // A set that started from a seed: the guard saved with its roots, per root.
  std::map<std::string,json> seed_match;
  if(binding.has_saved() && binding.description().details.contains("seed_guard")) {
    excited["units"]["seed_match"]={{"omega_seed","hartree"},{"omega_shift","hartree"},
                                    {"seed_overlap","dimensionless"},{"best_overlap","dimensionless"}};
    excited["definitions"]["seed_match"]=
        "the root against the seed it started from: seed_overlap with the seed root of the same identity "
        "(seed_slot), best_overlap with any seed root (best_seed_slot), both normalized in the response metric; "
        "omega_shift = excitation_energy - omega_seed";
    for(const auto& root:binding.description().details.at("seed_guard").at("roots")) {
      json match=json::object();
      for(const char* field:{"seed_slot","omega_seed","omega_shift","seed_overlap","best_seed_slot","best_overlap"})
        if(root.contains(field)) match[field]=root.at(field);
      seed_match[root.at("root_id").template get<std::string>()]=std::move(match);
    }
  }
  const auto item=[&](const std::string& root,const json& values,const json& missing) {
    json fields=values;
    fields["root"]=root;
    fields["availability"]=missing;
    fields["status"]=values.empty() ? "unavailable" : missing.empty() ? "complete" : "partial";
    if(const auto match=seed_match.find(root);match!=seed_match.end()) fields["seed_match"]=match->second;
    session.name_state(fields,set_id,root);
    return fields;
  };
  if(binding.complete()) {
    session.timing().begin("excited_properties",{{"task_id",binding.info().task_id}});
    set_response_protocol(ctx.world,ctx.L,thresh);
    auto coulop=poperatorT(CoulombOperatorPtr(ctx.world,ctx.gs.params().lo(),0.001*thresh));
    ctx.gs.prepare(ctx.world,0.001*thresh,coulop,ctx.fock_json);
    session.timing().begin("load");
    const auto state=binding.load([&](const auto& directory,const auto& description) {
      return response_shared::load_excited_functions<Type,ClosedShell>(ctx.world,directory,description,binding.controls().protocol);
    },true);
    session.timing().end();session.timing().begin("transition_properties");
    const auto analysis=analyze_es_bundle<Type,ClosedShell>(ctx.world,ctx.gs.orbitals_alpha(),state);
    std::vector<std::vector<double>> y_norms;
    if constexpr(std::is_same_v<Type,Full>) {
      for(std::size_t i=0;i<state.roots.size();++i) {
        if(!analysis[i].dipole_available) {y_norms.emplace_back();continue;}
        const auto& root=state.roots[i];
        std::vector<double> norms;norms.reserve(root.y_alpha.size());
        for(const auto& function:root.y_alpha) norms.push_back(function.norm2());
        y_norms.push_back(std::move(norms));
      }
    }
    if(ctx.print_level>=PrintLevel::Normal)
      print_es_analysis(ctx.world,analysis,std::is_same_v<Type,TDA>,y_norms,
                        madness::iteration::prints(ctx.verbosity,madness::iteration::needed::assembled_values));
    checkpoint::detail::on_root(ctx.world,[&] {
      std::ofstream output(ctx.calc_dir+"/es_analysis__"+protocol_key()+".json");
      output.exceptions(std::ios::failbit|std::ios::badbit);
      output<<es_analysis_to_json(analysis,std::is_same_v<Type,TDA>,y_norms,state.stable_index).dump(2)<<'\n';
      output.close();return std::string{};
    });
    const std::string not_computed="not computed by analysis diagnostic";
    for(std::size_t i=0;i<state.roots.size();++i) {
      json values{{"excitation_energy",analysis[i].omega}}, missing=json::object();
      if(analysis[i].dipole_available) {
        values["oscillator_strength"]=analysis[i].oscillator;
        json dipole=json::object();
        for(std::size_t a=0;a<3;++a) dipole[std::string(1,"xyz"[a])]=analysis[i].dipole[a];
        values["transition_dipole"]=dipole;
      } else {
        missing["oscillator_strength"]=not_computed;
        missing["transition_dipole"]=not_computed;
      }
      if(analysis[i].quadrupole_available) {
        json quadrupole=json::object();
        for(std::size_t a=0;a<3;++a) for(std::size_t b=0;b<3;++b)
          quadrupole[std::string(1,"xyz"[a])+"xyz"[b]]=analysis[i].quadrupole[a][b];
        values["transition_quadrupole"]=quadrupole;
      } else {
        missing["transition_quadrupole"]=not_computed;
      }
      excited["values"].push_back(item(make_root_id(state.stable_index[i]),values,missing));
    }
    session.timing().end();session.timing().end();
  } else {
    const auto reason=set_name+" has not met the requested criteria";
    json missing=json::object();
    for(const char* field:{"excitation_energy","oscillator_strength","transition_dipole","transition_quadrupole"})
      missing[field]=reason;
    for(const auto& id:binding.expected().states) excited["values"].push_back(item(id,json::object(),missing));
  }
  session.set_properties({{"excited_states",excited}});
}

// ---------------------------------------------------------------------------
// VBC quadratic-source build (beta-ii-b, closed-shell). Not iterative: load the
// two converged first-order states at this protocol, build the source, save.
// ---------------------------------------------------------------------------

namespace detail_exec {
// Exact qualified checkpoints only. Property readers never fall back to guesses.
inline std::optional<ResponseStateXY<ClosedShell>> qualified_fd_as_xy(
    ExecutorContext& ctx,const Perturbation& perturbation,double frequency) {
  auto& session=*ctx.shared_response;
  const auto found=session.states.find(fd_node_id(perturbation,frequency));
  if(found==session.states.end() || !found->second->complete()) return std::nullopt;
  auto& binding=*found->second;
  madness::checkpoint::records::require(binding.scientific_settings().at("frequency")==frequency &&
      binding.scientific_settings().at("perturbation")==perturbation.description(),
      "distinct requested frequencies share a response task identity");
  session.timing().begin("load",{{"task_id",binding.info().task_id}});
  ResponseStateXY<ClosedShell> functions;
  if(is_static_freq(frequency)) {
    auto x=binding.load<ResponseStateX<ClosedShell>>(true);
    functions.x_alpha=std::move(x.x_alpha);functions.y_alpha=functions.x_alpha;
  } else functions=binding.load<ResponseStateXY<ClosedShell>>(true);
  session.timing().end();return functions;
}
inline ResponseStateXY<ClosedShell> load_shared_vbc(ExecutorContext& ctx,
    response_shared::SecondOrderNode& binding) {
  auto& timing=ctx.shared_response->timing();
  timing.begin("load",{{"task_id",binding.id()}});
  auto source=binding.load([&](const auto& directory,const auto& description) {
    const madness::iteration::Protocol current{madness::FunctionDefaults<3>::get_thresh(),madness::FunctionDefaults<3>::get_k()};
    madness::checkpoint::records::require(description.protocol==current &&
        description.details.at("archives")==nlohmann::json::array({"response"}),
        "VBC node representation or payload differs");
    return ResponseStateXY<ClosedShell>::load(ctx.world,
        madness::native_archive_name((directory/"response").string()),IoBackend::Native);
  });
  timing.end();return source;
}
} // namespace detail_exec

inline NodeResult
solve_vbc_closed_shell(ExecutorContext &ctx, const CalcNode &node, double thresh) {
  using namespace madness;
  World &world = ctx.world;
  GroundState &gs = ctx.gs;

  auto& session=*ctx.shared_response;
  session.timing().begin("vbc",{{"task_id",node.id}});
  auto& B_state=session.at(fd_node_id(node.pert,node.freq));
  auto& C_state=session.at(fd_node_id(node.pert_c,node.freq_c));
  if(!B_state.complete() || !C_state.complete()) {
    session.timing().end();return {false,"",{}};
  }
  const auto definition=response_shared::vbc_source_settings({thresh,default_k_for_thresh(thresh)});
  auto* shared=&session.bind_second_order(node.id,vbc_node_name(node.pert,node.pert_c,node.freq,node.freq_c),definition,
      nlohmann::json::array({response_shared::lineage_entry("B",B_state),response_shared::lineage_entry("C",C_state)}));
  if(shared->reusable()) {
    session.record_vbc(node.id,*shared);session.timing().end();
    return {true,protocol_key(thresh,default_k_for_thresh(thresh)),{}};
  }
  session.timing().begin("prepare");

  const bool gs_verbose = ctx.log_prefix.empty();  // F2d-ii: quiet in subworlds
  set_response_protocol(world, ctx.L, thresh, /*override_k=*/-1, gs_verbose);
  const double t0 = FunctionDefaults<3>::get_thresh();
  {
    auto coulop = poperatorT(
        CoulombOperatorPtr(world, gs.params().lo(), 0.001 * t0));
    gs.prepare(world, 0.001 * t0, coulop, ctx.fock_json, gs_verbose);
  }

  // The two first-order states (X,Y) at this protocol, from their nodes.
  auto B = detail_exec::qualified_fd_as_xy(ctx,node.pert,node.freq);
  auto C = detail_exec::qualified_fd_as_xy(ctx,node.pert_c,node.freq_c);
  if (!B || !C) {
    if (world.rank() == 0)
      madness::print("[CALC] solve_vbc:", node.id,
                     "— a converged FD prerequisite is missing; skipping.");
    return NodeResult{/*converged=*/false, /*reached_protocol_key=*/"", {}};
  }

  session.timing().end();session.timing().begin("build_source");
  const double c_xc = gs.hf_exchange_coefficient();
  const double lo   = gs.params().lo();
  auto g0 = build_response_ground_state_closed_shell(world, gs, c_xc, lo);

  // Raw one-electron perturbation operators for B and C (dipole only for now).
  // Generic per-kind operator (dipole -> r, nuclear -> dV/dQ) so a VBC pair can
  // mix types (Raman = dipole x nuclear).
  auto VB_op = perturbation_operator(world, gs, node.pert);
  auto VC_op = perturbation_operator(world, gs, node.pert_c);

  // The (P,Q) builder sums the two photon orderings in one call.
  auto vbc_src = tpa::quadratic_source(world, g0, *B, *C, VB_op, VC_op);
  if (world.rank() == 0)
    madness::print(ctx.log_prefix, "[CALC] solve_vbc: source =",
                   "tpa::quadratic_source (P,Q)");
  session.timing().end();session.timing().begin("checkpoint");
  shared->save({"response"},[&](const std::filesystem::path& directory) {
    vbc_src.save(world,native_archive_name((directory/"response").string()),IoBackend::Native);
  });
  session.record_vbc(node.id,*shared);
  session.timing().end();session.timing().end();session.publish_timings();

  if (world.rank() == 0) {
    if (ctx.log_prefix.empty())
      madness::print("[CALC] solve_vbc: built", node.id, "at", protocol_key());
    else
      madness::print(ctx.log_prefix, "[CALC] solve_vbc: built", node.id, "at",
                     protocol_key());
  }
  return NodeResult{/*converged=*/true, protocol_key(), {}};
}

// ---------------------------------------------------------------------------
// FdResponseExecutor — dispatches the runtime (Type × Shell) for FD nodes.
// ---------------------------------------------------------------------------

class FdResponseExecutor : public ICalcExecutor {
public:
  explicit FdResponseExecutor(ExecutorContext ctx) : ctx_(ctx) {}

  /// Work this executor declines for the run's ground state, and why: a VBC
  /// node (a quadratic source) on an open-shell ground state. The same answer
  /// on every rank, so the scheduler's audit can name it.
  [[nodiscard]] std::optional<std::string> out_of_scope(const CalcNode &node) const override {
    if (node.kind == CalcKind::VBC && !ctx_.gs.is_spin_restricted())
      return std::string("the quadratic response source (VBC) and its contraction are implemented "
                         "for a closed-shell ground state only, and this ground state is open-shell");
    return std::nullopt;
  }

  NodeResult run_protocol(const WorkItem &item) override {
    const CalcNode &node = *item.node;
    // Each solve receives one threshold, but its role belongs to the full
    // requested ladder. This also applies in subworld executors.
    ctx_.final_protocol = std::none_of(
        node.protocols.begin(), node.protocols.end(),
        [&](double threshold) { return threshold < item.thresh; });

    if (node.kind == CalcKind::ES)
      return node.tda ? solve_shared_excited<TDA>(ctx_,node.n_roots,item.thresh) :
                        solve_shared_excited<Full>(ctx_,node.n_roots,item.thresh);
    if (const auto why = out_of_scope(node)) {
      if (ctx_.world.rank() == 0)
        madness::print("[CALC] run_protocol:", node.id, "is out of scope:", *why, "— skipping.");
      return NodeResult{/*converged=*/false, /*reached_protocol_key=*/"", {}};
    }
    if (node.kind == CalcKind::VBC) {
      if (ctx_.world.rank() == 0 && ctx_.print_level >= PrintLevel::Verbose)
        madness::print("[CALC] run_protocol: VBC", node.id,
                       " B=", node.pert.description(), "@", node.freq,
                       " C=", node.pert_c.description(), "@", node.freq_c,
                       " thresh=", item.thresh);
      return solve_vbc_closed_shell(ctx_, node, item.thresh);
    }
    if (node.kind == CalcKind::DerivedFD) {
      // DerivedFD should have been promoted to FD on expansion.
      if (ctx_.world.rank() == 0)
        madness::print("[CALC] run_protocol: node", node.id,
                       "kind DerivedFD should be promoted to FD — skipping");
      return NodeResult{/*converged=*/false, /*reached_protocol_key=*/"", {}};
    }

    const bool restricted = ctx_.gs.is_spin_restricted();
    const bool is_static   = detail_exec::is_static_freq(node.freq);

    if (ctx_.world.rank() == 0 && ctx_.print_level >= PrintLevel::Verbose) {
      madness::print("[CALC] run_protocol: node", node.id,
                     " pert=", node.pert.description(),
                     " freq=", node.freq, " thresh=", item.thresh,
                     " type=", (is_static ? "static" : "full"),
                     " shell=", (restricted ? "closed" : "open"),
                     " action=", node_action_name(item.action));
    }

    if (restricted && is_static)
      return solve_fd_protocol<Static, ClosedShell>(ctx_, node.pert, node.freq, item.thresh);
    if (restricted && !is_static)
      return solve_fd_protocol<Full, ClosedShell>(ctx_, node.pert, node.freq, item.thresh);
    if (!restricted && is_static)
      return solve_fd_protocol<Static, OpenShell>(ctx_, node.pert, node.freq, item.thresh);
    return solve_fd_protocol<Full, OpenShell>(ctx_, node.pert, node.freq, item.thresh);
  }

private:
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
  // Rank 0 only: the snapshot of what is saved (node_snapshot.hpp), read at
  // the start of every pass. Required.
  std::function<nlohmann::json()> saved_states;
  std::function<StatePlan(const CalcNode&)> plan_derived;
  std::function<void()> collect_subworlds;
  std::function<bool(const CalcNode&)> dependencies_ready;
  std::function<bool(const std::string&)> qualified_bundle;
  SeedStrategy seed = SeedStrategy::NearestConverged;
  int          fd_subworlds = 0;   // F2f: subworlds PER NODE (0=off, 1=node-aligned)
  int          fd_subworld_ranks = 0;   // ranks per subworld; >0 = universe split (may span nodes)
  // The session's plan for every state (D-016); schedule() requires one.
  StatePlans   state_plans;
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
// So the schedule must come out byte-identical on every rank: same snapshot
// in (rank 0 builds it from the node folders after a fence and broadcasts
// it), deterministic logic out (sorted ramp + std::map perturbation order + pure
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
  /// pass: each pass reads the snapshot of the node folders, expands any
  /// converged ES bundle into concrete DerivedFD nodes, re-schedules (so
  /// actions reflect current disk) and runs the first non-empty wave.
  /// Terminates when the schedule is empty: every rung is scheduled at most
  /// once per launch (`solved_`), and the node's own solver decides within
  /// that one run whether it converges, so no wave can repeat and no state
  /// can be scheduled again and again. On exit the plan is audited against
  /// the snapshot: planned VBC nodes that never converged at their top rung
  /// (out of scope for this ground state, or gated out by unconverged
  /// prerequisites) are reported via stop_reason 'complete_with_dropped_beta'
  /// + diag.dropped_beta, and recorded under run_summary/dropped_work in the
  /// metadata — a run with silently dropped beta/raman work is never labelled
  /// plain 'complete' (review fix).
  /// Returns the R1c scheduler trace (doc 16 L3) -> Output.diagnostics: the
  /// wave-by-wave reconcile actions (id/thresh/action per item), stop_reason,
  /// and pass count. Built identically on every rank (schedule() is
  /// deterministic), so the value is rank-consistent.
  /// F2 (doc 32 §5): inner per-subworld solve — build the ground state in `sub`
  /// from the archive, solve exactly `items` (its owned FD subset) at `thresh`,
  /// save to disk under the `node_index` metadata shard. Supplied by the
  /// orchestrator (which holds the archive path). Empty ⇒ the single-World path.
  using SubworldSolve =
      std::function<void(madness::World &sub, const std::vector<WorkItem> &items,
                         double thresh, int gid,
                         const std::string &log_prefix)>;

  nlohmann::json run(madness::World &world, ICalcExecutor &exec,
                     SubworldSolve fan_out = {}) {
    const std::vector<double> ramp = global_ramp();
    // Output only: the end-of-run audit is recorded in it. Nothing is read
    // from it to decide.
    const std::string meta_path = calc_dir_ + "/response_metadata.json";
    if (!policy_.saved_states)
      throw std::logic_error("the response scheduler needs the snapshot of the saved states");

    nlohmann::json diag;                       // R1c scheduler trace
    diag["schedule"] = nlohmann::json::array();
    int pass = 0;


    for (;;) {
      world.gop.fence();
      // What is saved, from the node folders, judged against this request.
      const nlohmann::json meta_json = detail_exec::snapshot_rank0(world, policy_.saved_states);
      // Grow the DAG from any ES bundle converged ON DISK (idempotent and
      // restart-safe: reads roots from the snapshot, not in-memory run state).
      expand_converged_es(world, meta_json);
      auto waves = schedule(dag_, ramp, meta_json, policy_.state_plans, solved_);
      // Native shared completion is judged against this invocation's request.
      // Legacy analysis flags can describe an earlier, looser request.
      if(policy_.dependencies_ready) for(auto& wave:waves)
        wave.erase(std::remove_if(wave.begin(),wave.end(),[&](const WorkItem& item) {
          return !policy_.dependencies_ready(*item.node);
        }),wave.end());
      waves.erase(std::remove_if(waves.begin(), waves.end(),
                                 [](const auto &w) { return w.empty(); }),
                  waves.end());
      if (waves.empty()) {
        // Review fix (confirmed MEDIUM — silent beta drop): an empty schedule
        // does NOT mean every planned node executed. Honest-climb can walk a
        // stubborn FD prerequisite up the whole ladder unconverged; the
        // dependent VBC node is then gated out of every wave (gated nodes are
        // invisible to schedule()), and the beta/raman it feeds is dropped
        // without a trace. Audit the plan against the snapshot: any VBC node
        // without a converged entry at its top rung never delivered.
        // (audit_dropped_work, calc_manager.hpp: undelivered VBC sources and
        // two-photon legs whose ES bundle never converged.)
        // Work this executor declines for this ground state, by node id,
        // with the reason (out_of_scope; decided on every rank alike).
        std::map<std::string, std::string> refused;
        for (const auto &n : dag_)
          if (const auto why = exec.out_of_scope(n)) refused[n.id] = *why;
        const nlohmann::json dropped = audit_dropped_work(dag_, meta_json, refused);
        // Dropped work ⇒ the beta/raman/2PA rows it feeds cannot be assembled.
        // Record the audit in the metadata (rank 0, through the layer) even
        // when empty, so a later completing run CLEARS a stale drop list.
        if (world.rank() == 0) {
          auto meta_out = ResponseMetadata::load_or_create(meta_path);
          meta_out.set_dropped_work(dropped);
          meta_out.save();
        }

        if (dropped.empty()) {
          if (world.rank() == 0)
            madness::print("[CALC] run: nothing left to schedule — done");
          diag["stop_reason"] = "complete";
        } else {
          if (world.rank() == 0) {
            madness::print("[CALC] run: WARNING —", (int)dropped.size(),
                           "planned VBC/beta node(s) were never solved "
                           "(beta/raman rows will be missing):");
            for (const auto &d : dropped)
              madness::print("    dropped:", d.value("id", std::string("?")),
                             "—", d.value("reason", std::string("?")));
          }
          diag["stop_reason"] = "complete_with_dropped_beta";
          diag["dropped_beta"] = dropped;
        }
        break;
      }

      const std::string sig = wave_signature(waves.front());

      // R1c: record this wave (id/thresh/action per item) + protocol markers.
      // A wave is one protocol level, so wthresh = the front item's threshold.
      const auto &wave = waves.front();
      const double wthresh = wave.front().thresh;
      const int    wk      = default_k_for_thresh(wthresh);
      int pidx = 0;
      for (size_t i = 0; i < ramp.size(); ++i)
        if (std::abs(ramp[i] - wthresh) <= 1e-12 * std::max(ramp[i], wthresh)) {
          pidx = static_cast<int>(i);
          break;
        }
      nlohmann::json wrec = {{"pass", pass}, {"protocol_index", pidx},
                             {"thresh", wthresh}, {"k", wk},
                             {"items", nlohmann::json::array()}};
      for (const auto &it : wave)
        wrec["items"].push_back({{"id", it.node->id},
                                 {"thresh", it.thresh},
                                 {"action", node_action_name(it.action)}});
      diag["schedule"].push_back(std::move(wrec));
      // A planned rung is scheduled once; the session's node owns completion.
      for (const auto &it : wave)
        solved_.insert(it.node->id + "|" + protocol_key_at(it.thresh));

      if (world.rank() == 0) {
        madness::print("[CALC] run: wave of", (int)wave.size(),
                       "protocol step(s): {", sig, "}");
        madness::print("PROTOCOL_START  pass=", pass, "  protocol_index=", pidx,
                       "  thresh=", wthresh, "  k=", wk,
                       "  items=", (int)wave.size());
      }
      const bool fan_enabled = fan_out && (policy_.fd_subworlds > 0 || policy_.fd_subworld_ranks > 0);
      if (!fan_enabled) {
        // Single-World reference path — byte-for-byte the pre-F2 behaviour
        // (every rank solves each wave item together, in wave order).
        for (const auto &item : wave)
          exec.run_protocol(item);
      } else {
        // F2 (doc 32 §5): split the wave into the fan-out-eligible subset (FD /
        // NuclearFD / VBC — all independent per item) and the remainder (ES, which
        // stays single-World). The two sub-phases run SEQUENTIALLY on all universe
        // ranks — never concurrently across communicators (would deadlock). VBC
        // items (F2g) read their converged FD prerequisites from the shared
        // calc_dir; the scheduler's dependency gate guarantees those archives are
        // already on disk from a prior pass, so there is no intra-wave race.
        std::vector<WorkItem> fan_items, rest;
        for (const auto &it : wave) {
          const CalcKind k = it.node->kind;
          ((k == CalcKind::FD || k == CalcKind::NuclearFD || k == CalcKind::VBC)
               ? fan_items
               : rest)
              .push_back(it);
        }

        if (!fan_items.empty()) {
          // F2f: groups_per_node = policy_.fd_subworlds (1 = node-aligned; larger
          // = sub-node / NUMA packing for small systems). G = total subworlds.
          // Guard the pool construction with the same collective-error discipline
          // as the fan_out() call below (review: it was unprotected). This
          // converts a SYMMETRIC/post-collective failure — World ctor, bad_alloc,
          // a Split all ranks fail — into a clean collective abort instead of a
          // hang. (An ASYMMETRIC throw inside make_subworld_pool's own internal
          // collectives — ranks_per_host gather/broadcast, Split, gop.fence — is
          // inherently unrecoverable in MPI: the non-throwing ranks are already
          // blocked inside those collectives and never reach the max() below.)
          // Right-size the pool to THIS wave (2026-09-10). With a fixed
          // groups_per_node = deck `subworlds` (8 = one rank per NUMA domain), a
          // wave of 3 static legs left 5 of 8 subworlds without an item; they
          // waited in the universe collective below for the whole wave (~16 min
          // per leg at thresh 1e-6) and ThreadPool::await killed them after 900 s
          // (closeout job 2162152). The number of items is known here, so use
          //   groups_per_node = clamp(ceil(n_items / n_nodes), 1, min(subworlds, ranks_per_node)):
          // 3 items on one 8-rank node -> 3 subworlds of 3/3/2 ranks (no idle
          // rank, each leg gets 2-3 ranks); 24 items on 3 nodes -> 24 x 1 rank;
          // 12 items on 3 nodes -> 12 x 2 ranks. Subworlds stay node-aligned, so
          // with fewer items than nodes some nodes still idle for the wave.
          //
          // Large-system regime (deck `subworld_ranks` R > 0): the per-node
          // packing is replaced by a UNIVERSE-level contiguous split into
          //   G = clamp(min(n_items, universe_size / R), 1, universe_size)
          // subworlds of ~universe_size/G ranks each, which may span nodes. R is
          // the SCF sizing (10-20 MOs per rank): valinomycin (300 MOs) runs one
          // SCF on 6-8 nodes x 8 ranks, so R = 48-64 gives each independent
          // state the same footprint as the ground-state solve. With fewer
          // items than universe/R, the blocks grow (all ranks stay busy).
          int gpn_wave = std::max(1, policy_.fd_subworlds);
          int G_universe = 0;   // > 0 selects the universe-level split
          {
            const auto rph = madness::ranks_per_host(world);   // collective
            const int n_nodes = std::max<int>(1, static_cast<int>(rph.size()));
            int rpn = world.size();
            for (const auto &kv : rph) rpn = std::min<int>(rpn, static_cast<int>(kv.second.size()));
            const int n_items = static_cast<int>(fan_items.size());
            if (policy_.fd_subworld_ranks > 0) {
              const int R = std::min(policy_.fd_subworld_ranks, world.size());
              G_universe = std::max(1, std::min(n_items, world.size() / R));
            } else {
              const int want = (n_items + n_nodes - 1) / n_nodes;   // ceil
              gpn_wave = std::max(1, std::min({policy_.fd_subworlds, rpn, want}));
            }
          }
          NodeSubworldInfo info;
          std::shared_ptr<madness::World> sub;
          std::string pool_err;
          try {
            sub = (G_universe > 0)
                      ? make_subworld_pool_universe(world, G_universe, &info)
                      : make_subworld_pool(world, gpn_wave, &info);
          } catch (const std::exception &e) { pool_err = e.what(); }
            catch (...) { pool_err = "unknown exception in make_subworld_pool"; }
          int pool_bad = pool_err.empty() ? 0 : 1;
          world.gop.max(pool_bad);
          if (pool_bad) {
            if (!pool_err.empty())
              madness::print("[SUBWORLD-POOL-ERROR] universe rank", world.rank(),
                             ":", pool_err);
            MADNESS_EXCEPTION("subworld pool construction failed (see "
                              "[SUBWORLD-POOL-ERROR]); aborting collectively", 0);
          }
          const int G = info.n_subworlds;
          if (G <= 1) {
            // One subworld total (single node, P=1) ⇒ no real partition (subworld
            // == universe). Literal G=0 path, no redundant GS reload (doc 32 §7).
            if (world.rank() == 0)
              madness::print("SUBWORLD_FANOUT skipped: n_subworlds=", G,
                             " (nodes=", info.n_nodes, " groups_per_node=",
                             info.groups_per_node, ") — universe path");
            sub.reset();
            for (const auto &it : fan_items) exec.run_protocol(it);
          } else {
            std::vector<WorkItem> mine;   // round-robin by global gid (deterministic)
            for (size_t i = 0; i < fan_items.size(); ++i)
              if (static_cast<int>(i) % G == info.gid)
                mine.push_back(fan_items[i]);
            // F2d: per-subworld log tag "[g{gid}/{G} {host}] " (doc 32 §5.6).
            const std::string lp = "[g" + std::to_string(info.gid) + "/" +
                                   std::to_string(G) + " " + info.hostname + "] ";
            if (world.rank() == 0)
              madness::print("SUBWORLD_FANOUT  pass=", pass, "  n_subworlds=", G,
                             "  groups_per_node=", info.groups_per_node,
                             "  ranks_per_subworld=", info.subworld_size,
                             "  nodes_per_subworld=", info.nodes_spanned,
                             "  fan_items=", (int)fan_items.size(),
                             (G_universe > 0 ? "  (universe split, deck subworld_ranks="
                                             : "  (right-sized from deck subworlds="),
                             (G_universe > 0 ? policy_.fd_subworld_ranks : policy_.fd_subworlds),
                             ")");
            // S1 pmap discipline: point the default pmap at the subworld so
            // everything built inside is subworld-local; restore BEFORE reset.
            // Exception safety: a subworld-local throw must NOT skip the pmap
            // restore or the universe fence below — the failing ranks would
            // reach the driver's top-level catch and finalize while every
            // other rank blocks in the fence forever. Catch, restore, fence,
            // then agree collectively and abort together. (A throw that leaves
            // SIBLING ranks of the same subworld inside a subworld collective
            // is inherently unrecoverable in MPI — this handles the
            // synchronized-failure cases: archive open, shard save, guards.)
            std::string fan_err;
            // Save the universe's CURRENT pmap object and restore THAT object
            // afterwards. set_default_pmap(world) builds a NEW map, so every
            // function created after the fan-out (seed bundles reloaded from
            // disk, ES intermediates, V0x/gamma products) was distributed
            // differently from the ground-state functions built before it
            // (orbitals, V_local, K0, Q); operations mixing the two silently
            // lost contributions: LiH sigma orbital <phi|V_local|phi> came out
            // -0.008 instead of -0.194 and Q left 1e-2 occupied overlaps, and
            // every seeded ES solve after a fan-out diverged (attempts 3-11,
            // 2026-09-10). With subworlds 0 the same solve converges.
            const auto universe_pmap = madness::FunctionDefaults<3>::get_pmap();
            madness::FunctionDefaults<3>::set_default_pmap(*sub);
            try {
              fan_out(*sub, mine, wthresh, info.gid, lp);
              sub->gop.fence();
            } catch (const std::exception &e) {
              fan_err = e.what();
            } catch (...) {
              fan_err = "unknown exception in subworld fan-out";
            }
            madness::FunctionDefaults<3>::set_pmap(universe_pmap);   // the ORIGINAL object, not a new map
            sub.reset();
            if (!fan_err.empty())
              madness::print("[FANOUT-ERROR] universe rank", world.rank(),
                             "(gid", info.gid, "):", fan_err);
            int fan_bad = fan_err.empty() ? 0 : 1;
            world.gop.max(fan_bad);
            if (fan_bad) {
              world.gop.fence();
              MADNESS_EXCEPTION(
                  "subworld fan-out failed on at least one subworld "
                  "(see [FANOUT-ERROR]); aborting collectively", 0);
            }
          }
          world.gop.fence();
          // F2b/F2g: collapse the per-group FD+VBC metadata shards into the
          // canonical file (universe rank 0; through the metadata layer). Only
          // when we actually sharded (G > 1). The merge itself can throw
          // (corrupt shard, ENOSPC) — that must be collective too, not a
          // rank-0-only unwind that strands everyone else at the fence.
          int merge_bad = 0;
          if (world.rank() == 0 && G > 1) {
            try {
              ResponseMetadata::merge_state_shards(calc_dir_, G);
            } catch (const std::exception &e) {
              merge_bad = 1;
              madness::print("[SHARD-MERGE-ERROR]", e.what());
            }
          }
          world.gop.max(merge_bad);
          if (merge_bad)
            MADNESS_EXCEPTION("per-wave shard merge failed on rank 0 "
                              "(see [SHARD-MERGE-ERROR])", 0);
          world.gop.fence();
          if(policy_.collect_subworlds) policy_.collect_subworlds();
        } else if ((policy_.fd_subworlds > 0 || policy_.fd_subworld_ranks > 0) && !rest.empty() &&
                   world.rank() == 0) {
          // Review MED: this wave is ES-only (ES is deliberately excluded from
          // fan-out — it stays single-World), so the requested `subworlds` had
          // no effect on it. Say so, or the knob looks silently ignored.
          madness::print("SUBWORLD_FANOUT skipped: wave is ES-only (",
                         (int)rest.size(),
                         "item(s)); ES runs single-World, `subworlds` applies "
                         "only to FD/VBC waves");
        }
        for (const auto &it : rest) exec.run_protocol(it);   // ES: universe
      }
      world.gop.fence();
      if (world.rank() == 0)
        madness::print("PROTOCOL_DONE  pass=", pass, "  protocol_index=", pidx,
                       "  thresh=", wthresh, "  items=", (int)wave.size());
      ++pass;
    }
    diag["passes"] = pass;
    return diag;
  }

private:
  /// Coarse->fine union of every node's protocol ladder.
  std::vector<double> global_ramp() const {
    std::vector<double> all;
    for (const auto &n : dag_)
      for (double t : n.protocols) all.push_back(t);
    return detail_calc::sorted_ramp(std::move(all));
  }

  /// Stable, rank-deterministic name of a wave for the log (sorted
  /// id@protocol list; protocol_key, not freq_key, so two rungs of one node
  /// read differently).
  static std::string wave_signature(const std::vector<WorkItem> &wave) {
    std::vector<std::string> parts;
    parts.reserve(wave.size());
    for (const auto &it : wave)
      parts.push_back(it.node->id + "@" +
                      protocol_key(it.thresh, default_k_for_thresh(it.thresh)));
    std::sort(parts.begin(), parts.end());
    std::string s;
    for (const auto &p : parts) { s += p; s += ';'; }
    return s;
  }

  /// Turn each symbolic DerivedFD "*" whose ES bundle has CONVERGED ON DISK into
  /// one promoted FD node per root (freq = root energy at the ES bundle's top
  /// protocol). Reads roots from the snapshot of the excited set's node
  /// folder, so it is restart-safe (works on a resumed run where the ES
  /// converged in a previous process) and idempotent (new nodes are deduped
  /// against existing dag_ ids; re-running is a no-op). Deterministic across
  /// ranks: same snapshot + same dag_ in, same nodes appended in the same order.
  void expand_converged_es(madness::World &world, const nlohmann::json &meta) {
    if (!meta.contains("excited_states")) return;
    std::map<std::string,std::pair<double,std::string>> have;
    for(const auto& n:dag_) if(n.kind==CalcKind::FD)
      have.emplace(n.id,std::make_pair(n.freq,n.pert.description()));

    std::vector<CalcNode> additions;
    for (const auto &sym : dag_) {
      if (!sym.is_symbolic()) continue;
      if (sym.prerequisites.empty()) continue;          // no ES dependency
      const std::string &es_id = sym.prerequisites.front();
      // Resolve the ES node to its top protocol (derived FDs start once the ES
      // bundle has reached its finest protocol).
      const CalcNode *es = nullptr;
      for (const auto &n : dag_) if (n.id == es_id) { es = &n; break; }
      if (!es || es->protocols.empty()) continue;
      if(policy_.qualified_bundle && !policy_.qualified_bundle(es_id)) continue;
      const std::string es_key = protocol_key_at(es->protocols.back());
      if (!meta["excited_states"].contains(es_key)) continue;
      const auto &b = meta["excited_states"][es_key];
      if (!detail_calc::es_entry_matches(b, *es)) continue;  // another bundle's roots (C5)
      if (!b.value("converged", false)) continue;       // ES not converged yet
      if (!b.contains("roots") || !b["roots"].is_array()) continue;

      for (size_t i = 0; i < b["roots"].size(); ++i) {
        const double w = b["roots"][i].value(
            "omega", std::numeric_limits<double>::quiet_NaN());
        if (!(w == w)) continue;                         // skip NaN omega
        const double fd_freq = sym.es_freq_factor * w;   // e.g. ωₙ/2 (off-pole)
        CalcNode c;
        // A derived FD at ω = factor·(root energy) is an ordinary FD point ->
        // kind FD, so the FD executor solves it and skip/restart/seed logic is
        // uniform. factor=0.5 puts it at the two-photon resonance ωₙ/2, off the
        // linear-response pole at ω = ωₙ.
        c.kind       = CalcKind::FD;
        c.pert       = sym.pert;
        c.freq       = fd_freq;
        c.protocols  = sym.protocols;
        c.es_root_id = b["roots"][i].value("root_id",make_es_root_label(static_cast<int>(i)));
        c.es_freq_factor=sym.es_freq_factor;
        c.prerequisites={es_id};
        c.id         = fd_node_id(c.pert, fd_freq);  // id keys on the node's TRUE freq (ωₙ/2), not ωₙ
        const auto physical=std::make_pair(c.freq,c.pert.description());
        if(const auto existing=have.find(c.id);existing!=have.end()) {
          madness::checkpoint::records::require(existing->second==physical,
              "distinct requested frequencies share a response task identity");
          continue;
        }
        if(policy_.plan_derived) policy_.state_plans[c.id]=policy_.plan_derived(c);
        have.emplace(c.id,physical);
        additions.push_back(std::move(c));
      }
    }
    if (!additions.empty()) {
      if (world.rank() == 0)
        madness::print("[CALC] expanded", (int)additions.size(),
                       "derived FD node(s) at ω = es_freq_factor·ωₙ from converged ES roots");
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
  SolvedRungs           solved_;   // rungs this run solved (planned states)
};

// ---------------------------------------------------------------------------
// beta property assembly (Tier-A). A post-run step: contract the converged A
// states against each VBC source (+ zeta terms) into the beta tensor, recorded
// under properties/beta. Runs at one protocol (typically the top of the ramp).
// ---------------------------------------------------------------------------

inline char beta_axis_name(int a) { return a == 0 ? 'x' : a == 1 ? 'y' : 'z'; }

// Why a β component is missing: the folders of the states it needs that have
// not met the final request (a VBC source not built counts).
inline std::string missing_beta_inputs(const response_shared::Session& session,const std::string& source_id,
                                       const std::vector<std::string>& driven_ids) {
  std::vector<std::string> names;
  for(const auto& id:driven_ids)
    if(session.knows(id) && !session.complete(id)) names.push_back(session.node_name(id));
  const auto source=session.second_order.find(source_id);
  if(source==session.second_order.end() || !source->second->reusable())
    names.push_back(session.knows(source_id) ? session.node_name(source_id) : source_id);
  std::string text;
  for(const auto& name:names) text+=(text.empty() ? "" : ", ")+name;
  return (text.empty() ? std::string("a required state") : text)+" has not met the requested criteria";
}

// Contract each VBC source node with the qualified driven states of this
// request.
inline void assemble_shared_beta(ExecutorContext& ctx,const ResponsePlan& plan,double thresh) {
  using namespace madness;
  using json=nlohmann::json;
  auto& session=*ctx.shared_response;
  auto& world=ctx.world;auto& gs=ctx.gs;
  session.timing().begin("quadratic_properties");
  const bool supported=gs.is_spin_restricted();
  std::optional<ResponseGroundState> ground;
  if(supported) {
    session.timing().begin("prepare");
    set_response_protocol(world,ctx.L,thresh);
    auto coulop=poperatorT(CoulombOperatorPtr(world,gs.params().lo(),0.001*thresh));
    gs.prepare(world,0.001*thresh,coulop,ctx.fock_json);
    ground=build_response_ground_state_closed_shell(world,gs,gs.hf_exchange_coefficient(),gs.params().lo());
    session.timing().end();
  }
  std::map<std::string,json> beta_tensors,raman_tensors;
  const auto missing_inputs=[&](const std::string& source_id,const auto& request,const std::string& A_id) {
    return missing_beta_inputs(session,source_id,{A_id,fd_node_id(request.pert_b,request.freq_b),
                                                  fd_node_id(request.pert_c,request.freq_c)});
  };
  json beta_rows=json::array(),raman_rows=json::array();
  for(const auto& request:expand_vbc_requests(plan.vbc,static_cast<int>(gs.molecule().natom()))) {
    const auto id=vbc_node_id(request.pert_b,request.pert_c,request.freq_b,request.freq_c);
    const bool nuclear=request.pert_c.kind==Perturbation::Kind::NuclearDisplacement;
    // One item per frequency pair holds the whole tensor. A polarizability
    // gradient's pair is (laser frequency, 0); its components are keyed by
    // nuclear coordinate, then by dipole pair.
    const auto tensor_key=json::array({request.freq_b,request.freq_c}).dump();
    auto& tensor=(nuclear ? raman_tensors : beta_tensors)[tensor_key];
    if(tensor.is_null()) {
      tensor={{"components",json::object()},{"availability",json::object()}};
      if(!nuclear) {tensor["frequency_b"]=request.freq_b;tensor["frequency_c"]=request.freq_c;}
      else {
        tensor["frequency"]=request.freq_b;
        if(request.freq_c!=0.0) tensor["frequency_c"]=request.freq_c;
      }
    }
    // A component of a polarizability gradient is written under its coordinate.
    const auto slot=[&](const char* key) -> json& {
      return nuclear ? tensor[key][request.pert_c.description()] : tensor[key];
    };
    session.name_state(tensor,fd_node_id(request.pert_b,request.freq_b));
    session.name_state(tensor,fd_node_id(request.pert_c,request.freq_c));
    session.name_state(tensor,id);
    for(int axis=0;axis<3;++axis)
      session.name_state(tensor,fd_node_id(Perturbation::dipole(axis),request.freq_b+request.freq_c));
    std::optional<ResponseStateXY<ClosedShell>> source,B,C;
    const auto selected=session.second_order.find(id);
    if(supported && selected!=session.second_order.end() && selected->second->reusable()) {
      B=detail_exec::qualified_fd_as_xy(ctx,request.pert_b,request.freq_b);
      C=detail_exec::qualified_fd_as_xy(ctx,request.pert_c,request.freq_c);
      if(B && C) {
        source=detail_exec::load_shared_vbc(ctx,*selected->second);
      }
    }
    for(int axis=0;axis<3;++axis) {
      const auto A_id=fd_node_id(Perturbation::dipole(axis),request.freq_b+request.freq_c);
      if(!session.states.contains(A_id)) continue;
      std::string component=std::string(1,"xyz"[axis])+"xyz"[request.pert_b.axis];
      if(!nuclear) component+="xyz"[request.pert_c.axis];
      if(!source || !B || !C || !session.at(A_id).complete()) {
        slot("availability")[component]=supported ? missing_inputs(id,request,A_id)
            : "open-shell quadratic source and contraction are not implemented";
        continue;
      }
      auto A=detail_exec::qualified_fd_as_xy(ctx,Perturbation::dipole(axis),request.freq_b+request.freq_c);
      session.timing().begin("contraction",{{"task_id",id},{"component",component}});
      auto VA=dipole_operator(world,axis);
      const double value=beta::beta_abc<ClosedShell>(world,*ground,*A,*source,*B,*C,VA);
      session.timing().end();
      slot("components")[component]=value;
      if(world.rank()==0 && madness::iteration::prints(ctx.verbosity,madness::iteration::needed::assembled_values))
        print(nuclear ? "[RAMAN]" : "[BETA]",component,
          "fB=",request.freq_b,"fC=",request.freq_c,"value=",value);
      // Preserve the analysis mirror consumed by the existing Raman analysis.
      double residual=0.0;
      for(const auto& state_id:{A_id,fd_node_id(request.pert_b,request.freq_b),fd_node_id(request.pert_c,request.freq_c)})
        residual=std::max(residual,session.at(state_id).measured_values().at("values").at("bsh_residual").get<double>());
      (nuclear ? raman_rows : beta_rows).push_back({{"A",std::string(1,"xyz"[axis])},
          {"B",request.pert_b.description()},{"C",request.pert_c.description()},
          {"freq_b",request.freq_b},{"freq_c",request.freq_c},{"beta",value},
          {"converged",true},{"max_bsh_residual",residual}});
    }
  }
  auto publish=[&](const char* name,const char* units,auto& tensors) {
    if(tensors.empty()) return;
    json values=json::array();
    for(auto& [key,tensor]:tensors) {
      tensor["status"]=tensor.at("availability").empty() ? "complete" : tensor.at("components").empty() ? "unavailable" : "partial";
      values.push_back(std::move(tensor));
    }
    // Ordered by frequency, as the keys' text does not order numbers.
    std::stable_sort(values.begin(),values.end(),[](const json& a,const json& b) {
      const auto key=[](const json& v) {
        return std::make_pair(v.value("frequency_b",v.value("frequency",0.0)),v.value("frequency_c",0.0));
      };
      return key(a)<key(b);
    });
    session.set_properties({{name,{{"units",units},{"value_type","real"},{"values",values}}}});
  };
  publish("hyperpolarizability","atomic_units",beta_tensors);
  publish("polarizability_gradient","atomic_units",raman_tensors);
  checkpoint::detail::on_root(world,[&] {
    auto meta=ResponseMetadata::load_or_create(ctx.calc_dir+"/response_metadata.json");
    const auto key=protocol_key(thresh,default_k_for_thresh(thresh));
    meta.replace_property("beta",key,beta_rows);meta.replace_property("raman",key,raman_rows);
    meta.save();return std::string{};
  });
  session.timing().end();
}

inline void assemble_beta(ExecutorContext &ctx, const ResponsePlan &plan,
                          double thresh) {
  if (plan.vbc.empty()) return;
  assemble_shared_beta(ctx,plan,thresh);
}

// ---------------------------------------------------------------------------
// assemble_tpa — Tier-B two-photon-absorption assembly (post-solve, off the
// critical path; mirrors assemble_beta). Needs a converged ES bundle + the
// derived dipole FD at omega_f/2 (both produced by the resonant/ES plan). For
// each root f: load X_f + the three mu_a responses at omega_f/2, contract via
// tpa::tpa_moment (the beta-residue with a homogeneous C-channel), record the
// S tensor + delta^|| under properties/tpa. ClosedShell/TDA only.
// The residue form is a CANDIDATE validated against refs/dalton_tpa.json.
// ---------------------------------------------------------------------------
/// Load the ES bundle of EsType and return per-root (omega, X_f as XY). Full
/// (RPA) carries the real de-excitation Y; TDA sets Y=0. Templated so the
/// EsType Storage (X-only vs X,Y) is handled at compile time.
template <typename EsType>
inline bool load_tpa_es_xy(madness::World &world,
                           std::vector<double> &omega_out,
                           std::vector<ResponseStateXY<ClosedShell>> &Xf_out,
                           response_shared::ExcitedSet& binding,
                           std::vector<std::string>& root_ids) {
  if(!binding.complete()) return false;
  const auto state=binding.load([&](const auto& directory,const auto& description) {
    return response_shared::load_excited_functions<EsType,ClosedShell>(world,directory,description,binding.controls().protocol);
  },true);
  const long nr = state.omega.dim(0);
  for (long f = 0; f < nr; ++f) {
    omega_out.push_back(state.omega(f));
    root_ids.push_back(make_root_id(state.stable_index[f]));
    ResponseStateXY<ClosedShell> Xf;
    Xf.x_alpha = state.roots[f].x_alpha;
    if constexpr (std::is_same_v<EsType, Full>) {
      Xf.y_alpha = state.roots[f].y_alpha;               // RPA de-excitation
    } else {
      Xf.y_alpha = madness::copy(world, Xf.x_alpha);      // TDA: y = 0
      madness::scale(world, Xf.y_alpha, 0.0);
    }
    Xf_out.push_back(std::move(Xf));
  }
  return true;
}

/// The definitions of the two-photon quantities results.json holds per root
/// (tpa::observables, madness/chem/DerivedQuantities.h). The cross-section
/// convention is DALTON's, written so a reader can convert to another one.
inline nlohmann::json tpa_definitions() {
  return {
      {"Df","(1/30) sum_ij S_ii S_jj"},
      {"Dg","(1/30) sum_ij S_ij S_ij; equals (1/30) sum_ij S_ij S_ji since S is symmetric"},
      {"D_linear","2 Df + 4 Dg: the strength for two photons of the same linear polarization "
                  "(Monson and McClain 1970; F = G = H = 2)"},
      {"D_circular","-2 Df + 6 Dg: the strength for two photons of the same circular polarization "
                    "(F = -2, G = H = 3)"},
      {"R","(-Df + 3 Dg) / (Df + 2 Dg), the circular-to-linear ratio; 0 when the denominator is 0"},
      {"sigma_linear_gm","sigma = N pi^2 alpha a0^5 omega^2 D_linear g / c in GM (1 GM = 1e-50 cm^4 s per photon), "
                         "with N = 8 (prefactor 8 pi^2), omega = E/2 the photon energy (E the excitation energy), "
                         "and lineshape g = 1/Gamma at resonance, a constant, not a normalized Lorentzian; "
                         "linewidth Gamma = 0.1 eV, a full width at half maximum. This is DALTON's AU_TO_GM"},
      {"sigma_circular_gm","as sigma_linear_gm with D_circular"}};
}

inline void assemble_tpa(ExecutorContext &ctx, const ResponsePlan &plan,
                         double thresh) {
  using namespace madness;
  World &world = ctx.world;
  GroundState &gs = ctx.gs;
  if (plan.es.empty()) return;
  if (!gs.is_spin_restricted()) {
    if (world.rank() == 0)
      print("[TPA] open-shell TPA not implemented — SKIPPED (ClosedShell only).");
    return;
  }

  nlohmann::json shared_values=nlohmann::json::array();
  auto finish_shared=[&] {
    nlohmann::json property{
        {"value_type","real"},{"scientific_status","preliminary"},
        {"units",{{"transition_tensor","atomic_units"},{"energy","hartree"}}},
        {"values",shared_values}};
    // The definitions travel with the numbers, once any root has them.
    if(std::any_of(shared_values.begin(),shared_values.end(),[](const auto& v){return v.contains("derived");})) {
      property["definitions"]=tpa_definitions();
      property["derived_units"]={{"Df","atomic_units"},{"Dg","atomic_units"},{"D_linear","atomic_units"},
          {"D_circular","atomic_units"},{"R","dimensionless"},{"sigma_linear_gm","GM"},{"sigma_circular_gm","GM"}};
    }
    ctx.shared_response->set_properties({{"two_photon_absorption",property}});
    ctx.shared_response->timing().end();
  };
  const auto& request=plan.es.front();
  const auto set_id=es_node_id(request.tda,request.n_roots);
  // A root's value comes from the set and the three dipole states at half its
  // excitation energy (none known until the set has converged).
  // The set counts through the root's own metrics.
  auto named=[&](nlohmann::json value,std::optional<double> half) {
    ctx.shared_response->name_state(value,set_id,value.at("root").get<std::string>());
    if(half) for(int a=0;a<3;++a) ctx.shared_response->name_state(value,fd_node_id(Perturbation::dipole(a),*half));
    return value;
  };
  auto unavailable=[&](const std::string& id,const std::string& reason,std::optional<double> half=std::nullopt) {
    shared_values.push_back(named({{"root",id},{"status","unavailable"},
        {"components",nlohmann::json::object()},{"availability",{{"two_photon_absorption",reason}}}},half));
  };
  // The folders of the dipole states at `half` that have not met the request.
  auto unmet_half=[&](double half) {
    std::string text;
    for(int a=0;a<3;++a) {
      const auto id=fd_node_id(Perturbation::dipole(a),half);
      const auto& session=*ctx.shared_response;
      if(!session.complete(id))
        text+=(text.empty() ? "" : ", ")+(session.knows(id) ? session.node_name(id) : id);
    }
    return (text.empty() ? std::string("a dipole state at half the root energy") : text)+" has not met the requested criteria";
  };
  ctx.shared_response->timing().begin("two_photon_properties");
  auto& shared_bundle=ctx.shared_response->excited_at(set_id);
  if(!shared_bundle.complete()) {
    for(const auto& id:shared_bundle.expected().states) unavailable(id,ctx.shared_response->node_name(set_id)+" has not met the requested criteria");
    finish_shared();return;
  }
  set_response_protocol(world, ctx.L, thresh);
  const double t0 = FunctionDefaults<3>::get_thresh();
  {
    auto coulop = poperatorT(
        CoulombOperatorPtr(world, gs.params().lo(), 0.001 * t0));
    gs.prepare(world, 0.001 * t0, coulop, ctx.fock_json);
  }
  const double c_xc = gs.hf_exchange_coefficient();
  const double lo   = gs.params().lo();
  auto g0 = build_response_ground_state_closed_shell(world, gs, c_xc, lo);
  const std::string key = protocol_key();

  // The ES method of the requested bundle (tda vs full/RPA). Dalton's TPA
  // reference is full RPA.
  const int es_full=request.tda ? 0 : 1;
  std::vector<std::string> root_ids;
  std::vector<double> omegas;
  std::vector<ResponseStateXY<ClosedShell>> Xfs;
  const double t_es0 = madness::wall_time();
  const bool loaded =
      es_full ? load_tpa_es_xy<Full>(world, omegas, Xfs, shared_bundle, root_ids)
              : load_tpa_es_xy<TDA>(world, omegas, Xfs, shared_bundle, root_ids);
  if (world.rank() == 0) {
    printf("[TPA timing] ES bundle load: %.1f s\n", madness::wall_time() - t_es0);
    fflush(stdout);
  }
  if (!loaded) {
    if (world.rank() == 0)
      print("[TPA] no ES bundle under", ctx.calc_dir, "— SKIPPED");
    finish_shared();return;
  }
  const long nroots = static_cast<long>(omegas.size());
  if (world.rank() == 0)
    print("[TPA] ES method:", es_full ? "Full (RPA)" : "TDA");

  std::array<real_function_3d, 3> mu_op{dipole_operator(world, 0),
                                        dipole_operator(world, 1),
                                        dipole_operator(world, 2)};

  if (world.rank() == 0)
    print("\n=== TPA assembly  protocol_key=", key, "  n_roots=", nroots, " ===");

  std::vector<nlohmann::json> rows;
  // rank-0 table data (S tensor + observables per root) for the Dalton-style print.
  std::vector<int>                     tbl_f;
  std::vector<double>                  tbl_w;
  std::vector<madness::Tensor<double>> tbl_S;
  std::vector<tpa::Observables>        tbl_o;
  for (long f = 0; f < nroots; ++f) {
    const double wf = omegas[f];
    const double wf_half = 0.5 * wf;
    const ResponseStateXY<ClosedShell> &Xf = Xfs[f];   // real Y for RPA, 0 for TDA

    std::array<ResponseStateXY<ClosedShell>, 3> mu_resp;
    bool ok = true;
    const double t_fd0 = madness::wall_time();
    for (int a = 0; a < 3; ++a) {
      auto r = detail_exec::qualified_fd_as_xy(ctx,Perturbation::dipole(a),wf_half);
      if (!r) { ok = false; break; }
      mu_resp[a] = std::move(*r);
    }
    if (world.rank() == 0)
      printf("[TPA timing] root %ld/%ld: 3 FD loads @ w=%.5f: %.1f s\n",
             f + 1, nroots, wf_half, madness::wall_time() - t_fd0);
             fflush(stdout);
    if (!ok) {
      if (world.rank() == 0)
        print("[TPA] root", f, " omega_f=", wf,
              " — missing dipole FD @", wf_half, " — skip");
      unavailable(root_ids[f],unmet_half(wf_half),wf_half);
      continue;
    }

    // ---- cross-code diagnostics (PrintLevel >= Normal):
    // vector norms, transition dipoles/oscillators, and alpha(omega_f/2) from
    // the SAME loaded states the contraction uses — line-for-line comparable
    // with the patched-DALTON QRSMONORM + QRLRVE output (TPA_SCOPING §5n/§5o).
    if (madness::iteration::prints(ctx.verbosity, madness::iteration::needed::assembled_values)) {
      const double xx = inner(Xf.x_alpha, Xf.x_alpha);
      const double yy = inner(Xf.y_alpha, Xf.y_alpha);
      if (world.rank() == 0)
        printf("  [TPA diag] root %ld  omega=%.8f (%.3f eV)\n"
               "    Xf norms: |x|=%.10f  |y|=%.10f  x2-y2=%.10f "
               "(DALTON C(EXCI): Z2-Y2=0.5 => ours/DALTON = sqrt2)\n",
               f, wf, madness::derived::excitation_energy_ev(wf), std::sqrt(xx), std::sqrt(yy),
               xx - yy);
      std::array<vecfuncT, 3> mu_amo;
      for (int a = 0; a < 3; ++a) mu_amo[a] = mul(world, mu_op[a], g0.amo, true);
      // transition dipole + oscillator (Parker-normalized: M = -sqrt2 * t)
      for (int a = 0; a < 3; ++a) {
        const double t = inner(mu_amo[a], Xf.x_alpha) + inner(mu_amo[a], Xf.y_alpha);
        if (world.rank() == 0 && std::abs(t) > 1e-6)
          printf("    transition dipole %c: t=%+.8f  M=-sqrt2*t=%+.8f  "
                 "osc=(2/3)w*2t^2=%.6f\n",
                 "xyz"[a], t, -std::sqrt(2.0) * t,
                 (2.0 / 3.0) * wf * 2.0 * t * t);
      }
      // alpha(omega_f/2) full 3x3 from the loaded FD states (validated formula)
      for (int a = 0; a < 3; ++a) {
        const double nx = inner(mu_resp[a].x_alpha, mu_resp[a].x_alpha);
        const double ny = inner(mu_resp[a].y_alpha, mu_resp[a].y_alpha);
        if (world.rank() == 0)
          printf("    N^%c(w/2) norms: |x|=%.10f  |y|=%.10f\n", "xyz"[a],
                 std::sqrt(nx), std::sqrt(ny));
      }
      if (world.rank() == 0) printf("    alpha(w=%.6f):", wf_half);
      for (int a = 0; a < 3; ++a) {
        const double al = -2.0 * (inner(mu_resp[a].x_alpha, mu_amo[a]) +
                                  inner(mu_resp[a].y_alpha, mu_amo[a]));
        if (world.rank() == 0) printf("  %c%c=%.6f", "xyz"[a], "xyz"[a], al);
      }
      if (world.rank() == 0) printf("   (DALTON: QRLRVE <<A;A>> at same w)\n");
    }

    madness::Tensor<double> S;
    ctx.shared_response->timing().begin("contraction",{{"state",root_ids[f]}});
    // PRODUCTION composition (TPA_SCOPING §5m + §5q, verified 6/6 vs DALTON
    // d-aug-QZ 2026-07-22 via verify_e3_k6.py):
    //   S = C_N * ( S_1e + S_E3corr ),   C_N = sqrt(2)
    // S_1e: mu-operator terms of V^{bc} only (tpa_moment_residue_1e; cheap).
    // S_E3corr: the two-electron part from the Parker-style state-free (P,Q)
    // source (symmetrized builder; it returns the sum of both orderings).
    // The printed text, prefactor included, is unchanged from before the
    // other compositions were removed.
    if (world.rank() == 0) {
      print("[TPA] contraction: single-residue, S = sqrt2*(1e + E3corr),"
            " 2e composition =", "Parker (P,Q) spec [production]",
            " prefactor =", 1.0);
      fflush(stdout);
    }
    const auto S_1e =
        tpa::tpa_moment_residue_1e(world, g0, Xf, mu_resp, mu_op, 1.0);
    madness::Tensor<double> S_e3c(3L, 3L);
    const double max_asym = 0.0;
    for (int b3 = 0; b3 < 3; ++b3) {
      for (int c3 = b3; c3 < 3; ++c3) {
        if (world.rank() == 0) {
          printf("  [TPA e3corr] root %ld pair %c%c%s (%s) ...\n", f,
                 "xyz"[b3], "xyz"[c3],
                 (b3 != c3 ? " (both orderings)" : ""), "Parker spec");
          fflush(stdout);
        }
        // 0.5*<f|P_sym,Q_sym> is one ordering in the units of the outer
        // sqrt2*(1e+E3corr).
        auto pq = source_spec::assemble_source(
            world, g0,
            tpa::tpa_pq_spec_sym(world, g0,
                                 mu_resp[static_cast<size_t>(b3)],
                                 mu_resp[static_cast<size_t>(c3)]));
        const double sp = inner(world, Xf.x_alpha, pq[0]).sum() +
                          inner(world, Xf.y_alpha, pq[1]).sum();
        S_e3c(b3, c3) = S_e3c(c3, b3) = 0.5 * sp;
      }
    }
    S = madness::copy(S_1e);
    S += S_e3c;
    S.scale(std::sqrt(2.0));
    if (world.rank() == 0 && madness::iteration::prints(ctx.verbosity, madness::iteration::needed::assembled_values))
      printf("  [TPA] root %ld assembled: S = sqrt2*(1e + E3corr)   "
             "b<->c max asym = %.2e\n", f, max_asym);
    const auto obs = tpa::observables(S, wf);
    {
      ctx.shared_response->timing().end();
      nlohmann::json components=nlohmann::json::object();
      for(int a=0;a<3;++a) for(int b=0;b<3;++b) components[std::string(1,"xyz"[a])+"xyz"[b]]=S(a,b);
      shared_values.push_back(named({{"root",root_ids[f]},{"energy",wf},{"frequency",wf_half},{"status","complete"},
          {"components",components},{"availability",nlohmann::json::object()},
          {"derived",{{"Df",obs.Df},{"Dg",obs.Dg},{"D_linear",obs.D_linear},{"D_circular",obs.D_circular},
              {"R",obs.R},{"sigma_linear_gm",obs.sigma_linear_gm},{"sigma_circular_gm",obs.sigma_circular_gm}}}},wf_half));
    }


    if (world.rank() == 0) {
      nlohmann::json Sj = nlohmann::json::array();
      for (int a = 0; a < 3; ++a) {
        nlohmann::json ra = nlohmann::json::array();
        for (int b = 0; b < 3; ++b) ra.push_back(S(a, b));
        Sj.push_back(ra);
      }
      rows.push_back({{"es_root_id", static_cast<int>(f)},
                      {"omega", wf},
                      {"omega_ev", madness::derived::excitation_energy_ev(wf)},
                      {"source_files", nlohmann::json::array()}, // sources are not saved
                      {"writer_nproc", static_cast<int>(world.size())},
                      {"S", Sj},
                      {"Df", obs.Df}, {"Dg", obs.Dg},
                      {"D_linear", obs.D_linear}, {"D_circular", obs.D_circular},
                      {"R", obs.R},
                      {"sigma_linear_gm", obs.sigma_linear_gm},
                      {"sigma_circular_gm", obs.sigma_circular_gm},
                      {"delta_parallel", obs.D_linear}});
      {
        nlohmann::json legs = nlohmann::json::array();
        bool all_conv = true; double max_res = 0.0;
        for (int a = 0; a < 3; ++a) {
          auto leg = nlohmann::json{
              {"converged",true},{"bsh_residual",ctx.shared_response->at(fd_node_id(Perturbation::dipole(a),wf_half)).measured_values().at("values").at("bsh_residual")}};
          all_conv = all_conv && leg.value("converged", false);
          max_res  = std::max(max_res, leg.value("bsh_residual", 0.0));
          legs.push_back(std::move(leg));
        }
        rows.back()["row_accuracy"]     = std::move(legs);
        rows.back()["converged"]        = all_conv;
        rows.back()["max_bsh_residual"] = max_res;
      }
      tbl_f.push_back(static_cast<int>(f)); tbl_w.push_back(wf);
      tbl_S.push_back(S); tbl_o.push_back(obs);
    }
    world.gop.fence();
  }

  // Dalton-style output (matches rspvec.F QRSMO): tensor S table + the
  // Monson-McClain Df/Dg/D(lin,circ)/sigma/R summary. No point-group symmetry
  // in MRA, so one manifold (no Sym column). Diffable against Dalton .out.
  if (world.rank() == 0 && !tbl_f.empty() &&
      madness::iteration::prints(ctx.verbosity, madness::iteration::needed::assembled_values)) {
    const double H2EV = madness::derived::excitation_hartree_to_ev;
    printf("\n                  +--------------------------------+\n");
    printf("                  | Two-photon transition tensor S |\n");
    printf("                  +--------------------------------+\n");
    printf("     -----------------------------------------------------------------\n");
    printf("      No  Energy(eV)     Sxx     Syy     Szz     Sxy     Sxz     Syz\n");
    printf("     -----------------------------------------------------------------\n");
    for (size_t r = 0; r < tbl_f.size(); ++r)
      printf("     %3d   %8.2f  %8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n",
             tbl_f[r] + 1, tbl_w[r] * H2EV, tbl_S[r](0, 0), tbl_S[r](1, 1),
             tbl_S[r](2, 2), tbl_S[r](0, 1), tbl_S[r](0, 2), tbl_S[r](1, 2));
    printf("     -----------------------------------------------------------------\n");
    printf("\n     D = 2*Df+4*Dg (Linear);  D = -2*Df+6*Dg (Circular)\n");
    printf("     Df = sum_ij S_ii*S_jj /30;   Dg = sum_ij S_ij^2 /30\n");
    printf("     sigma = D*(E/2)^2*AU_TO_GM (GM, 0.1 eV FWHM);  R = (-Df+3Dg)/(Df+2Dg)\n");
    printf("\n                   +-----------------------------------+\n");
    printf("                   |   Two-photon absorption summary   |\n");
    printf("                   +-----------------------------------+\n");
    printf("      No Energy(eV) Polarization      Df         Dg          D        sigma      R\n");
    printf("     ---------------------------------------------------------------------------------\n");
    for (size_t r = 0; r < tbl_f.size(); ++r) {
      const auto &o = tbl_o[r];
      printf("     %3d %8.2f  Linear    %10.3E %10.3E %10.3E %10.3E %6.2f\n",
             tbl_f[r] + 1, tbl_w[r] * H2EV, o.Df, o.Dg, o.D_linear, o.sigma_linear_gm, o.R);
      printf("     %3d %8.2f  Circular  %10.3E %10.3E %10.3E %10.3E %6.2f\n",
             tbl_f[r] + 1, tbl_w[r] * H2EV, o.Df, o.Dg, o.D_circular, o.sigma_circular_gm, o.R);
    }
    printf("     ---------------------------------------------------------------------------------\n");
    fflush(stdout);
  }

  checkpoint::detail::on_root(world,[&] {
    auto meta=ResponseMetadata::load_or_create(ctx.calc_dir+"/response_metadata.json");
    meta.replace_property("tpa",key,rows);meta.save();return std::string{};
  });
  finish_shared();
}

// ---------------------------------------------------------------------------
// Polarizability assembly (Tier-A, closed-shell). After the FD states converge,
// contract each converged response against the dipole perturbation vectors:
//
//   alpha_ij(omega) = -2 * ( <x_i | v_j> + <y_i | v_j> )
//
// where x_i,y_i are the converged FD response (direction i, frequency omega) and
// v_j = Q(mu_j * phi0) is the dipole perturbation source (direction j). For the
// static case y = x, so this reduces to -4<x|v> -- the
// same factors as molresponse_v2 PropertyManager::compute_alpha (-2 dynamic /
// -4 static). Printed as [ALPHA ...] and recorded under properties/alpha.
// ---------------------------------------------------------------------------
inline void assemble_alpha(ExecutorContext &ctx, const ResponsePlan &plan,
                           double thresh) {
  using namespace madness;
  World &world = ctx.world;
  GroundState &gs = ctx.gs;

  // Review finding (confirmed HIGH): the assembly readers below are hardcoded
  // to the ClosedShell archive layout; an open-shell state deserializes
  // MISALIGNED (crash or garbage alpha). The open-shell SOLVES are fine and
  // saved — refuse the assembly loudly instead of fabricating numbers.
  if (!gs.is_spin_restricted()) {
    if (world.rank() == 0)
      print("[ASSEMBLE] open-shell property assembly is not implemented yet — "
            "response states are solved and saved, but alpha assembly is "
            "SKIPPED (ClosedShell-only readers; see REVIEW_FINDINGS).");
    {
      std::set<int> axes;
      std::set<double> frequencies;
      for (const auto& request:plan.fd) if (request.pert.kind==Perturbation::Kind::Dipole) {
        axes.insert(request.pert.axis); frequencies.insert(request.freq);
      }
      nlohmann::json tensors=nlohmann::json::array();
      for (double frequency:frequencies) {
        nlohmann::json availability=nlohmann::json::object();
        for (int row:axes) for (int column:axes)
          availability[std::string{char('x'+row),char('x'+column)}]=
              "open-shell polarizability contraction is not implemented";
        nlohmann::json tensor={{"frequency",frequency},{"status","unavailable"},
            {"components",nlohmann::json::object()},{"availability",availability}};
        for (int axis:axes) ctx.shared_response->name_state(tensor,fd_node_id(Perturbation::dipole(axis),frequency));
        tensors.push_back(std::move(tensor));
      }
      ctx.shared_response->set_properties({{"polarizability",{{"units","bohr^3"},{"value_type","real"},{"values",tensors}}}});
    }
    return;
  }

  // Dipole directions + frequencies present in the plan.
  std::set<int>    axes;
  std::set<double> freqs;
  for (const auto &r : plan.fd)
    if (r.pert.kind == Perturbation::Kind::Dipole) {
      axes.insert(r.pert.axis);
      freqs.insert(r.freq);
    }
  if (axes.empty() || freqs.empty()) return;
  const std::vector<int>    ax(axes.begin(),  axes.end());
  const std::vector<double> fr(freqs.begin(), freqs.end());

  set_response_protocol(world, ctx.L, thresh);
  const double t0 = FunctionDefaults<3>::get_thresh();
  {
    auto coulop = poperatorT(
        CoulombOperatorPtr(world, gs.params().lo(), 0.001 * t0));
    gs.prepare(world, 0.001 * t0, coulop, ctx.fock_json);
  }
  const std::string key = protocol_key();

  // Dipole perturbation source per direction (built once at this protocol).
  std::map<int, vector_real_function_3d> v;
  for (int a : ax) v[a] = dipole_perturbation(world, gs, a);

  if (world.rank() == 0)
    madness::print("=== polarizability (alpha) assembly  protocol_key=", key,
                   " ===");

  {
    nlohmann::json tensors=nlohmann::json::array(),qualified=nlohmann::json::array();
    for (double frequency:fr) {
      nlohmann::json components=nlohmann::json::object(),availability=nlohmann::json::object();
      nlohmann::json matrix=nlohmann::json::array();
      double max_residual=0.0;
      for (int axis:ax) {
        auto& binding=ctx.shared_response->at(fd_node_id(Perturbation::dipole(axis),frequency));
        nlohmann::json row=nlohmann::json::array();
        if (!binding.complete()) {
          for (int column:ax) availability[std::string{char('x'+axis),char('x'+column)}]=binding.name()+" has not met the requested criteria";
          continue;
        }
        iteration::ElapsedTimer<> load_timer;
        ResponseStateXY<ClosedShell> functions;
        if (detail_exec::is_static_freq(frequency)) {
          auto x=binding.load<ResponseStateX<ClosedShell>>(true);
          functions.x_alpha=std::move(x.x_alpha);functions.y_alpha=functions.x_alpha;
        } else functions=binding.load<ResponseStateXY<ClosedShell>>(true);
        const double load_seconds=load_timer.seconds();
        ctx.shared_response->timing().append({{"name","load"},{"task_id",binding.info().task_id},{"seconds",load_seconds}});
        const auto frame=binding.saved_metrics();
        if (world.rank()==0) {
          // Reusable states bypass the solver node and first load here.
          if (binding.reusable() && iteration::prints(ctx.verbosity,iteration::needed::restart_summary))
            iteration::write_restart_summary(std::cout,"response",frame,binding.controls().protocol,
                iteration::assess(binding.controls(),binding.expected(),frame),true,response_shared::metric_labels());
          if (iteration::prints(ctx.verbosity,iteration::needed::checkpoint_timing))
            iteration::write_checkpoint_timing(std::cout,"response",binding.info().task_id,load_seconds);
        }
        max_residual=std::max(max_residual,*frame.states.front().metrics.at(iteration::Metric::bsh_residual).value);
        for (int column:ax) {
          const double value=-2.0*(inner(functions.x_alpha,v[column])+inner(functions.y_alpha,v[column]));
          components[std::string{char('x'+axis),char('x'+column)}]=value;row.push_back(value);
        }
        matrix.push_back(row);
      }
      nlohmann::json tensor={{"frequency",frequency},{"status",availability.empty() ? "complete" : components.empty() ? "unavailable" : "partial"},
                             {"components",components},{"availability",availability}};
      for (int axis:ax) ctx.shared_response->name_state(tensor,fd_node_id(Perturbation::dipole(axis),frequency));
      tensors.push_back(std::move(tensor));
      // Keep the legacy analysis mirror for qualified tensors only.
      if (availability.empty()) {
        std::string directions;for (int axis:ax) directions+=char('x'+axis);
        qualified.push_back(
            {{"omega",frequency},{"directions",directions},{"alpha",matrix},{"converged",true},{"max_bsh_residual",max_residual}});
      }
    }
    replace_alpha_properties(world,ctx.calc_dir+"/response_metadata.json",key,qualified);
    ctx.shared_response->set_properties({{"polarizability",{{"units","bohr^3"},{"value_type","real"},{"values",tensors}}}});
  }
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_CALC_CALC_EXECUTOR_HPP
