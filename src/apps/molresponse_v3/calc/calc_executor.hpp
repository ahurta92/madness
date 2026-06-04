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
#include "../kernels/beta.hpp"
#include "../kernels/vbc.hpp"
#include "../solvers/es_save_load.hpp"
#include "../solvers/es_solver.hpp"
#include "../solvers/fd_solver.hpp"
#include "../solvers/iterate_protocol.hpp"
#include "../solvers/response_metadata.hpp"
#include "../solvers/response_state.hpp"
#include "../solvers/vbc_save_load.hpp"

#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <optional>
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
  // Cross-type seed toggle (optional): a derived FD starts from its converged
  // ES-root vector rather than the perturbation source. Default OFF — measured
  // neutral-to-slightly-worse at ωₙ/2 (off-resonance, H2); the ES root is the
  // right guess only for a future AT-resonance (factor 1.0) derived FD. FUTURE:
  // make the seed selectable per node, and allow a MIXTURE of roots to target
  // in-between frequencies (not just a single root's vector).
  bool              seed_derived_from_es_root = false;
  // Excited-state (Full / TDA-warmup) solve settings — defaults for the Full
  // closed-shell path (random guess, 10 warmup iters, 2x oversample, KAIN).
  ESGuessMode       es_guess              = ESGuessMode::SolidHarmonics;  // sweep-validated default
  int               es_tda_warmup_iters   = 10;
  double            es_warmup_oversample  = 2.0;
  int               es_kain_maxsub        = 8;
  double            es_maxrotn            = 0.5;
  // Delay KAIN onset in the MAIN ES solve by this many iters (pure BSH +
  // step-restriction first, so roots stabilize before KAIN starts recording
  // history). 0 = KAIN from iter 1 (previous behaviour). Distinct from
  // es_tda_warmup_iters, which is the separate oversampled-warmup pre-pass.
  int               es_main_kain_delay    = 0;
  // Full-deflation locking of converged ES roots (workstream B). ON by default
  // (sweep-validated: required to converge near-degenerate h2o/c2h4 clusters).
  bool              es_lock_converged     = true;
  // Cache the (promoted) TDA warmup guess to a side bundle es_warmup__<pkey> and
  // reuse it on later runs, so iterating on the Full main solve skips the warmup.
  // Opt-in; the main solve never overwrites the cache (stable fixed start point).
  bool              es_warmup_cache       = false;
};

// ---------------------------------------------------------------------------
namespace detail_exec {

inline bool is_static_freq(double freq) { return std::abs(freq) < 1e-12; }

/// Parse a derived-FD provenance label "es_root_NNNN" -> NNNN; -1 otherwise.
inline int es_root_index(const std::string &label) {
  const std::string pre = "es_root_";
  if (label.rfind(pre, 0) != 0) return -1;
  try { return std::stoi(label.substr(pre.size())); } catch (...) { return -1; }
}

/// Load a converged first-order state as an (X,Y) pair for the VBC build.
/// Dynamic (freq != 0): the Full FD state directly. Static (freq == 0): the
/// Static FD state with Y = X (the static limit). nullopt if not on disk.
template <class Shell>
inline std::optional<ResponseStateXY<Shell>>
load_fd_as_xy(madness::World &world, const std::string &calc_dir,
              const Perturbation &pert, double freq) {
  if (is_static_freq(freq)) {
    auto r = try_load_fd_state<Static, Shell>(world, calc_dir, pert, freq);
    if (!r) return std::nullopt;
    ResponseStateXY<Shell> xy;
    xy.x_alpha = madness::copy(world, r->state.responses[0].x_alpha);
    xy.y_alpha = madness::copy(world, r->state.responses[0].x_alpha);  // static: Y = X
    if constexpr (std::is_same_v<Shell, OpenShell>) {
      xy.x_beta = madness::copy(world, r->state.responses[0].x_beta);
      xy.y_beta = madness::copy(world, r->state.responses[0].x_beta);
    }
    return xy;
  }
  auto r = try_load_fd_state<Full, Shell>(world, calc_dir, pert, freq);
  if (!r) return std::nullopt;
  return r->state.responses[0];
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

/// Solve (pert, freq) at a single protocol `thresh`. `action` is the reconcile
/// verdict: Fresh starts from the perturbation guess (never loads, so a
/// diverged archive can't poison the solve); Restart/Resume load the nearest
/// converged-or-partial seed via try_load_fd_state and re-project it in the
/// first prepare(). Saves the result (+ metrics) through save_fd_state.
template <typename Type, typename Shell>
NodeResult solve_fd_protocol(ExecutorContext &ctx, const Perturbation &pert,
                         double freq, double thresh, NodeAction action,
                         const std::string &es_root_id = std::string()) {
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

  // Initial-guess precedence (the source carrier IS a valid Storage):
  //   1. disk FD partial at this (pert, freq) — restart (skipped if Fresh);
  //   2. else the converged ES-root vector if this is a derived FD and the
  //      cross-type seed is enabled (optional, ctx.seed_derived_from_es_root);
  //   3. else the perturbation source.
  typename Solver::State s0;
  s0.responses.resize(1);
  s0.responses[0] = tgt.responses[0].source;
  const char *seed_kind = "source";

  bool seeded = false;
  if (action != NodeAction::Fresh) {
    auto loaded =
        try_load_fd_state<Type, Shell>(world, ctx.calc_dir, pert, freq);
    if (loaded) {
      s0 = std::move(loaded->state);
      seeded = true;
      seed_kind = "fd_restart";
    }
  }
  // Cross-type seed: a derived FD (es_root_id set) starts from its converged
  // ES-root vector. The root is ResponseStateX (X only); for FD Full we seed
  // x = y = X. Only Full/ClosedShell derived FDs use it, and only when enabled.
  if constexpr (std::is_same_v<Type, Full> &&
                std::is_same_v<Shell, ClosedShell>) {
    if (!seeded && ctx.seed_derived_from_es_root && !es_root_id.empty()) {
      const int idx = detail_exec::es_root_index(es_root_id);
      if (idx >= 0) {
        auto esb = try_load_es_bundle<TDA, ClosedShell>(world, ctx.calc_dir);
        if (esb && idx < static_cast<int>(esb->state.roots.size())) {
          s0.responses[0].x_alpha =
              madness::copy(world, esb->state.roots[idx].x_alpha);
          s0.responses[0].y_alpha =
              madness::copy(world, esb->state.roots[idx].x_alpha);
          seeded = true;
          seed_kind = "es_root";
        }
      }
    }
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
                               /*converged=*/converged_now(st, solv),
                               /*seed=*/seed_kind);
  };

  solvers::IterateProtocolPolicy pp;
  pp.max_iters_per_step = ctx.max_iters;
  const std::vector<double> one_protocol = {thresh};
  auto sf = solvers::iterate_protocol(solver, std::move(s0), one_protocol,
                                      prepare, post_step, pp);

  NodeResult r;
  r.converged = converged_now(sf, solver);
  r.reached_protocol_key = protocol_key();  // active defaults reflect this protocol step
  if (world.rank() == 0)
    madness::print("[CALC] fd solve: pert=", pert.description(), " freq=", freq,
                   " thresh=", thresh, " seed=", seed_kind, " iters=", sf.iter,
                   " converged=", r.converged);
  return r;
}

// ---------------------------------------------------------------------------
// Excited-state bundle solve (15b-i: TDA closed-shell only).
// ---------------------------------------------------------------------------

/// Solve a TDA closed-shell excited-state bundle of `n_roots` at one protocol
/// step. Mirrors solve_fd_protocol: set protocol -> prepare GS -> build problem
/// -> guess (Fresh = fresh solid-harmonics trial; else nearest on-disk bundle
/// via the restart precedence) -> iterate_protocol over the single step (guarded
/// set_gs/reproject prepare + save_es_roots post-step) -> return the converged
/// excitation energies in NodeResult::es_root_freqs (which drive DerivedFD
/// expansion). The bundle saves to a per-protocol subdir `es__<key>` under
/// calc_dir so lower-protocol restart works.
inline NodeResult solve_es_tda_closed_shell(ExecutorContext &ctx, int n_roots,
                                            double thresh, NodeAction action) {
  using namespace madness;
  using Solver = ESSolver<TDA, ClosedShell>;
  World &world = ctx.world;
  GroundState &gs = ctx.gs;

  set_response_protocol(world, ctx.L, thresh);
  const double t0 = FunctionDefaults<3>::get_thresh();
  {
    auto coulop = poperatorT(
        CoulombOperatorPtr(world, gs.params().lo(), 0.001 * t0));
    gs.prepare(world, 0.001 * t0, coulop, ctx.fock_json);
  }

  const double c_xc = gs.hf_exchange_coefficient();
  const double lo   = gs.params().lo();
  auto problem = build_es_problem_tda<ClosedShell>(world, gs, n_roots, c_xc, lo);

  // Warmup policy (KAIN on, kain_maxsub/maxrotn from ctx); the main solve
  // disables the warmup window — the fresh guess is an oversampled TDA warmup
  // (>= 2x roots by default), same recipe as the Full path.
  ConvergencePolicy warm_policy = ctx.policy;
  warm_policy.kain                     = true;
  warm_policy.kain_maxsub              = ctx.es_kain_maxsub;
  warm_policy.maxrotn                  = ctx.es_maxrotn;
  warm_policy.tda_warmup_iters         = ctx.es_tda_warmup_iters;
  warm_policy.warmup_oversample_factor = ctx.es_warmup_oversample;
  ConvergencePolicy main_policy = warm_policy;
  main_policy.tda_warmup_iters = ctx.es_main_kain_delay;
  main_policy.lock_converged   = ctx.es_lock_converged;  // warmup never locks

  Solver::State s0;
  bool seeded = false;
  if (action != NodeAction::Fresh) {
    auto loaded = try_load_es_bundle<TDA, ClosedShell>(world, ctx.calc_dir);
    if (loaded) { s0 = std::move(loaded->state); seeded = true; }
  }
  if (!seeded) {
    const long n_warm = std::max<long>(
        n_roots, static_cast<long>(
                     std::ceil(ctx.es_warmup_oversample *
                               static_cast<double>(n_roots))));
    s0 = run_oversampled_tda_warmup<ClosedShell>(
        world, gs, n_roots, n_warm, ctx.es_tda_warmup_iters, warm_policy,
        c_xc, lo, ctx.print_level, ctx.es_guess);
  }

  Solver solver(world, std::move(problem), main_policy, ctx.print_level);

  // Single-step guard (as in solve_fd_protocol): re-prepare the ground state
  // only on an actual protocol change; otherwise just re-project the roots.
  double prepared_t = t0;
  auto prepare = [&](double th, Solver &solv, Solver::State &st) {
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
      solv.set_gs(build_response_ground_state_closed_shell(
          world, gs, gs.hf_exchange_coefficient(), gs.params().lo()));
      prepared_t = new_t;
    }
    const int    k  = FunctionDefaults<3>::get_k();
    const double tt = FunctionDefaults<3>::get_thresh();
    for (auto &root : st.roots)
      for (auto &fn : root.x_alpha) fn = madness::project(fn, k, tt);
    for (auto &rho : st.rho_alpha_prev) rho = madness::project(rho, k, tt);
  };

  // Report/save convergence with the SAME criterion the iteration stops on
  // (ESSolver::converged = energy+density, not the jittering BSH amplitude).
  auto converged_now = [](const Solver::State &st, const Solver &sv) {
    return !st.diverged && sv.converged(st);
  };

  auto post_step = [&](double, Solver &solv, Solver::State &st) {
    const std::string dir = ctx.calc_dir + "/es__" + protocol_key();
    save_es_roots<TDA, ClosedShell>(world, st, dir,
                                    /*converged=*/converged_now(st, solv));
  };

  solvers::IterateProtocolPolicy pp;
  pp.max_iters_per_step = ctx.max_iters;
  const std::vector<double> one_protocol = {thresh};
  auto sf = solvers::iterate_protocol(solver, std::move(s0), one_protocol,
                                      prepare, post_step, pp);

  NodeResult r;
  r.converged = converged_now(sf, solver);
  r.reached_protocol_key = protocol_key();
  for (long i = 0; i < sf.omega.size(); ++i)
    r.es_root_freqs.push_back(sf.omega(i));
  return r;
}

// ---------------------------------------------------------------------------
// Excited-state bundle solve — Full (paired X,Y) closed-shell (15b-ii).
// ---------------------------------------------------------------------------

/// Solve a Full (paired X,Y) closed-shell ES bundle of `n_roots` at one protocol
/// step via the direct ESSolver<Full, ClosedShell>. Fresh guess: a random-guess
/// TDA warmup (es_tda_warmup_iters, oversampled by es_warmup_oversample)
/// promoted to (X, Y=0) — the Full solver couples X and Y from iter 1. Restart:
/// load the (X,Y) bundle directly (no promote). Saves a proper "full" (X,Y)
/// bundle. Convergence is residual-based (same fields as TDA).
inline NodeResult solve_es_full_closed_shell(ExecutorContext &ctx, int n_roots,
                                             double thresh, NodeAction action) {
  using namespace madness;
  using Solver = ESSolver<Full, ClosedShell>;
  World &world = ctx.world;
  GroundState &gs = ctx.gs;

  set_response_protocol(world, ctx.L, thresh);
  const double t0 = FunctionDefaults<3>::get_thresh();
  {
    auto coulop = poperatorT(
        CoulombOperatorPtr(world, gs.params().lo(), 0.001 * t0));
    gs.prepare(world, 0.001 * t0, coulop, ctx.fock_json);
  }

  const double c_xc = gs.hf_exchange_coefficient();
  const double lo   = gs.params().lo();
  auto problem = build_es_problem_full<ClosedShell>(world, gs, n_roots, c_xc, lo);

  // Warmup policy (KAIN on, kain_maxsub / maxrotn from ctx); the main solve
  // disables the TDA warmup window (it was done explicitly for the fresh guess).
  ConvergencePolicy warm_policy = ctx.policy;
  warm_policy.kain                     = true;
  warm_policy.kain_maxsub              = ctx.es_kain_maxsub;
  warm_policy.maxrotn                  = ctx.es_maxrotn;
  warm_policy.tda_warmup_iters         = ctx.es_tda_warmup_iters;
  warm_policy.warmup_oversample_factor = ctx.es_warmup_oversample;
  ConvergencePolicy main_policy = warm_policy;
  main_policy.tda_warmup_iters = ctx.es_main_kain_delay;
  main_policy.lock_converged   = ctx.es_lock_converged;  // warmup never locks

  Solver::State s0;
  bool seeded = false;
  if (action != NodeAction::Fresh) {
    auto loaded = try_load_es_bundle<Full, ClosedShell>(world, ctx.calc_dir);
    if (loaded) { s0 = std::move(loaded->state); seeded = true; }
  }
  if (!seeded) {
    const std::string warm_cache =
        ctx.calc_dir + "/es_warmup__" + protocol_key();
    bool used_cache = false;
    if (ctx.es_warmup_cache) {
      int have = 0;
      if (world.rank() == 0)
        have = std::filesystem::exists(warm_cache + "/roots.json") ? 1 : 0;
      world.gop.broadcast(have, 0);
      if (have) {
        try {
          s0 = load_es_roots<Full, ClosedShell>(world, warm_cache);
          used_cache = true;
          if (world.rank() == 0)
            madness::print("[CALC] solve_es_full: reusing cached TDA warmup guess"
                           " from", warm_cache, "(skipping warmup)");
        } catch (const std::exception &e) {
          if (world.rank() == 0)
            madness::print("[CALC] solve_es_full: cached warmup unusable (",
                           e.what(), ") — running a fresh warmup");
        }
      }
    }
    if (!used_cache) {
      const long n_warm = std::max<long>(
          n_roots, static_cast<long>(
                       std::ceil(ctx.es_warmup_oversample *
                                 static_cast<double>(n_roots))));
      auto tda = run_oversampled_tda_warmup<ClosedShell>(
          world, gs, n_roots, n_warm, ctx.es_tda_warmup_iters, warm_policy,
          c_xc, lo, ctx.print_level, ctx.es_guess);
      s0 = promote_tda_to_full_closed_shell(world, tda);
      if (ctx.es_warmup_cache) {
        save_es_roots<Full, ClosedShell>(world, s0, warm_cache,
                                         /*converged=*/false,
                                         /*register_aggregate=*/false);
        if (world.rank() == 0)
          madness::print("[CALC] solve_es_full: cached TDA warmup guess ->",
                         warm_cache);
      }
    }
  }

  Solver solver(world, std::move(problem), main_policy, ctx.print_level);

  double prepared_t = t0;
  auto prepare = [&](double th, Solver &solv, Solver::State &st) {
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
      solv.set_gs(build_response_ground_state_closed_shell(world, gs, c_xc, lo));
      prepared_t = new_t;
    }
    const int    k  = FunctionDefaults<3>::get_k();
    const double tt = FunctionDefaults<3>::get_thresh();
    for (auto &root : st.roots)
      for (auto *blk : root.blocks())
        for (auto &fn : *blk) fn = madness::project(fn, k, tt);
    for (auto &rho : st.rho_alpha_prev) rho = madness::project(rho, k, tt);
  };

  // Report/save convergence with the SAME criterion the iteration stops on
  // (ESSolver::converged = energy+density, not the jittering BSH amplitude).
  auto converged_now = [](const Solver::State &st, const Solver &sv) {
    return !st.diverged && sv.converged(st);
  };

  auto post_step = [&](double, Solver &solv, Solver::State &st) {
    const std::string dir = ctx.calc_dir + "/es__" + protocol_key();
    save_es_roots<Full, ClosedShell>(world, st, dir,
                                     /*converged=*/converged_now(st, solv));
  };

  solvers::IterateProtocolPolicy pp;
  pp.max_iters_per_step = ctx.max_iters;
  const std::vector<double> one_protocol = {thresh};
  auto sf = solvers::iterate_protocol(solver, std::move(s0), one_protocol,
                                      prepare, post_step, pp);

  NodeResult r;
  r.converged = converged_now(sf, solver);
  r.reached_protocol_key = protocol_key();
  for (long i = 0; i < sf.omega.size(); ++i)
    r.es_root_freqs.push_back(sf.omega(i));
  return r;
}

// ---------------------------------------------------------------------------
// VBC quadratic-source build (beta-ii-b, closed-shell). Not iterative: load the
// two converged first-order states at this protocol, build the source, save.
// ---------------------------------------------------------------------------

inline NodeResult
solve_vbc_closed_shell(ExecutorContext &ctx, const CalcNode &node, double thresh) {
  using namespace madness;
  World &world = ctx.world;
  GroundState &gs = ctx.gs;

  set_response_protocol(world, ctx.L, thresh);
  const double t0 = FunctionDefaults<3>::get_thresh();
  {
    auto coulop = poperatorT(
        CoulombOperatorPtr(world, gs.params().lo(), 0.001 * t0));
    gs.prepare(world, 0.001 * t0, coulop, ctx.fock_json);
  }

  // The two converged first-order states (X,Y) at this protocol. The
  // prerequisites_converged gate guarantees they are converged before we run.
  auto B = detail_exec::load_fd_as_xy<ClosedShell>(world, ctx.calc_dir, node.pert,   node.freq);
  auto C = detail_exec::load_fd_as_xy<ClosedShell>(world, ctx.calc_dir, node.pert_c, node.freq_c);
  if (!B || !C) {
    if (world.rank() == 0)
      madness::print("[CALC] solve_vbc:", node.id,
                     "— a converged FD prerequisite is missing; skipping.");
    return NodeResult{/*converged=*/false, /*reached_protocol_key=*/"", {}};
  }

  const double c_xc = gs.hf_exchange_coefficient();
  const double lo   = gs.params().lo();
  auto g0 = build_response_ground_state_closed_shell(world, gs, c_xc, lo);

  // Raw one-electron perturbation operators for B and C (dipole only for now).
  auto VB_op = dipole_operator(world, node.pert.axis);
  auto VC_op = dipole_operator(world, node.pert_c.axis);

  auto vbc_src = vbc::compute_vbc<ClosedShell>(world, g0, *B, *C, VB_op, VC_op);
  save_vbc_state<ClosedShell>(world, vbc_src, ctx.calc_dir, node.id, /*converged=*/true);

  if (world.rank() == 0)
    madness::print("[CALC] solve_vbc: built", node.id, "at", protocol_key());
  return NodeResult{/*converged=*/true, protocol_key(), {}};
}

// ---------------------------------------------------------------------------
// FdResponseExecutor — dispatches the runtime (Type × Shell) for FD nodes.
// ---------------------------------------------------------------------------

class FdResponseExecutor : public ICalcExecutor {
public:
  explicit FdResponseExecutor(ExecutorContext ctx) : ctx_(ctx) {}

  NodeResult run_protocol(const WorkItem &item) override {
    const CalcNode &node = *item.node;

    if (node.kind == CalcKind::ES) {
      const bool restricted = ctx_.gs.is_spin_restricted();
      if (restricted) {
        if (ctx_.world.rank() == 0)
          madness::print("[CALC] run_protocol: ES bundle", node.id,
                         " n_roots=", node.n_roots, " thresh=", item.thresh,
                         (node.tda ? " (TDA closed-shell)"
                                   : " (Full closed-shell)"),
                         " action=", action_name(item.action));
        return node.tda
                   ? solve_es_tda_closed_shell(ctx_, node.n_roots, item.thresh,
                                               item.action)
                   : solve_es_full_closed_shell(ctx_, node.n_roots, item.thresh,
                                                item.action);
      }
      // Open-shell excited states are out of scope (closed-shell only; future).
      if (ctx_.world.rank() == 0)
        madness::print("[CALC] run_protocol: ES bundle", node.id,
                       "on an OPEN-SHELL ground state — open-shell excited "
                       "states are OUT OF SCOPE (closed-shell only; future "
                       "research). Skipping.");
      return NodeResult{/*converged=*/false, /*reached_protocol_key=*/"", {}};
    }
    if (node.kind == CalcKind::VBC) {
      if (ctx_.gs.is_spin_restricted()) {
        if (ctx_.world.rank() == 0)
          madness::print("[CALC] run_protocol: VBC", node.id,
                         " B=", node.pert.description(), "@", node.freq,
                         " C=", node.pert_c.description(), "@", node.freq_c,
                         " thresh=", item.thresh);
        return solve_vbc_closed_shell(ctx_, node, item.thresh);
      }
      if (ctx_.world.rank() == 0)
        madness::print("[CALC] run_protocol: VBC", node.id,
                       "on an open-shell ground state is out of scope. Skipping.");
      return NodeResult{/*converged=*/false, /*reached_protocol_key=*/"", {}};
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
                                                 item.thresh, item.action, node.es_root_id);
    if (restricted && !is_static)
      return solve_fd_protocol<Full, ClosedShell>(ctx_, node.pert, node.freq,
                                               item.thresh, item.action, node.es_root_id);
    if (!restricted && is_static)
      return solve_fd_protocol<Static, OpenShell>(ctx_, node.pert, node.freq,
                                               item.thresh, item.action, node.es_root_id);
    return solve_fd_protocol<Full, OpenShell>(ctx_, node.pert, node.freq,
                                          item.thresh, item.action, node.es_root_id);
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
      // Grow the DAG from any ES bundle converged ON DISK (idempotent and
      // restart-safe: reads roots from metadata, not in-memory run state).
      expand_converged_es(world, meta.json());
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
      for (const auto &item : waves.front())
        exec.run_protocol(item);
      world.gop.fence();
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

  /// Turn each symbolic DerivedFD "*" whose ES bundle has CONVERGED ON DISK into
  /// one promoted FD node per root (freq = root energy at the ES bundle's top
  /// protocol). Reads roots straight from response_metadata.json so it is
  /// restart-safe (works on a resumed run where the ES converged in a previous
  /// process) and idempotent (new nodes are deduped against existing dag_ ids;
  /// re-running is a no-op). Deterministic across ranks: same metadata + same
  /// dag_ in, same nodes appended in the same order.
  void expand_converged_es(madness::World &world, const nlohmann::json &meta) {
    if (!meta.contains("excited_states")) return;
    std::set<std::string> have;
    for (const auto &n : dag_) have.insert(n.id);

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
      const std::string es_key = protocol_key_at(es->protocols.back());
      if (!meta["excited_states"].contains(es_key)) continue;
      const auto &b = meta["excited_states"][es_key];
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
        c.es_root_id = make_es_root_label(static_cast<int>(i));  // provenance
        c.id         = fd_node_id(c.pert, w);
        if (have.count(c.id)) continue;                  // dedup -> idempotent
        have.insert(c.id);
        additions.push_back(std::move(c));
      }
    }
    if (!additions.empty()) {
      if (world.rank() == 0)
        madness::print("[CALC] expanded", (int)additions.size(),
                       "derived FD node(s) at ω = es_freq_factor·ωₙ from converged ES roots (metadata)");
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
};

// ---------------------------------------------------------------------------
// beta property assembly (Tier-A). A post-run step: contract the converged A
// states against each VBC source (+ zeta terms) into the beta tensor, recorded
// under properties/beta. Runs at one protocol (typically the top of the ramp).
// ---------------------------------------------------------------------------

inline char beta_axis_name(int a) { return a == 0 ? 'x' : a == 1 ? 'y' : 'z'; }

inline void assemble_beta(ExecutorContext &ctx, const ResponsePlan &plan,
                          double thresh) {
  using namespace madness;
  World &world = ctx.world;
  GroundState &gs = ctx.gs;
  if (plan.vbc.empty()) return;

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

  if (world.rank() == 0)
    madness::print("\n=== beta assembly  protocol_key=", key, " ===");

  for (const auto &vr : plan.vbc) {
    const double ws = vr.freq_b + vr.freq_c;   // omega_sigma = omega_B + omega_C
    const std::string vbc_id =
        vbc_node_id(vr.pert_b, vr.pert_c, vr.freq_b, vr.freq_c);

    auto vbc = load_vbc<ClosedShell>(world, ctx.calc_dir, vbc_id);
    auto B = detail_exec::load_fd_as_xy<ClosedShell>(world, ctx.calc_dir, vr.pert_b, vr.freq_b);
    auto C = detail_exec::load_fd_as_xy<ClosedShell>(world, ctx.calc_dir, vr.pert_c, vr.freq_c);
    if (!vbc || !B || !C) {
      if (world.rank() == 0)
        madness::print("[BETA] skip", vbc_id, "— missing VBC or FD input");
      continue;
    }

    for (int a = 0; a < 3; ++a) {
      const Perturbation pA = Perturbation::dipole(a);
      auto xA = detail_exec::load_fd_as_xy<ClosedShell>(world, ctx.calc_dir, pA, ws);
      if (!xA) {
        if (world.rank() == 0)
          madness::print("[BETA] skip A=", beta_axis_name(a), "— missing FD",
                         pA.description(), "@", ws);
        continue;
      }
      auto VA_op = dipole_operator(world, a);
      const double b = beta::beta_abc<ClosedShell>(world, g0, *xA, *vbc, *B, *C, VA_op);

      if (world.rank() == 0) {
        madness::print("[BETA]  A=", beta_axis_name(a),
                       "  B=", vr.pert_b.description(),
                       "  C=", vr.pert_c.description(),
                       "  fB=", vr.freq_b, "  fC=", vr.freq_c,
                       "  beta=", b);
        auto meta = ResponseMetadata::load_or_create(
            ctx.calc_dir + "/response_metadata.json");
        meta.add_property("beta", key,
                          nlohmann::json{{"A", std::string(1, beta_axis_name(a))},
                                         {"B", vr.pert_b.description()},
                                         {"C", vr.pert_c.description()},
                                         {"freq_b", vr.freq_b},
                                         {"freq_c", vr.freq_c},
                                         {"beta", b}});
        meta.save();
      }
      world.gop.fence();
    }
  }
}

// ---------------------------------------------------------------------------
// Polarizability assembly (Tier-A, closed-shell). After the FD states converge,
// contract each converged response against the dipole perturbation vectors:
//
//   alpha_ij(omega) = -2 * ( <x_i | v_j> + <y_i | v_j> )
//
// where x_i,y_i are the converged FD response (direction i, frequency omega) and
// v_j = Q(mu_j * phi0) is the dipole perturbation source (direction j). For the
// static case load_fd_as_xy supplies y = x, so this reduces to -4<x|v> -- the
// same factors as molresponse_v2 PropertyManager::compute_alpha (-2 dynamic /
// -4 static). Printed as [ALPHA ...] and recorded under properties/alpha.
// ---------------------------------------------------------------------------
inline void assemble_alpha(ExecutorContext &ctx, const ResponsePlan &plan,
                           double thresh) {
  using namespace madness;
  World &world = ctx.world;
  GroundState &gs = ctx.gs;

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

  for (double w : fr) {
    Tensor<double> alpha(static_cast<long>(ax.size()),
                         static_cast<long>(ax.size()));
    bool complete = true;
    for (size_t i = 0; i < ax.size() && complete; ++i) {
      auto Xi = detail_exec::load_fd_as_xy<ClosedShell>(
          world, ctx.calc_dir, Perturbation::dipole(ax[i]), w);
      if (!Xi) {
        if (world.rank() == 0)
          madness::print("[ALPHA] skip omega=", w, " — missing FD dipole_",
                         beta_axis_name(ax[i]));
        complete = false;
        break;
      }
      for (size_t j = 0; j < ax.size(); ++j) {
        const double aij =
            -2.0 * (inner(Xi->x_alpha, v[ax[j]]) + inner(Xi->y_alpha, v[ax[j]]));
        alpha(static_cast<long>(i), static_cast<long>(j)) = aij;
      }
    }
    if (!complete) continue;

    if (world.rank() == 0) {
      madness::print("[ALPHA] omega=", w, "  directions=",
                     [&]{ std::string d; for (int a : ax) d += beta_axis_name(a); return d; }());
      for (size_t i = 0; i < ax.size(); ++i)
        for (size_t j = 0; j < ax.size(); ++j)
          madness::print("[ALPHA]  alpha_", beta_axis_name(ax[i]),
                         beta_axis_name(ax[j]), " (omega=", w, ") =",
                         alpha(static_cast<long>(i), static_cast<long>(j)));
      madness::print("[ALPHA] tensor(omega=", w, ") =\n", alpha);

      auto meta = ResponseMetadata::load_or_create(
          ctx.calc_dir + "/response_metadata.json");
      nlohmann::json row;
      row["omega"] = w;
      std::string dirs; for (int a : ax) dirs += beta_axis_name(a);
      row["directions"] = dirs;
      nlohmann::json mat = nlohmann::json::array();
      for (size_t i = 0; i < ax.size(); ++i) {
        nlohmann::json r = nlohmann::json::array();
        for (size_t j = 0; j < ax.size(); ++j)
          r.push_back(alpha(static_cast<long>(i), static_cast<long>(j)));
        mat.push_back(r);
      }
      row["alpha"] = mat;
      meta.add_property("alpha", key, row);
      meta.save();
    }
    world.gop.fence();
  }
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_CALC_CALC_EXECUTOR_HPP
