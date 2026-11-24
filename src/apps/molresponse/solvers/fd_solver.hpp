#ifndef MOLRESPONSE_V3_SOLVERS_FD_SOLVER_HPP
#define MOLRESPONSE_V3_SOLVERS_FD_SOLVER_HPP

// =========================================================================
// FDSolver<Type, Shell> — frequency-dependent response, batched over
// (perturbation, frequency) responses at a single protocol step.
//
// step() shape vs ESSolver:
//
//   1. per-channel: V0x, E0x, gamma, rho            (no T0x, no E0x_full)
//   2. (skip — no Lambda assembly)
//   3. (skip — omega fixed, no subspace eigenproblem)
//   4. (skip — no rotation)
//   5. Theta = V0x - E0x + gamma + perturbation_source
//   6. BSH apply at fixed omega for this channel → new_x [+ new_y]
//      residual = || (x,y)_old - (x,y)_new ||
//   7. explosion guard + banner
//
// Same `ConvergencePolicy`, same `iterate` / `iterate_protocol`,
// same Storage-typed kernel signatures as ESSolver. The differences
// are concentrated here: no T0x / Lambda / sygv / rotation, BSH ops
// could be cached per response (deferred — built inside K::bsh_apply
// for now, same as ES).
//
// TODO(streaming): step() currently materializes V0x, E0x, gamma as
// separate vectors of Storage and assembles theta from them — peak
// is ~6M Storage instances in flight (in + out + V0x + E0x + gamma
// + theta). Refactor to stream the assembly per response: theta starts
// as V0x, then E0x is computed as a temporary and subtracted in
// place, then gamma is added in place, then perturbation_source is
// added. Peak drops to in + out + theta + ~1 temporary ≈ 2M + 3.
// This also sets up future per-response MacroTaskQ parallelism: each
// worker owns 1 response's buffers, not the whole bundle.
//
// TODO(VBC): for cubic response properties (β, dα/dQ for Raman) we
// need a Kernels<Type, Shell>::compute_second_order_vector kernel
// mirroring molresponse_v2/VBCMacrotask.hpp::compute_vbc_i. Pattern:
// given two first-order response vectors B and C, build ζ_BC =
// Σ y_B·x_C, apply compute_g (a 3-input J+K-like operator), apply
// V_pert·c_x, project onto virtuals, combine. Same shape but on 2×
// the storage. Plan the Kernels API so this plugs in cleanly.
// =========================================================================

#include "function_kinetics.hpp"
#include "function_metric_diagnostics.hpp"
#include "fd_explicit_accuracy.hpp"
#include "function_residuals.hpp"
#include "../kernels/exchange_ctx.hpp"  // exch::assemble_theta_tensor (--fd-tensor gate)
#include "../kernels/full.hpp"     // Kernels<Full, ClosedShell>
#include "../kernels/static.hpp"   // Kernels<Static, ClosedShell>
#include "../kernels/tags.hpp"
#include "convergence_policy.hpp"
#include "fd_problem.hpp"
#include "iterate.hpp"
#include "response_state.hpp"
#include "response_subspace_kain.hpp"
#include <madness/chem/MetricHelpers.h>

#include <madness/mra/mra.h>

#include <algorithm>
#include <fstream>
#include <string>
#include <type_traits>
#include <limits>
#include <utility>
#include <vector>

namespace molresponse_v3 {

template <typename Type, typename Shell>
class FDSolver {
public:
  using K        = Kernels<Type, Shell>;
  using Storage  = StorageOf_t<Type, Shell>;
  using Channel  = FDPerturbationOf_t<Type, Shell>;
  using vecfuncT = std::vector<madness::real_function_3d>;

  struct State {
    std::vector<Storage>                    responses;
    // Scalars only. History belongs to this solve/rung/reference and is reset
    // on continuation; the function archives remain the source of state.
    std::vector<FunctionKineticSample> kinetic_history;
    std::vector<FunctionResidualSample> bsh_history;
    std::vector<madness::FunctionResidualBlocks> function_residuals;
    std::vector<FunctionKinetics> function_kinetics;
    std::vector<double>                     last_bsh_residual;
    std::vector<madness::real_function_3d>  rho_alpha_prev;
    std::vector<double>                     last_density_residual;
    /// Per-channel ||theta||, the BSH source norm at the last iteration. This is
    /// the constant in the error bound the convergence guide states,
    /// |d alpha| <= |c| ||theta|| ||r_B||, which until now nothing in the code
    /// evaluated -- the bound was documented and unevaluable. Recorded per
    /// channel so the budget can be checked at run time.
    std::vector<double>                     last_theta_norm;
    /// Per-channel source overlap sum_i <v phi_i | x_i> (+ y block for Full):
    /// the diagonal property estimate of this channel, tracked per iteration
    /// so the property's own convergence can be read off the log (2026-09-11).
    std::vector<double>                     last_property;
    /// Per-channel ||x|| (after KAIN) and ||rho1||: the scales of the relative
    /// FD gate, ||dx|| <= atol + rtol ||x|| and ||drho|| <= atol + rtol ||rho1||
    /// (review C1). Empty for a state loaded from disk, which then falls back
    /// to the absolute gate.
    std::vector<double>                     last_x_norm;
    std::vector<double>                     last_rho_norm;
    /// Normalised distance per iteration (newest last), ONE TRACK PER GATED
    /// QUANTITY: [0] = max bsh/bsh_target, [1] = max drho/density_target.
    /// A stall needs every track flat (ConvergencePolicy::plateau); reset when
    /// the thresh changes.
    std::vector<std::vector<double>>        gate_history;
    double                                  gate_thresh = -1.0;
    int                                     iter = 0;
    bool                                    diverged = false;
    /// Set by step() when the residual has plateaued above the targets
    /// (ConvergencePolicy::stall_*). converged() returns true on it so the
    /// loop exits; the executor records it and applies best-effort acceptance.
    bool                                    stalled = false;
    /// How the last solve's rung ended (B5, rung_exit), set by iterate<>.
    madness::convergence::StopReason stop
        = madness::convergence::StopReason::in_progress;
    bool                                    rung_exhausted = false;
    nlohmann::json explicit_accuracy = nlohmann::json::array();
    nlohmann::json saved_explicit_accuracy = nlohmann::json::array();
    nlohmann::json progress_criteria_history = nlohmann::json::array();
    nlohmann::json configured_progress_history = nlohmann::json::array();
    nlohmann::json configured_progress = nlohmann::json::object();
  };

  FDSolver(madness::World &world, FDProblem<Type, Shell> target,
           ConvergencePolicy policy,
           PrintLevel print_level = PrintLevel::Normal,
           std::string log_prefix = {})
      : world_(world), target_(std::move(target)),
        policy_(policy), kain_(world, policy),
        print_level_(print_level), log_prefix_(std::move(log_prefix)) {
    record_function_metrics_ = !policy_.explicit_accuracy_targets.empty();
    MADNESS_CHECK_THROW(policy_.explicit_progress_configuration.empty() || target_.n_responses() == 1,
        "configured FD progress requires one named channel per solve");
    refresh_convergence_targets();
    print_header();
  }

  // -- accessors / config -------------------------------------------------
  madness::World&     world() const             { return world_; }
  PrintLevel          print_level() const       { return print_level_; }
  void                set_print_level(PrintLevel p) { print_level_ = p; }
  ConvergencePolicy   convergence_policy() const{ return policy_; }
  void                set_convergence_policy(ConvergencePolicy p) {
    policy_ = p;
    kain_.set_policy(p);
    refresh_convergence_targets();
  }
  void                set_target(FDProblem<Type, Shell> t) {
    target_ = std::move(t);
    // KAIN's stored iterates live at the old k/thresh; drop them on
    // protocol-driven target swaps.
    reset_shared();
  }
  const FDProblem<Type, Shell>& target() const { return target_; }

  /// Enable per-iter convergence logging to a CSV file. One row per
  /// (iter, channel) appended at the end of each step(). Pass empty
  /// string to disable.
  void set_log_path(const std::string &path) {
    log_path_ = path;
    log_header_written_ = false;
  }

  void refresh_convergence_targets() {
    targets_ = policy_.effective_for_thresh(
        madness::FunctionDefaults<3>::get_thresh());
  }
  const ConvergencePolicy::Targets& targets() const { return targets_; }
  void record_function_metrics(bool enabled) { record_function_metrics_ = enabled || !policy_.explicit_accuracy_targets.empty(); }

  /// <z|(A-w)|z> for the operator this solver actually solves, formed directly.
  ///
  /// Read off bsh_apply rather than assumed: the fixed point satisfies
  ///
  ///     (T - eps_i - w) x_i + (V0 x)_i - (E0 x)_i + gamma_i(x) = -v_i
  ///
  /// and theta(z) = V0 z - E0 z + gamma(z) + v collects every term but the
  /// kinetic and diagonal-energy ones, so
  ///
  ///     <z|(A-w)|z> = sum_i [ <z_i|T|z_i> - (eps_i + w) <z_i|z_i> ]
  ///                   + <z|theta(z)> - <z|v>
  ///
  /// No Green's operator and no projector enter, which is the point: the two
  /// earlier attempts at this both went wrong by inferring the quadratic form
  /// from the BSH map (first ignoring its level shift, then assuming Q passed
  /// through the inner product because theta was in the Q space -- compute_V0x
  /// returns V*x unprojected, so it is not). <z_i|T|z_i> comes from
  /// madness::Kinetic, which uses the gradient trick and never applies a
  /// Laplacian to an MRA function.
  ///
  /// Collective. Static only: at finite w the operator is symmetric in the
  /// paired metric, not in this inner product.
  double quadratic_form(const Storage &z, int r) {
    auto rho   = K::compute_density(world_, target_.gs, z);
    auto theta = assemble_theta_for(z, r, rho);
    madness::Kinetic<double, 3> T(world_);
    const auto zb = z.blocks();
    double acc = 0.0;
    for (std::size_t b = 0; b < zb.size(); ++b) {
      const auto &zv = *zb[b];
      if (zv.empty()) continue;
      std::vector<madness::real_function_3d> v(zv.begin(), zv.end());
      madness::Tensor<double> tm = T(v, v);            // collective
      madness::Tensor<double> ov = madness::inner(world_, v, v);
      const double w = target_.responses[r].omega;
      for (std::size_t i = 0; i < v.size(); ++i)
        acc += tm(static_cast<long>(i), static_cast<long>(i))
             - (target_.gs.aeps(static_cast<long>(i)) + w) * ov(static_cast<long>(i));
    }
    return acc + storage_inner(world_, z, theta)
               - source_overlap(world_, z, target_.responses[r]);
  }

  /// Stationary (Hylleraas) estimate of channel `r`'s source overlap <v|x>.
  ///
  /// The equation is (A-w)x = -v, so the functional stationary at the solution
  /// is S[z] = 2<v|z> + <z|(A-w)|z>, with S[x] = <v|x>. The naive estimate
  /// <v|z> is linear in z and therefore first order in the response error; S is
  /// stationary, so the SAME vector gives an error quadratic in it -- exactly
  /// -<dx|(A-w)|dx>, and exactly, because S is quadratic and the expansion
  /// terminates.
  ///
  /// Evaluated at the converged vector itself: no extra BSH cycle, one kernel
  /// build for theta(z) and one gradient pass for the kinetic term.
  double stationary_property(const State &s, int r) {
    if (r < 0 || r >= static_cast<int>(s.responses.size()))
      return std::numeric_limits<double>::quiet_NaN();
    if constexpr (!std::is_same_v<Type, Static>) {
      return std::numeric_limits<double>::quiet_NaN();
    } else {
      const Storage &z = s.responses[r];
      return 2.0 * source_overlap(world_, z, target_.responses[r])
           + quadratic_form(z, r);
    }
  }

  /// Is S actually stationary at the converged vector? Perturb it along a fixed
  /// direction by eps and 2*eps and compare the two deviations: a stationary
  /// functional gives a ratio near 4, a first-order one near 2.
  ///
  /// This is the check the two earlier attempts skipped. It needs no reference
  /// value, no sweep and no curve fit -- it interrogates the functional itself,
  /// on one converged state, so a wrong quadratic form is caught before any
  /// compute is spent on measuring its slope.
  /// NB the size of `eps`. The signal here is eps^2 <d|(A-w)|d>, recovered by
  /// cancellation from terms of order ||theta|| ~ 10; at eps = 1e-3 that signal
  /// is ~1e-6, which is the MRA representation noise at thresh = 1e-6, and the
  /// ratio then measures the noise rather than the functional (the first run of
  /// this check returned ~8, i.e. a cubic law, which a quadratic functional
  /// cannot produce). eps must be large enough that eps^2 clears the floor.
  double stationarity_ratio(const State &s, int r, double eps = 1e-1) {
    if constexpr (!std::is_same_v<Type, Static>) {
      return std::numeric_limits<double>::quiet_NaN();
    } else {
      const Storage &x = s.responses[r];
      const double S0 = stationary_property(s, r);
      // direction: the perturbation source, Q-projected like x and independent
      // of the error already in x
      auto dir = target_.responses[r].source;
      const double dn = storage_norm(world_, dir);
      if (!(dn > 0.0)) return std::numeric_limits<double>::quiet_NaN();
      auto probe = [&](double a) {
        State t = s;
        t.responses[r] = x;
        t.responses[r].axpy(world_, a / dn, dir);
        return stationary_property(t, r);
      };
      const double d1 = std::abs(probe(eps)      - S0);
      const double d2 = std::abs(probe(2.0 * eps) - S0);
      return (d1 > 0.0) ? d2 / d1 : std::numeric_limits<double>::quiet_NaN();
    }
  }

  /// The LINEAR coefficient of S about the converged vector, extracted by a
  /// central difference so no fitting is involved:
  ///
  ///     S(+eps) - S(-eps) = 2 eps alpha,   alpha = 2 <d | (A-w) x + v>
  ///
  /// alpha is twice the projected residual contracted with the probe direction
  /// d (which lies in the Q space, being the perturbation source). If the
  /// quadratic form matched the equation the solver solves, alpha would BE the
  /// residual overlap and would fall with it. Reported beside ||r_B|| so the
  /// two can be compared directly at several convergence levels.
  double linear_coefficient(const State &s, int r, double eps = 1e-1) {
    if constexpr (!std::is_same_v<Type, Static>) {
      return std::numeric_limits<double>::quiet_NaN();
    } else {
      const Storage &x = s.responses[r];
      auto dir = target_.responses[r].source;
      const double dn = storage_norm(world_, dir);
      if (!(dn > 0.0)) return std::numeric_limits<double>::quiet_NaN();
      auto probe = [&](double a) {
        State t = s;
        t.responses[r] = x;
        t.responses[r].axpy(world_, a / dn, dir);
        return stationary_property(t, r);
      };
      return (probe(eps) - probe(-eps)) / (2.0 * eps);
    }
  }

  /// One parseable line per channel: the naive and stationary estimates, the
  /// residual, ||theta||, the documented error bound evaluated, and the
  /// stationarity self-check. No-op unless the executor asks (--hylleraas).
  void report_property_accuracy(const State &s, const std::string &pert,
                                double freq, double c_conv = -2.0) {
    for (int r = 0; r < target_.n_responses(); ++r) {
      const double naive = (r < static_cast<int>(s.last_property.size()))
                               ? s.last_property[r] : 0.0;
      const double res   = (r < static_cast<int>(s.last_bsh_residual.size()))
                               ? s.last_bsh_residual[r] : 0.0;
      const double tn    = (r < static_cast<int>(s.last_theta_norm.size()))
                               ? s.last_theta_norm[r] : 0.0;
      const double S     = stationary_property(s, r);      // collective
      // two probe sizes: a correct quadratic form gives 4 at BOTH, and the
      // agreement between them is what says the probe is above the noise floor
      const double ratio  = stationarity_ratio(s, r, 1e-1);   // collective
      const double ratio2 = stationarity_ratio(s, r, 5e-2);   // collective
      const double alin   = linear_coefficient(s, r);         // collective
      plog("[PROPACC] pert=", pert, " freq=", freq,
           " thresh=", madness::FunctionDefaults<3>::get_thresh(),
           " k=", madness::FunctionDefaults<3>::get_k(),
           " ch=", r, " iters=", s.iter,
           " res=", res, " theta_norm=", tn,
           " bound=", std::abs(c_conv) * tn * res,
           " naive=", naive, " S=", S, " diff=", S - naive,
           " stat_ratio=", ratio, " stat_ratio_half=", ratio2,
           " lin_coeff=", alin);
    }
  }

private:
  // F2d (doc 32 §5.6): gated print. Centralizes the Normal/rank-0 guard and the
  // per-subworld tag. EMPTY prefix ⇒ print(...) verbatim (G=0 byte-identical);
  // non-empty ⇒ the "[g{gid}/{G} {host}] " tag leads each line.
  template <class... Ts> void plog(Ts &&...a) const {
    if (print_level_ < PrintLevel::Normal || world_.rank() != 0) return;
    if (log_prefix_.empty()) print(std::forward<Ts>(a)...);
    else                     print(log_prefix_, std::forward<Ts>(a)...);
  }

  void print_header() const {
    plog("");
    plog("FDSolver<", type_name(), ", ClosedShell>  n_responses =",
          target_.n_responses(),
          " thresh =", madness::FunctionDefaults<3>::get_thresh(),
          " c_xc =", target_.gs.c_xc);
    plog("  policy: dconv_user =", policy_.dconv_user,
          " bsh_target =", targets_.bsh_residual,
          " density_target =", targets_.density_residual);
    // D-096: the gate named once; each row then says how far each quantity
    // is from it (its value over the gate).
    plog("  gate: ||dx|| <= bsh_target + rtol ||x||;"
         " ||dρ|| <= density_target + rtol ||ρ1||; rtol =",
         targets_.rtol,
         "; iteration 1 has no previous density, so its dρ is unavailable"
         " (first sample)");
    plog("  responses:");
    for (int c = 0; c < target_.n_responses(); ++c) {
      plog("    [", c, "]  omega =", target_.responses[c].omega);
    }
    plog("");
  }

  static constexpr const char* type_name() {
    if constexpr (std::is_same_v<Type, Static>) return "Static";
    if constexpr (std::is_same_v<Type, Full>)   return "Full";
    if constexpr (std::is_same_v<Type, TDA>)    return "TDA";
    return "Unknown";
  }

  // D-096: one row per iteration at Normal, each gated quantity with its
  // distance from the gate (value over gate, the tracks of gate_history);
  // iteration 1 has no previous density, so its dρ is "unavailable (first
  // sample)", never zero (the record says the same). Channel detail at
  // Verbose.
  void print_iter_banner(const State &out) const {
    if (print_level_ < PrintLevel::Normal || world_.rank() != 0) return;
    double max_res  = 0.0;
    for (double r : out.last_bsh_residual) max_res = std::max(max_res, r);
    double max_drho = 0.0;
    for (double r : out.last_density_residual) max_drho = std::max(max_drho, r);
    const bool first_sample = out.iter <= 1;
    auto track = [&](std::size_t t) {
      return out.gate_history.size() > t && !out.gate_history[t].empty()
                 ? out.gate_history[t].back()
                 : 0.0;
    };
    if (first_sample)
      plog("iter", out.iter, "  bsh =", max_res, " (", track(0),
           "of gate)   dρ = unavailable (first sample)");
    else
      plog("iter", out.iter, "  bsh =", max_res, " (", track(0),
           "of gate)   dρ =", max_drho, " (", track(1), "of gate)");
    if (print_level_ >= PrintLevel::Verbose) {
      for (size_t c = 0; c < out.last_bsh_residual.size(); ++c) {
        double dr = (c < out.last_density_residual.size())
                        ? out.last_density_residual[c] : 0.0;
        double pr = (c < out.last_property.size()) ? out.last_property[c] : 0.0;
        if (first_sample)
          plog("  ch", c, "  omega =", target_.responses[c].omega,
               "  res =", out.last_bsh_residual[c],
               "  dρ = unavailable (first sample)  <v|x> =", pr);
        else
          plog("  ch", c, "  omega =", target_.responses[c].omega,
               "  res =", out.last_bsh_residual[c], "  dρ =", dr,
               "  <v|x> =", pr);
      }
    }
  }

  /// Source overlap sum over blocks of sum_i <source_i | x_i>: the channel's
  /// diagonal property estimate (alpha_vv up to the convention factor). Cheap
  /// (n_occ inner products per block) and collective.
  static double source_overlap(madness::World &world, const Storage &x,
                               const Channel &ch) {
    double s = 0.0;
    const auto xb = x.blocks();
    const auto vb = ch.source.blocks();
    for (std::size_t b = 0; b < xb.size(); ++b) {
      const auto &xv = *xb[b];
      const auto &vv = *vb[b];
      const std::size_t n = std::min(xv.size(), vv.size());
      if (n == 0) continue;
      madness::Tensor<double> ov = madness::inner(world,
          std::vector<madness::real_function_3d>(xv.begin(), xv.begin() + n),
          std::vector<madness::real_function_3d>(vv.begin(), vv.begin() + n));
      s += ov.sum();
    }
    return s;
  }

  /// theta = V0z - E0z + gamma(z) + V_p for an ARBITRARY z, mirroring the
  /// assembly step() performs inline. The stationary estimate must see the same
  /// operator the iteration solves.
  Storage assemble_theta_for(const Storage &z, int r,
                             const madness::real_function_3d &rho) {
    const double thr = madness::FunctionDefaults<3>::get_thresh();
    Storage theta = K::compute_V0x(world_, target_.gs, z);
    {
      auto E0x = K::compute_E0x(world_, target_.gs, z);
      theta.axpy(world_, -1.0, E0x);
    }
    {
      auto gamma = K::compute_gamma(world_, target_.gs, z, rho);
      theta.axpy(world_, +1.0, gamma);
    }
    add_perturbation_source(world_, theta, target_.responses[r]);
    theta.truncate_all(world_, thr);
    return theta;
  }

  /// Inner product of two Storages, block by block. source_overlap is this
  /// against the channel's perturbation; the stationary property estimate and
  /// the error-bound constant need the general form.
  static double storage_inner(madness::World &world, const Storage &a,
                              const Storage &b) {
    double s = 0.0;
    const auto ab = a.blocks();
    const auto bb = b.blocks();
    for (std::size_t k = 0; k < ab.size() && k < bb.size(); ++k) {
      const auto &av = *ab[k];
      const auto &bv = *bb[k];
      const std::size_t n = std::min(av.size(), bv.size());
      if (n == 0) continue;
      madness::Tensor<double> ov = madness::inner(world,
          std::vector<madness::real_function_3d>(av.begin(), av.begin() + n),
          std::vector<madness::real_function_3d>(bv.begin(), bv.begin() + n));
      s += ov.sum();
    }
    return s;
  }

  static double storage_norm(madness::World &world, const Storage &a) {
    return std::sqrt(std::max(0.0, storage_inner(world, a, a)));
  }

  /// Plateau bookkeeping: append this iteration's normalised distance for EACH
  /// gated quantity (bsh, density) to its own track (reset on a thresh change)
  /// and set out.stalled per policy. Tracking them separately is what keeps a
  /// slow-moving density from masking a residual that is still falling.
  void update_stall(State &out, const State &in) const {
    const double thr = madness::FunctionDefaults<3>::get_thresh();
    out.gate_history = (in.gate_thresh == thr)
                           ? in.gate_history
                           : std::vector<std::vector<double>>(2);
    out.gate_thresh  = thr;
    if (out.gate_history.size() < 2) out.gate_history.resize(2);
    double g_bsh = 0.0, g_rho = 0.0;
    for (std::size_t c = 0; c < out.last_bsh_residual.size(); ++c) {
      g_bsh = std::max(g_bsh, out.last_bsh_residual[c] /
                                  fd_gate(targets_.bsh_residual, targets_.rtol, x_norm(out, c)));
      if (c < out.last_density_residual.size())
        g_rho = std::max(g_rho, out.last_density_residual[c] /
                                    fd_gate(targets_.density_residual, targets_.rtol, rho_norm(out, c)));
    }
    // iter 1 has no previous density, so drho is 0 and would read as "target
    // met"; mark it as no measurement so the track cannot be compared against it.
    out.gate_history[0].push_back(g_bsh);
    out.gate_history[1].push_back(out.iter > 1 ? g_rho
                                   : std::numeric_limits<double>::infinity());
    // The engine's plateau detector belongs to the engine check (B8): it
    // applies when the engine decides. When the written targets decide, only
    // the configured progress rule (progress.*) can declare a stall.
    out.stalled = policy_.decider == madness::convergence::Decider::engine_check
                  && !out.diverged && policy_.plateau(out.gate_history);
    if (out.stalled && print_level_ >= PrintLevel::Normal && world_.rank() == 0) {
      const std::size_t n = out.gate_history[0].size();
      const std::size_t w = static_cast<std::size_t>(policy_.stall_window);
      // D-133: "not improving" covers a quantity that rose; say which.
      plog("[STALL] iter", out.iter, ": no gated quantity improved by",
           100.0 * policy_.stall_ratio, "% over", policy_.stall_window,
           "iters — bsh",
           stall_change(out.gate_history[0][n - 1],
                        out.gate_history[0][n - 1 - w])
               + ", density",
           stall_change(out.gate_history[1][n - 1],
                        out.gate_history[1][n - 1 - w]),
           "(normalised to target); stopping.");
    }
  }

  /// Add the perturbation source to theta in place. Differs across
  /// Storage shapes — Storage = ResponseStateX adds to x_alpha only;
  /// ResponseStateXY adds to both x_alpha and y_alpha. For symmetric
  /// (real, time-independent) dipole perturbations the X and Y source
  /// are equal, so caller fills source.x_alpha == source.y_alpha for
  /// Full; for Static there's only one component anyway.
  static void add_perturbation_source(madness::World &world,
                                       Storage &theta,
                                       const Channel &ch) {
    // X-channel α
    madness::gaxpy(world, 1.0, theta.x_alpha, 1.0, ch.source.x_alpha);
    madness::truncate(world, theta.x_alpha);
    // Y-channel α (Full only)
    if constexpr (std::is_same_v<Type, Full>) {
      madness::gaxpy(world, 1.0, theta.y_alpha, 1.0, ch.source.y_alpha);
      madness::truncate(world, theta.y_alpha);
    }
    // β-spin (OpenShell only): X-channel β, and Y-channel β for Full.
    if constexpr (std::is_same_v<Shell, OpenShell>) {
      madness::gaxpy(world, 1.0, theta.x_beta, 1.0, ch.source.x_beta);
      madness::truncate(world, theta.x_beta);
      if constexpr (std::is_same_v<Type, Full>) {
        madness::gaxpy(world, 1.0, theta.y_beta, 1.0, ch.source.y_beta);
        madness::truncate(world, theta.y_beta);
      }
    }
  }

  // Shared by iteration and saved-state diagnostics. This is the raw map:
  // it never touches acceleration history, iteration counters or stop flags.
  Storage raw_fixed_point(const Storage& response, const madness::real_function_3d& rho,
                          int r, double& theta_norm) {
    const double thr = madness::FunctionDefaults<3>::get_thresh();
    // Inc-2: build the φ-only g0 exchange tensor ONCE per protocol (cached on the
    // ground state) so assemble_theta_tensor doesn't rebuild it every iteration.
    // Only for the ClosedShell Static/Full tensor path (matches the θ gate below).
    if constexpr (std::is_same_v<typename K::State, ResponseStateX<ClosedShell>> ||
                  std::is_same_v<typename K::State, ResponseStateXY<ClosedShell>>) {
      if (policy_.exchange_tensor && target_.gs.g0_alpha.empty())
        target_.gs.g0_alpha =
            exch::build_g0(world_, target_.gs, thr * 0.1, policy_.exchange_tile);
    }

    // θ = V0x − E0x + γ. Gate 1 (policy_.exchange_tensor, --fd-tensor): the
    // fused tensor-exchange assembly (exch::, ClosedShell Static/Full only);
    // gate 0 (default) = the per-op reference path, streamed via Storage::axpy,
    // UNTOUCHED. doc 28 §4 Inc 1; A/B-to-thresh between the two gates.
    typename K::State theta;
    bool tensor_done = false;
    if constexpr (std::is_same_v<typename K::State, ResponseStateX<ClosedShell>> ||
                  std::is_same_v<typename K::State, ResponseStateXY<ClosedShell>>) {
      if (policy_.exchange_tensor) {
        theta = exch::assemble_theta_tensor(world_, target_.gs,
                                            response, rho,
                                            policy_.exchange_tile);
        tensor_done = true;
      }
    }
    if (!tensor_done) {
      theta = K::compute_V0x(world_, target_.gs, response);
      {
        auto E0x = K::compute_E0x(world_, target_.gs, response);
        theta.axpy(world_, -1.0, E0x);                        // theta -= E0x
      }
      {
        auto gamma = K::compute_gamma(world_, target_.gs,
                                      response, rho);
        theta.axpy(world_, +1.0, gamma);                      // theta += γ
      }
    }
    add_perturbation_source(world_, theta, target_.responses[r]); // theta += V_p
    theta.truncate_all(world_, thr);
    // The constant in the error bound docs/guides/convergence.md states,
    // |d alpha| <= |c| ||theta|| ||r_B||. Until now nothing evaluated it, so
    // the bound was documented and uncheckable; theta is already in hand here.
    theta_norm = storage_norm(world_, theta);

    // BSH apply
    return K::bsh_apply(world_, target_.gs, response,
                              theta, target_.responses[r].omega);
  }

public:
  // D-096: the solve's end is said once, by the executor's identified line
  // (the state, the rung, the exit, the iterations); this block is detail,
  // at Verbose.
  void print_final(const State &s, bool converged) const {
    if (print_level_ < PrintLevel::Verbose || world_.rank() != 0)
      return;
    plog("");
    // Diverged-first: converged() returns true on diverged so iterate<>
    // would otherwise mislabel the exit as "Converged".
    if (s.diverged)       plog("Stopped at iter", s.iter,
                                "(diverged — residual exceeded "
                                "explosion guard).");
    else if (s.stalled && !policy_.explicit_progress_configuration.empty())
      plog("Stopped at iter", s.iter, "(configured progress stalled; requested accuracy unmet).");
    else if (s.stalled)   plog("Stopped at iter", s.iter,
                                "(stalled — residual plateau above the "
                                "targets for", policy_.stall_window,
                                "iters; not converged).");
    else if (converged)   plog("Converged in", s.iter, "iters.");
    else                  plog("Stopped at iter", s.iter,
                                "(max iters reached, not converged).");
    if (!s.last_bsh_residual.empty()) {
      plog("  residuals   =");
      for (size_t c = 0; c < s.last_bsh_residual.size(); ++c) {
        double dr = (c < s.last_density_residual.size())
                        ? s.last_density_residual[c] : 0.0;
        double pr = (c < s.last_property.size()) ? s.last_property[c] : 0.0;
        plog("    ch", c, "  omega =", target_.responses[c].omega,
              "  bsh =", s.last_bsh_residual[c], "  dρ =", dr,
              "  <v|x> =", pr);
      }
    }
    plog("");
  }

  // Shared measurements describe the current functions. Only advance() can
  // consume this evaluation's proposal and touch KAIN/comparison history.
  void reset_shared() {
    kain_.reset();
    shared_proposal_.clear(); shared_density_.clear(); shared_previous_density_.clear();
    shared_kinetic_.clear(); shared_previous_kinetic_.clear();
    shared_previous_info_.reset(); shared_state_ = nullptr;
  }

  madness::iteration::IterationMetrics evaluate(
      State& state, const madness::iteration::IterationInfo& info) {
    namespace it = madness::iteration;
    MADNESS_CHECK_THROW(state.responses.size() == target_.responses.size() &&
        !state.responses.empty() && info.iteration == static_cast<std::size_t>(state.iter),
        "FD evaluation state does not match its iteration context");
    shared_state_ = nullptr;
    shared_proposal_.clear(); shared_density_.clear(); shared_kinetic_.clear();
    state.last_bsh_residual.clear(); state.last_density_residual.clear();
    state.last_x_norm.clear(); state.last_rho_norm.clear();
    state.diverged = false;
    it::IterationMetrics result{info, {}};
    for (std::size_t channel = 0; channel < state.responses.size(); ++channel) {
      const auto& response = state.responses[channel];
      auto rho = K::compute_density(world_, target_.gs, response);
      double theta_norm = 0.0;
      auto proposal = raw_fixed_point(response, rho, static_cast<int>(channel), theta_norm);
      const auto residuals = measure_function_residuals(world_, response, proposal);
      std::vector<it::ComponentEvidence> components;
      // flatten() has the storage's canonical block order. Build each block's
      // norms explicitly so map ordering cannot swap X/Y/spin scales.
      const auto add = [&](const char* block, const auto& functions) {
        const auto& values = residuals.at(block);
        for (std::size_t i = 0; i < functions.size(); ++i)
          components.push_back({std::string(block) + ":" + std::to_string(i),
                                values.at(i), functions[i].norm2()});
      };
      add("x_alpha", response.x_alpha);
      if constexpr (rs::detail::has_y_alpha<Storage>::value) add("y_alpha", response.y_alpha);
      if constexpr (rs::detail::has_x_beta<Storage>::value) add("x_beta", response.x_beta);
      if constexpr (rs::detail::has_y_beta<Storage>::value) add("y_beta", response.y_beta);
      it::StateMetrics values;
      values.id = state.responses.size() == 1 ? "response" : "response_" + std::to_string(channel);
      if (shared_previous_info_) {
        const auto& old = *shared_previous_info_;
        values.comparison = it::Comparison{old.task_id, old.protocol, old.iteration, values.id};
      }
      auto residual = it::maximum_residual(std::move(components));
      state.last_bsh_residual.push_back(residual.value.value_or(0.0));
      state.last_x_norm.push_back(storage_norm(world_, response));
      const double density_norm = rho.norm2();
      it::FieldNorms field{density_norm, {}, {}};
      if (values.comparison) {
        field.previous = shared_previous_density_.at(channel).norm2();
        field.change = (rho - shared_previous_density_.at(channel)).norm2();
      }
      auto density = it::maximum_field_change({field}, values.comparison);
      double kinetic = 0.0;
      for (const auto& [block, functions] : measure_function_kinetics(world_, response))
        for (double value : functions) kinetic += value;
      const std::optional<double> previous = values.comparison ?
          std::optional<double>(shared_previous_kinetic_.at(channel)) : std::nullopt;
      values.metrics = {{it::Metric::bsh_residual, residual},
          {it::Metric::density_change, density},
          {it::Metric::energy_change, it::MetricValue::not_applicable("driven_frequency_is_input")},
          {it::Metric::kinetic_change, it::scalar_change(kinetic, previous, values.comparison)}};
      state.diverged = state.diverged || residual.status == it::Availability::invalid ||
          !std::isfinite(kinetic) || fd_exploded(residual.value.value_or(0.0),
              std::max(state.last_x_norm.back(), storage_norm(world_, proposal)), policy_.explosion_guard);
      state.last_density_residual.push_back(density.value.value_or(0.0));
      state.last_rho_norm.push_back(density_norm);
      shared_proposal_.push_back(std::move(proposal));
      shared_density_.push_back(std::move(rho));
      shared_kinetic_.push_back(kinetic);
      result.states.push_back(std::move(values));
    }
    shared_info_ = info; shared_state_ = &state;
    return result;
  }

  void advance(State& state) {
    MADNESS_CHECK_THROW(shared_state_ == &state && !state.diverged &&
        shared_info_.iteration == static_cast<std::size_t>(state.iter),
        "FD update requires its current evaluated proposal");
    const int diagnostics = print_level_ >= PrintLevel::Verbose ? 1 : 0;
    kain_.apply(state.responses, shared_proposal_, diagnostics);
    state.responses = std::move(shared_proposal_);
    shared_previous_density_ = std::move(shared_density_);
    shared_previous_kinetic_ = std::move(shared_kinetic_);
    shared_previous_info_ = shared_info_;
    ++state.iter;
    shared_state_ = nullptr;
  }

  /// One outer iteration. Streamed-theta layout: for each response
  /// in turn, build theta in place (V0x → += gamma → -= E0x → += source)
  /// then BSH-apply, freeing the intermediates before moving to the next
  /// response. Peak memory ≈ `in + out + 1 theta + 1 temp + 1 new_x`
  /// instead of the materialize-V0x/E0x/gamma-for-all-responses-up-front
  /// layout's `~5M Storage` peak.
  State step(State in) {
    State out;
    out.iter      = in.iter + 1;
    out.responses = in.responses;
    out.progress_criteria_history = std::move(in.progress_criteria_history);
    out.configured_progress_history = std::move(in.configured_progress_history);
    if (record_function_metrics_) {
      out.kinetic_history = std::move(in.kinetic_history);
      out.bsh_history = std::move(in.bsh_history);
      if (out.kinetic_history.empty()) {
        FunctionKineticSample initial{in.iter, {}};
        for (const auto& response : in.responses)
          initial.channels.push_back(measure_function_kinetics(world_, response));
        out.kinetic_history.push_back(std::move(initial));
      }
    }
    const int M   = target_.n_responses();
    const double thr = madness::FunctionDefaults<3>::get_thresh();

    out.last_density_residual.assign(M, 0.0);
    out.last_bsh_residual.assign(M, 0.0);
    out.last_property.assign(M, 0.0);
    out.last_theta_norm.assign(M, 0.0);
    out.last_x_norm.assign(M, 0.0);
    out.last_rho_norm.assign(M, 0.0);
    out.rho_alpha_prev.resize(M);

    for (int r = 0; r < M; ++r) {
      // ρ (one density function per response; kept for next iter's Δρ check)
      auto rho = K::compute_density(world_, target_.gs, in.responses[r]);
      out.last_rho_norm[r] = rho.norm2();
      if (!in.rho_alpha_prev.empty()) {
        auto drho = rho - in.rho_alpha_prev[r];
        out.last_density_residual[r] = drho.norm2();
      }

      auto x_new = raw_fixed_point(in.responses[r], rho, r, out.last_theta_norm[r]);
      out.last_bsh_residual[r] =
          K::compute_residual_norm(world_, in.responses[r], x_new);

      // Capture the raw fixed-point difference before KAIN changes the output.
      if (record_function_metrics_)
        out.function_residuals.push_back(measure_function_residuals(world_, in.responses[r], x_new));
      out.responses[r]      = std::move(x_new);
      out.rho_alpha_prev[r] = std::move(rho);
    }

    // ---- KAIN + step restriction (over flat in/out responses) ------------
    // diag_level pulled from print_level_; gated uniformly across ranks
    // inside ResponseSubspaceKain (collective work is unconditional).
    const int kain_diag =
        (print_level_ >= PrintLevel::Verbose) ? 1 : 0;
    kain_.apply(in.responses, out.responses, kain_diag);

    // ||x|| after KAIN (the iterate that is saved), then the explosion guard,
    // relative to it.
    for (int r = 0; r < M; ++r)
      out.last_x_norm[r] = storage_norm(world_, out.responses[r]);
    if (record_function_metrics_)
      out.bsh_history.push_back({in.iter, out.function_residuals, out.last_bsh_residual,
          out.last_x_norm, targets_.bsh_residual, targets_.rtol,
          out.last_density_residual, out.last_rho_norm, targets_.density_residual,
          !in.rho_alpha_prev.empty()});
    for (int r = 0; r < M; ++r) {
      if (fd_exploded(out.last_bsh_residual[r], out.last_x_norm[r],
                      policy_.explosion_guard)) { out.diverged = true; break; }
    }

    // Property trace (after KAIN: this is the iterate that will be saved).
    for (int r = 0; r < M; ++r)
      out.last_property[r] = source_overlap(world_, out.responses[r],
                                            target_.responses[r]);
    if (record_function_metrics_) {
      FunctionKineticSample sample{out.iter, {}};
      for (const auto& response : out.responses)
        sample.channels.push_back(measure_function_kinetics(world_, response));
      out.kinetic_history.push_back(std::move(sample));
    }
    update_stall(out, in);
    if (!policy_.explicit_accuracy_targets.empty()) {
      out.explicit_accuracy = evaluate_requested_accuracy(out);
      for (const auto& criterion : out.explicit_accuracy)
        MADNESS_CHECK_THROW(criterion.at("status") != "invalid", "invalid explicit FD accuracy evidence");
      if (!policy_.explicit_progress_configuration.empty()) {
        out.progress_criteria_history.push_back(out.explicit_accuracy.at(0).at("criteria"));
        out.configured_progress = fd_progress_decision(out.progress_criteria_history,
            policy_.explicit_progress_configuration, madness::FunctionDefaults<3>::get_thresh(),
            out.diverged, policy_.diagnostic_iterations_per_rung,
            /*fixed_intermediate_work=*/true);
        MADNESS_CHECK_THROW(out.diverged || out.configured_progress.at("action") != "stop_invalid",
            "invalid configured FD progress evidence");
        out.configured_progress_history.push_back(out.configured_progress);
        if (out.configured_progress.at("action_applied").template get<bool>()) out.stalled = true;
      }
      const bool input_satisfied = std::all_of(out.explicit_accuracy.begin(), out.explicit_accuracy.end(),
          [](const auto& criterion) { return criterion.at("status") == "satisfied"; });
      // The saved-state check is part of the targets' decision (criteria mode);
      // the engine check doesn't need it.
      if (policy_.decider == madness::convergence::Decider::requested_criteria
          && !out.diverged && !out.stalled && input_satisfied) {
        const auto saved = evaluate_saved_state(out);
        out.saved_explicit_accuracy = evaluate_requested_accuracy(out, &saved);
        for (const auto& criterion : out.saved_explicit_accuracy)
          MADNESS_CHECK_THROW(criterion.at("status") != "invalid", "invalid saved FD accuracy evidence");
      }
    }

    print_iter_banner(out);
    append_convergence_log(out);
    return out;
  }

  /// Evaluate the serialized output without taking an extra solver step.
  /// rho_alpha_prev is the last genuine step input, including a one-step run.
  State evaluate_saved_state(const State& saved) {
    State evidence;
    evidence.iter = saved.iter;
    for (int r = 0; r < target_.n_responses(); ++r) {
      auto response = saved.responses[r].copy(world_);
      if (record_function_metrics_)
        evidence.function_kinetics.push_back(measure_function_kinetics(world_, response));
      auto rho = K::compute_density(world_, target_.gs, response);
      evidence.last_rho_norm.push_back(rho.norm2());
      if (!saved.rho_alpha_prev.empty())
        evidence.last_density_residual.push_back((rho - saved.rho_alpha_prev[r]).norm2());
      double theta_norm = 0.0;
      auto next = raw_fixed_point(response, rho, r, theta_norm);
      evidence.last_bsh_residual.push_back(K::compute_residual_norm(world_, response, next));
      if (record_function_metrics_)
        evidence.function_residuals.push_back(measure_function_residuals(world_, response, next));
      evidence.last_x_norm.push_back(storage_norm(world_, response));
    }
    return evidence;
  }

  /// Apply only requested criteria to the genuine input or reassessed output.
  nlohmann::json evaluate_requested_accuracy(const State& state, const State* validation = nullptr) const {
    using json = nlohmann::json;
    json results = json::array();
    const auto& measured = validation ? *validation : state;
    for (int channel = 0; channel < target_.n_responses(); ++channel) {
      const auto index = static_cast<std::size_t>(channel);
      const auto scalar = [&](const std::vector<double>& values, bool available, const char* units) {
        const bool present = index < values.size();
        const bool valid = present && std::isfinite(values[index]) && values[index] >= 0.0;
        return json{{"value", available && valid ? json(values[index]) : json(nullptr)},
                    {"availability", present && !valid ? "invalid" : available && valid ? "measured" : "unavailable"},
                    {"units", units}};
      };
      json metrics;
      metrics["bsh_residual"] = scalar(measured.last_bsh_residual, true, "response_amplitude_L2_atomic_units");
      if (index < measured.function_residuals.size())
        metrics["bsh_residual"]["per_function"] = madness::function_residual_record(measured.function_residuals[index]);
      metrics["density_change"] = scalar(measured.last_density_residual, validation || state.iter > 1,
          "kernel_induced_density_L2_atomic_units");
      metrics["density_change"]["reduction"] = "kernel_induced_density_L2_norm";
      const double threshold = madness::FunctionDefaults<3>::get_thresh();
      const int k = madness::FunctionDefaults<3>::get_k();
      const auto history_record = [&](std::size_t sample) {
        const auto& item = state.kinetic_history.at(sample);
        return function_kinetic_record(item.channels.at(index), item.iteration, "", threshold, k, false,
            item.iteration == 0 ? "initial_solver_iterate" : "solver_iterate_after_acceleration");
      };
      const auto count = state.kinetic_history.size();
      if (validation && index < measured.function_kinetics.size()) {
        const auto current = function_kinetic_record(measured.function_kinetics[index], state.iter, "", threshold, k,
            true, "saved_response_before_archive");
        metrics["delta_kinetic_energy"] = kinetic_change_record(current, count >= 2 ? history_record(count - 2) : json(nullptr));
      } else if (!validation && count >= 2) {
        metrics["delta_kinetic_energy"] = kinetic_change_record(history_record(count - 2),
            count >= 3 ? history_record(count - 3) : json(nullptr));
      }
      results.push_back(madness::evaluate_explicit_criteria(fd_accuracy_measurements(metrics), policy_.explicit_accuracy_targets));
    }
    return results;
  }

  /// Every leg within its relative gate (review C1): ||dx|| <= atol + rtol ||x||
  /// and ||drho|| <= atol + rtol ||rho1|| (fd_leg_within). No iteration rule and
  /// no diverged/stalled check. This residual-only diagnostic is not a strict
  /// convergence verdict; callers must also respect the solver's eligibility.
  bool within_targets(const State &s) const {
    if (s.last_bsh_residual.empty()) return false;
    for (std::size_t c = 0; c < s.last_bsh_residual.size(); ++c) {
      const double drho = c < s.last_density_residual.size() ? s.last_density_residual[c] : 0.0;
      if (!fd_leg_within(s.last_bsh_residual[c], x_norm(s, c), drho, rho_norm(s, c), targets_))
        return false;
    }
    return true;
  }

  /// Converged: diverged or stalled (the loop exits; the executor reports
  /// them), or every leg converged by fd_leg_converged: at least two
  /// iterations, so the density change has been measured, and both relative
  /// gates hold.
  bool legacy_converged(const State &s) const {
    if (s.diverged) return true;
    if (s.stalled)  return true;   // exit; executor reads State::stalled
    return engine_check_passed(s);
  }

  /// The engine check alone (B8): every leg converged by fd_leg_converged.
  /// Divergence and stalls are not verdicts on convergence; they are the
  /// caller's validity check and rung_exit's stall exit.
  [[nodiscard]] bool engine_check_passed(const State& s) const {
    if (s.iter < policy_.min_iters_before_conv) return false;
    if (s.last_bsh_residual.empty()) return false;
    for (std::size_t c = 0; c < s.last_bsh_residual.size(); ++c) {
      const double drho = c < s.last_density_residual.size() ? s.last_density_residual[c] : 0.0;
      if (!fd_leg_converged(s.iter, s.last_bsh_residual[c], x_norm(s, c), drho,
                            rho_norm(s, c), targets_))
        return false;
    }
    return true;
  }

  bool converged(const State& state) const {
    if (state.diverged || state.stalled) return true; // loop exit, never eligibility
    return madness::convergence::converged(policy_.decider,
                                           legacy_converged(state),
                                           requested_criteria_met(state));
  }

  /// B5: how this solve's rung ends after `iterations` iterations of this
  /// attempt, by the shared rung_exit; this solver only supplies its
  /// verdicts. iterate<> checks divergence before asking.
  [[nodiscard]] madness::convergence::StopReason
  stop_reason(const State& state, int iterations, int cap,
              int forced_iterations, bool final_protocol = true) const {
    return madness::convergence::rung_exit(
        {policy_.decider, engine_check_passed(state),
         requested_criteria_met(state), state.stalled, iterations, cap,
         forced_iterations, final_protocol});
  }

  /// Every written target met on the solver's input and on the saved state
  /// (B8; the latter is evaluated in criteria mode only).
  [[nodiscard]] bool requested_criteria_met(const State& state) const {
    if (policy_.explicit_accuracy_targets.empty()
        || state.explicit_accuracy.size() != state.responses.size()
        || state.saved_explicit_accuracy.size() != state.responses.size())
      return false;
    for (const auto* evidence : {&state.explicit_accuracy, &state.saved_explicit_accuracy})
      for (const auto& criterion : *evidence)
        if (criterion.at("status") != "satisfied") return false;
    return true;
  }

  /// ||x|| / ||rho1|| of leg c, or 0 when not recorded (a state loaded from
  /// disk): the gate then falls back to the absolute target.
  static double x_norm(const State &s, std::size_t c) {
    return c < s.last_x_norm.size() ? s.last_x_norm[c] : 0.0;
  }
  static double rho_norm(const State &s, std::size_t c) {
    return c < s.last_rho_norm.size() ? s.last_rho_norm[c] : 0.0;
  }

  void save(const State &s, const std::string &path_prefix) const {
    for (int c = 0; c < target_.n_responses(); ++c) {
      s.responses[c].save(world_, path_prefix + ".channel_"
                                   + std::to_string(c));
    }
  }

private:
  /// Append one row per channel at the end of step(). Schema:
  ///   iter,protocol_thresh,state,omega,bsh_residual,density_residual,diverged,
  ///   property,gate,stalled
  /// `state` here is the channel index; `omega` is the fixed channel
  /// frequency (not changing iter-to-iter like ES). `property` is the source
  /// overlap <v|x> of the channel, `gate` the normalised gate distance
  /// (plateau detector). Header is written once per file. Only rank 0 writes.
  void append_convergence_log(const State &s) {
    if (log_path_.empty() || world_.rank() != 0) return;
    std::ofstream out(log_path_, std::ios::app);
    if (!out) return;
    if (!log_header_written_) {
      out << "iter,protocol_thresh,state,omega,bsh_residual,"
          << "density_residual,diverged,property,gate,stalled\n";
      log_header_written_ = true;
    }
    const double pthr = madness::FunctionDefaults<3>::get_thresh();
    const int M = target_.n_responses();
    out.precision(12);
    for (int c = 0; c < M; ++c) {
      const double bsh  = (static_cast<size_t>(c) < s.last_bsh_residual.size())
                              ? s.last_bsh_residual[c] : 0.0;
      const double drho = (static_cast<size_t>(c) < s.last_density_residual.size())
                              ? s.last_density_residual[c] : 0.0;
      const double prop = (static_cast<size_t>(c) < s.last_property.size())
                              ? s.last_property[c] : 0.0;
      const double gate = (s.gate_history.empty() || s.gate_history[0].empty())
                              ? 0.0 : s.gate_history[0].back();
      out << s.iter << ',' << pthr << ',' << c << ','
          << target_.responses[c].omega << ',' << bsh << ','
          << drho << ',' << (s.diverged ? 1 : 0) << ','
          << prop << ',' << gate << ',' << (s.stalled ? 1 : 0) << '\n';
    }
  }

  State* shared_state_ = nullptr;
  std::vector<Storage> shared_proposal_;
  std::vector<madness::real_function_3d> shared_density_, shared_previous_density_;
  std::vector<double> shared_kinetic_, shared_previous_kinetic_;
  madness::iteration::IterationInfo shared_info_;
  std::optional<madness::iteration::IterationInfo> shared_previous_info_;

  madness::World            &world_;
  FDProblem<Type, Shell>      target_;
  ConvergencePolicy          policy_;
  ConvergencePolicy::Targets targets_{};
  ResponseSubspaceKain<Storage>   kain_;
  PrintLevel                 print_level_ = PrintLevel::Normal;
  std::string                log_prefix_;   // F2d: per-subworld tag ("" = none)
  std::string                log_path_;
  bool                       log_header_written_ = false;
  bool                       record_function_metrics_ = false;
};

// Call-site shape (cf. ESSolver):
//
//   FDProblem<Static, ClosedShell> tgt{...};                // gs + responses
//   FDSolver<Static, ClosedShell> solver(world, tgt,
//                                        ConvergencePolicy{1e-4});
//   typename decltype(solver)::State s0 = ...;             // initial guess
//   auto sf = solvers::iterate_protocol(
//       solver, s0, {1e-4, 1e-6}, prepare_fn,
//       solvers::IterateProtocolPolicy{25});
//   solver.save(sf, "fd_solve");

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_FD_SOLVER_HPP
