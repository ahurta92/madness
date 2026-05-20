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
// could be cached per channel (deferred — built inside K::bsh_apply
// for now, same as ES).
// =========================================================================

#include "../kernels/full.hpp"     // Kernels<Full, ClosedShell>
#include "../kernels/static.hpp"   // Kernels<Static, ClosedShell>
#include "../kernels/tags.hpp"
#include "convergence_policy.hpp"
#include "fd_problem.hpp"
#include "iterate.hpp"
#include "response_state.hpp"
#include "response_subspace_kain.hpp"

#include <madness/mra/mra.h>

#include <algorithm>
#include <fstream>
#include <string>
#include <type_traits>
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
    std::vector<double>                     last_bsh_residual;
    std::vector<madness::real_function_3d>  rho_alpha_prev;
    std::vector<double>                     last_density_residual;
    int                                     iter = 0;
    bool                                    diverged = false;
  };

  FDSolver(madness::World &world, FDProblem<Type, Shell> target,
           ConvergencePolicy policy,
           PrintLevel print_level = PrintLevel::Normal)
      : world_(world), target_(std::move(target)),
        policy_(policy), kain_(world, policy),
        print_level_(print_level) {
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
    kain_.reset();
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

private:
  void print_header() const {
    if (print_level_ < PrintLevel::Normal || world_.rank() != 0) return;
    print("");
    print("FDSolver<", type_name(), ", ClosedShell>  n_responses =",
          target_.n_responses(),
          " thresh =", madness::FunctionDefaults<3>::get_thresh(),
          " c_xc =", target_.gs.c_xc);
    print("  policy: dconv_user =", policy_.dconv_user,
          " bsh_target =", targets_.bsh_residual,
          " density_target =", targets_.density_residual);
    print("  responses:");
    for (int c = 0; c < target_.n_responses(); ++c) {
      print("    [", c, "]  omega =", target_.responses[c].omega);
    }
    print("");
  }

  static constexpr const char* type_name() {
    if constexpr (std::is_same_v<Type, Static>) return "Static";
    if constexpr (std::is_same_v<Type, Full>)   return "Full";
    if constexpr (std::is_same_v<Type, TDA>)    return "TDA";
    return "Unknown";
  }

  void print_iter_banner(const State &out) const {
    if (print_level_ < PrintLevel::Normal || world_.rank() != 0) return;
    double max_res  = 0.0;
    for (double r : out.last_bsh_residual) max_res = std::max(max_res, r);
    double max_drho = 0.0;
    for (double r : out.last_density_residual) max_drho = std::max(max_drho, r);
    print("iter", out.iter,
          "  max_res =", max_res, "  max_dρ =", max_drho);
    if (print_level_ >= PrintLevel::Verbose) {
      for (size_t c = 0; c < out.last_bsh_residual.size(); ++c) {
        double dr = (c < out.last_density_residual.size())
                        ? out.last_density_residual[c] : 0.0;
        print("  ch", c, "  omega =", target_.responses[c].omega,
              "  res =", out.last_bsh_residual[c],
              "  dρ =", dr);
      }
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

public:
  void print_final(const State &s, bool converged) const {
    if (print_level_ < PrintLevel::Normal || world_.rank() != 0) return;
    print("");
    // Diverged-first: converged() returns true on diverged so iterate<>
    // would otherwise mislabel the exit as "Converged".
    if (s.diverged)       print("Stopped at iter", s.iter,
                                "(diverged — residual exceeded "
                                "explosion guard).");
    else if (converged)   print("Converged in", s.iter, "iters.");
    else                  print("Stopped at iter", s.iter,
                                "(max iters reached, not converged).");
    if (!s.last_bsh_residual.empty()) {
      print("  residuals   =");
      for (size_t c = 0; c < s.last_bsh_residual.size(); ++c) {
        double dr = (c < s.last_density_residual.size())
                        ? s.last_density_residual[c] : 0.0;
        print("    ch", c, "  omega =", target_.responses[c].omega,
              "  bsh =", s.last_bsh_residual[c], "  dρ =", dr);
      }
    }
    print("");
  }

  /// One outer iteration. Skips stages 2-4 of the ESSolver shape.
  State step(State in) {
    State out;
    out.iter      = in.iter + 1;
    out.responses  = in.responses;
    const int M   = target_.n_responses();

    // ---- 1. per-channel building blocks ----------------------------------
    std::vector<Storage> V0x(M), E0x(M), gamma(M);
    std::vector<madness::real_function_3d> rho_alpha(M);
    out.last_density_residual.assign(M, 0.0);
    for (int c = 0; c < M; ++c) {
      rho_alpha[c] = K::compute_density(world_, target_.gs, out.responses[c]);
      gamma[c]     = K::compute_gamma(world_, target_.gs, out.responses[c],
                                      rho_alpha[c]);
      V0x[c]       = K::compute_V0x(world_, target_.gs, out.responses[c]);
      E0x[c]       = K::compute_E0x(world_, target_.gs, out.responses[c]);
      if (!in.rho_alpha_prev.empty()) {
        auto drho = rho_alpha[c] - in.rho_alpha_prev[c];
        out.last_density_residual[c] = drho.norm2();
      }
    }
    out.rho_alpha_prev = std::move(rho_alpha);

    // ---- 5. Theta = V0x - E0x + gamma + perturbation_source --------------
    std::vector<Storage> theta(M);
    for (int c = 0; c < M; ++c) {
      theta[c] = K::compute_theta(world_, V0x[c], E0x[c], gamma[c]);
      add_perturbation_source(world_, theta[c], target_.responses[c]);
    }

    // ---- 6. BSH apply + residual -----------------------------------------
    // Each channel's raw BSH update goes into out.responses[c], and the
    // per-channel BSH residual is the raw ||x_old − x_BSH|| (it's
    // recorded BEFORE KAIN damping so convergence tracking measures
    // the dynamics, not the damped step).
    out.last_bsh_residual.assign(M, 0.0);
    for (int c = 0; c < M; ++c) {
      auto x_new = K::bsh_apply(world_, target_.gs, out.responses[c],
                                theta[c], target_.responses[c].omega);
      out.last_bsh_residual[c] =
          K::compute_residual_norm(world_, out.responses[c], x_new);
      out.responses[c] = std::move(x_new);
    }

    // ---- 7. KAIN acceleration + step restriction --------------------------
    // Operates on the flat concat of all responses. Reuses the
    // residual ||x_old − x_BSH|| (computed by KAIN internally on the
    // flat vector). policy_.kain=false disables and falls back to
    // the raw BSH iteration.
    kain_.apply(in.responses, out.responses);

    // Explosion guard.
    for (double r : out.last_bsh_residual) {
      if (r > policy_.explosion_guard) { out.diverged = true; break; }
    }

    print_iter_banner(out);
    append_convergence_log(out);
    return out;
  }

  /// Converged: per-channel BSH residual < bsh_target AND
  ///            per-channel Δρ < density_target  (iter >= 2)
  ///            OR diverged.
  bool converged(const State &s) const {
    if (s.diverged) return true;
    if (s.iter < policy_.min_iters_before_conv) return false;
    if (s.last_bsh_residual.empty()) return false;

    double mx_bsh = 0.0;
    for (double r : s.last_bsh_residual) mx_bsh = std::max(mx_bsh, r);
    if (mx_bsh >= targets_.bsh_residual) return false;

    if (s.iter <= 1) return true;
    double mx_drho = 0.0;
    for (double r : s.last_density_residual) mx_drho = std::max(mx_drho, r);
    return mx_drho < targets_.density_residual;
  }

  void save(const State &s, const std::string &path_prefix) const {
    for (int c = 0; c < target_.n_responses(); ++c) {
      s.responses[c].save(world_, path_prefix + ".channel_"
                                   + std::to_string(c));
    }
  }

private:
  /// Append one row per channel at the end of step(). Schema:
  ///   iter,protocol_thresh,state,omega,bsh_residual,density_residual,diverged
  /// `state` here is the channel index; `omega` is the fixed channel
  /// frequency (not changing iter-to-iter like ES). Header is written
  /// once per file. Only rank 0 writes.
  void append_convergence_log(const State &s) {
    if (log_path_.empty() || world_.rank() != 0) return;
    std::ofstream out(log_path_, std::ios::app);
    if (!out) return;
    if (!log_header_written_) {
      out << "iter,protocol_thresh,state,omega,bsh_residual,"
          << "density_residual,diverged\n";
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
      out << s.iter << ',' << pthr << ',' << c << ','
          << target_.responses[c].omega << ',' << bsh << ','
          << drho << ',' << (s.diverged ? 1 : 0) << '\n';
    }
  }

  madness::World            &world_;
  FDProblem<Type, Shell>      target_;
  ConvergencePolicy          policy_;
  ConvergencePolicy::Targets targets_{};
  ResponseSubspaceKain<Storage>   kain_;
  PrintLevel                 print_level_ = PrintLevel::Normal;
  std::string                log_path_;
  bool                       log_header_written_ = false;
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
