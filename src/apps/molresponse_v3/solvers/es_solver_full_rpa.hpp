#ifndef MOLRESPONSE_V3_SOLVERS_ES_SOLVER_FULL_RPA_HPP
#define MOLRESPONSE_V3_SOLVERS_ES_SOLVER_FULL_RPA_HPP

// =========================================================================
// ESSolverFullRPA<ClosedShell>
//
// Closed-shell excited-state solver via the symmetric Casida reduction:
//
//   [ A  B ] [X]       [ I  0 ] [X]
//   [ B  A ] [Y]  = ω  [ 0 −I ] [Y]
//   (X+Y =: u,  X−Y =: v)   ⇒   (A−B)(A+B) u = ω² u
//
// The iteration variable is `u ∈ ResponseStateX<ClosedShell>`, half the
// storage of the direct (X,Y) solver. The MRA-friendly preconditioner
// for the resulting (M − ω²) Hessian factorizes as a CHAINED BSH:
//
//     1/((ε_a − ε_i)² − ω²)  =  G(+ω) · G(−ω)
//
// where G(±ω) is the standard SCF BSH Green's function. Per Davidson
// step we apply both BSH operators in sequence to the residual.
//
// Algorithm (per outer iter):
//
//   0. Top-of-iter Q + Gram-Schmidt orthonormalize the u-bundle.
//   1. Au_s   = K::apply_AplusB(u_s),   BAu_s = K::apply_AminusB(Au_s).
//   2. Subspace M_sub = ⟨u | BAu⟩,  S = ⟨u | u⟩,  sygv → ω²_s, U.
//      Rotate u and BAu by U (sygv preserves the symmetric Davidson
//      shape on the rotated basis).
//   3. Residual r_s = BAu_s − ω²_s · u_s  (Davidson residual).
//   4. Preconditioned step:
//        δu_s = −2 · BSH(+ω_s) · BSH(−ω_s) · r_s,   Q-projected.
//      Then u_new_s = u_s − δu_s, truncate.
//   5. KAIN accelerate u (history is independent across slots, same
//      shape as the TDA solver's u-bundle history).
//
// Final (X, Y) recovery (call recover_xy() at convergence — NOT part
// of step()):
//   v_s = (1/ω_s) · K::apply_AplusB(u_s),
//   X_s = ½(u_s + v_s),  Y_s = ½(u_s − v_s).
//
// What this solver does NOT do (deliberately):
//   * No FDSolver-style streaming-theta dichotomy. There is no Λ to
//     stream; M·u is built in one apply pair.
//   * No Y storage during iteration. Y materializes only at the end.
// =========================================================================

#include "../kernels/assembly.hpp"
#include "../kernels/common_ops.hpp"            // make_bsh_operators, bsh_shift
#include "../kernels/full.hpp"                  // Kernels<Full, ClosedShell>
#include "../kernels/response_space_ops.hpp"
#include "../kernels/tags.hpp"
#include "convergence_policy.hpp"
#include "iterate.hpp"
#include "response_state.hpp"
#include "response_subspace_kain.hpp"

#include <madness/mra/mra.h>
#include <madness/tensor/tensor.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

namespace molresponse_v3 {

/// Read-only problem definition for ESSolverFullRPA<ClosedShell>.
/// Symmetric with `ESProblem<TDA, *>` — same fields, distinct type to
/// keep ctor calls type-safe.
template <typename Shell>
struct ESProblemFullRPA {
  ResponseGroundState gs;
  int                 n_roots = 0;
};

template <typename Shell>
class ESSolverFullRPA {
  static_assert(std::is_same_v<Shell, ClosedShell>,
                "ESSolverFullRPA is closed-shell only — open-shell RPA "
                "via the symmetric reduction is not implemented.");

public:
  using K       = Kernels<Full, ClosedShell>;
  using Storage = ResponseStateX<ClosedShell>;     // u (and δu)
  using vecfuncT = std::vector<madness::real_function_3d>;

  /// Iteration state. Mirrors `ESSolver::State` but residuals are the
  /// Davidson residual ‖BAu − ω²·u‖ (not BSH-step residuals — there's
  /// no BSH-as-fixed-point step in this solver).
  struct State {
    std::vector<Storage>     roots;          // u_s
    madness::Tensor<double>  omega;          // ω_s = √(max(ω²_s, 0))
    std::vector<double>      last_residual;  // ‖BAu_s − ω²·u_s‖
    std::vector<double>      last_u_change;  // ‖u_new − u_old‖ per slot
    int                      iter     = 0;
    bool                     diverged = false;
  };

  ESSolverFullRPA(madness::World &world,
                  ESProblemFullRPA<Shell> problem,
                  ConvergencePolicy policy,
                  PrintLevel print_level = PrintLevel::Normal)
      : world_(world), gs_(std::move(problem.gs)),
        n_roots_(problem.n_roots),
        policy_(policy), kain_(world, policy),
        print_level_(print_level) {
    refresh_convergence_targets();
    print_header();
  }

  // -- accessors / config -------------------------------------------------
  madness::World&    world() const             { return world_; }
  PrintLevel         print_level() const       { return print_level_; }
  void               set_print_level(PrintLevel p) { print_level_ = p; }
  ConvergencePolicy  convergence_policy() const{ return policy_; }
  void               set_convergence_policy(ConvergencePolicy p) {
    policy_ = p;
    kain_.set_policy(p);
    refresh_convergence_targets();
  }
  void set_gs(ResponseGroundState gs) { gs_ = std::move(gs); kain_.reset(); }
  const ResponseGroundState& gs()      const { return gs_; }
  int                        n_roots() const { return n_roots_; }
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
    print("ESSolverFullRPA<ClosedShell>  num_roots =", n_roots_,
          " thresh =", madness::FunctionDefaults<3>::get_thresh(),
          " c_xc =", gs_.c_xc);
    print("  policy: dconv_user =", policy_.dconv_user,
          " residual_target =", targets_.bsh_residual,
          " explosion_guard =", policy_.explosion_guard);
    print("");
  }

  void print_iter_banner(const State &out) const {
    if (print_level_ < PrintLevel::Normal || world_.rank() != 0) return;
    double mx_r = 0.0, mx_du = 0.0;
    for (double r : out.last_residual) mx_r  = std::max(mx_r,  r);
    for (double r : out.last_u_change) mx_du = std::max(mx_du, r);
    print("iter", out.iter, "  omega =", out.omega,
          "  max_res =", mx_r, "  max_|du| =", mx_du);
    if (print_level_ >= PrintLevel::Verbose) {
      for (size_t s = 0; s < out.last_residual.size(); ++s) {
        const double du =
            (s < out.last_u_change.size()) ? out.last_u_change[s] : 0.0;
        print("  root", s, "  res =", out.last_residual[s],
              "  |du| =", du,
              "  omega =", out.omega(static_cast<long>(s)));
      }
    }
  }

public:
  /// Called by iterate<> after the loop exits.
  void print_final(const State &s, bool converged) const {
    if (print_level_ < PrintLevel::Normal || world_.rank() != 0) return;
    print("");
    if (s.diverged)     print("Stopped at iter", s.iter,
                              "(diverged — residual exceeded explosion guard).");
    else if (converged) print("Converged in", s.iter, "iters.");
    else                print("Stopped at iter", s.iter,
                              "(max iters reached, not converged).");
    print("  omega_final =", s.omega);
    if (!s.last_residual.empty()) {
      print("  residuals   =");
      for (size_t i = 0; i < s.last_residual.size(); ++i) {
        const double du = (i < s.last_u_change.size())
                              ? s.last_u_change[i] : 0.0;
        print("    root", i, ":  res =", s.last_residual[i],
              "  |du| =", du);
      }
    }
    print("");
  }

  /// One outer Davidson iter.
  State step(State in) {
    State out;
    out.iter  = in.iter + 1;
    out.roots = in.roots;
    const int M = n_roots_;

    // ---- 0. top-of-iter Q + GS -------------------------------------------
    rs::project(out.roots, gs_.Qa, gs_.Qb);
    rs::orthonormalize(world_, out.roots);

    // ---- 1. Au_s = (A+B)·u_s,  BAu_s = (A−B)·Au_s ------------------------
    std::vector<Storage> Au(M), BAu(M);
    for (int s = 0; s < M; ++s) {
      Au[s]  = K::apply_AplusB (world_, gs_, out.roots[s]);
      BAu[s] = K::apply_AminusB(world_, gs_, Au[s]);
    }

    // ---- 2. Subspace M_sub = ⟨u | BAu⟩,  S = ⟨u | u⟩ ---------------------
    //
    // BAu is the action of M = (A−B)(A+B) on u. Davidson subspace
    // matrix is then ⟨u_i | M | u_j⟩ = ⟨u_i | BAu_j⟩. The metric S is
    // positive-definite because the u_s are real vecfuncs in a flat
    // concat. sygv (via rs::diagonalize) returns ω² and U with omega
    // ascending and slot-identity preserved (so KAIN history stays
    // aligned across rotation).
    auto Msub  = rs::inner(out.roots, BAu);
    auto Smat  = rs::inner(out.roots, out.roots);
    auto diag  = rs::diagonalize(Msub, Smat,
                                 /*thresh_degenerate=*/-1.0,
                                 policy_.cluster_unmix_factor);
    auto &om2  = diag.omega;   // ω² (eigenvalues of M; should be ≥ 0)
    auto &U    = diag.U;

    // Take ω = √(max(ω², 0)). Clip negatives to small positive for the
    // BSH apply (which can't handle ω < 0 sensibly here). Negative ω²
    // signals a triplet instability or a bad seed; log and let
    // convergence stall if it persists.
    out.omega = madness::Tensor<double>(M);
    int n_clipped = 0;
    for (int s = 0; s < M; ++s) {
      const double w2 = om2(s);
      if (w2 < 0.0) ++n_clipped;
      out.omega(s) = std::sqrt(std::max(w2, 1.0e-6));
    }
    if (n_clipped > 0 && print_level_ >= PrintLevel::Normal
        && world_.rank() == 0) {
      print("[RPA] WARN: iter", out.iter, "  ", n_clipped,
            "negative ω² eigenvalue(s) clipped — may indicate triplet "
            "instability or bad seed.");
    }

    // ---- 3. Rotate u and BAu by U (subspace rotation, slot-preserving). --
    rs::transform(world_, out.roots, U);
    rs::transform(world_, BAu,       U);

    // ---- 4. Davidson residual + chained-BSH preconditioned step ----------
    //
    // For each root s:
    //   r_s   = BAu_s − ω²_s · u_s             (residual)
    //   δu_s  = G(+ω_s) · G(−ω_s) · (−2 r_s)   (chained BSH preconditioner)
    //   u_s  ← Q( u_s − δu_s ),  truncate
    //
    // Sign convention: BSH applies the kernel (T − ε − ω)⁻¹ with the
    // existing factor convention; we propagate r through both green's
    // functions and let the −2 absorb the kinetic factor exactly like
    // the linear TDA bsh_apply does. The signs here are consistent
    // with the linear case in the ω → 0 limit (where (A+B)(A−B) →
    // (T0+V0)² and the chained Green's function squared matches the
    // standard preconditioner). If the iteration diverges, FLIP the
    // sign of δu — it's the first thing to try.
    std::vector<Storage> x_pre_kain = out.roots;   // KAIN snapshot
    out.last_residual.assign(M, 0.0);
    out.last_u_change.assign(M, 0.0);
    for (int s = 0; s < M; ++s) {
      const double w   = out.omega(s);
      const double w2s = om2(s);

      // r = BAu − ω²·u
      auto r = madness::copy(world_, BAu[s].x_alpha);
      madness::gaxpy(world_, 1.0, r, -w2s, out.roots[s].x_alpha);
      out.last_residual[s] = madness::norm2(world_, r);

      // Chained-BSH preconditioner: δu = G(+ω) G(−ω) (−2 r).
      madness::scale(world_, r, -2.0);
      auto bsh_plus  = common_ops::make_bsh_operators(world_, gs_.aeps,
                                                       +w, gs_.lo);
      auto bsh_minus = common_ops::make_bsh_operators(world_, gs_.aeps,
                                                       -w, gs_.lo);
      auto step1 = apply(world_, bsh_minus, r);
      auto du    = apply(world_, bsh_plus,  step1);
      du = gs_.Qa(du);
      madness::truncate(world_, du);

      // u_new = u − δu
      auto u_new = madness::copy(world_, out.roots[s].x_alpha);
      madness::gaxpy(world_, 1.0, u_new, -1.0, du);
      madness::truncate(world_, u_new);

      // |du| as a step-size diagnostic.
      out.last_u_change[s] = madness::norm2(world_, du);
      out.roots[s].x_alpha = std::move(u_new);
    }

    // ---- 5. KAIN acceleration on u --------------------------------------
    const int kain_diag =
        (print_level_ >= PrintLevel::Verbose) ? 1 : 0;
    kain_.apply(x_pre_kain, out.roots, kain_diag);

    // Explosion guard
    for (double rv : out.last_residual)
      if (rv > policy_.explosion_guard) { out.diverged = true; break; }

    print_iter_banner(out);
    append_convergence_log(out);
    return out;
  }

  /// Converged when:
  ///   iter ≥ min_iters_before_conv  AND
  ///   max_s ‖BAu_s − ω²_s · u_s‖ < targets_.bsh_residual.
  /// (The "density residual" notion doesn't have a direct analog in
  /// this solver — `last_u_change` is informative but not gated.)
  bool converged(const State &s) const {
    if (s.diverged) return true;
    if (s.iter < policy_.min_iters_before_conv) return false;
    if (s.last_residual.empty()) return false;
    double mx = 0.0;
    for (double r : s.last_residual) mx = std::max(mx, r);
    return mx < targets_.bsh_residual;
  }

  void save(const State &s, const std::string &path_prefix) const {
    for (int b = 0; b < n_roots_; ++b)
      s.roots[b].save(world_, path_prefix + ".u_root_" + std::to_string(b));
  }

  void sort_state_by_omega(State &s) const {
    const long M = s.omega.dim(0);
    if (M <= 1) return;
    std::vector<long> perm(M);
    std::iota(perm.begin(), perm.end(), 0L);
    std::sort(perm.begin(), perm.end(),
              [&](long a, long b) { return s.omega(a) < s.omega(b); });
    bool ident = true;
    for (long i = 0; i < M; ++i)
      if (perm[i] != i) { ident = false; break; }
    if (ident) return;

    auto omega_new = madness::copy(s.omega);
    std::vector<Storage> roots_new(M);
    std::vector<double>  res_new(M, 0.0), du_new(M, 0.0);
    for (long i = 0; i < M; ++i) {
      const long src = perm[i];
      omega_new(i)  = s.omega(src);
      roots_new[i]  = std::move(s.roots[src]);
      if (src < static_cast<long>(s.last_residual.size()))
        res_new[i] = s.last_residual[src];
      if (src < static_cast<long>(s.last_u_change.size()))
        du_new[i] = s.last_u_change[src];
    }
    s.omega         = omega_new;
    s.roots         = std::move(roots_new);
    s.last_residual = std::move(res_new);
    s.last_u_change = std::move(du_new);
  }

  /// Recover the (X, Y) pair for each root from converged u_s.
  ///
  ///   v_s = (1/ω_s) · (A+B)·u_s,  X_s = ½(u_s + v_s),  Y_s = ½(u_s − v_s).
  ///
  /// Returns a ResponseStateXY<ClosedShell> per root, packaged the same
  /// way an ESSolver<Full, ClosedShell>::State carries them so callers
  /// can apply existing (X,Y)-aware post-processing (density assembly,
  /// oscillator strengths, save_es_roots, etc.).
  std::vector<ResponseStateXY<ClosedShell>>
  recover_xy(const State &s) const {
    const int M = static_cast<int>(s.roots.size());
    std::vector<ResponseStateXY<ClosedShell>> out(M);
    for (int r = 0; r < M; ++r) {
      const double w = s.omega(r);
      MADNESS_CHECK(w > 1.0e-9);
      auto Au = K::apply_AplusB(world_, gs_, s.roots[r]);
      auto v  = madness::copy(world_, Au.x_alpha);
      madness::scale(world_, v, 1.0 / w);

      // X = ½(u + v),  Y = ½(u − v)
      auto X = madness::copy(world_, s.roots[r].x_alpha);
      auto Y = madness::copy(world_, s.roots[r].x_alpha);
      madness::gaxpy(world_, 1.0, X,  1.0, v);
      madness::gaxpy(world_, 1.0, Y, -1.0, v);
      madness::scale(world_, X, 0.5);
      madness::scale(world_, Y, 0.5);
      madness::truncate(world_, X);
      madness::truncate(world_, Y);
      out[r].x_alpha = std::move(X);
      out[r].y_alpha = std::move(Y);
    }
    return out;
  }

private:
  void append_convergence_log(const State &s) {
    if (log_path_.empty() || world_.rank() != 0) return;
    std::ofstream out(log_path_, std::ios::app);
    if (!out) return;
    if (!log_header_written_) {
      out << "iter,protocol_thresh,state,omega,davidson_residual,"
          << "u_change,diverged\n";
      log_header_written_ = true;
    }
    const double pthr = madness::FunctionDefaults<3>::get_thresh();
    const long   M    = s.omega.dim(0);
    out.precision(12);
    for (long st = 0; st < M; ++st) {
      const double du = (static_cast<size_t>(st) < s.last_u_change.size())
                            ? s.last_u_change[st] : 0.0;
      const double r  = (static_cast<size_t>(st) < s.last_residual.size())
                            ? s.last_residual[st] : 0.0;
      out << s.iter << ',' << pthr << ',' << st << ','
          << s.omega(st) << ',' << r << ',' << du << ','
          << (s.diverged ? 1 : 0) << '\n';
    }
  }

  madness::World            &world_;
  ResponseGroundState        gs_;
  int                        n_roots_ = 0;
  ConvergencePolicy          policy_;
  ConvergencePolicy::Targets targets_{};
  ResponseSubspaceKain<Storage> kain_;
  PrintLevel                 print_level_ = PrintLevel::Normal;
  std::string                log_path_;
  bool                       log_header_written_ = false;
};

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_ES_SOLVER_FULL_RPA_HPP
