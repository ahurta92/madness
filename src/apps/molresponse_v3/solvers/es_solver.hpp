#ifndef MOLRESPONSE_V3_SOLVERS_ES_SOLVER_HPP
#define MOLRESPONSE_V3_SOLVERS_ES_SOLVER_HPP

// =========================================================================
// ESSolver<Type, Shell> — one skeleton, six instantiations.
//
// step() body is identical across (Type, Shell). Every per-root
// kernel call has the signature K::foo(world_, gs_, state[,
// intermediate]); the Kernels<Type, Shell> specialization decides
// which fields it reads. Bundle-level ops (subspace matrix,
// diagonalize, rotate) are shared free functions in
// kernels/bundle_helpers.hpp.
//
// FDSolver<Type, Shell> will share the SAME Kernels<T,S> bodies for
// per-root quantities (gamma, V0·x, residual) — the FD/ES axis is in
// the outer iteration shape only.
//
// Convergence policy and per-iter density tracking live here so the
// FD solver can reuse [[ConvergencePolicy]] without duplication.
// =========================================================================

#include "../kernels/response_space_ops.hpp"
#include "../kernels/tags.hpp"
#include "../kernels/tda.hpp"   // Kernels<TDA, *>; add kernels/full.hpp etc.
#include "convergence_policy.hpp"
#include "iterate.hpp"
#include "response_state.hpp"
#include "response_subspace_kain.hpp"

#include <madness/mra/mra.h>
#include <madness/tensor/tensor.h>

#include <algorithm>
#include <fstream>
#include <string>
#include <type_traits>
#include <vector>

namespace molresponse_v3 {

/// Read-only problem definition for ESSolver<Type, Shell>. Symmetric
/// with FDProblem<Type, Shell> in fd_problem.hpp:
///
///   gs       — the prepared ground state (zeroth-order data) the
///              kernels read from. Same object FDSolver uses.
///   n_roots  — solver-level count of excitation roots to track.
///
/// `n_roots` is a solver concern, not a ground-state property, so it
/// lives on the problem alongside gs rather than inside the gs struct.
template <typename Type, typename Shell>
struct ESProblem {
  ResponseGroundState gs;
  int                 n_roots = 0;
};

template <typename Type, typename Shell>
class ESSolver {
public:
  using K       = Kernels<Type, Shell>;
  using Storage = StorageOf_t<Type, Shell>;
  using vecfuncT = std::vector<madness::real_function_3d>;

  struct State {
    std::vector<Storage>                    roots;
    madness::Tensor<double>                 omega;
    std::vector<double>                     last_bsh_residual;
    /// Per-root previous-iter density (for Δρ tracking). Empty on
    /// iter 0, populated thereafter; index s tracks the same root
    /// slot as roots[s] before rotation.
    std::vector<madness::real_function_3d>  rho_alpha_prev;
    std::vector<double>                     last_density_residual;
    int                                     iter = 0;
    /// Set by step() when the explosion guard trips; iterate<>
    /// terminates on this so we don't burn iters on a runaway state.
    bool                                    diverged = false;
  };

  /// Preferred ctor — caller supplies a problem and a policy.
  ESSolver(madness::World &world, ESProblem<Type, Shell> problem,
           ConvergencePolicy policy,
           PrintLevel print_level = PrintLevel::Normal)
      : world_(world), gs_(std::move(problem.gs)),
        n_roots_(problem.n_roots),
        policy_(policy), kain_(world, policy),
        print_level_(print_level) {
    refresh_convergence_targets();
    print_header();
  }

  /// Future variant — low-memory Lambda assembly:
  ///   The current step() materializes T0x, V0x, E0x_full, gamma as
  ///   four separate response_spaces before summing them into Lambda.
  ///   For larger systems where holding three or four response_spaces
  ///   simultaneously is prohibitive, an alternative is to fold each
  ///   term in place: build Lambda by += T0x, += V0x, -= E0x_full, +=
  ///   gamma, freeing each piece as soon as it has been added. Doing
  ///   so would change compute_lambda's signature from "take all
  ///   pieces" to a streaming-accumulator form. Leaving the current
  ///   shape for clarity; revisit when memory pressure on a target
  ///   system demands it.

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

  /// Swap in a fresh ground state — used by iterate_protocol's
  /// prepare callback after the ground-state side has been re-prepped
  /// at a new (k, thresh). Carries everything protocol-sensitive that
  /// isn't part of State (phi0, V_local, coulop, Q, fock). KAIN's
  /// stored iterates are reset because they live at the old wavelet
  /// basis. n_roots is unchanged.
  void set_gs(ResponseGroundState gs) {
    gs_ = std::move(gs);
    kain_.reset();
  }
  const ResponseGroundState& gs()      const { return gs_; }
  int                        n_roots() const { return n_roots_; }

  /// Enable per-iter convergence logging to a CSV file. One row per
  /// (iter, state) appended at the end of each step(). Header is
  /// written on the first row of a fresh file. Pass empty string to
  /// disable.
  void set_log_path(const std::string &path) {
    log_path_ = path;
    log_header_written_ = false;
  }

  /// Called by iterate_protocol after FunctionDefaults<3>::thresh has
  /// changed. Re-derives effective BSH / density targets from the
  /// fresh thresh.
  void refresh_convergence_targets() {
    targets_ = policy_.effective_for_thresh(
        madness::FunctionDefaults<3>::get_thresh());
  }
  const ConvergencePolicy::Targets& targets() const { return targets_; }

private:
  void print_header() const {
    if (print_level_ < PrintLevel::Normal || world_.rank() != 0) return;
    print("");
    print("ESSolver<TDA, ClosedShell>  num_roots =", n_roots_,
          " thresh =", madness::FunctionDefaults<3>::get_thresh(),
          " c_xc =", gs_.c_xc);
    print("  policy: dconv_user =", policy_.dconv_user,
          " bsh_target =", targets_.bsh_residual,
          " density_target =", targets_.density_residual,
          " explosion_guard =", policy_.explosion_guard);
    print("");
  }

  void print_iter_banner(const State &out) const {
    if (print_level_ < PrintLevel::Normal || world_.rank() != 0) return;
    double max_res  = 0.0;
    for (double r : out.last_bsh_residual) max_res = std::max(max_res, r);
    double max_drho = 0.0;
    for (double r : out.last_density_residual) max_drho = std::max(max_drho, r);
    print("iter", out.iter, "  omega =", out.omega,
          "  max_res =", max_res, "  max_dρ =", max_drho);

    if (print_level_ >= PrintLevel::Verbose) {
      for (size_t s = 0; s < out.last_bsh_residual.size(); ++s) {
        double dr = (s < out.last_density_residual.size())
                        ? out.last_density_residual[s] : 0.0;
        print("  root", s, "  res =", out.last_bsh_residual[s],
              "  dρ =", dr,
              "  omega =", out.omega(static_cast<long>(s)));
      }
    }
  }

  // Debug-tier per-iter dumps: A, S, U, omega tensors + selected
  // <X|.|X> inner-product matrices.
  void print_debug_iter(const madness::Tensor<double> &A,
                        const madness::Tensor<double> &S_mat,
                        const madness::Tensor<double> &omega_new,
                        const madness::Tensor<double> &U) const {
    if (print_level_ < PrintLevel::Debug || world_.rank() != 0) return;
    print("[DEBUG] A (subspace):");  print(A);
    print("[DEBUG] S (overlap):");    print(S_mat);
    print("[DEBUG] omega:");          print(omega_new);
    print("[DEBUG] U:");              print(U);
  }

  /// Slot-identity diagnostic — emitted per iter at Verbose+. Shows
  /// whether dominance-sort swapped slots (perm != identity) and
  /// whether the post-fixup U(i,i) stayed near ±1 (slot identity
  /// preserved) or dropped (slot picked up character from another).
  /// Pair this with the per-state residual columns to correlate
  /// slot-flip events with KAIN history corruption.
  void print_rot_slots(int iter,
                       const rs::DiagonalizeResult &r) const {
    if (print_level_ < PrintLevel::Verbose || world_.rank() != 0) return;
    const long M = static_cast<long>(r.dominance_perm.size());
    bool ident = true;
    for (long i = 0; i < M; ++i)
      if (r.dominance_perm[i] != i) { ident = false; break; }
    printf("[ROT-SLOTS] iter=%d  perm=[", iter);
    for (long i = 0; i < M; ++i)
      printf(" %ld", r.dominance_perm[i]);
    printf(" ]  U_diag=[");
    for (long i = 0; i < r.U_diag.dim(0); ++i)
      printf(" %+.3f", r.U_diag(i));
    printf(" ]  omega_asc=[");
    for (long i = 0; i < r.omega_ascending.dim(0); ++i)
      printf(" %+.5f", r.omega_ascending(i));
    printf(" ]%s\n", ident ? "  identity" : "  REORDERED");
    fflush(stdout);
  }

  void print_debug_norms(int s, const State &in,
                         const std::vector<madness::real_function_3d> &V0x_s,
                         const std::vector<madness::real_function_3d> &T0x_s,
                         const std::vector<madness::real_function_3d> &gamma_s,
                         const std::vector<madness::real_function_3d> &theta_s) const {
    if (print_level_ < PrintLevel::Debug) return;
    // Collective: every rank must call norm2.
    double nx = madness::norm2(world_, in.roots[s].x_alpha);
    double nv = madness::norm2(world_, V0x_s);
    double nt = madness::norm2(world_, T0x_s);
    double ng = madness::norm2(world_, gamma_s);
    double nh = madness::norm2(world_, theta_s);
    if (world_.rank() != 0) return;
    printf("[NORMS] iter=%d root=%d  |x|=%.3e  |V0x|=%.3e  |T0x|=%.3e  "
           "|gamma|=%.3e  |theta|=%.3e\n",
           in.iter, s, nx, nv, nt, ng, nh);
    fflush(stdout);
  }

public:
  /// Called by iterate<> after the loop exits. Reports converged or
  /// max-iter-reached banner.
  void print_final(const State &s, bool converged) const {
    if (print_level_ < PrintLevel::Normal || world_.rank() != 0) return;
    print("");
    // iterate<> sets `converged=true` whenever solver.converged()
    // returns true. converged() also returns true on diverged so the
    // loop exits — so distinguish here by checking s.diverged first.
    if (s.diverged)       print("Stopped at iter", s.iter,
                                "(diverged — residual exceeded "
                                "explosion guard).");
    else if (converged)   print("Converged in", s.iter, "iters.");
    else                  print("Stopped at iter", s.iter,
                                "(max iters reached, not converged).");
    print("  omega_final =", s.omega);
    if (!s.last_bsh_residual.empty()) {
      print("  residuals   =");
      for (size_t i = 0; i < s.last_bsh_residual.size(); ++i) {
        double dr = (i < s.last_density_residual.size())
                        ? s.last_density_residual[i] : 0.0;
        print("    root", i, ":  bsh =", s.last_bsh_residual[i],
              "  dρ =", dr);
      }
    }
    print("");
  }

  /// One outer iteration. Dispatches on `policy_.stream_theta`:
  ///   - false (default): step_rotate_pieces() — keeps V0x/E0x/gamma
  ///     for all roots in memory, rotates them by U, assembles Theta
  ///     from the rotated pieces. ~8M Storage peak. Faster (no
  ///     recompute) but heavier.
  ///   - true: step_recompute_pieces() — drops V0x/E0x/gamma after
  ///     building Lambda, rotates only X, then recomputes V0x/E0x/
  ///     gamma against the rotated X to stream Theta in place.
  ///     ~3-4M Storage peak. ~2× the per-iter kernel work but big
  ///     memory savings.
  State step(State in) {
    return policy_.stream_theta
        ? step_recompute_pieces(std::move(in))
        : step_rotate_pieces(std::move(in));
  }

private:
  /// Top-of-iter Q-projection + orthonormalization (mode picked by
  /// `policy_.orthonormalize_mode`). Shared by both step() variants.
  void apply_top_of_iter_discipline(State &out, int M) {
    // Q-project per spin.
    for (int s = 0; s < M; ++s) {
      out.roots[s].x_alpha = gs_.Qa(out.roots[s].x_alpha);
      if constexpr (std::is_same_v<Shell, OpenShell>) {
        out.roots[s].x_beta = gs_.Qb(out.roots[s].x_beta);
      }
    }

    // Flatten — α (+ β for OpenShell) concat per state. The flat
    // inner product gives the OpenShell bundle metric for free.
    response_space flats(M);
    for (int s = 0; s < M; ++s) flats[s] = out.roots[s].flatten();

    if (policy_.orthonormalize_mode ==
        ConvergencePolicy::OrthonormalizeMode::Lowdin) {
      // diag_level is derived from print_level_, which is uniform
      // across ranks, so all ranks enter or skip the diagnostic
      // collectives in lock-step inside lowdin_orthonormalize.
      const int diag_level =
          (print_level_ >= PrintLevel::Verbose) ? 1 : 0;
      const int dropped = rs::lowdin_orthonormalize(
          world_, flats, /*thresh=*/-1.0, diag_level);
      if (dropped > 0 && print_level_ >= PrintLevel::Normal &&
          world_.rank() == 0) {
        print("[LOWDIN] iter", out.iter, ": dropped", dropped,
              "near-singular direction(s) (bundle partially collapsed)");
      }
    } else {
      // Gram-Schmidt, 2-pass + normalize. Order-dependent: slot 0
      // anchored, higher slots projected against lower.
      for (int pass = 0; pass < 2; ++pass) {
        for (int i = 0; i < M; ++i) {
          for (int j = 0; j < i; ++j) {
            const double c = madness::inner(flats[j], flats[i]);
            madness::gaxpy(world_, 1.0, flats[i], -c, flats[j]);
          }
          const double n2 = madness::inner(flats[i], flats[i]);
          const double nrm = std::sqrt(std::abs(n2));
          if (nrm > 1.0e-12)
            madness::scale(world_, flats[i], 1.0 / nrm);
        }
      }
    }

    for (int s = 0; s < M; ++s) out.roots[s].from_flat(flats[s]);
  }

  /// Pack a per-root Storage into a flat vecfuncT (α + β concat).
  /// Used by both variants for the subspace rs::* ops.
  vecfuncT pack_storage(const Storage &s) const {
    vecfuncT v = s.x_alpha;
    if constexpr (std::is_same_v<Shell, OpenShell>) {
      v.insert(v.end(), s.x_beta.begin(), s.x_beta.end());
    }
    return v;
  }

  /// Unpack a flat vecfuncT back into a Storage's α (+ β) components.
  void unpack_storage(Storage &s, const vecfuncT &v,
                     std::size_t n_alpha, std::size_t n_beta) const {
    if constexpr (std::is_same_v<Shell, OpenShell>) {
      s.x_alpha = vecfuncT(v.begin(), v.begin() + n_alpha);
      s.x_beta  = vecfuncT(v.begin() + n_alpha,
                           v.begin() + n_alpha + n_beta);
    } else {
      (void)n_beta;
      (void)n_alpha;
      s.x_alpha = v;
    }
  }

  /// Final stage shared by both variants — KAIN apply, explosion
  /// guard, banner, log. Mutates out in place.
  void finalize_iter(const State &in, State &out,
                     const std::vector<Storage> &x_pre_bsh) {
    (void)in;
    const int kain_diag =
        (print_level_ >= PrintLevel::Verbose) ? 1 : 0;
    kain_.apply(x_pre_bsh, out.roots, kain_diag);
    for (double r : out.last_bsh_residual) {
      if (r > policy_.explosion_guard) { out.diverged = true; break; }
    }
    print_iter_banner(out);
    append_convergence_log(out);
  }

public:
  /// "Keep pieces, rotate them, assemble Theta from rotated pieces."
  /// Default variant — fastest, highest memory.
  /// Algorithm:
  ///   1. Per root, build V0x, T0x, E0x_full, E0x (no-diag), gamma.
  ///   2. Assemble Lambda = T0x + V0x - E0x_full + gamma.
  ///   3. Build A = <X | Lambda>, S = <X | X>, diagonalize → omega, U.
  ///   4. Rotate {X, V0x, E0x, gamma} by U  (T0x and E0x_full are
  ///      Lambda-only and can be discarded).
  ///   5. Assemble Theta = V0x - E0x + gamma  (from rotated pieces).
  ///   6. BSH apply → new X; residual = || X_old - X_new ||.
  State step_rotate_pieces(State in) {
    State out;
    out.iter  = in.iter + 1;
    out.roots = in.roots;
    const int M = n_roots_;

    // ---- 0. top-of-iter discipline ----------------------------------------
    // Strip ground-orbital contamination that accumulated during the
    // previous iter's BSH + KAIN steps; without this, the residual
    // ground-state component is amplified by BSH into the
    // ω = -ε_core ghost eigenvector. Then re-orthonormalize the
    // bundle so the subspace step sees a clean basis. Mirrors legacy
    // ExcitedResponse.cpp:2483-2517 (top-of-iter Q + Gram-Schmidt or
    // Löwdin, depending on policy.orthonormalize_mode).
    apply_top_of_iter_discipline(out, M);

    // ---- 1. per-root building blocks --------------------------------------
    // Intermediates are Storage-typed: TDA-ClosedShell carries only
    // x_alpha; Full-ClosedShell will carry x_alpha + y_alpha; OpenShell
    // adds the *_beta members. Kernel signatures don't change between
    // (Type, Shell); the wrapping is invisible to step().
    std::vector<Storage> V0x(M), T0x(M), E0x_full(M), E0x(M), gamma(M);
    std::vector<madness::real_function_3d> rho_alpha(M);
    out.last_density_residual.assign(M, 0.0);
    for (int s = 0; s < M; ++s) {
      rho_alpha[s]  = K::compute_density(world_, gs_, out.roots[s]);
      gamma[s]      = K::compute_gamma(world_, gs_, out.roots[s],
                                       rho_alpha[s]);
      V0x[s]        = K::compute_V0x(world_, gs_, out.roots[s]);
      T0x[s]        = K::compute_T0x(world_, gs_, out.roots[s]);
      E0x_full[s]   = K::compute_E0x_full(world_, gs_, out.roots[s]);
      E0x[s]        = K::compute_E0x(world_, gs_, out.roots[s]);

      // Δρ vs previous iter (per slot — see note in State).
      if (!in.rho_alpha_prev.empty()) {
        auto drho = rho_alpha[s] - in.rho_alpha_prev[s];
        out.last_density_residual[s] = drho.norm2();
      }
    }
    out.rho_alpha_prev = std::move(rho_alpha);

    // ---- 2. assemble Lambda per root --------------------------------------
    std::vector<Storage> lambda(M);
    for (int s = 0; s < M; ++s) {
      lambda[s] = K::compute_lambda(world_, T0x[s], V0x[s],
                                    E0x_full[s], gamma[s]);
    }

    // ---- 3. subspace A, S, diagonalize ------------------------------------
    // Subspace ops live on response_space (= vector<vecfuncT>). For
    // OpenShell ES TDA, the subspace matrix elements sum over both
    // spin blocks:
    //   A_ij = <x_α_i | Λ_α_j> + <x_β_i | Λ_β_j>
    // We get this for free by CONCATENATING α and β into a single
    // vecfuncT per state — rs::inner sums over all orbitals already.
    // Same rotation U applies to both spin blocks (it mixes ROOTS,
    // not spins), so transform() on the concat is identical to
    // transforming each spin separately.
    constexpr bool open_shell = std::is_same_v<Shell, OpenShell>;
    std::size_t n_alpha = 0, n_beta = 0;
    if constexpr (open_shell) {
      n_alpha = out.roots[0].x_alpha.size();
      n_beta  = out.roots[0].x_beta.size();
    }
    auto pack = [&](const Storage &s) -> vecfuncT {
      vecfuncT v = s.x_alpha;
      if constexpr (open_shell) {
        v.insert(v.end(), s.x_beta.begin(), s.x_beta.end());
      }
      return v;
    };
    auto unpack = [&](Storage &s, const vecfuncT &v) {
      if constexpr (open_shell) {
        s.x_alpha = vecfuncT(v.begin(), v.begin() + n_alpha);
        s.x_beta  = vecfuncT(v.begin() + n_alpha,
                             v.begin() + n_alpha + n_beta);
      } else {
        s.x_alpha = v;
      }
    };

    response_space X_rs(M), Lambda_rs(M);
    for (int s = 0; s < M; ++s) {
      X_rs[s]      = pack(out.roots[s]);
      Lambda_rs[s] = pack(lambda[s]);
    }
    auto A     = rs::inner(X_rs, Lambda_rs);
    auto S_mat = rs::inner(X_rs, X_rs);
    auto diag_result = rs::diagonalize(A, S_mat,
                                       /*thresh_degenerate=*/-1.0,
                                       policy_.cluster_unmix_factor);
    auto &omega_new = diag_result.omega;
    auto &U         = diag_result.U;
    out.omega = omega_new;
    print_debug_iter(A, S_mat, omega_new, U);
    print_rot_slots(out.iter, diag_result);

    // ---- 4. rotate the response_spaces needed for Theta -------------------
    response_space V0x_rs(M), E0x_rs(M), gamma_rs(M);
    for (int s = 0; s < M; ++s) {
      V0x_rs[s]   = pack(V0x[s]);
      E0x_rs[s]   = pack(E0x[s]);
      gamma_rs[s] = pack(gamma[s]);
    }
    rs::transform(world_, X_rs,    U);
    rs::transform(world_, V0x_rs,  U);
    rs::transform(world_, E0x_rs,  U);
    rs::transform(world_, gamma_rs, U);
    for (int s = 0; s < M; ++s) {
      unpack(out.roots[s], X_rs[s]);
      unpack(V0x[s],       V0x_rs[s]);
      unpack(E0x[s],       E0x_rs[s]);
      unpack(gamma[s],     gamma_rs[s]);
    }

    // ---- 5. assemble Theta from rotated pieces ---------------------------
    std::vector<Storage> theta(M);
    for (int s = 0; s < M; ++s) {
      theta[s] = K::compute_theta(world_, V0x[s], E0x[s], gamma[s]);
    }

    // ---- 6. BSH apply + residual ------------------------------------------
    // out.roots[s] currently holds the rotated state (pre-BSH). BSH
    // updates each root in place; we capture the raw BSH residual
    // BEFORE KAIN damping (it's what tracks convergence of the
    // underlying dynamics, not the damped step size).
    //
    // Save the rotated pre-BSH state so KAIN's "x_old" is the
    // self-consistent input to the BSH step, not the pre-rotation
    // input. Without this, KAIN's history would mix pre/post rotation
    // and produce garbage iterates.
    std::vector<Storage> x_pre_bsh = out.roots;
    out.last_bsh_residual.assign(M, 0.0);
    for (int s = 0; s < M; ++s) {
      auto x_new = K::bsh_apply(world_, gs_, out.roots[s],
                                theta[s], omega_new(s));
      out.last_bsh_residual[s] =
          K::compute_residual_norm(world_, out.roots[s], x_new);
      print_debug_norms(s, in, V0x[s].x_alpha, T0x[s].x_alpha,
                        gamma[s].x_alpha, theta[s].x_alpha);
      out.roots[s] = std::move(x_new);
    }

    // ---- 7. KAIN acceleration + step restriction --------------------------
    // Feed (rotated x_pre_bsh) → (raw BSH output in out.roots).
    // policy_.kain=false disables.
    {
      const int kain_diag =
          (print_level_ >= PrintLevel::Verbose) ? 1 : 0;
      kain_.apply(x_pre_bsh, out.roots, kain_diag);
    }

    // Explosion guard — if any per-root BSH residual blows up past
    // policy_.explosion_guard, mark the state diverged and let iterate<>
    // stop. Legacy iterate_excited.cc:164 uses the same heuristic.
    for (double r : out.last_bsh_residual) {
      if (r > policy_.explosion_guard) { out.diverged = true; break; }
    }

    print_iter_banner(out);
    append_convergence_log(out);
    return out;
  }

  /// "Stream Lambda, drop pieces, recompute V0x/E0x/gamma after
  /// rotation to assemble Theta in place." Memory-conscious variant.
  /// Algorithm:
  ///   0. Top-of-iter Q + GS (shared).
  ///   1+2. Per root, stream Lambda = T0x + V0x − E0x_full + gamma,
  ///        freeing each kernel temporary after it's folded in.
  ///   3. Subspace A = <X|Λ>, S = <X|X>, diagonalize → omega, U.
  ///   4. Rotate X only (drop Lambda).
  ///   5. RECOMPUTE V0x, E0x_nodi, gamma against rotated X, stream
  ///      Theta = V0x − E0x + gamma in place.
  ///   6. BSH apply → new X.
  ///   7. KAIN (shared).
  ///
  /// Peak memory: ~3-4M Storage (vs ~8M for step_rotate_pieces).
  /// Cost: ~2× the per-iter kernel work (each V0x/E0x/gamma computed
  /// twice — once for Lambda, once for Theta).
  State step_recompute_pieces(State in) {
    State out;
    out.iter  = in.iter + 1;
    out.roots = in.roots;
    const int M = n_roots_;
    const double thr = madness::FunctionDefaults<3>::get_thresh();

    // ---- 0. top-of-iter discipline ----------------------------------------
    apply_top_of_iter_discipline(out, M);

    // ---- 1 + 2. per-root Lambda assembly (streamed) -----------------------
    std::vector<Storage> lambda(M);
    std::vector<madness::real_function_3d> rho_alpha(M);
    out.last_density_residual.assign(M, 0.0);
    for (int s = 0; s < M; ++s) {
      rho_alpha[s] = K::compute_density(world_, gs_, out.roots[s]);
      if (!in.rho_alpha_prev.empty()) {
        auto drho = rho_alpha[s] - in.rho_alpha_prev[s];
        out.last_density_residual[s] = drho.norm2();
      }

      // Stream Lambda = T0x + V0x − E0x_full + gamma. Each temporary
      // is scoped: built, axpy'd into lambda[s], then freed.
      lambda[s] = K::compute_T0x(world_, gs_, out.roots[s]);
      {
        auto V = K::compute_V0x(world_, gs_, out.roots[s]);
        lambda[s].axpy(world_, +1.0, V);
      }
      {
        auto Efull = K::compute_E0x_full(world_, gs_, out.roots[s]);
        lambda[s].axpy(world_, -1.0, Efull);
      }
      {
        auto G = K::compute_gamma(world_, gs_, out.roots[s], rho_alpha[s]);
        lambda[s].axpy(world_, +1.0, G);
      }
      lambda[s].truncate_all(world_, thr);
    }
    out.rho_alpha_prev = std::move(rho_alpha);

    // ---- 3. subspace A, S, diagonalize ------------------------------------
    constexpr bool open_shell = std::is_same_v<Shell, OpenShell>;
    std::size_t n_alpha = 0, n_beta = 0;
    if constexpr (open_shell) {
      n_alpha = out.roots[0].x_alpha.size();
      n_beta  = out.roots[0].x_beta.size();
    }
    response_space X_rs(M), Lambda_rs(M);
    for (int s = 0; s < M; ++s) {
      X_rs[s]      = pack_storage(out.roots[s]);
      Lambda_rs[s] = pack_storage(lambda[s]);
    }
    auto A     = rs::inner(X_rs, Lambda_rs);
    auto S_mat = rs::inner(X_rs, X_rs);
    auto diag_result = rs::diagonalize(A, S_mat,
                                       /*thresh_degenerate=*/-1.0,
                                       policy_.cluster_unmix_factor);
    auto &omega_new = diag_result.omega;
    auto &U         = diag_result.U;
    out.omega = omega_new;
    print_debug_iter(A, S_mat, omega_new, U);
    print_rot_slots(out.iter, diag_result);

    // ---- 4. rotate X only; drop Lambda ------------------------------------
    rs::transform(world_, X_rs, U);
    for (int s = 0; s < M; ++s)
      unpack_storage(out.roots[s], X_rs[s], n_alpha, n_beta);
    lambda.clear();  // free Lambda — not needed past here

    // ---- 5 + 6. Recompute V0x/E0x/gamma against rotated X, stream Theta
    //              + BSH apply per root --------------------------------------
    // KAIN's "x_old" is the rotated state (pre-BSH), same as
    // step_rotate_pieces. Snapshot before BSH overwrites out.roots.
    std::vector<Storage> x_pre_bsh = out.roots;
    out.last_bsh_residual.assign(M, 0.0);
    for (int s = 0; s < M; ++s) {
      // ρ on rotated X (different from the pre-rotation ρ we used for
      // Lambda — gamma here needs the rotated coupling).
      auto rho_rot = K::compute_density(world_, gs_, out.roots[s]);

      // Stream Theta = V0x − E0x_nodi + gamma, scoped temporaries.
      auto theta = K::compute_V0x(world_, gs_, out.roots[s]);
      {
        auto E = K::compute_E0x(world_, gs_, out.roots[s]);
        theta.axpy(world_, -1.0, E);
      }
      {
        auto G = K::compute_gamma(world_, gs_, out.roots[s], rho_rot);
        theta.axpy(world_, +1.0, G);
      }
      theta.truncate_all(world_, thr);

      // BSH
      auto x_new = K::bsh_apply(world_, gs_, out.roots[s],
                                theta, omega_new(s));
      out.last_bsh_residual[s] =
          K::compute_residual_norm(world_, out.roots[s], x_new);
      out.roots[s] = std::move(x_new);
    }

    // ---- 7. KAIN + guard + banner + log -----------------------------------
    finalize_iter(in, out, x_pre_bsh);
    return out;
  }

  /// Converged when:
  ///   iter >= min_iters_before_conv  AND
  ///   max_s BSH-residual[s] < targets_.bsh_residual  AND
  ///   max_s density-residual[s] < targets_.density_residual
  ///                              (skipped on iter 1: rho_prev empty)
  ///   OR  state.diverged is set (we exit either way).
  bool converged(const State &s) const {
    if (s.diverged) return true;  // exit; print_final reports diverged
    if (s.iter < policy_.min_iters_before_conv) return false;
    if (s.last_bsh_residual.empty()) return false;

    double mx_bsh = 0.0;
    for (double r : s.last_bsh_residual) mx_bsh = std::max(mx_bsh, r);
    if (mx_bsh >= targets_.bsh_residual) return false;

    // First iter: rho_prev was empty going in, density residuals stay
    // zero — accept BSH-only convergence in that one case.
    if (s.iter <= 1) return true;

    double mx_drho = 0.0;
    for (double r : s.last_density_residual) mx_drho = std::max(mx_drho, r);
    return mx_drho < targets_.density_residual;
  }

  void save(const State &s, const std::string &path_prefix) const {
    for (int b = 0; b < n_roots_; ++b) {
      s.roots[b].save(world_, path_prefix + ".root_" + std::to_string(b));
    }
  }

  /// Permute the state so omega is ascending. The diagonal-dominance
  /// sort inside rs::diagonalize preserves slot identity from the
  /// initial guess (so per-slot densities don't ping-pong across
  /// near-degeneracies), but the initial guess slot order is
  /// arbitrary. Call this once after convergence to get a canonical
  /// ordering for output / archive comparison against legacy.
  void sort_state_by_omega(State &s) const {
    const long M = s.omega.dim(0);
    if (M <= 1) return;

    std::vector<long> perm(M);
    std::iota(perm.begin(), perm.end(), 0L);
    std::sort(perm.begin(), perm.end(),
              [&](long a, long b) { return s.omega(a) < s.omega(b); });

    // No-op fast path.
    bool already_sorted = true;
    for (long i = 0; i < M; ++i) {
      if (perm[i] != i) { already_sorted = false; break; }
    }
    if (already_sorted) return;

    auto omega_new = madness::copy(s.omega);
    std::vector<Storage> roots_new(M);
    std::vector<double>  bsh_new(M, 0.0), drho_new(M, 0.0);
    std::vector<madness::real_function_3d> rho_prev_new(s.rho_alpha_prev.size());
    for (long i = 0; i < M; ++i) {
      const long src = perm[i];
      omega_new(i)  = s.omega(src);
      roots_new[i]  = std::move(s.roots[src]);
      if (src < static_cast<long>(s.last_bsh_residual.size()))
        bsh_new[i] = s.last_bsh_residual[src];
      if (src < static_cast<long>(s.last_density_residual.size()))
        drho_new[i] = s.last_density_residual[src];
      if (src < static_cast<long>(s.rho_alpha_prev.size()))
        rho_prev_new[i] = s.rho_alpha_prev[src];
    }
    s.omega                 = omega_new;
    s.roots                 = std::move(roots_new);
    s.last_bsh_residual     = std::move(bsh_new);
    s.last_density_residual = std::move(drho_new);
    s.rho_alpha_prev        = std::move(rho_prev_new);

    if (print_level_ >= PrintLevel::Normal && world_.rank() == 0) {
      print("  [sort_state_by_omega] permuted slots to ascending omega:");
      for (long i = 0; i < M; ++i) {
        print("    new slot", i, "<- old slot", perm[i],
              "  omega =", s.omega(i));
      }
    }
  }

private:
  /// Append one row per state at the end of step() to a CSV file at
  /// `log_path_`. Schema:
  ///   iter,protocol_thresh,state,omega,bsh_residual,density_residual,diverged
  /// Header is written exactly once per file (tracked via flag). Only
  /// rank 0 writes; CSV is robust against concurrent reads.
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
    const long M = s.omega.dim(0);
    out.precision(12);
    for (long st = 0; st < M; ++st) {
      const double drho = (static_cast<size_t>(st) < s.last_density_residual.size())
                              ? s.last_density_residual[st] : 0.0;
      const double bsh  = (static_cast<size_t>(st) < s.last_bsh_residual.size())
                              ? s.last_bsh_residual[st] : 0.0;
      out << s.iter << ',' << pthr << ',' << st << ','
          << s.omega(st) << ',' << bsh << ','
          << drho << ',' << (s.diverged ? 1 : 0) << '\n';
    }
  }

  madness::World            &world_;
  ResponseGroundState        gs_;
  int                        n_roots_ = 0;
  ConvergencePolicy          policy_;
  ConvergencePolicy::Targets targets_{};
  ResponseSubspaceKain<Storage>   kain_;
  PrintLevel                 print_level_ = PrintLevel::Normal;
  std::string                log_path_;
  bool                       log_header_written_ = false;
};

// Build helpers (ResponseGroundState + initial guess) live in
// build_response_ground_state.hpp so this header doesn't depend on
// GroundState or the legacy ESSolverGuess.

} // namespace molresponse_v3

// Call-site shape — identical across (Type, Shell):
//
//   auto problem = build_es_problem_tda<ClosedShell>(world, gs, n_roots);
//   ESSolver<TDA, ClosedShell> solver(world, problem,
//                                     ConvergencePolicy{1e-4});
//   typename ESSolver<TDA, ClosedShell>::State s0 =
//       build_initial_guess(...);
//   auto sf = solvers::iterate(solver, s0, {/*max_iters=*/25});
//   solver.save(sf, save_prefix);
//
// With protocol ladder:
//   sf = solvers::iterate_protocol(
//          solver, s0,
//          /*thresholds=*/{1e-4, 1e-6, 1e-7},
//          /*set_protocol=*/[&](double th){
//             set_response_protocol(world, L, th, override_k);
//          },
//          /*max_iters_per_step=*/25);

#endif // MOLRESPONSE_V3_SOLVERS_ES_SOLVER_HPP
