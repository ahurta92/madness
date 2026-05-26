#ifndef MOLRESPONSE_V3_KERNELS_RESPONSE_SPACE_OPS_HPP
#define MOLRESPONSE_V3_KERNELS_RESPONSE_SPACE_OPS_HPP

// =========================================================================
// Type/Shell-agnostic operations on a response_space.
//
// A response_space here is the same object as in the original
// molresponse codebase and Hurtado thesis: a collection of vectors of
// orbital functions, one per excited state (or per perturbation),
// arranged as an (n_states × n_orbitals) "matrix". Each row is one
// state's response orbitals.
//
//   using response_space = std::vector<vector_real_function_3d>;
//
// rs::inner mirrors the efficient form of molresponse's
// response_space_inner (src/apps/molresponse/response_functions.h):
// transpose to (n_orbitals × n_states) once, then accumulate
// matrix_inner per orbital index. For m states and n orbitals that's
// n batched matrix_inners of size m×m rather than m² scalar inners.
// =========================================================================

#include <madness/mra/mra.h>
#include <madness/tensor/tensor.h>
#include <madness/tensor/tensor_lapack.h>

#include <utility>
#include <vector>

namespace molresponse_v3 {

/// A response_space: n_states × n_orbitals matrix of 3D functions.
/// Each row is one state's response orbitals. Thesis / molresponse
/// vocabulary.
using response_space = std::vector<madness::vector_real_function_3d>;

namespace rs {

/// Inner-product matrix between two response_spaces.
///
///   result(i, j) = sum_p <a[i][p] | b[j][p]>
///
/// Efficient form (ported from molresponse::response_space_inner):
/// transpose each operand to (n_orbitals × n_states), then one
/// matrix_inner per orbital index. n_orbitals batched calls instead
/// of m² scalar inners.
inline madness::Tensor<double>
inner(const response_space &a, const response_space &b) {
  MADNESS_CHECK(!a.empty());
  MADNESS_CHECK(!b.empty());
  MADNESS_CHECK(a[0].size() == b[0].size());

  const long m_a   = static_cast<long>(a.size());
  const long m_b   = static_cast<long>(b.size());
  const long n_orb = static_cast<long>(a[0].size());

  madness::World &world = a[0][0].world();

  // matrix_inner requires compressed form.
  for (const auto &vi : a) compress(world, vi, false);
  for (const auto &vi : b) compress(world, vi, false);
  world.gop.fence();

  // Transposed views: aT[p][i] = a[i][p].
  // Function<double,3> has shared_ptr semantics; this aliases rather
  // than copies the underlying impls.
  std::vector<madness::vector_real_function_3d> aT(n_orb);
  std::vector<madness::vector_real_function_3d> bT(n_orb);
  for (long p = 0; p < n_orb; ++p) {
    aT[p].resize(m_a);
    bT[p].resize(m_b);
    for (long i = 0; i < m_a; ++i) aT[p][i] = a[i][p];
    for (long j = 0; j < m_b; ++j) bT[p][j] = b[j][p];
  }

  madness::Tensor<double> result(m_a, m_b);
  for (long p = 0; p < n_orb; ++p) {
    result += matrix_inner(world, aT[p], bT[p]);
    world.gop.fence();
  }
  return result;
}

// Forward declaration — `lowdin_orthonormalize` calls `transform`.
inline void transform(madness::World &world, response_space &X,
                      const madness::Tensor<double> &U);

/// Symmetric Löwdin orthonormalization of a bundle of M states.
///
/// Given X = [x₁, …, x_M] (each xᵢ a vecfuncT — packed α or α+β
/// concat for OpenShell), this produces
///
///     X' = X · S^{-1/2}      with S = ⟨X|X⟩
///
/// so that ⟨X'|X'⟩ = I exactly. Among orthonormalizations of X, this
/// one minimises ‖X' − X‖_F — equivalently, the unique orthonormal
/// basis closest to X. Unlike Gram-Schmidt, no slot is privileged:
/// all slots are perturbed symmetrically toward the orthonormal
/// limit. This is what we want at top-of-iter when slot identity
/// must survive across iters for KAIN's per-slot history.
///
/// Algorithm (S real symmetric ⇒ V·D·Vᵀ via syev):
///
///   1. S = ⟨X|X⟩;  symmetrize S := 0.5·(S + Sᵀ) (cheap, robust).
///   2. syev(S) → V·D·Vᵀ. Eigenvalues d_i in ascending order.
///   3. Build S^{-1/2} = V · diag(d_i^{-1/2}) · Vᵀ, with directions
///      where d_i ≤ d_floor (default 10·thresh) DROPPED (set to 0).
///      Catches the near-linear-dependent case: if two states have
///      collapsed onto each other, S has a small eigenvalue and
///      1/√d would amplify noise; dropping is the safer choice.
///   4. X' = X · S^{-1/2} via rs::transform.
///
/// Returns the number of directions dropped (0 in the well-conditioned
/// case). A non-zero return is a signal that the bundle has partially
/// collapsed; caller may want to inject random perturbations or stop.
///
/// `thresh < 0` (default) → use FunctionDefaults<3>::get_thresh().
///
/// `diag_level`: 0 = silent; 1 = print a one-line [LOWDIN] per call
/// with S-eigenvalues, condition number, n_dropped, max per-slot
/// movement ‖X'ᵢ − Xᵢ‖. The diagnostic block performs additional
/// collective work (snapshot copy + per-slot norm2 + per-slot sub)
/// — it must be entered or skipped UNIFORMLY across all ranks. The
/// caller is responsible for passing the same diag_level on every
/// rank (which is the case when it comes from a uniformly-set
/// print_level).
///
/// **Collective-call discipline.** All madness numerical operations
/// (`rs::inner`, `rs::transform`, `madness::copy`, `madness::sub`,
/// `madness::norm2`, `madness::scale`) are collective and must be
/// called by every rank in lock-step. They are placed at function
/// scope, never inside `if (rank == 0)` blocks. Only `printf`/`print`
/// are rank-gated. Tensor-only ops (`syev`, `inner` on Tensors,
/// `transpose`, element-wise loops on a replicated Tensor) are local
/// and run identically on every rank by virtue of their inputs being
/// replicated.
inline int
lowdin_orthonormalize(madness::World &world, response_space &X,
                      double thresh    = -1.0,
                      int diag_level   = 0) {
  if (thresh < 0.0)
    thresh = madness::FunctionDefaults<3>::get_thresh();

  const long M = static_cast<long>(X.size());
  if (M == 0) return 0;
  if (M == 1) {
    // Single-state shortcut: just normalize. inner+scale are collective.
    const double n2  = madness::inner(X[0], X[0]);
    const double nrm = std::sqrt(std::abs(n2));
    if (nrm > 1.0e-12) madness::scale(world, X[0], 1.0 / nrm);
    if (diag_level >= 1 && world.rank() == 0) {
      printf("[LOWDIN] M=1  norm=%.6e  (single-state normalize)\n", nrm);
      fflush(stdout);
    }
    return 0;
  }

  // ---- 0. Optional pre-snapshot (uniform across ranks) ------------------
  // Allocated only when diag_level >= 1 — but the BRANCH on diag_level
  // is on a value that is identical on every rank, so all ranks
  // either build X_before or don't. No collective is gated on rank.
  response_space X_before;
  if (diag_level >= 1) {
    X_before.resize(M);
    for (long i = 0; i < M; ++i) X_before[i] = madness::copy(world, X[i]);
  }

  // ---- 1. Overlap, symmetrized ------------------------------------------
  auto S = rs::inner(X, X);
  S = 0.5 * (S + madness::transpose(S));

  // ---- 2. Symmetric eigen-decomposition  S = V · diag(D) · Vᵀ ----------
  madness::Tensor<double> V;
  madness::Tensor<double> D;
  madness::syev(S, V, D);

  // ---- 3. S^{-1/2}, regularized -----------------------------------------
  const double d_floor = 10.0 * thresh;
  madness::Tensor<double> Dhalf_inv(M);
  int n_dropped = 0;
  for (long i = 0; i < M; ++i) {
    if (D(i) <= d_floor) {
      // includes any spuriously-negative eigenvalues from noise
      Dhalf_inv(i) = 0.0;
      ++n_dropped;
    } else {
      Dhalf_inv(i) = 1.0 / std::sqrt(D(i));
    }
  }

  // S^{-1/2}(j,i) = Σ_k V(j,k) · Dhalf_inv(k) · V(i,k)
  // Compute as V_scaled · Vᵀ where V_scaled(j,k) = V(j,k) · Dhalf_inv(k).
  madness::Tensor<double> V_scaled = madness::copy(V);
  for (long j = 0; j < M; ++j)
    for (long k = 0; k < M; ++k)
      V_scaled(j, k) *= Dhalf_inv(k);
  auto S_invsqrt = madness::inner(V_scaled, madness::transpose(V));

  // ---- 4. X' = X · S^{-1/2} via rs::transform --------------------------
  // rs::transform convention is X_new[j] = Σ_i U(i,j) · X_old[i], which
  // matches X' = X · M with M = S^{-1/2}.
  rs::transform(world, X, S_invsqrt);

  // ---- 5. Diagnostic (collective work uniform across ranks; printf rank-0) -
  if (diag_level >= 1) {
    // Per-slot movement: ‖X'ᵢ − Xᵢ‖. sub+norm2 are collective and run on
    // every rank for every i (no rank gate). We compute all M norms
    // first, then a single rank-0 printf consumes them.
    std::vector<double> movement(M, 0.0);
    for (long i = 0; i < M; ++i) {
      auto diff   = madness::sub(world, X[i], X_before[i]);
      movement[i] = madness::norm2(world, diff);
    }

    // Verify orthonormality post-Löwdin: <X'|X'> should be ≈ I.
    // This re-runs rs::inner on the rotated bundle — purely
    // diagnostic, costs one more m×m inner product per iter.
    auto S_after = rs::inner(X, X);
    double off_diag_max = 0.0;
    double diag_dev_max = 0.0;
    for (long i = 0; i < M; ++i) {
      diag_dev_max = std::max(diag_dev_max,
                              std::fabs(S_after(i, i) - 1.0));
      for (long j = 0; j < M; ++j)
        if (i != j)
          off_diag_max = std::max(off_diag_max,
                                  std::fabs(S_after(i, j)));
    }

    if (world.rank() == 0) {
      // Condition number from extreme eigenvalues (D is ascending).
      const long i_min = 0, i_max = M - 1;
      const double d_min = D(i_min);
      const double d_max = D(i_max);
      const double cond  = (d_min > 1.0e-30) ? d_max / d_min
                                              : std::numeric_limits<double>::infinity();
      double m_max = 0.0;
      for (double m : movement) m_max = std::max(m_max, m);

      printf("[LOWDIN]  S_evals=[");
      for (long i = 0; i < M; ++i) printf(" %+.3e", D(i));
      printf(" ]  cond=%.2e  d_floor=%.2e  dropped=%d"
             "  ||X'-X||_max=%.3e  |S'-I|_diag=%.2e  |S'-I|_off=%.2e\n",
             cond, d_floor, n_dropped, m_max,
             diag_dev_max, off_diag_max);
      fflush(stdout);
    }
  }

  return n_dropped;
}

/// Output of `diagonalize` — includes evals + U plus diagnostic
/// data the caller can log to track slot identity across iters:
///
///   omega_ascending — eigenvalues in ascending order, indexed by the
///                     ω-rank (NOT by output slot).
///   dominance_perm  — `dominance_perm[i] = j` means output slot i was
///                     filled by the ω-rank-j eigenvector during the
///                     diagonal-dominance swap loop. Identity perm
///                     means the sygvp ascending order matched the
///                     input slot ordering — slots are "well-aligned".
///   U_diag          — diag(U) AFTER all fixups, indexed by output
///                     slot. Values near ±1 mean the rotation kept slot
///                     identity; values near 0 mean the slot picked up
///                     character from another slot via the rotation.
struct DiagonalizeResult {
  madness::Tensor<double> omega;          ///< eigenvalues, output-slot indexed
  madness::Tensor<double> U;              ///< rotation, output-slot columns
  madness::Tensor<double> omega_ascending; ///< pre-dominance ω order
  std::vector<long>       dominance_perm; ///< i ↦ source ω-rank of slot i
  madness::Tensor<double> U_diag;         ///< diag(U) post-fixups, per slot
};

/// Solve A U = S U diag(omega). Returns evals + U + diagnostic data
/// (DiagonalizeResult).
///
/// Beyond `sygvp`, this performs the same identity-preserving fixups
/// that the legacy molresponse::get_fock_transformation does — without
/// them, per-slot response densities oscillate near degeneracies (H2's
/// ω₂ ≈ ω₃ case):
///
///   0. Symmetrize: A := 0.5·(A + Aᵀ), S := 0.5·(S + Sᵀ). Numerical
///      noise (and the asymmetry in <X|Λ> when Λ isn't computed
///      self-adjointly) can push A off-symmetric by O(thresh); sygvp
///      assumes symmetric. Mirrors ExcitedResponse.cpp:835,846.
///   1. SVD of S — drop columns whose singular value < 10·thresh_degen
///      (subspace reduction; rare in practice, included for robustness).
///   2. Generalized eigenproblem on the (possibly reduced) A, S.
///   3. Diagonal-dominance sort: repeatedly swap pairs of columns so the
///      diagonal of U holds the largest entries. This keeps slot s
///      tracking the eigenvalue closest to slot s's input axis, so
///      iter-to-iter slot identity is continuous.
///   4. Phase fix: flip the sign of any column with U(i,i) < 0 so the
///      diagonal is non-negative — fixes the canonical phase choice.
///   5. Cluster unmixing: within each near-degenerate eigenvalue
///      cluster (|ω_i − ω_j| < cluster_factor · thresh_degenerate ·
///      max(|ω_i|, 1)), replace the cluster's U block with the polar-
///      decomposition orthogonal factor of that block — picks the
///      rotation closest to identity within the degenerate subspace,
///      so eigenvectors don't rotate freely between iterations.
///   6. Transform U back to the original size if step 1 reduced it.
///
/// `thresh_degenerate < 0` (default) → use FunctionDefaults<3>::thresh.
/// `cluster_factor`: TDA convention is 100 (loose), Full/RPA is 10
///   (tighter — matches legacy molresponse Full path, see
///   ExcitedResponse.cpp:1073 / TDDFT.cc:3097). The caller decides.
inline DiagonalizeResult
diagonalize(const madness::Tensor<double> &A_in,
            const madness::Tensor<double> &S_in,
            double thresh_degenerate = -1.0,
            double cluster_factor    = 100.0) {
  using madness::Slice;
  using madness::_;

  if (thresh_degenerate < 0.0) {
    thresh_degenerate = madness::FunctionDefaults<3>::get_thresh();
  }

  madness::Tensor<double> A = madness::copy(A_in);
  madness::Tensor<double> S = madness::copy(S_in);

  // ---- 0. Symmetrize (cheap, harmless if already symmetric) ------------
  A = 0.5 * (A + madness::transpose(A));
  S = 0.5 * (S + madness::transpose(S));

  // ---- 1. SVD of overlap; identify singular columns --------------------
  madness::Tensor<double> l_vecs, s_vals, r_vecs;
  {
    auto S_copy = madness::copy(S);
    madness::svd(S_copy, l_vecs, s_vals, r_vecs);
  }

  std::size_t num_sv = 0;
  for (long i = 0; i < s_vals.dim(0); ++i) {
    if (s_vals(i) < 10.0 * thresh_degenerate) ++num_sv;
  }
  const std::size_t size_l = static_cast<std::size_t>(s_vals.dim(0));
  const std::size_t size_s = size_l - num_sv;

  // ---- 2a. Reduce subspace if needed ------------------------------------
  if (num_sv > 0) {
    S = madness::Tensor<double>(size_s, size_s);
    for (std::size_t i = 0; i < size_s; ++i) S(i, i) = s_vals(i);

    auto l_vecs_s = madness::copy(l_vecs(_, Slice(0, static_cast<long>(size_s) - 1)));
    madness::Tensor<double> work(size_l, size_s);
    madness::mxm(size_l, size_s, size_l, work.ptr(), A.ptr(), l_vecs_s.ptr());
    A = madness::Tensor<double>(size_s, size_s);
    auto l_vecs_t = madness::transpose(l_vecs);
    madness::mxm(size_s, size_s, size_l, A.ptr(), l_vecs_t.ptr(), work.ptr());
  }

  // ---- 2b. Generalized eigenproblem -------------------------------------
  madness::Tensor<double> U_small, evals;
  madness::sygvp(madness::World::get_default(), A, S, 1, U_small, evals);

  const long nmo = A.dim(0);

  // Snapshot ascending-order evals BEFORE the dominance swap so the
  // caller can compare "what sygvp returned" vs "what slots ended up
  // holding". omega_ascending(k) = k-th smallest ω.
  madness::Tensor<double> omega_ascending = madness::copy(evals);

  // Track the dominance permutation: perm[i] = which sygvp column
  // ended up in slot i after the swap loop. Identity perm means the
  // sygvp ascending order matched the input slot ordering.
  std::vector<long> perm(nmo);
  for (long i = 0; i < nmo; ++i) perm[i] = i;

  // ---- 3. Diagonal-dominance sort (preserve slot identity) -------------
  bool switched = true;
  while (switched) {
    switched = false;
    for (long i = 0; i < nmo; ++i) {
      for (long j = i + 1; j < nmo; ++j) {
        const double sold = U_small(i, i) * U_small(i, i)
                          + U_small(j, j) * U_small(j, j);
        const double snew = U_small(i, j) * U_small(i, j)
                          + U_small(j, i) * U_small(j, i);
        if (snew > sold) {
          auto tmp = madness::copy(U_small(_, i));
          U_small(_, i) = U_small(_, j);
          U_small(_, j) = tmp;
          std::swap(evals[i], evals[j]);
          std::swap(perm[i], perm[j]);
          switched = true;
        }
      }
    }
  }

  // ---- 4. Phase fix (column-wise sign so U(i,i) ≥ 0) -------------------
  for (long i = 0; i < nmo; ++i) {
    if (U_small(i, i) < 0.0) U_small(_, i).scale(-1.0);
  }

  // ---- 5. Polar-decomposition unmixing within degenerate clusters -----
  // cluster_factor=100 matches legacy TDA path; 10 matches legacy Full.
  long ilo = 0;
  while (ilo < nmo - 1) {
    long ihi = ilo;
    while (std::fabs(evals[ilo] - evals[ihi + 1]) <
           thresh_degenerate * cluster_factor *
               std::max(std::fabs(evals[ilo]), 1.0)) {
      ++ihi;
      if (ihi == nmo - 1) break;
    }
    const long nclus = ihi - ilo + 1;
    if (nclus > 1) {
      auto q = madness::copy(U_small(Slice(ilo, ihi), Slice(ilo, ihi)));
      madness::Tensor<double> W(nclus, nclus), VH(nclus, nclus), sigma(nclus);
      madness::svd(q, W, sigma, VH);
      q = madness::transpose(madness::inner(W, VH));
      U_small(_, Slice(ilo, ihi)) =
          madness::inner(U_small(_, Slice(ilo, ihi)), q);
    }
    ilo = ihi + 1;
  }

  // ---- 6. Transform U back to original size if reduced -----------------
  madness::Tensor<double> U;
  if (num_sv > 0) {
    madness::Tensor<double> temp_U(size_l, size_l);
    temp_U(Slice(0, static_cast<long>(size_s) - 1),
           Slice(0, static_cast<long>(size_s) - 1)) = madness::copy(U_small);
    for (std::size_t i = size_s; i < size_l; ++i) temp_U(i, i) = 1.0;
    U = madness::Tensor<double>(size_l, size_l);
    madness::mxm(size_l, size_l, size_l,
                 U.ptr(), l_vecs.ptr(), temp_U.ptr());
  } else {
    U = std::move(U_small);
  }

  // Slot-wise diag(U) post-fixups, for diagnostic. After subspace
  // expansion (step 6) any slots beyond size_s are pure identity.
  madness::Tensor<double> U_diag(U.dim(0));
  for (long i = 0; i < U.dim(0); ++i) U_diag(i) = U(i, i);

  DiagonalizeResult out;
  out.omega           = std::move(evals);
  out.U               = std::move(U);
  out.omega_ascending = std::move(omega_ascending);
  out.dominance_perm  = std::move(perm);
  out.U_diag          = std::move(U_diag);
  return out;
}

/// Rotate a response_space in place: X_new[j] = sum_i U(i, j) * X_old[i].
///
/// Implementation: transpose to (n_orb × n_states), dispatch to the
/// native `madness::transform` (vmra.h) once per orbital, untranspose.
/// The native routine uses tensor-train sparsity via FunctionImpl's
/// vtransform pathway — much faster than n×m manual gaxpys on the
/// nested response_space layout, especially as n_orb grows.
inline void transform(madness::World &world, response_space &X,
                      const madness::Tensor<double> &U) {
  const long n = static_cast<long>(X.size());
  if (n == 0) return;
  const long n_orb = static_cast<long>(X[0].size());
  MADNESS_CHECK(U.dim(0) == n && U.dim(1) == n);

  // Transposed view: XT[p][i] = X[i][p].
  // Function<double,3> has shared_ptr semantics — this aliases impls
  // rather than copying tree data.
  std::vector<madness::vector_real_function_3d> XT(n_orb);
  for (long p = 0; p < n_orb; ++p) {
    XT[p].resize(n);
    for (long i = 0; i < n; ++i) XT[p][i] = X[i][p];
  }

  // Native transform per orbital: YT[p] = XT[p] · U (sum over states).
  std::vector<madness::vector_real_function_3d> YT(n_orb);
  for (long p = 0; p < n_orb; ++p) {
    YT[p] = madness::transform(world, XT[p], U, /*fence=*/false);
  }
  world.gop.fence();

  // Untranspose into Y[i][p] = YT[p][i].
  response_space Y(n);
  for (long i = 0; i < n; ++i) Y[i].resize(n_orb);
  for (long p = 0; p < n_orb; ++p)
    for (long i = 0; i < n; ++i)
      Y[i][p] = YT[p][i];

  const double thresh = madness::FunctionDefaults<3>::get_thresh();
  for (auto &v : Y) truncate(world, v, thresh);
  X = std::move(Y);
}

} // namespace rs
} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_KERNELS_RESPONSE_SPACE_OPS_HPP
