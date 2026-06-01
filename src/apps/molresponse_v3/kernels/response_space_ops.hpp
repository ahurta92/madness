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

// ----------------------------------------------------------------------
// State-templated overloads. They flatten each State into the canonical
// component ordering (`flatten()` already does α (+β) concat per shell)
// then dispatch to the response_space versions above. The OpenShell
// bundle metric ⟨x_α|x_α⟩ + ⟨x_β|x_β⟩ and the slot-preserving rotation
// fall out of `flatten()` for free, so solvers no longer need their own
// pack/unpack scaffolding.
//
// Overload resolution: the non-template response_space variants win when
// the argument IS a response_space (`std::vector<vector_real_function_3d>`);
// the templates kick in only when the argument is a
// `std::vector<ResponseStateX<…>>` etc.
//
// Header-only — caller is responsible for having included the relevant
// `ResponseStateX/XY<…>` definitions before instantiation. No include
// of response_state.hpp here, to keep kernels/ free of a hard dependency
// on solvers/.
// ----------------------------------------------------------------------

template <typename State>
inline madness::Tensor<double>
inner(const std::vector<State> &a, const std::vector<State> &b) {
  response_space ra(a.size()), rb(b.size());
  for (std::size_t i = 0; i < a.size(); ++i) ra[i] = a[i].flatten();
  for (std::size_t j = 0; j < b.size(); ++j) rb[j] = b[j].flatten();
  return inner(ra, rb);
}

template <typename State>
inline void
transform(madness::World &world, std::vector<State> &X,
          const madness::Tensor<double> &U) {
  response_space rX(X.size());
  for (std::size_t i = 0; i < X.size(); ++i) rX[i] = X[i].flatten();
  transform(world, rX, U);
  for (std::size_t i = 0; i < X.size(); ++i) X[i].from_flat(rX[i]);
}

} // namespace rs
} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_KERNELS_RESPONSE_SPACE_OPS_HPP
