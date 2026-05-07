#ifndef MOLRESPONSE_V3_ESSOLVERGUESS_HPP
#define MOLRESPONSE_V3_ESSOLVERGUESS_HPP

#include "GroundState.hpp"
#include "ResponseFunctions.hpp"
#include "ResponseKernel.hpp" // for orthonormalize_bundle, ResponseType

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>

#include <vector>

namespace molresponse_v3 {

using namespace madness;

// -------------------------------------------------------------------------
// Simple Gaussian functor used to build a localized envelope around atoms.
// Independent copy of the same form used in molresponse_legacy/TDDFT.h:66
// and molresponse/ExcitedResponse.cpp:122 — kept local so this header is
// self-contained.
// -------------------------------------------------------------------------
template <std::size_t NDIM>
class ESGaussianGuess : public FunctionFunctorInterface<double, NDIM> {
  typedef Vector<double, NDIM> coordT;

public:
  ESGaussianGuess(const coordT &origin, double exponent,
                  std::vector<int> ijk = std::vector<int>(NDIM))
      : origin(origin), exponent(exponent), ijk(std::move(ijk)) {}

  coordT origin;
  double exponent;
  std::vector<int> ijk;

  double operator()(const coordT &xyz) const {
    double r2 = 0.0, prefac = 1.0;
    for (std::size_t i = 0; i < NDIM; i++) {
      const double d = xyz[i] - origin[i];
      r2 += d * d;
      prefac *= std::pow(xyz[i], ijk[i]);
    }
    return prefac * std::exp(-exponent * r2);
  }
};

/// Build a sum-of-Gaussians envelope localized at every atom of the
/// molecule. Width controlled by `exponent` (smaller = more diffuse).
inline real_function_3d build_atom_envelope(World &world,
                                            const Molecule &molecule,
                                            double exponent = 0.01) {
  real_function_3d envelope = real_factory_3d(world);
  for (const auto &atom : molecule.get_atoms()) {
    Vector<double, 3> origin{atom.x, atom.y, atom.z};
    real_function_3d g = real_factory_3d(world).functor(real_functor_3d(
        new ESGaussianGuess<3>(origin, exponent, std::vector<int>{0, 0, 0})));
    envelope = envelope + g;
  }
  return envelope;
}

/// In-place: add random noise (per leaf) to each function in the vector.
/// Mirrors molresponse_legacy/TDDFT.cc:3720-3750 (`add_randomness`).
inline void add_random_noise(World &world, vector_real_function_3d &v,
                             double magnitude = 1.0e3) {
  auto noise = [magnitude](const Key<3> & /*key*/, Tensor<double> &x) {
    Tensor<double> y(x.size());
    y.fillrandom();
    y.scale(magnitude);
    x = x + y;
  };
  for (auto &f : v)
    f.unaryop(noise);
}

/// Initial guess for **TDA closed-shell (RHF)** excited states.
///
/// Algorithm (mirrors `molresponse_legacy/TDDFT.cc::create_random_guess` /
/// `molresponse/ExcitedResponse::make_random_trial`):
///   1. Create N×n_occ zero functions.
///   2. Add per-leaf random noise.
///   3. Multiply by sum-of-Gaussians envelope localized at each atom.
///   4. Q-project to remove occupied-space components.
///   5. Orthonormalize across the N roots.
inline std::vector<RealResponseState>
make_initial_guess_tda_rhf(World &world, GroundState &gs, long num_roots,
                           double envelope_exponent = 0.01,
                           double random_magnitude = 1.0e3) {

  MADNESS_CHECK(gs.is_spin_restricted());
  const long n_alpha = gs.num_alpha();

  // Atom-localized Gaussian envelope (one shared envelope across roots —
  // randomness comes from add_random_noise).
  auto envelope = build_atom_envelope(world, gs.molecule(), envelope_exponent);

  std::vector<RealResponseState> X(num_roots);
  for (long s = 0; s < num_roots; s++) {
    X[s] = RealResponseState::allocate(world, n_alpha, /*n_beta=*/0,
                                       /*include_y=*/false);
    // Step 1+2: noise into each x_alpha[i]
    add_random_noise(world, X[s].x_alpha, random_magnitude);
    // Step 3: localize at atoms
    for (auto &f : X[s].x_alpha)
      f = envelope * f;
  }

  // Step 4: Q-project each state's x_alpha against ground orbitals
  const auto &Q = gs.Q_alpha();
  for (auto &state : X) {
    state.x_alpha = Q(state.x_alpha);
  }

  // Step 5: Gram-Schmidt across roots (TDA metric is identity)
  orthonormalize_bundle(world, ResponseType::TDA, X);

  return X;
}

// Cartesian moment functor: r → x^i · y^j · z^k.
// Direct port of the local copy in molresponse_legacy/TDDFT.cc (the
// comment there says SCF.cc's version "wasn't linking right" so it was
// duplicated). Kept self-contained inside this header.
class BS_MomentFunctor : public FunctionFunctorInterface<double, 3> {
private:
  const int i, j, k;

public:
  BS_MomentFunctor(int i, int j, int k) : i(i), j(j), k(k) {}
  explicit BS_MomentFunctor(const std::vector<int> &x)
      : i(x[0]), j(x[1]), k(x[2]) {}
  double operator()(const Vector<double, 3> &r) const override {
    double xi = 1.0, yj = 1.0, zk = 1.0;
    for (int p = 0; p < i; ++p) xi *= r[0];
    for (int p = 0; p < j; ++p) yj *= r[1];
    for (int p = 0; p < k; ++p) zk *= r[2];
    return xi * yj * zk;
  }
};

inline double kronecker(std::size_t l, std::size_t n) {
  return (l == n) ? 1.0 : 0.0;
}

/// Solid spherical harmonics up to angular momentum L = `n`, keyed by {l, m}.
/// Direct port of `TDDFT::solid_harmonics` (legacy TDDFT.cc:413-478).
/// `inline` so it can live in a header without ODR violations.
inline std::map<std::vector<int>, real_function_3d>
solid_harmonics(World &world, int n) {
  std::map<std::vector<int>, real_function_3d> result;

  // Basic x, y, z, constant, and zero.
  real_function_3d x = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{1, 0, 0})));
  real_function_3d y = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{0, 1, 0})));
  real_function_3d z = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{0, 0, 1})));
  real_function_3d c = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{0, 0, 0})));
  real_function_3d zero = real_factory_3d(world);

  // Seed (assume n >= 1).
  result[std::vector<int>{0, 0}] = copy(c);
  result[std::vector<int>{0, -1}] = zero;
  result[std::vector<int>{0, 1}] = zero;
  result[std::vector<int>{-1, 0}] = zero;

  // Generate the solid harmonics recursively (verbatim from legacy).
  for (int l = 0; l < n; l++) {
    result[std::vector<int>{l + 1, l + 1}] =
        sqrt(pow(2, kronecker(l, 0) * (2 * l) / (2 * l + 1))) *
        (x * result[std::vector<int>{l, l}] -
         (1 - kronecker(l, 0) * y * result[std::vector<int>{l, -l}]));
    result[std::vector<int>{l + 1, -l - 1}] =
        sqrt(pow(2, kronecker(l, 0) * (2 * l) / (2 * l + 1))) *
        (y * result[std::vector<int>{l, l}] +
         (1 - kronecker(l, 0) * x * result[std::vector<int>{l, -l}]));

    // Recursion needs out-of-range neighbors as zero placeholders.
    result[std::vector<int>{l + 1, l + 2}] = zero;
    result[std::vector<int>{l + 1, -l - 2}] = zero;

    for (int m = -l; m < l + 1; m++) {
      result[std::vector<int>{l + 1, m}] =
          1.0 / std::sqrt((l + m + 1) * (l - m + 1)) *
          ((2 * l + 1) * z * result[std::vector<int>{l, m}] -
           sqrt((l + m) * (l - m)) * (x * x + y * y + z * z) *
               result[std::vector<int>{l - 1, m}]);
    }
  }

  // Drop zeros and the constant {0,0}.
  for (auto it = result.begin(); it != result.end();) {
    if (it->second.norm2() == 0)
      it = result.erase(it);
    else
      ++it;
  }
  result.erase(std::vector<int>{0, 0});

  return result;
}

/// Initial guess for **TDA closed-shell (RHF)** excited states using
/// solid-harmonic-times-orbital trial functions.
///
/// Mirrors `TDDFT::create_trial_functions` (legacy TDDFT.cc:535-595):
///   1. Build solid harmonics up to L = max(2, ceil(sqrt(num_roots/n_occ) - 1))
///      so that solids.size() * n_occ >= num_roots.
///   2. For each (orbital, harmonic) pair, build a *single-excitation* trial
///      state: only one orbital position is non-zero and equals
///      `harmonic * orbital[mirror_index]`. This yields a richer subspace
///      than putting the same harmonic on every orbital.
///   3. Q-project each state to remove ground-state components.
///   4. Truncate.
///   5. Gram-Schmidt across the bundle.
///
/// Returns `num_roots` states (or fewer if not enough harmonic-orbital
/// pairs are available, which only happens for very large num_roots).
inline std::vector<RealResponseState>
create_solid_harmonics_guess(World &world, GroundState &gs, long num_roots) {
  MADNESS_CHECK(gs.is_spin_restricted());
  const long n_occ = gs.num_alpha();
  MADNESS_CHECK(n_occ >= 1);

  // Same heuristic as legacy create_trial_functions:
  // ensure solids.size() * n_occ > num_roots.
  const int max_l = static_cast<int>(std::max(
      2.0,
      std::ceil(std::sqrt(static_cast<double>(num_roots) /
                          static_cast<double>(n_occ)) -
                1.0)));
  auto solids = solid_harmonics(world, max_l);
  if (world.rank() == 0) {
    print("ES guess: created ", solids.size(),
          " solid harmonics (L<=", max_l, ") for ", num_roots,
          " trial states across ", n_occ, " occupied orbitals");
  }

  const auto &orbitals = gs.orbitals_alpha();
  std::vector<RealResponseState> X;
  X.reserve(num_roots);

  long count = 0;
  for (long i = 0; i < n_occ && count < num_roots; ++i) {
    for (const auto &kv : solids) {
      if (count >= num_roots) break;
      const real_function_3d &harmonic = kv.second;
      // Single-excitation trial: one occupied position non-zero.
      RealResponseState state = RealResponseState::allocate(
          world, n_occ, /*n_beta=*/0, /*include_y=*/false);
      const long pos = count % n_occ;
      const long orb_idx = n_occ - (count % n_occ) - 1;
      state.x_alpha[pos] = harmonic * orbitals[orb_idx];
      X.push_back(std::move(state));
      ++count;
    }
  }

  // Q-project each state.
  const auto &Q = gs.Q_alpha();
  for (auto &state : X) {
    state.x_alpha = Q(state.x_alpha);
  }

  // Truncate.
  const double thresh = FunctionDefaults<3>::get_thresh();
  for (auto &state : X) {
    truncate(world, state.x_alpha, thresh);
  }

  // Gram-Schmidt across roots (TDA metric is identity).
  orthonormalize_bundle(world, ResponseType::TDA, X);

  return X;
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_ESSOLVERGUESS_HPP
