#ifndef MOLRESPONSE_V3_ESSOLVERGUESS_HPP
#define MOLRESPONSE_V3_ESSOLVERGUESS_HPP

#include "GroundState.hpp"
#include "ResponseFunctions.hpp"
#include "ResponseKernel.hpp"   // for orthonormalize_bundle, ResponseType

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
    ESGaussianGuess(const coordT& origin, double exponent,
                    std::vector<int> ijk = std::vector<int>(NDIM))
        : origin(origin), exponent(exponent), ijk(std::move(ijk)) {}

    coordT origin;
    double exponent;
    std::vector<int> ijk;

    double operator()(const coordT& xyz) const {
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
inline real_function_3d build_atom_envelope(World& world,
                                             const Molecule& molecule,
                                             double exponent = 0.01) {
    real_function_3d envelope = real_factory_3d(world);
    for (const auto& atom : molecule.get_atoms()) {
        Vector<double, 3> origin{atom.x, atom.y, atom.z};
        real_function_3d g = real_factory_3d(world).functor(real_functor_3d(
            new ESGaussianGuess<3>(origin, exponent,
                                    std::vector<int>{0, 0, 0})));
        envelope = envelope + g;
    }
    return envelope;
}

/// In-place: add random noise (per leaf) to each function in the vector.
/// Mirrors molresponse_legacy/TDDFT.cc:3720-3750 (`add_randomness`).
inline void add_random_noise(World& world, vector_real_function_3d& v,
                              double magnitude = 1.0e3) {
    auto noise = [magnitude](const Key<3>& /*key*/, Tensor<double>& x) {
        Tensor<double> y(x.size());
        y.fillrandom();
        y.scale(magnitude);
        x = x + y;
    };
    for (auto& f : v) f.unaryop(noise);
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
inline std::vector<RealResponseState> make_initial_guess_tda_rhf(
    World& world,
    GroundState& gs,
    long num_roots,
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
        for (auto& f : X[s].x_alpha) f = envelope * f;
    }

    // Step 4: Q-project each state's x_alpha against ground orbitals
    const auto& Q = gs.Q_alpha();
    for (auto& state : X) {
        state.x_alpha = Q(state.x_alpha);
    }

    // Step 5: Gram-Schmidt across roots (TDA metric is identity)
    orthonormalize_bundle(world, ResponseType::TDA, X);

    return X;
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_ESSOLVERGUESS_HPP
