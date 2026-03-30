#ifndef MOLRESPONSE_V3_PERTURBATIONS_HPP
#define MOLRESPONSE_V3_PERTURBATIONS_HPP

#include "GroundState.hpp"
#include <madness/chem/SCF.h>
#include <madness/chem/SCFOperators.h>
#include <madness/chem/molecular_functors.h>
#include <madness/mra/mra.h>

namespace molresponse_v3 {

using namespace madness;

/// Construct the dipole perturbation right-hand side for one direction.
///
/// Returns Q * (mu_dir * phi_i) for each occupied orbital phi_i.
/// Uses MomentFunctor from SCF.h for the coordinate function.
///
/// @param world  MADNESS world
/// @param gs     ground state (provides orbitals and Q projector)
/// @param axis   0=x, 1=y, 2=z
/// @return       Q-projected perturbation vector (vecfuncT)
inline vector_real_function_3d
dipole_perturbation(World& world, const GroundState& gs, int axis) {
    MADNESS_CHECK(axis >= 0 && axis <= 2);

    std::vector<int> dir(3, 0);
    dir[axis] = 1;
    real_function_3d mu = real_factory_3d(world).functor(
        real_functor_3d(new MomentFunctor(dir)));
    mu.truncate(FunctionDefaults<3>::get_thresh());

    auto rhs = mul(world, mu, gs.orbitals_alpha());
    rhs = gs.Q()(rhs);
    truncate(world, rhs, FunctionDefaults<3>::get_thresh());
    return rhs;
}

/// Construct the dipole perturbation for beta orbitals (open-shell).
inline vector_real_function_3d
dipole_perturbation_beta(World& world, const GroundState& gs, int axis) {
    MADNESS_CHECK(axis >= 0 && axis <= 2);
    MADNESS_CHECK(!gs.is_spin_restricted());

    std::vector<int> dir(3, 0);
    dir[axis] = 1;
    real_function_3d mu = real_factory_3d(world).functor(
        real_functor_3d(new MomentFunctor(dir)));
    mu.truncate(FunctionDefaults<3>::get_thresh());

    // For unrestricted, beta Q projector would be built from beta orbitals.
    // For now, apply the same Q (alpha) — proper beta Q is future work.
    auto rhs = mul(world, mu, gs.orbitals_beta());
    rhs = gs.Q()(rhs);
    truncate(world, rhs, FunctionDefaults<3>::get_thresh());
    return rhs;
}

/// Construct the nuclear displacement perturbation right-hand side.
///
/// Uses DNuclear from SCFOperators.h which handles nuclear correlation
/// factors and core potentials correctly.
///
/// @param world  MADNESS world
/// @param gs     ground state
/// @param atom   atom index
/// @param axis   0=x, 1=y, 2=z
/// @return       Q-projected perturbation vector
inline vector_real_function_3d
nuclear_perturbation(World& world, const GroundState& gs, int atom, int axis) {
    DNuclear<double, 3> dV(world, &gs.scf(), atom, axis);
    auto rhs = dV(gs.orbitals_alpha());
    rhs = gs.Q()(rhs);
    truncate(world, rhs, FunctionDefaults<3>::get_thresh());
    return rhs;
}

/// Construct the magnetic (Lz) perturbation right-hand side.
///
/// Uses Lz from SCFOperators.h. Returns complex functions since
/// Lz = -i(x d/dy - y d/dx) maps real orbitals to complex.
///
/// NOTE: The Lz<T,NDIM> template in SCFOperators.h has a known issue
/// when T=double — its operator() returns real but internal code expects
/// complex. Must be instantiated as Lz<std::complex<double>, 3> with
/// complex input orbitals, or the SCFOperators template needs fixing.
/// For now, this is left as future work for magnetic response.
///
// inline std::vector<Function<std::complex<double>, 3>>
// magnetic_perturbation_lz(World& world, const GroundState& gs);

/// Construct the raw dipole operator function (without applying to orbitals).
/// Useful for property computation: alpha_ij = <x_i | mu | phi_j>
inline real_function_3d
dipole_operator(World& world, int axis) {
    MADNESS_CHECK(axis >= 0 && axis <= 2);
    std::vector<int> dir(3, 0);
    dir[axis] = 1;
    real_function_3d mu = real_factory_3d(world).functor(
        real_functor_3d(new MomentFunctor(dir)));
    mu.truncate(FunctionDefaults<3>::get_thresh());
    return mu;
}

/// Construct all three dipole operators (x, y, z).
inline std::array<real_function_3d, 3>
dipole_operators(World& world) {
    return {dipole_operator(world, 0),
            dipole_operator(world, 1),
            dipole_operator(world, 2)};
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_PERTURBATIONS_HPP
