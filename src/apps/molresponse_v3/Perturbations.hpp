#ifndef MOLRESPONSE_V3_PERTURBATIONS_HPP
#define MOLRESPONSE_V3_PERTURBATIONS_HPP

#include "GroundState.hpp"
#include <madness/chem/SCF.h>
#include <madness/chem/SCFOperators.h>
#include <madness/chem/molecular_functors.h>
#include <madness/mra/mra.h>

namespace molresponse_v3 {

using namespace madness;

/// Stable identity for an FD perturbation source (doc 13). `description()`
/// is the canonical string used in archive filenames and as the
/// `<pert>` key in response_metadata.json — it must round-trip across runs
/// and contain only filename-safe characters.
///
/// This is identity only; the existing dipole_perturbation / magnetic /
/// nuclear free functions below still build the RHS. Callers dispatch on
/// `kind` and pass `axis`/`atom` to whichever builder applies.
struct Perturbation {
    enum class Kind { Dipole, NuclearDisplacement, Magnetic };

    Kind kind = Kind::Dipole;
    int  axis = 0;        // 0=x, 1=y, 2=z
    int  atom = -1;       // for NuclearDisplacement; -1 = N/A

    /// Canonical filename/metadata-key fragment:
    ///   dipole_x   magnetic_y   nuc_3_z
    std::string description() const {
        static const char *ax[] = {"x", "y", "z"};
        const char *a = (axis >= 0 && axis <= 2) ? ax[axis] : "?";
        switch (kind) {
            case Kind::Dipole:               return std::string("dipole_") + a;
            case Kind::Magnetic:             return std::string("magnetic_") + a;
            case Kind::NuclearDisplacement:
                return std::string("nuc_") + std::to_string(atom) + "_" + a;
        }
        return "unknown";
    }

    static Perturbation dipole(int axis)            { return {Kind::Dipole,   axis, -1}; }
    static Perturbation magnetic(int axis)          { return {Kind::Magnetic, axis, -1}; }
    static Perturbation nuclear(int atom, int axis) { return {Kind::NuclearDisplacement, axis, atom}; }
};

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

    auto rhs = mul(world, mu, gs.orbitals_beta());
    rhs = gs.Q_beta()(rhs);
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
