#pragma once

// The quantities people quote, computed from the tensors and values a run
// produces: pure functions of numbers, so each formula is written once and
// unit-tested (test_derived_quantities.cc). results.json carries what they
// return under each item's `derived` (DerivedResults.h); two-photon absorption
// calls them while it assembles. Nothing here multiplies by a convention
// factor: a quantity is linear or quadratic in the tensor as the code produces
// it.
#include <madness/constants.h>

#include <array>
#include <cmath>

namespace madness::derived {

using Vector3 = std::array<double, 3>;
/// t[i][j]: the component named by axes i, j (x = 0, y = 1, z = 2).
using Tensor2 = std::array<std::array<double, 3>, 3>;
/// b[i][j][k]: the hyperpolarizability component named `ijk`.
using Tensor3 = std::array<Tensor2, 3>;

inline double norm(const Vector3& v) { return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); }

// ---------------------------------------------------------------- ground state

/// A dipole moment in atomic units (e a0) as Debye.
inline double debye_from_atomic_units(double dipole) {
    return dipole * constants::atomic_unit_of_electric_dipole_moment / constants::Debye;
}

// -------------------------------------------------------------- polarizability

/// (xx + yy + zz) / 3.
inline double polarizability_isotropic(double xx, double yy, double zz) { return (xx + yy + zz) / 3.0; }

/// sqrt( 1/2 [ (xx-yy)^2 + (yy-zz)^2 + (zz-xx)^2 ] + 3 [ xy^2 + xz^2 + yz^2 ] ).
inline double polarizability_anisotropy(double xx, double yy, double zz, double xy, double xz, double yz) {
    const double diagonal = (xx - yy) * (xx - yy) + (yy - zz) * (yy - zz) + (zz - xx) * (zz - xx);
    return std::sqrt(0.5 * diagonal + 3.0 * (xy * xy + xz * xz + yz * yz));
}

// --------------------------------------------------------- hyperpolarizability

/// beta_i = (1/3) sum_j ( beta_ijj + beta_jij + beta_jji ).
inline Vector3 hyperpolarizability_vector(const Tensor3& b) {
    Vector3 v{};
    for (int i = 0; i < 3; ++i) {
        double sum = 0.0;
        for (int j = 0; j < 3; ++j) sum += b[i][j][j] + b[j][i][j] + b[j][j][i];
        v[i] = sum / 3.0;
    }
    return v;
}

/// (3/5) sum_i beta_i mu_i / |mu|. The caller decides that the dipole is not
/// zero; this divides by its norm.
inline double hyperpolarizability_parallel(const Vector3& beta_vector, const Vector3& dipole) {
    const double projection = beta_vector[0] * dipole[0] + beta_vector[1] * dipole[1] + beta_vector[2] * dipole[2];
    return 0.6 * projection / norm(dipole);
}

/// The two orientation averages of hyper-Rayleigh scattering.
struct HyperRayleighAverages {
    double zzz = 0.0; ///< <beta_ZZZ^2>
    double zxx = 0.0; ///< <beta_ZXX^2>
};

/// <beta_ZZZ^2> and <beta_ZXX^2>: the squares of the laboratory components
/// beta_ZZZ = sum R_Zi R_Zj R_Zk beta_ijk and beta_ZXX = sum R_Zi R_Xj R_Xk
/// beta_ijk, averaged over all orientations R of the molecule. This is the
/// definition of the averages in R. Bersohn, Y.-H. Pao and H. L. Frisch,
/// J. Chem. Phys. 45, 3184 (1966), and S. J. Cyvin, J. E. Rauch and
/// J. C. Decius, J. Chem. Phys. 43, 4083 (1965), with no symmetry assumed
/// among the indices of beta (no Kleinman symmetry).
///
/// The average is evaluated from the definition with a quadrature that is
/// exact for it, rather than from the expanded sums of products of components
/// those papers give, so there is no expansion to transcribe. In ZYZ Euler
/// angles the integrand is a polynomial of degree 6 in the entries of R. Its
/// angle-averaged part over alpha and gamma has Fourier orders |m| <= 6, so 8
/// equally spaced points in each are exact; what is left is a polynomial of
/// degree at most 6 in cos(beta), integrated exactly by 4-point Gauss-Legendre.
/// The unit test checks the result against a far finer grid.
inline HyperRayleighAverages hyper_rayleigh_averages(const Tensor3& b) {
    const double s = std::sqrt(6.0 / 5.0);
    const double inner = std::sqrt((3.0 - 2.0 * s) / 7.0), outer = std::sqrt((3.0 + 2.0 * s) / 7.0);
    const double w_inner = (18.0 + std::sqrt(30.0)) / 36.0, w_outer = (18.0 - std::sqrt(30.0)) / 36.0;
    const std::array<double, 4> nodes{-outer, -inner, inner, outer};
    const std::array<double, 4> weights{w_outer, w_inner, w_inner, w_outer};
    constexpr int n = 8;
    const double pi = constants::pi;
    HyperRayleighAverages sum;
    double total = 0.0;
    for (int a = 0; a < n; ++a)
        for (int c = 0; c < n; ++c)
            for (std::size_t g = 0; g < nodes.size(); ++g) {
                const double alpha = 2.0 * pi * a / n, gamma = 2.0 * pi * c / n;
                const double ca = std::cos(alpha), sa = std::sin(alpha), cg = std::cos(gamma), sg = std::sin(gamma);
                const double cb = nodes[g], sb = std::sqrt(1.0 - cb * cb);
                // Rows X and Z of R = Rz(alpha) Ry(beta) Rz(gamma).
                const Vector3 X{ca * cb * cg - sa * sg, -ca * cb * sg - sa * cg, ca * sb};
                const Vector3 Z{-sb * cg, sb * sg, cb};
                double zzz = 0.0, zxx = 0.0;
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                        for (int k = 0; k < 3; ++k) {
                            zzz += Z[i] * Z[j] * Z[k] * b[i][j][k];
                            zxx += Z[i] * X[j] * X[k] * b[i][j][k];
                        }
                sum.zzz += weights[g] * zzz * zzz;
                sum.zxx += weights[g] * zxx * zxx;
                total += weights[g];
            }
    sum.zzz /= total;
    sum.zxx /= total;
    return sum;
}

/// sqrt( <beta_ZZZ^2> + <beta_ZXX^2> ).
inline double hyper_rayleigh(const HyperRayleighAverages& averages) { return std::sqrt(averages.zzz + averages.zxx); }

/// Below this hrs (atomic units) the tensor is numerical noise, and the
/// depolarization ratio, a ratio of two noise numbers, is not reported.
inline constexpr double hyper_rayleigh_floor = 1e-8;

/// <beta_ZZZ^2> / <beta_ZXX^2>.
inline double hyper_rayleigh_depolarization(const HyperRayleighAverages& averages) {
    return averages.zzz / averages.zxx;
}

// ------------------------------------------------------------- excited states

/// The hartree-to-eV factor of excitation energies: CODATA, from constants.h.
inline constexpr double excitation_hartree_to_ev = constants::hartree_electron_volt_relationship;

/// An excitation energy in hartree as eV.
inline double excitation_energy_ev(double energy) { return energy * excitation_hartree_to_ev; }

/// The vacuum wavelength in nm of a photon of `energy` hartree: 1e7 / (E in cm^-1).
inline double wavelength_nm(double energy) { return 1e7 / (energy * constants::au2invcm); }

// ------------------------------------------------------ two-photon absorption

/// The orientation-averaged two-photon strengths of a transition tensor S
/// (P. R. Monson and W. M. McClain, J. Chem. Phys. 53, 29 (1970)), as DALTON's
/// QRSMO reports them (rspvec.F:2803-2821):
///   Df = (1/30) sum_ij S_ii S_jj,  Dg = (1/30) sum_ij S_ij S_ij,
///   D_linear = 2 Df + 4 Dg,  D_circular = -2 Df + 6 Dg,
///   R = (-Df + 3 Dg) / (Df + 2 Dg) (0 when the denominator is 0).
/// D = F Df + G Dg + H Dh with Dh = (1/30) sum_ij S_ij S_ji, which equals Dg
/// for the symmetric S the code builds, so H is folded into G.
struct TwoPhotonStrengths {
    double Df = 0.0, Dg = 0.0, D_linear = 0.0, D_circular = 0.0, R = 0.0;
};

inline TwoPhotonStrengths two_photon_strengths(const Tensor2& S) {
    double df = 0.0, dg = 0.0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            df += S[i][i] * S[j][j];
            dg += S[i][j] * S[i][j];
        }
    TwoPhotonStrengths o;
    o.Df = df / 30.0;
    o.Dg = dg / 30.0;
    o.D_linear = 2.0 * o.Df + 4.0 * o.Dg;
    o.D_circular = -2.0 * o.Df + 6.0 * o.Dg;
    const double den = o.Df + 2.0 * o.Dg;
    o.R = (den != 0.0) ? (-o.Df + 3.0 * o.Dg) / den : 0.0;
    return o;
}

/// The two-photon cross section's conversion from a.u. to GM, DALTON's
/// AU_TO_GM (rspvec.F:1795-6): 8 pi^2 alpha a0[pm]^5 / (c[cm/s] Gamma[au]),
/// with Gamma = 0.1 eV taken as a full width at half maximum, about 2.170.
/// In sigma = N pi^2 alpha a0^5 omega^2 delta g / c this is N = 8 and a
/// lineshape value g = 1/Gamma at resonance (not a normalized Lorentzian,
/// whose peak for that full width would be 2 / (pi Gamma)). The constants are
/// DALTON's, kept as they are; they differ from constants.h in the last digits.
inline double two_photon_au_to_gm() {
    constexpr double PI = 3.14159265358979324;
    constexpr double ALPHA = 7.2973525693e-3;  // fine-structure constant
    constexpr double A0_PM = 52.9177210903;    // Bohr radius in pm
    constexpr double C_CMS = 2.99792458e10;    // speed of light, cm/s
    constexpr double FWHM_AU = 0.0036749326;   // 0.1 eV in a.u.
    return 8.0 * PI * PI * ALPHA * std::pow(A0_PM, 5.0) / (C_CMS * FWHM_AU);
}

/// The two-photon cross section in GM of `strength` (a.u.) for a transition
/// of `excitation_energy` hartree: strength (E/2)^2 AU_TO_GM, both photons at
/// half the excitation energy.
inline double two_photon_cross_section_gm(double strength, double excitation_energy) {
    const double ph2 = (0.5 * excitation_energy) * (0.5 * excitation_energy);
    return strength * ph2 * two_photon_au_to_gm();
}

} // namespace madness::derived
