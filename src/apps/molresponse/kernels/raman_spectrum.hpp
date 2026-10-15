#ifndef MOLRESPONSE_V3_KERNELS_RAMAN_SPECTRUM_HPP
#define MOLRESPONSE_V3_KERNELS_RAMAN_SPECTRUM_HPP

// ===========================================================================
// raman_spectrum — the Tier-B (chemist-facing) step of vibrational Raman:
// Cartesian polarizability derivatives d(alpha_AB)/dQ over the 3N nuclear
// coordinates, projected on normal modes, reduced to the rotational
// invariants and intensities. Pure tensor algebra: no World, no I/O, so the
// formulas can be unit-tested on synthetic tensors.
//
// For one mode q with derivative tensor a_ij = d(alpha_ij)/dq (a.u.):
//   mean polarizability derivative  a' = tr(a)/3
//   alpha'^2 = a'^2
//   beta'^2  = 1/2 sum_ij ( 3 a_ij^2 - a_ii a_jj )     (anisotropy; symmetric a)
//   activity (linear polarization, 90 degrees)  S = 45 alpha'^2 + 4 beta'^2
//   depolarization ratio                        rho = 3 beta'^2 / S
// Normal-mode derivative: a_q = sum_k (d alpha / d x_k) L_kq, where L is the
// Cartesian displacement per unit mass-weighted normal coordinate
// (L = M^{-1/2} * eigenvector, Molecule::massweights()). a_q then carries
// bohr^4 * m_e^-1 ... in atomic units; multiplying by AU_TO_A4_AMU gives the
// customary Angstrom^4 / amu. These definitions and the constant are the ones
// the removed molresponse_v2 engine used (MolresponseLib.hpp compute_raman),
// so its archived LiH / H2O numbers remain comparable.
// ===========================================================================

#include <madness/tensor/tensor.h>

#include <cmath>
#include <vector>

namespace molresponse_v3::raman {

/// a.u. (bohr^4 / m_e) -> Angstrom^4 / amu: (0.529177210903)^4 * 1822.888486.
constexpr double AU_TO_A4_AMU = 142.9435756;
/// Hartree -> cm^-1 (same value as madness::constants::au2invcm).
constexpr double AU_TO_INVCM = 219474.63136314;

struct Invariants {
  double alpha2 = 0.0;   ///< alpha'^2 (a.u.)
  double beta2  = 0.0;   ///< beta'^2  (a.u.)
  double activity = 0.0; ///< 45 alpha'^2 + 4 beta'^2 (a.u.)
  double depolarization = 0.0; ///< 3 beta'^2 / activity (dimensionless)
};

/// Rotational invariants of one 3x3 polarizability-derivative tensor.
inline Invariants invariants(const madness::Tensor<double> &a) {
  Invariants r;
  double tr = 0.0;
  for (int i = 0; i < 3; ++i) tr += a(i, i);
  const double mean = tr / 3.0;
  r.alpha2 = mean * mean;
  double b2 = 0.0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      b2 += 0.5 * (3.0 * a(i, j) * a(i, j) - a(i, i) * a(j, j));
  r.beta2 = b2;
  r.activity = 45.0 * r.alpha2 + 4.0 * r.beta2;
  r.depolarization = (r.activity != 0.0) ? 3.0 * r.beta2 / r.activity : 0.0;
  return r;
}

/// Project Cartesian derivatives D (9 x 3N; row = 3*A + B for the dipole
/// axes A,B; column = 3*atom + axis) on the mode matrix L (3N x n_modes):
/// returns 9 x n_modes.
inline madness::Tensor<double>
project_on_modes(const madness::Tensor<double> &D,
                 const madness::Tensor<double> &L) {
  return madness::inner(D, L);
}

/// Column `q` of a projected (9 x n_modes) table as a 3x3 tensor a(A,B).
inline madness::Tensor<double>
mode_tensor(const madness::Tensor<double> &Dq, long q) {
  madness::Tensor<double> a(3L, 3L);
  for (int A = 0; A < 3; ++A)
    for (int B = 0; B < 3; ++B) a(A, B) = Dq(3 * A + B, q);
  return a;
}

/// Indices of the vibrational modes: |nu| above `cutoff_cm1` (the projected
/// translations/rotations sit at ~0 after compute_frequencies(project_tr)).
inline std::vector<long>
vibrational_modes(const madness::Tensor<double> &freq_au, double cutoff_cm1) {
  std::vector<long> m;
  for (long i = 0; i < freq_au.dim(0); ++i)
    if (std::abs(freq_au(i) * AU_TO_INVCM) > cutoff_cm1) m.push_back(i);
  return m;
}

} // namespace molresponse_v3::raman

#endif // MOLRESPONSE_V3_KERNELS_RAMAN_SPECTRUM_HPP
