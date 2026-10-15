// ===========================================================================
// test_raman_spectrum.cpp — the Tier-B Raman formulas on synthetic tensors
// (kernels/raman_spectrum.hpp). Pure C++, no MPI, no World.
// ===========================================================================

#include "../kernels/raman_spectrum.hpp"

#include <cmath>
#include <cstdio>

using namespace molresponse_v3::raman;
using madness::Tensor;

namespace {
int failed = 0;
#define EXPECT(cond, label)                                                \
  do {                                                                     \
    if (cond) { std::printf("  [PASS]  %s\n", label); }                    \
    else      { std::printf("  [FAIL]  %s\n", label); ++failed; }          \
  } while (0)
bool near(double a, double b, double tol = 1e-12) { return std::abs(a - b) < tol; }
}  // namespace

int main() {
  std::printf("=== invariants: isotropic derivative ===\n");
  {
    Tensor<double> a(3L, 3L);
    for (int i = 0; i < 3; ++i) a(i, i) = 2.0;
    const auto r = invariants(a);
    EXPECT(near(r.alpha2, 4.0), "alpha'^2 = (tr/3)^2 = 4");
    EXPECT(near(r.beta2, 0.0), "beta'^2 = 0 for an isotropic tensor");
    EXPECT(near(r.activity, 180.0), "activity = 45 alpha'^2");
    EXPECT(near(r.depolarization, 0.0), "depolarization ratio 0 (fully polarized)");
  }
  std::printf("=== invariants: traceless anisotropic derivative ===\n");
  {
    // diag(1, -1, 0): mean 0; beta'^2 = 1/2 sum (3 a_ij^2 - a_ii a_jj)
    //   = 1/2 [ 3(1+1) - (1*1 + 1*(-1) + (-1)*1 + (-1)(-1)) ] = 1/2 [6 - 0] = 3
    Tensor<double> a(3L, 3L);
    a(0, 0) = 1.0; a(1, 1) = -1.0;
    const auto r = invariants(a);
    EXPECT(near(r.alpha2, 0.0), "alpha'^2 = 0 for a traceless tensor");
    EXPECT(near(r.beta2, 3.0), "beta'^2 = 3");
    EXPECT(near(r.activity, 12.0), "activity = 4 beta'^2");
    EXPECT(near(r.depolarization, 0.75), "depolarization ratio 3/4 (depolarized limit)");
  }
  std::printf("=== invariants: off-diagonal symmetric ===\n");
  {
    Tensor<double> a(3L, 3L);
    a(0, 1) = a(1, 0) = 1.5;   // beta'^2 = 1/2 * 3 * (1.5^2 * 2) = 6.75
    const auto r = invariants(a);
    EXPECT(near(r.beta2, 6.75), "beta'^2 from one off-diagonal pair");
  }
  std::printf("=== projection and mode extraction ===\n");
  {
    // 2 atoms (6 Cartesian coordinates), one mode along z of atom 0 minus z of atom 1.
    Tensor<double> D(9L, 6L);
    // d(alpha_zz)/dz_0 = 2, d(alpha_zz)/dz_1 = -2 ; row 3*2+2 = 8
    D(8, 2) = 2.0; D(8, 5) = -2.0;
    Tensor<double> L(6L, 1L);
    L(2, 0) = 0.5; L(5, 0) = -0.5;
    const auto Dq = project_on_modes(D, L);
    EXPECT(Dq.dim(0) == 9 && Dq.dim(1) == 1, "projected table is 9 x n_modes");
    EXPECT(near(Dq(8, 0), 2.0), "a_zz along the mode = 2*0.5 + (-2)(-0.5) = 2");
    const auto a = mode_tensor(Dq, 0);
    EXPECT(near(a(2, 2), 2.0) && near(a(0, 0), 0.0), "mode tensor row/col mapping (3A+B)");
    const auto r = invariants(a);
    // mean 2/3, alpha'^2 = 4/9; beta'^2 = 1/2 (3*4 - 4) = 4
    EXPECT(near(r.alpha2, 4.0 / 9.0) && near(r.beta2, 4.0), "invariants of the mode tensor");
  }
  std::printf("=== vibrational mode selection ===\n");
  {
    Tensor<double> f(6L);
    f(0L) = 1e-9; f(1L) = -2e-9; f(2L) = 0.0; f(3L) = 3e-9; f(4L) = 1e-9;
    f(5L) = 1427.94 / AU_TO_INVCM;   // ~6.506e-3 Eh
    const auto m = vibrational_modes(f, 1.0);
    EXPECT(m.size() == 1 && m[0] == 5, "one vibrational mode above 1 cm^-1");
    EXPECT(near(f(5L) * AU_TO_INVCM, 1427.94, 1e-9), "cm^-1 round trip");
  }
  std::printf("=== unit constant ===\n");
  {
    const double a0 = 0.529177210903, amu = 1822.888486209;
    EXPECT(near(AU_TO_A4_AMU, std::pow(a0, 4) * amu, 1e-3),
           "AU_TO_A4_AMU = a0^4 * (amu/m_e) to 1e-3");
  }
  std::printf("\n%s: %d failure(s)\n", failed == 0 ? "ALL PASS" : "FAILED", failed);
  return failed == 0 ? 0 : 1;
}
