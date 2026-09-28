// =========================================================================
// The relative FD convergence gate (review finding C1).
//
// The FD gate used to compare absolute norms: ||dx|| < 5 dconv and
// ||drho|| < 5 dconv, whatever the size of the response. The relative accuracy
// it imposed therefore scaled as 1/||x||: a Li leg (alpha ~ 164) or a nuclear
// leg (~1e3 a dipole leg) was held ~1000x tighter than a He dipole leg, and a
// leg whose first step already came in under 5 dconv was accepted at
// iteration 1 with no density check. Now
//     ||dx||   <= atol   + rtol ||x||
//     ||drho|| <= atol_r + rtol ||rho1||
// with atol = today's 5 max(thresh, dconv) and rtol = 5 dconv, a leg needs
// two iterations before it can pass, and the explosion guard and the
// step-restriction cap are relative.
//
// Pure C++; no World needed.
// =========================================================================

#include "../solvers/convergence_policy.hpp"

#include <cstdio>

using namespace molresponse_v3;

namespace {

int failed = 0;

void expect(bool cond, const char *label) {
  std::printf("  [%s]  %s\n", cond ? "PASS" : "FAIL", label);
  if (!cond) ++failed;
}

} // namespace

int main() {
  ConvergencePolicy p;
  p.dconv_user = 1e-4;
  const auto t = p.effective_for_thresh(1e-4);   // dconv = 1e-4

  std::printf("=== targets ===\n");
  expect(t.bsh_residual == 5e-4 && t.density_residual == 5e-4,
         "atol is today's absolute gate, 5 dconv");
  expect(t.rtol == 5e-4, "rtol = 5 dconv");
  expect(fd_gate(t.bsh_residual, t.rtol, 0.0) == 5e-4,
         "a zero-norm leg keeps exactly the absolute gate");
  expect(fd_gate(t.bsh_residual, t.rtol, 1000.0) > 0.5,
         "a leg with ||x|| = 1e3 gets a gate of atol + 0.5");

  std::printf("=== the absolute gate's two failure modes ===\n");
  // Large leg: ||x|| = 1e3, stopped at ||dx|| = 4e-3 (relative 4e-6).
  expect(!(4e-3 < t.bsh_residual),
         "absolute gate: a 1e3-norm leg at relative 4e-6 fails (the old behaviour)");
  expect(fd_leg_within(4e-3, 1000.0, 4e-3, 1000.0, t),
         "relative gate: the same leg passes");
  // Small leg at iteration 1: ||dx|| under atol, density never measured.
  expect(!fd_leg_converged(1, 1e-4, 0.1, 0.0, 0.1, t),
         "iteration 1 never passes, even under the gate (drho not yet measured)");
  expect(fd_leg_converged(2, 1e-4, 0.1, 1e-4, 0.1, t),
         "iteration 2 passes when both gates hold");

  std::printf("=== what must still fail ===\n");
  expect(!fd_leg_converged(3, 8e-4, 0.5, 1e-5, 0.5, t),
         "a unit-scale leg just over atol + rtol||x|| fails on dx");
  expect(!fd_leg_converged(3, 1e-5, 0.5, 8e-4, 0.5, t),
         "... and on drho");
  expect(!fd_leg_converged(3, 0.6, 1000.0, 1e-5, 1000.0, t),
         "a large leg at relative 6e-4 (> rtol) fails");
  expect(fd_leg_converged(3, 0.2, 1000.0, 1e-5, 1000.0, t),
         "... at relative 2e-4 (the Li floor) it passes");

  std::printf("=== explosion guard ===\n");
  expect(!fd_exploded(500.0, 1.0, 1e3), "unit leg: 500 is not an explosion");
  expect(fd_exploded(2e3, 1.0, 1e3), "unit leg: 2e3 is an explosion (as before)");
  expect(!fd_exploded(2e3, 10.0, 1e3), "a 10-norm leg at 2e3 is not (relative guard)");
  expect(fd_exploded(2e4, 10.0, 1e3), "... at 2e4 it is");
  expect(fd_exploded(2e3, 0.1, 1e3), "a small leg keeps the absolute floor of the guard");

  std::printf("=== step-restriction cap ===\n");
  expect(step_cap(0.5, 0.0) == 0.5 && step_cap(0.5, 1.0) == 0.5,
         "a normalized root keeps the absolute cap");
  expect(step_cap(0.5, 36.0) == 18.0, "a 36-norm leg (Li) gets 0.5 * 36");

  std::printf("\n%s: %d failure(s)\n", failed == 0 ? "ALL PASS" : "FAILED",
              failed);
  return failed == 0 ? 0 : 1;
}
