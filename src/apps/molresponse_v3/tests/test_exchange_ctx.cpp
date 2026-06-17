// test_exchange_ctx — Inc 1 of the FD exchange tensor layer (docs 26/27/28).
//
// Proves the new tensor-layer θ (build_ctx + assemble_theta, which builds the
// shared Tx/g0 convolution tensors once and contracts them) matches the REFERENCE
// θ (compute_V0x − compute_E0x + compute_gamma) for Static ClosedShell. Both use
// the SAME gs/state/ρ/E0x — so the only difference is the exchange path (the
// Exchange operator vs the explicit Tx-build + dot-contraction). Same math, a
// different truncation/accumulation order → expect A/B-to-thresh (NOT bit-identical).
//
//   mpirun -np <N> ./test_exchange_ctx [n_occ]
// PASS iff ||θ_new − θ_ref|| / ||θ_ref|| < 1e-3 (truncation-level; report actual).

#include "../kernels/exchange_ctx.hpp"
#include "../kernels/static.hpp"          // Kernels<Static,ClosedShell> reference

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>         // CoulombOperatorPtr
#include <madness/mra/vmra.h>
#include <madness/world/MADworld.h>

#include <cmath>
#include <cstdlib>
#include <vector>

using namespace madness;
using namespace molresponse_v3;
using K      = Kernels<Static, ClosedShell>;
using StateX = ResponseStateX<ClosedShell>;

namespace {
class GaussFunctor : public FunctionFunctorInterface<double, 3> {
  coord_3d c_; double a_;
public:
  GaussFunctor(const coord_3d &c, double a) : c_(c), a_(a) {}
  double operator()(const coord_3d &r) const override {
    const double x = r[0]-c_[0], y = r[1]-c_[1], z = r[2]-c_[2];
    return std::exp(-a_ * (x*x + y*y + z*z));
  }
};

// Synthetic-but-valid Static ClosedShell GS. K0=null → the reference V0x uses the
// Exchange operator directly (apply_exchange), which is exactly what we A/B against.
ResponseGroundState make_gs(World &world, int n_occ) {
  ResponseGroundState g0;
  g0.amo.resize(n_occ);
  for (int i = 0; i < n_occ; ++i) {
    coord_3d c{{-3.0 + 1.3*i, 0.4*((i%3)-1), -0.2*(i%2)}};
    g0.amo[i] = real_factory_3d(world).functor(
        real_functor_3d(new GaussFunctor(c, 0.7 + 0.05*i)));
  }
  truncate(world, g0.amo);
  g0.V_local_alpha = real_factory_3d(world).functor(
      real_functor_3d(new GaussFunctor(coord_3d{{0.,0.,0.}}, 0.3)));
  g0.V_local_alpha.scale(-1.0); g0.V_local_alpha.truncate();
  Tensor<double> Fnd(n_occ, n_occ);                 // off-diagonal Fock (zero diag)
  for (int i=0;i<n_occ;++i) for (int j=0;j<n_occ;++j)
    Fnd(i,j) = (i==j) ? 0.0 : -0.05/(1.0+std::abs(i-j));
  g0.focka_no_diag = Fnd;
  g0.aeps = Tensor<double>(n_occ);
  for (int i=0;i<n_occ;++i) g0.aeps(i) = -(0.6 + 0.1*i);
  g0.c_xc = 1.0; g0.lo = 1.0e-4;
  g0.coulop = poperatorT(CoulombOperatorPtr(
      world, g0.lo, FunctionDefaults<3>::get_thresh()*1.0e-1));
  g0.Qa = QProjector<double,3>(g0.amo);
  g0.K0_alpha = nullptr;                            // → direct Exchange in the reference
  return g0;
}
} // namespace

int main(int argc, char **argv) {
  World &universe = initialize(argc, argv);
  startup(universe, argc, argv, true);

  // Scope guard: all MADNESS Functions/operators/QProjector below MUST be
  // destroyed BEFORE finalize(). gs.Qa (a QProjector) holds orbital copies; if
  // gs leaves main-scope AFTER finalize(), those Function dtors throw inside a
  // destructor -> std::terminate -> SIGABRT (the rc=134 seen right after PASS).
  bool ok = false;
  {
  FunctionDefaults<3>::set_k(8);
  FunctionDefaults<3>::set_thresh(1e-5);
  FunctionDefaults<3>::set_cubic_cell(-12.0, 12.0);
  const double vtol = FunctionDefaults<3>::get_thresh() * 0.1;

  const int n_occ = (argc > 1) ? std::atoi(argv[1]) : 4;

  auto gs = make_gs(universe, n_occ);

  // response state x (n_occ gaussians)
  StateX state;
  state.x_alpha.resize(n_occ);
  for (int p = 0; p < n_occ; ++p) {
    coord_3d c{{-2.5 + 1.1*p, 0.2*((p%3)-1), 0.1}};
    state.x_alpha[p] = real_factory_3d(universe).functor(
        real_functor_3d(new GaussFunctor(c, 0.8 + 0.04*p)));
  }
  truncate(universe, state.x_alpha);
  universe.gop.fence();

  auto rho = K::compute_density(universe, gs, state);

  // --- REFERENCE θ = V0x − E0x + γ  (op-level kernels) ---
  StateX theta_ref = K::compute_V0x(universe, gs, state);
  { auto E0x = K::compute_E0x(universe, gs, state);
    theta_ref.axpy(universe, -1.0, E0x); }
  { auto gam = K::compute_gamma(universe, gs, state, rho);
    theta_ref.axpy(universe, +1.0, gam); }
  theta_ref.truncate_all(universe, FunctionDefaults<3>::get_thresh());

  // --- NEW θ via the tensor layer (build Tx/g0 once, contract) ---
  auto g0  = exch::build_g0(universe, gs, vtol);
  auto ctx = exch::build_ctx_static_cs(universe, gs, state, rho, vtol);
  auto theta_new = exch::assemble_theta_static_cs(universe, gs, state, g0, ctx);
  universe.gop.fence();

  // --- compare (relative L2 over the n_occ orbital functions) ---
  double num = 0.0, den = 0.0;
  for (int p = 0; p < n_occ; ++p) {
    auto d = theta_ref.x_alpha[p] - theta_new.x_alpha[p];
    num += d.norm2() * d.norm2();
    den += theta_ref.x_alpha[p].norm2() * theta_ref.x_alpha[p].norm2();
  }
  const double rel = (den > 0.0) ? std::sqrt(num/den) : 0.0;
  ok = (rel < 1e-3);

  if (universe.rank() == 0) {
    print("\n=== FD exchange tensor-layer θ equivalence (Static CS, n_occ=", n_occ, ") ===");
    print("  θ_ref = compute_V0x − compute_E0x + compute_gamma  (Exchange operator)");
    print("  θ_new = assemble_theta(build_ctx) [shared Tx + cached g0, dot-contracted]");
    print("  ||θ_new − θ_ref|| / ||θ_ref|| =", rel, " (tol 1e-3; same math, diff truncation)");
    print("\nEXCHANGE_CTX_TEST:", ok ? "PASS" : "FAIL");
  }

  universe.gop.fence();
  }  // gs/state/theta_ref/theta_new/ctx/g0/rho destroyed here, before finalize()

  finalize();
  return ok ? 0 : 1;
}
