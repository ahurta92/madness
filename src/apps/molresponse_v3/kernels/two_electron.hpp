#ifndef MOLRESPONSE_V3_KERNELS_TWO_ELECTRON_HPP
#define MOLRESPONSE_V3_KERNELS_TWO_ELECTRON_HPP

// ===========================================================================
// Shared two-electron apply core (doc 15 / kernel unification).
//
// Every Kernels<Type,Shell> implements the SAME operation per output channel:
//
//     out_channel = Q( J[rho]*apply_to  -  c_xc * Sum_pairs K(bra, ket)(apply_to) )
//
// The ONLY things that differ between (Type, Shell) are (a) the density factor
// + which response channels build rho (handled by each kernel's two-state
// `compute_density(S1, S2)`), and (b) the list of exchange (bra, ket) pairs per
// output channel (the per-type pairing). This header factors out the common
// boilerplate -- the Coulomb multiply, the c_xc*K accumulation, projection and
// truncation -- so each kernel's `apply_g` / `compute_gamma` is just "build J,
// then name the exchange pairs."
//
// The unified picture: gamma involves THREE response states. S1, S2 build the
// density AND supply the exchange bra/kets; S3 is the state we apply to. Linear
// response is the special case S1 = chi, S2 = S3 = Phi (ground state dressed as
// {phi, ...}); the VBC quadratic source is the case where S1/S2/S3 are the
// B/C-derived response states. compute_gamma is thus apply_g(chi, Phi, Phi).
//
// Exchange bra/ket order is the Dalton-validated convention
// apply_exchange(bra, ket, apply_to, lo) taken verbatim from the original
// per-type compute_gamma bodies (and pinned by tests/test_kernel_equivalence).
// ===========================================================================

#include "common_ops.hpp"   // common_ops::apply_exchange, poperatorT

#include <madness/chem/projector.h>
#include <madness/mra/mra.h>

#include <initializer_list>
#include <vector>

namespace molresponse_v3::two_electron {

using vecfuncT = std::vector<madness::real_function_3d>;

/// One exchange pairing K(bra, ket) to be applied to a gamma channel. `bra` and
/// `ket` form the density-like pair; the kernel applies the resulting exchange
/// operator to the third state's channel (`apply_to` in apply_channel). The
/// references must outlive the apply_channel call (they always do -- callers
/// pass kernel-local lvalues inside a braced-init-list argument).
struct ExPair {
  const vecfuncT &bra;
  const vecfuncT &ket;
};

/// Shared per-output-channel core of every Kernels<Type,Shell>::apply_g and
/// compute_gamma:
///     out = Q( J * apply_to  -  c_xc * Sum_pairs K(pair.bra, pair.ket)(apply_to) ).
/// J is the already-applied Coulomb potential coulop(rho); the caller builds it
/// once per gamma and passes it to each channel.
inline vecfuncT
apply_channel(madness::World &world,
              const madness::real_function_3d &J,
              const vecfuncT &apply_to,
              std::initializer_list<ExPair> pairs,
              const madness::QProjector<double, 3> &Q,
              double c_xc, double lo) {
  using namespace madness;
  auto out = mul(world, J, apply_to, true);
  if (c_xc > 0.0) {
    for (const auto &p : pairs) {
      auto k = common_ops::apply_exchange(world, p.bra, p.ket, apply_to, lo);
      gaxpy(world, 1.0, out, -c_xc, k);
    }
  }
  out = Q(out);
  truncate(world, out);
  return out;
}

} // namespace molresponse_v3::two_electron

#endif // MOLRESPONSE_V3_KERNELS_TWO_ELECTRON_HPP
