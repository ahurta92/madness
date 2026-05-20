#ifndef MOLRESPONSE_V3_SOLVERS_RESPONSE_SUBSPACE_KAIN_HPP
#define MOLRESPONSE_V3_SOLVERS_RESPONSE_SUBSPACE_KAIN_HPP

// =========================================================================
// ResponseSubspaceKain<Storage> — hand-rolled per-state KAIN +
// per-orbital step restriction, modelled on SCF::update_subspace
// (src/madness/chem/SCF.cc:1916).
//
// Direct port of the SCF algorithm rather than XNonlinearSolver, so we
// can inspect the coefficient vector `c` and recover when KAIN's
// subspace becomes ill-conditioned. Each state (root for ES, channel
// for FD) gets its own independent history and Q matrix — different
// states converge at different rates and a shared history blends
// inconsistent iterates.
//
// Residual sign convention (matches SCF::update_subspace):
//
//     r = v − F(v)         (positive when v hasn't converged)
//     Q[i,j] = <v_i | r_j>
//     new_v = Σ_m c_m · (v_m − r_m) = Σ_m c_m · F(v_m)     ← KAIN(Q, rcond)
//
// where F is the BSH fixed-point map. The XNonlinearSolver-based
// KainAccelerator used the OPPOSITE residual sign (r = F(v) − v); we
// flip here so the SCF formula `c·(v − r)` works verbatim.
//
// Blowup safeguard (simplified vs SCF's `goto restart`):
//
//   while: c = KAIN(Q, rcond)
//          if |c|_max < 3:           accept
//          elif rcond < 0.01:        rcond *= 100, retry
//          else:                     clear history, return raw BSH
//                                    (skip KAIN this iter — simpler than
//                                    SCF's recursive restart)
//
// Step restriction (per-orbital, mirrors SCF::do_step_restriction at
// SCF.cc:2058) is applied AFTER KAIN, also AFTER raw BSH when KAIN is
// disabled or skipped. For each function p in the flat state:
//
//     if ||v[p] − new_v[p]|| > maxrotn:
//         new_v[p] := s · new_v[p] + (1−s) · v[p],  s = maxrotn / norm
//
// =========================================================================

#include "../kernels/tags.hpp"
#include "convergence_policy.hpp"
#include "response_state.hpp"

#include <madness/mra/mra.h>
#include <madness/tensor/solvers.h>   // madness::KAIN

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace molresponse_v3 {

template <typename Storage>
class ResponseSubspaceKain {
public:
  using vecfuncT = std::vector<madness::real_function_3d>;

  ResponseSubspaceKain(madness::World &world, const ConvergencePolicy &policy)
      : world_(world), policy_(policy) {}

  /// Apply per-state KAIN (if enabled) then per-orbital step restriction.
  ///   in[s]  = x_old for state s
  ///   out[s] = raw BSH output for state s; overwritten in place
  void apply(const std::vector<Storage> &in,
             std::vector<Storage>       &out) {
    if (states_.size() != in.size()) {
      states_.clear();
      states_.resize(in.size());
    }

    if (policy_.kain) {
      if (first_call_) {
        first_call_ = false;
        // Iter-0: skip KAIN extrapolation; history empty. Step
        // restriction still applies below.
      } else {
        for (std::size_t s = 0; s < in.size(); ++s) {
          apply_one_state(s, in[s], out[s]);
        }
      }
    }

    apply_step_restriction(in, out);
  }

  /// Drop all per-state KAIN histories. Call between protocol steps.
  void reset() {
    states_.clear();
    first_call_ = true;
  }

  void set_policy(const ConvergencePolicy &p) { policy_ = p; }

private:
  struct PerStateHistory {
    std::vector<std::pair<vecfuncT, vecfuncT>> entries; // (v_m, r_m)
    madness::Tensor<double>                    Q;        // m × m
  };

  /// One state's KAIN update. Mirrors SCF::update_subspace lines
  /// 1957-2028 step-for-step.
  void apply_one_state(std::size_t s_idx, const Storage &s_in,
                       Storage &s_out) {
    auto v = s_in.flatten();
    auto x_bsh = s_out.flatten();
    // SCF convention: residual r = v − F(v) (NOT F(v) − v)
    auto r = madness::sub(world_, v, x_bsh);

    auto &state = states_[s_idx];
    state.entries.push_back({v, r});
    const long m = static_cast<long>(state.entries.size());

    // ---- build new row + column of Q (SCF:1962-1979) -------------------
    madness::Tensor<double> ms(m), sm(m);
    for (long k = 0; k < m; ++k) {
      ms[k] = madness::inner(v, state.entries[k].second);
      sm[k] = madness::inner(state.entries[k].first, r);
    }
    madness::Tensor<double> newQ(m, m);
    if (m > 1) {
      newQ(madness::Slice(0, -2), madness::Slice(0, -2)) = state.Q;
    }
    newQ(m - 1, madness::_) = ms;
    newQ(madness::_, m - 1) = sm;
    state.Q = newQ;

    // ---- regularized solve with rcond escalation (SCF:1983-2004) -------
    double rcond = 1.0e-12;
    madness::Tensor<double> c;
    while (true) {
      c = madness::KAIN(state.Q, rcond);
      if (c.absmax() < 3.0) break;
      if (rcond < 0.01) {
        rcond *= 100.0;
        continue;
      }
      // Bail-out: clear history, leave s_out as raw BSH. KAIN will
      // start re-learning from the next iter. (Simpler than SCF's
      // recursive restart; the user opted for "clear and skip".)
      state.entries.clear();
      state.Q = madness::Tensor<double>();
      return;
    }

    // ---- linear combination (SCF:2013-2028) ----------------------------
    // new_v = Σ_k c_k · (v_k − r_k) = Σ_k c_k · F(v_k)
    auto new_v = madness::zero_functions<double, 3>(world_, v.size());
    for (long k = 0; k < m; ++k) {
      const auto &vm = state.entries[k].first;
      const auto &rm = state.entries[k].second;
      madness::gaxpy(world_, 1.0, new_v,  c(k), vm);
      madness::gaxpy(world_, 1.0, new_v, -c(k), rm);
    }

    // ---- sliding window (SCF:2030-2035) --------------------------------
    if (m >= policy_.kain_maxsub) {
      state.entries.erase(state.entries.begin());
      state.Q = state.Q(madness::Slice(1, -1), madness::Slice(1, -1));
    }

    s_out.from_flat(new_v);
  }

  /// Per-orbital step restriction (SCF::do_step_restriction, lines
  /// 2058-2084). Iterates per function in the flat state; each gets
  /// its own scale factor. A single bad orbital doesn't dominate the
  /// per-state norm and over-restrict the well-behaved ones.
  void apply_step_restriction(const std::vector<Storage> &in,
                              std::vector<Storage>       &out) {
    const double maxrotn = policy_.maxrotn;
    if (maxrotn <= 0.0) return;

    for (std::size_t s = 0; s < in.size(); ++s) {
      auto v_old = in[s].flatten();
      auto v_new = out[s].flatten();
      // Per-orbital diff norms via batched norm2s (collective-safe).
      auto diff  = madness::sub(world_, v_new, v_old);
      auto norms = madness::norm2s(world_, diff);
      bool changed = false;
      for (std::size_t p = 0; p < v_new.size(); ++p) {
        if (norms[p] > maxrotn) {
          const double scale = maxrotn / norms[p];
          // v_new[p] := scale·v_new[p] + (1−scale)·v_old[p]
          v_new[p].gaxpy(scale, v_old[p], 1.0 - scale, false);
          changed = true;
        }
      }
      if (changed) {
        world_.gop.fence();
        out[s].from_flat(v_new);
      }
    }
  }

  madness::World               &world_;
  ConvergencePolicy             policy_;
  std::vector<PerStateHistory>  states_;
  bool                          first_call_ = true;
};

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_RESPONSE_SUBSPACE_KAIN_HPP
