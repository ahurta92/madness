#pragma once

#include "ResponseVector.hpp"
#include <madness/mra/nonlinsol.h>
#include <madness/mra/vmra.h>
#include <madness/world/world.h>
#include <vector>

// ============================================================================
// ResponseBundle<R> — typed collection of M response vectors
// ============================================================================
//
// Analogous to the legacy response_space class in molresponse/response_functions.h,
// but built on typed ResponseVector structs (StaticRestrictedResponse, etc.)
// rather than raw vector_real_function_3d.
//
// Provides:
//   - Element access:     bundle[i]   → R&
//   - Bundle-wide flat:   to_flat()   → M * flat_size functions (concatenated)
//                         from_flat() ← scatter back into per-state flats
//   - KAIN-compatible arithmetic: inner(), +=, +, -, *
//   - Sync / flatten:     flatten_all(), sync_all()
//
// Bundle flat layout (M * flat_size functions total):
//   [ state[0].flat | state[1].flat | ... | state[M-1].flat ]
//
// Used by iterate_excited(World, ResponseBundle<R>&, ...) for bundle-level
// KAIN acceleration.  The bundle KAIN treats all M states as a single vector
// for the Jacobian, which is preferable to M independent per-state solvers
// because it captures inter-state couplings introduced by the subspace rotation.
// ============================================================================

template <typename R>
class ResponseBundle {
public:
    ResponseBundle() = default;
    explicit ResponseBundle(std::vector<R> states) : states_(std::move(states)) {}

    // --- Access ---
    R& operator[](size_t i)             { return states_.at(i); }
    const R& operator[](size_t i) const { return states_.at(i); }

    size_t size() const noexcept { return states_.size(); }

    size_t num_orbitals() const noexcept {
        return states_.empty() ? 0 : states_[0].num_orbitals();
    }

    /// Number of flat slots per state (N * num_channels).
    size_t flat_size() const noexcept {
        return states_.empty() ? 0 : response_all(states_[0]).size();
    }

    // Raw vector access — for passing to existing kernel functions.
    std::vector<R>&       states()       { return states_; }
    const std::vector<R>& states() const { return states_; }

    // Range-for support.
    auto begin()       { return states_.begin(); }
    auto end()         { return states_.end();   }
    auto begin() const { return states_.begin(); }
    auto end()   const { return states_.end();   }

    // --- Sync / flatten ---
    void flatten_all() { for (auto& s : states_) s.flatten(); }
    void sync_all()    { for (auto& s : states_) s.sync(); }

    // --- Bundle-wide flat ---

    /// Concatenate all per-state flats into a single vector (M * flat_size).
    madness::vector_real_function_3d to_flat() const {
        madness::vector_real_function_3d result;
        const size_t fsz = flat_size();
        result.reserve(states_.size() * fsz);
        for (const auto& s : states_) {
            const auto &all = response_all(s);
            result.insert(result.end(), all.begin(), all.end());
        }
        return result;
    }

    /// Scatter a bundle-wide flat back into per-state flats and sync.
    /// Precondition: f.size() == size() * flat_size()
    void from_flat(const madness::vector_real_function_3d& f) {
        const size_t n = flat_size();
        MADNESS_ASSERT(f.size() == states_.size() * n);
        for (size_t i = 0; i < states_.size(); ++i)
            assign_all_and_sync(states_[i],
                madness::vector_real_function_3d(f.begin() + i * n,
                                                 f.begin() + (i + 1) * n));
    }

    // --- KAIN-compatible arithmetic ---
    //
    // All operators work element-wise on the per-state flat vectors, then
    // call sync() to propagate changes back to the typed channel members.

    ResponseBundle& operator+=(const ResponseBundle& b) {
        MADNESS_ASSERT(states_.size() == b.states_.size());
        for (size_t i = 0; i < states_.size(); ++i) {
            auto &lhs_all = response_all(states_[i]);
            const auto &rhs_all = response_all(b.states_[i]);
            madness::World& w = lhs_all[0].world();
            madness::gaxpy(w, 1.0, lhs_all, 1.0, rhs_all);
            states_[i].sync();
        }
        return *this;
    }

    friend ResponseBundle operator+(ResponseBundle a, const ResponseBundle& b) {
        return a += b;
    }

    friend ResponseBundle operator-(ResponseBundle a, const ResponseBundle& b) {
        MADNESS_ASSERT(a.states_.size() == b.states_.size());
        for (size_t i = 0; i < a.states_.size(); ++i) {
            auto &lhs_all = response_all(a.states_[i]);
            const auto &rhs_all = response_all(b.states_[i]);
            madness::World& w = lhs_all[0].world();
            madness::gaxpy(w, 1.0, lhs_all, -1.0, rhs_all);
            a.states_[i].sync();
        }
        return a;
    }

    friend ResponseBundle operator*(ResponseBundle a, double s) {
        for (auto& st : a.states_) {
            auto &all = response_all(st);
            madness::World& w = all[0].world();
            madness::scale(w, all, s);
            st.sync();
        }
        return a;
    }

    friend ResponseBundle operator*(double s, ResponseBundle a) { return a * s; }

private:
    std::vector<R> states_;
};

// ============================================================================
// inner() — scalar inner product for XNonlinearSolver (found via ADL)
// ============================================================================
//
// XNonlinearSolver calls inner(u, r) to build the KAIN subspace matrix.
// Since ResponseBundle is in the global namespace, ADL from within
// madness::XNonlinearSolver will find this free function.
//
// Returns Σ_i <a[i].flat | b[i].flat>  (sum of per-state Euclidean inner products).

template <typename R>
double inner(const ResponseBundle<R>& a, const ResponseBundle<R>& b) {
    MADNESS_ASSERT(a.size() == b.size());
    if (a.size() == 0) return 0.0;
    madness::World& world = response_all(a[0])[0].world();
    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
        result += madness::inner(world, response_all(a[i]), response_all(b[i])).sum();
    return result;
}

// ============================================================================
// Allocator + KAIN type alias
// ============================================================================

/// Allocator for XNonlinearSolver<ResponseBundle<R>, double, ...>.
///
/// Creates a zero bundle of M states, each with N * R::num_channels() flat slots.
/// The typed channel members (x_alpha, y_alpha, ...) are also initialized to
/// zero via sync() after the flat vector is filled.
template <typename R>
struct ResponseBundleAllocator {
    madness::World& world;
    size_t M;  ///< Number of states in the bundle
    size_t N;  ///< Number of orbitals per state

    ResponseBundle<R> operator()() const {
        const size_t flat_sz = N * R::num_channels();
        std::vector<R> states;
        states.reserve(M);
        for (size_t i = 0; i < M; ++i) {
            // R(N) creates properly sized typed channel vectors (null functions),
            // then calls flatten() so flat.size() == flat_sz.
            // We overwrite flat with real zero functions and sync back.
            R s(N);
            assign_all_and_sync(
                s,
                madness::zero_functions_compressed<double, 3>(
                    world, static_cast<int>(flat_sz)));
            states.push_back(std::move(s));
        }
        return ResponseBundle<R>(std::move(states));
    }
};

/// KAIN solver type for bundle-level excited-state acceleration.
template <typename R>
using bundle_solver = madness::XNonlinearSolver<ResponseBundle<R>, double,
                                                ResponseBundleAllocator<R>>;

// ============================================================================
// Factory function
// ============================================================================

/// Create a zero ResponseBundle of M states (N orbitals each) for ResponseType R.
template <typename R>
ResponseBundle<R> make_response_bundle(madness::World& world, size_t M, size_t N) {
    return ResponseBundleAllocator<R>{world, M, N}();
}
