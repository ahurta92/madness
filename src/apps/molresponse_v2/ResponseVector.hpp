#pragma once
#include <madness/mra/vmra.h>
#include <madness/world/world.h>
#include "ResponseState.hpp"
#include <string>
#include <type_traits>
#include <variant> // C++ 17

using namespace madness;

// ============================================================================
// ResponseVector — Theory and Code Mapping
// ============================================================================
//
// TDDFT/TDHF linear response in the density-matrix formulation
// ------------------------------------------------------------
// The time-dependent one-particle density matrix obeys
//
//   i ∂γ/∂t = [F, γ]
//
// where F = h + g + v(t) is the time-dependent Fock operator and
// v(t) = Σ_C λ_C [ v^{C†}(ω_C) e^{+iω_C t} + v^C(ω_C) e^{-iω_C t} ]
// is an external perturbation (e.g. dipole operator along x, y, or z).
//
// Expanding γ = γ^0 + Σ_C λ_C γ^C + ... and collecting first-order terms
// gives the frequency-domain response equation:
//
//   ω_C γ^C = [F^0, γ^C] + [F^C, γ^0]
//
// First-order response density matrix
// ------------------------------------
// The first-order density matrix takes the form
//
//   γ^C(r,r') = Σ_i [ x_i^C(r) φ_i†(r')  +  φ_i(r) y_i^C(r') ]
//
// where φ_i are the N ground-state occupied orbitals and
//   x_i^C(r) — forward  response orbital (positive-frequency channel)
//   y_i^C(r) — backward response orbital (negative-frequency channel)
//
// These map directly to the data members of each ResponseVector concrete type:
//   x_i^C  ↔  x_alpha[i]   (closed-shell) or x_alpha[i], x_beta[i] (open-shell)
//   y_i^C  ↔  y_alpha[i]   (closed-shell) or y_alpha[i], y_beta[i] (open-shell)
//
// The response density used in property assembly is:
//   ρ¹(r) = Σ_i [ x_i^C(r) φ_i*(r)  +  y_i^{C*}(r) φ_i(r) ]
//   → see compute_response_density() in ResponseKernels.hpp
//
// Tamm-Dancoff Approximation (TDA / "static")
// -------------------------------------------
// Setting y ≡ 0 gives the TDA (Tamm-Dancoff Approximation):
//   γ^C ≈ Σ_i x_i^C(r) φ_i†(r')
// The y-channel is not stored.
//   TDA excited-state solves:   use TDARestrictedResponse (closed-shell)
//   ω=0 frequency response:     use StaticRestrictedResponse (closed-shell)
//   Open-shell TDA:             use StaticUnrestrictedResponse
//
// NOTE: TDARestrictedResponse and StaticRestrictedResponse have the same flat
// layout and the same x-only kernel, but are semantically distinct types:
//   - StaticRestrictedResponse: the ω=0 physical limit of full TDDFT
//   - TDARestrictedResponse: the B=0 approximation for excited states
//
// BSH integral-equation form
// --------------------------
// In MRA the differential equations are ill-conditioned; we solve instead:
//
//   x_p^C = -2 Ĝ(k_p^x) * [ V^0 x_p^C - Σ_{i≠p} ε_{ip} x_i^C + g_p'[γ^C]φ_p + V_p^C ]
//   y_p^C = -2 Ĝ(k_p^y) * [ V^0 y_p^C - Σ_{i≠p} ε_{ip} y_i^C + g_p'[γ^{C†}]φ_p + V_p^{C†} ]
//
// where Ĝ(k) is the Bound-State Helmholtz (BSH) Green's function and
//   k_p^x = sqrt(-2(ε_p + ω_C))   [x-channel BSH exponent]
//   k_p^y = sqrt(-2(ε_p - ω_C))   [y-channel BSH exponent]
//   ε_{ip} = ∫ φ_i†(r) F^0 φ_p(r) dr   [off-diagonal Fock, orbital localization]
//   V_p^C = Q̂ v^C φ_p,   Q̂ = 1 - γ^0  [projection onto unoccupied space]
//   g_p'[γ^C] φ_p  [response exchange-correlation kernel, see ResponseKernels.hpp]
//
// → k_p^x and k_p^y are the mu parameters in make_bsh_operators() /
//   make_excited_bsh_operators() in ResponseKernels.hpp.
//
// Symplectic metric (excited-state subspace rotation)
// ---------------------------------------------------
// For excited-state calculations the response states are rotated in a
// subspace that respects the generalized eigenvalue problem.
// The metric inner product is:
//
//   ⟨Φ|Φ'⟩ = ⟨x|x'⟩ − ⟨y|y'⟩     (y enters with MINUS sign)
//
// The minus sign arises from the symplectic structure of the TDDFT response
// matrix (Casida equation).  It is NOT arbitrary.
// → see response_metric_inner() in ResponseKernels.hpp
//
// Spin structure
// --------------
//   Restricted (closed-shell): x_alpha = x_beta by symmetry → only alpha stored.
//   Unrestricted (open-shell): independent alpha and beta spin channels.
//
// ============================================================================
// Flat vector layout — channel ordering within each concrete type
// ============================================================================
//
// Each struct stores a `flat` vector that concatenates all channels.
// This is the form consumed by KAIN, BSH application, and archive I/O.
//
//   StaticRestrictedResponse     (N slots, ω=0 frequency response, closed-shell):
//     flat = [ x_alpha[0..N-1] ]
//             channel 0 @ offset 0
//
//   TDARestrictedResponse        (N slots, TDA excited-state, closed-shell):
//     flat = [ x_alpha[0..N-1] ]
//             channel 0 @ offset 0
//     NOTE: same layout as StaticRestrictedResponse but a distinct type.
//
//   DynamicRestrictedResponse    (2N slots, full TDDFT closed-shell):
//     flat = [ x_alpha[0..N-1]  |  y_alpha[0..N-1] ]
//             channel 0 @ 0        channel 1 @ N
//
//   StaticUnrestrictedResponse   (2N slots, TDA open-shell):
//     flat = [ x_alpha[0..N-1]  |  x_beta[0..N-1] ]
//             channel 0 @ 0        channel 1 @ N
//
//   DynamicUnrestrictedResponse  (4N slots, full TDDFT open-shell):
//     flat = [ x_alpha[0..N-1]  |  y_alpha[0..N-1]  |  x_beta[0..N-1]  |  y_beta[0..N-1] ]
//             channel 0 @ 0        channel 1 @ N        channel 2 @ 2N     channel 3 @ 3N
//
// WARNING: This layout is fixed by the binary archive format (ResponseIO.hpp).
// Do NOT reorder channels without writing archive migration code.
//
// ============================================================================
// Sync/flatten invariant
// ============================================================================
//
// Each struct maintains TWO synchronized views of the same data:
//   - Typed channel members: x_alpha, y_alpha, x_beta, y_beta
//   - Concatenated flat vector: flat
//
// Contract:
//   After modifying flat directly (KAIN update, BSH apply, archive load):
//       call sync() to propagate flat → typed channels.
//   After modifying typed channels directly (guess init, pack_guess_states):
//       call flatten() to propagate typed channels → flat.
//   Never call both in sequence without an intervening modification —
//       one will silently overwrite the other's work.
//
// ============================================================================

struct StaticRestrictedResponse;
struct TDARestrictedResponse;
struct DynamicRestrictedResponse;
struct StaticUnrestrictedResponse;
struct DynamicUnrestrictedResponse;


using ResponseVector =
std::variant<StaticRestrictedResponse, DynamicRestrictedResponse,
    StaticUnrestrictedResponse, DynamicUnrestrictedResponse>;

/// Closed-shell Tamm-Dancoff (TDA) response state.
///
/// Represents the first-order density matrix in the TDA approximation (y ≡ 0):
///   γ^C(r,r') = Σ_i x_i^C(r) φ_i†(r')
///
/// BSH exponent per orbital p:  k_p = sqrt(-2(ε_p + ω_C))
/// Flat layout: [ x_alpha[0..N-1] ]   (N slots)
/// alpha_factor = -4.0
///   (factor 2 from two spin channels × factor 2 from restricted normalization)
struct StaticRestrictedResponse {
    vector_real_function_3d x_alpha;
    vector_real_function_3d flat;

    StaticRestrictedResponse() = default;

    explicit StaticRestrictedResponse(const size_t &n_orb) : x_alpha(n_orb) {
        flatten();
    }

    /// Number of ground-state occupied orbitals N.
    [[nodiscard]] size_t num_orbitals() const noexcept { return x_alpha.size(); }
    /// Number of logical channels stored in flat (1 for TDA restricted).
    [[nodiscard]] static constexpr size_t num_channels() noexcept { return 1; }
    /// Total flat slots: num_orbitals() * num_channels().
    [[nodiscard]] size_t num_flat_slots() const noexcept { return num_orbitals() * num_channels(); }
    /// Start index in flat for logical channel k (all channels are N-wide).
    [[nodiscard]] size_t channel_offset(size_t k) const noexcept { return k * num_orbitals(); }
    /// Prefactor for polarizability assembly: α_{AB} = alpha_factor() * ⟨x_A|V_B⟩.
    /// Derivation: ρ¹ = 2 Σ x_i φ_i* (spin factor 2); α = -2*2*⟨x|v⟩ = -4⟨x|v⟩.
    [[nodiscard]] static constexpr double alpha_factor() noexcept { return -4.0; }

    /// Propagate flat → typed channels.  Call after any direct modification of flat.
    void sync() {
        MADNESS_ASSERT(flat.size() == num_flat_slots());
        for (size_t i = 0; i < x_alpha.size(); ++i)
            x_alpha[i] = flat[channel_offset(0) + i];
    }

    /// Propagate typed channels → flat.  Call after modifying x_alpha directly.
    void flatten() { flat = x_alpha; }
};

/// Closed-shell Tamm-Dancoff Approximation (TDA) excited-state response.
///
/// Represents the excited-state density matrix in the TDA (B=0) approximation:
///   γ^C(r,r') = Σ_i x_i^C(r) φ_i†(r')    (y ≡ 0 by approximation, not by ω→0)
///
/// Distinguished from StaticRestrictedResponse (ω=0 frequency response) in that:
///   - TDA:    y=0 is an approximation; kernel uses A-matrix only (x-only density)
///   - Static: ω=0 is a physical limit; correct kernel includes B-matrix (x+y density)
///
/// Both types share the same flat layout and the same x-only kernel implementation,
/// but are kept as separate types so code paths are semantically clear.
///
/// Flat layout: [ x_alpha[0..N-1] ]   (N slots)
/// alpha_factor = -4.0
struct TDARestrictedResponse {
    vector_real_function_3d x_alpha;
    vector_real_function_3d flat;

    TDARestrictedResponse() = default;

    explicit TDARestrictedResponse(const size_t &n_orb) : x_alpha(n_orb) {
        flatten();
    }

    /// Number of ground-state occupied orbitals N.
    [[nodiscard]] size_t num_orbitals() const noexcept { return x_alpha.size(); }
    /// Number of logical channels stored in flat (1 for TDA restricted).
    [[nodiscard]] static constexpr size_t num_channels() noexcept { return 1; }
    /// Total flat slots: num_orbitals() * num_channels().
    [[nodiscard]] size_t num_flat_slots() const noexcept { return num_orbitals() * num_channels(); }
    /// Start index in flat for logical channel k (all channels are N-wide).
    [[nodiscard]] size_t channel_offset(size_t k) const noexcept { return k * num_orbitals(); }
    /// Prefactor for polarizability assembly: α_{AB} = alpha_factor() * ⟨x_A|V_B⟩.
    /// Derivation: ρ¹ = 2 Σ x_i φ_i* (spin factor 2); α = -2*2*⟨x|v⟩ = -4⟨x|v⟩.
    [[nodiscard]] static constexpr double alpha_factor() noexcept { return -4.0; }

    /// Propagate flat → typed channels.  Call after any direct modification of flat.
    void sync() {
        MADNESS_ASSERT(flat.size() == num_flat_slots());
        for (size_t i = 0; i < x_alpha.size(); ++i)
            x_alpha[i] = flat[channel_offset(0) + i];
    }

    /// Propagate typed channels → flat.  Call after modifying x_alpha directly.
    void flatten() { flat = x_alpha; }
};

/// Closed-shell full TDDFT/TDHF response state (frequency-dependent).
///
/// Represents the first-order density matrix including both x and y channels:
///   γ^C(r,r') = Σ_i [ x_i^C(r) φ_i†(r')  +  φ_i(r) y_i^C(r') ]
///
/// BSH exponents per orbital p:
///   x-channel: k_p^x = sqrt(-2(ε_p + ω_C))
///   y-channel: k_p^y = sqrt(-2(ε_p - ω_C))
///
/// Symplectic metric: ⟨Φ|Φ'⟩ = ⟨x|x'⟩ − ⟨y|y'⟩   (y enters with MINUS sign)
/// Flat layout: [ x_alpha[0..N-1] | y_alpha[0..N-1] ]   (2N slots)
/// alpha_factor = -2.0
struct DynamicRestrictedResponse {
    vector_real_function_3d x_alpha;
    vector_real_function_3d y_alpha;
    vector_real_function_3d flat;

    DynamicRestrictedResponse() = default;

    explicit DynamicRestrictedResponse(const size_t &n_orb)
        : x_alpha(n_orb), y_alpha(n_orb) {
        flatten();
    }

    /// Number of ground-state occupied orbitals N.
    [[nodiscard]] size_t num_orbitals() const noexcept { return x_alpha.size(); }
    /// Number of logical channels: 2 (x_alpha and y_alpha).
    [[nodiscard]] static constexpr size_t num_channels() noexcept { return 2; }
    /// Total flat slots: num_orbitals() * num_channels().
    [[nodiscard]] size_t num_flat_slots() const noexcept { return num_orbitals() * num_channels(); }
    /// Start index in flat for logical channel k (x=0, y=1).
    [[nodiscard]] size_t channel_offset(size_t k) const noexcept { return k * num_orbitals(); }
    /// Prefactor for polarizability assembly: α_{AB} = alpha_factor() * ⟨x_A|V_B⟩.
    [[nodiscard]] static constexpr double alpha_factor() noexcept { return -2.0; }

    /// Propagate flat → typed channels.  Call after any direct modification of flat.
    void sync() {
        MADNESS_ASSERT(flat.size() == num_flat_slots());
        for (size_t i = 0; i < x_alpha.size(); ++i) {
            x_alpha[i] = flat[channel_offset(0) + i];
            y_alpha[i] = flat[channel_offset(1) + i];
        }
    }

    /// Propagate typed channels → flat.  Call after modifying x_alpha or y_alpha directly.
    void flatten() {
        flat = x_alpha;
        flat.insert(flat.end(), y_alpha.begin(), y_alpha.end());
    }
};

/// Open-shell Tamm-Dancoff (TDA) response state.
///
/// Represents the TDA density matrix for unrestricted (open-shell) systems,
/// with independent alpha and beta spin channels:
///   γ^C = Σ_i [ x_i^{C,α}(r) φ_i^{α†}(r')  +  x_i^{C,β}(r) φ_i^{β†}(r') ]
///
/// Flat layout: [ x_alpha[0..N-1] | x_beta[0..N-1] ]   (2N slots)
/// alpha_factor = -2.0
///
/// NOTE: Unrestricted support in the solver is not yet fully implemented.
/// Most kernel functions assert on unrestricted types.
struct StaticUnrestrictedResponse {
    vector_real_function_3d x_alpha;
    vector_real_function_3d x_beta;
    vector_real_function_3d flat;

    StaticUnrestrictedResponse() = default;

    explicit StaticUnrestrictedResponse(const size_t &n_orb)
        : x_alpha(n_orb), x_beta(n_orb) {
        flatten();
    }

    /// Number of ground-state occupied orbitals N.
    [[nodiscard]] size_t num_orbitals() const noexcept { return x_alpha.size(); }
    /// Number of logical channels: 2 (x_alpha and x_beta).
    [[nodiscard]] static constexpr size_t num_channels() noexcept { return 2; }
    /// Total flat slots: num_orbitals() * num_channels().
    [[nodiscard]] size_t num_flat_slots() const noexcept { return num_orbitals() * num_channels(); }
    /// Start index in flat for logical channel k (x_alpha=0, x_beta=1).
    [[nodiscard]] size_t channel_offset(size_t k) const noexcept { return k * num_orbitals(); }
    /// Prefactor for polarizability assembly: α_{AB} = alpha_factor() * ⟨x_A|V_B⟩.
    [[nodiscard]] static constexpr double alpha_factor() noexcept { return -2.0; }

    /// Propagate flat → typed channels.  Call after any direct modification of flat.
    void sync() {
        MADNESS_ASSERT(flat.size() == num_flat_slots());
        for (size_t i = 0; i < x_alpha.size(); ++i) {
            x_alpha[i] = flat[channel_offset(0) + i];
            x_beta[i]  = flat[channel_offset(1) + i];
        }
    }

    /// Propagate typed channels → flat.  Call after modifying x_alpha or x_beta directly.
    void flatten() {
        flat = x_alpha;
        flat.insert(flat.end(), x_beta.begin(), x_beta.end());
    }
};

/// Open-shell full TDDFT/TDHF response state (frequency-dependent).
///
/// Represents the first-order density matrix for unrestricted (open-shell) systems,
/// with all four independent spin-frequency channels:
///   γ^C = Σ_i [ x_i^{C,α} φ_i^{α†}  +  φ_i^α y_i^{C,α}
///              + x_i^{C,β} φ_i^{β†}  +  φ_i^β y_i^{C,β} ]
///
/// Flat layout: [ x_alpha[0..N-1] | y_alpha[0..N-1] | x_beta[0..N-1] | y_beta[0..N-1] ]
///               channel 0 @ 0       channel 1 @ N      channel 2 @ 2N   channel 3 @ 3N
/// alpha_factor = -2.0
///
/// NOTE: Unrestricted support in the solver is not yet fully implemented.
/// Most kernel functions assert on unrestricted types.
struct DynamicUnrestrictedResponse {
    vector_real_function_3d x_alpha;
    vector_real_function_3d y_alpha;
    vector_real_function_3d x_beta;
    vector_real_function_3d y_beta;
    vector_real_function_3d flat;

    DynamicUnrestrictedResponse() = default;

    explicit DynamicUnrestrictedResponse(const size_t &n_orb)
        : x_alpha(n_orb), y_alpha(n_orb), x_beta(n_orb), y_beta(n_orb) {
        flatten();
        // NOTE: the commented-out line below was a prior bug (y_beta appended twice).
        // flatten() already includes all four channels correctly.
        // flat.insert(flat.end(), y_beta.begin(), y_beta.end());
    }

    /// Number of ground-state occupied orbitals N.
    [[nodiscard]] size_t num_orbitals() const noexcept { return x_alpha.size(); }
    /// Number of logical channels: 4 (x_alpha, y_alpha, x_beta, y_beta).
    [[nodiscard]] static constexpr size_t num_channels() noexcept { return 4; }
    /// Total flat slots: num_orbitals() * num_channels().
    [[nodiscard]] size_t num_flat_slots() const noexcept { return num_orbitals() * num_channels(); }
    /// Start index in flat for logical channel k (x_α=0, y_α=1, x_β=2, y_β=3).
    [[nodiscard]] size_t channel_offset(size_t k) const noexcept { return k * num_orbitals(); }
    /// Prefactor for polarizability assembly: α_{AB} = alpha_factor() * ⟨x_A|V_B⟩.
    [[nodiscard]] static constexpr double alpha_factor() noexcept { return -2.0; }

    /// Propagate typed channels → flat.  Call after modifying any typed channel directly.
    void flatten() {
        flat = x_alpha;
        flat.insert(flat.end(), y_alpha.begin(), y_alpha.end());
        flat.insert(flat.end(), x_beta.begin(), x_beta.end());
        flat.insert(flat.end(), y_beta.begin(), y_beta.end());
    }

    /// Propagate flat → typed channels.  Call after any direct modification of flat.
    void sync() {
        MADNESS_ASSERT(flat.size() == num_flat_slots());
        for (size_t i = 0; i < x_alpha.size(); ++i) {
            x_alpha[i] = flat[channel_offset(0) + i];
            y_alpha[i] = flat[channel_offset(1) + i];
            x_beta[i]  = flat[channel_offset(2) + i];
            y_beta[i]  = flat[channel_offset(3) + i];
        }
    }
};


/// Create a ResponseVector given the number of orbitals and two boolean flags.
/// @param n_orb     Number of ground-state occupied orbitals.
/// @param is_static True → TDA (no y-channel); false → full TDDFT.
/// @param is_unrestricted True → open-shell (separate alpha/beta).
inline ResponseVector make_response_vector(size_t n_orb, bool is_static,
                                           bool is_unrestricted) {
    if (!is_unrestricted && is_static) {
        return StaticRestrictedResponse(n_orb);
    } else if (!is_unrestricted && !is_static) {
        return DynamicRestrictedResponse(n_orb);
    } else if (is_unrestricted && is_static) {
        return StaticUnrestrictedResponse(n_orb);
    } else {
        return DynamicUnrestrictedResponse(n_orb);
    }
}

inline std::string response_variant_name(bool is_static,
                                         bool is_unrestricted) {
    if (!is_unrestricted && is_static) {
        return "static_restricted";
    } else if (!is_unrestricted && !is_static) {
        return "dynamic_restricted";
    } else if (is_unrestricted && is_static) {
        return "static_unrestricted";
    } else {
        return "dynamic_unrestricted";
    }
}

inline std::string response_variant_name(const ResponseVector &vec) {
    return std::visit(
        [](const auto &v) {
            using T = std::decay_t<decltype(v)>;
            if constexpr (std::is_same_v<T, StaticRestrictedResponse>) {
                return std::string("static_restricted");
            } else if constexpr (std::is_same_v<T, DynamicRestrictedResponse>) {
                return std::string("dynamic_restricted");
            } else if constexpr (std::is_same_v<T, StaticUnrestrictedResponse>) {
                return std::string("static_unrestricted");
            } else {
                return std::string("dynamic_unrestricted");
            }
        },
        vec);
}

inline bool response_variant_is_static(const std::string &variant) {
    return variant == "static_restricted" || variant == "static_unrestricted";
}

inline bool response_variant_is_unrestricted(const std::string &variant) {
    return variant == "static_unrestricted" ||
           variant == "dynamic_unrestricted";
}

inline bool response_variant_is_valid(const std::string &variant) {
    return variant == "static_restricted" ||
           variant == "dynamic_restricted" ||
           variant == "static_unrestricted" ||
           variant == "dynamic_unrestricted";
}

inline ResponseVector make_response_vector_from_variant(
    size_t num_orbitals, const std::string &variant) {
    return make_response_vector(num_orbitals, response_variant_is_static(variant),
                                response_variant_is_unrestricted(variant));
}

inline size_t response_num_orbitals(const ResponseVector &vec) {
    return std::visit(
        [](const auto &v) -> size_t { return v.x_alpha.size(); }, vec);
}

template <typename ResponseType>
constexpr bool response_has_y_channel_v =
    std::is_same_v<ResponseType, DynamicRestrictedResponse> ||
    std::is_same_v<ResponseType, DynamicUnrestrictedResponse>;

template <typename ResponseType>
constexpr bool response_is_unrestricted_v =
    std::is_same_v<ResponseType, StaticUnrestrictedResponse> ||
    std::is_same_v<ResponseType, DynamicUnrestrictedResponse>;

inline bool response_has_y_channel(const ResponseVector &vec) {
    return std::visit(
        [](const auto &v) {
            using T = std::decay_t<decltype(v)>;
            return response_has_y_channel_v<T>;
        },
        vec);
}

inline bool response_is_unrestricted(const ResponseVector &vec) {
    return std::visit(
        [](const auto &v) {
            using T = std::decay_t<decltype(v)>;
            return response_is_unrestricted_v<T>;
        },
        vec);
}

inline void flatten_response(ResponseVector &vec) {
    std::visit([](auto &v) { v.flatten(); }, vec);
}

inline void sync_response(ResponseVector &vec) {
    std::visit([](auto &v) { v.sync(); }, vec);
}

inline vector_real_function_3d &get_flat(ResponseVector &vec) {
    return std::visit(
        [](auto &v) -> vector_real_function_3d & { return v.flat; }, vec);
}

/// Get the flat vector from a (const) ResponseVector.
inline const vector_real_function_3d &get_flat(const ResponseVector &vec) {
    return std::visit(
        [](const auto &v) -> const vector_real_function_3d & { return v.flat; },
        vec);
}

/// Runtime dispatch of alpha_factor() for a ResponseVector variant.
/// Returns -4.0 for StaticRestricted (TDA, closed-shell) and -2.0 for all others.
/// See each struct's alpha_factor() for the physical derivation.
inline double response_alpha_factor(const ResponseVector &vec) {
    return std::visit([](const auto &v) { return v.alpha_factor(); }, vec);
}

/// Total number of flat slots for a ResponseVector: N * num_channels().
inline size_t response_num_flat_slots(const ResponseVector &vec) {
    return std::visit([](const auto &v) { return v.num_flat_slots(); }, vec);
}

/// Number of logical channels for a ResponseVector (1, 2, or 4).
inline size_t response_num_channels(const ResponseVector &vec) {
    return std::visit([](const auto &v) { return v.num_channels(); }, vec);
}

/// Atomically replace flat and sync typed channels.
///
/// Use this instead of the two-liner:
///   r.flat = expr;
///   r.sync();
///
/// This eliminates the footgun of assigning flat without re-syncing the typed
/// channel members (x_alpha, y_alpha, etc.).  The two-liner is safe only if
/// nothing reads the typed channels between the two lines; this helper makes
/// the atomicity explicit and auditable.
///
/// Do NOT use this when ensure_initialized_flat() must run between the
/// assignment and the sync — in that case keep the three-line pattern explicit:
///   r.flat = expr;
///   ensure_initialized_flat(world, r.flat);
///   r.sync();
template <typename ResponseType>
inline void assign_flat_and_sync(ResponseType &r, vector_real_function_3d f) {
    r.flat = std::move(f);
    r.sync();
}

/// Variant-dispatched assign_flat_and_sync for a ResponseVector.
inline void assign_flat_and_sync(ResponseVector &r, vector_real_function_3d f) {
    std::visit([&](auto &v) { assign_flat_and_sync(v, std::move(f)); }, r);
}
