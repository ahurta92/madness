#pragma once

#include <nlohmann/json.hpp>

namespace madness {
class SCF;

/// Effective scientific inputs for ordinary and frozen SCF. No functions or
/// comparison history are copied. Protocol/layout and convergence requests have
/// their own shared records; known output, seed and iteration controls are omitted.
/// Unknown parameters participate conservatively. Call after configuring SCF.
[[nodiscard]] nlohmann::json scf_scientific_settings(const SCF& scf);

/// Physical model for ordinary restart compatibility. Removes only known
/// numerical precision controls; full snapshots still qualify saved metrics.
[[nodiscard]] nlohmann::json scf_model_settings(const nlohmann::json& settings);

/// Some external operator inputs are not exposed as resolved data. Their
/// measurements remain useful, but cannot qualify a later solve without checks.
[[nodiscard]] bool scf_settings_allow_reuse(const nlohmann::json& settings);
} // namespace madness
