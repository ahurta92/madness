#pragma once

#include "../kernels/response_space_ops.hpp"
#include <nlohmann/json.hpp>
#include <cmath>
#include <map>
#include <string>
#include <vector>

namespace molresponse_v3 {
// Keys identify stored blocks; indices identify functions in the reference's
// orbital ordering. Static response has no fabricated Y or beta contribution.
using FunctionKinetics = std::map<std::string, std::vector<double>>;
struct FunctionKineticSample {
  int iteration;
  std::vector<FunctionKinetics> channels;
};

template<class Storage>
FunctionKinetics measure_function_kinetics(madness::World& world, const Storage& state) {
  FunctionKinetics values{{"x_alpha", rs::kinetic_diagonal(world, state.x_alpha)}};
  if constexpr (rs::detail::has_y_alpha<Storage>::value)
    values["y_alpha"] = rs::kinetic_diagonal(world, state.y_alpha);
  if constexpr (rs::detail::has_x_beta<Storage>::value)
    values["x_beta"] = rs::kinetic_diagonal(world, state.x_beta);
  if constexpr (rs::detail::has_y_beta<Storage>::value)
    values["y_beta"] = rs::kinetic_diagonal(world, state.y_beta);
  return values;
}

// Pure reporting: all collective numerical work has already completed.
inline nlohmann::json function_kinetic_record(const FunctionKinetics& values,
    int iteration, const std::string& attempt, double threshold, int k,
    bool saved, const char* phase) {
  using json = nlohmann::json;
  json blocks = json::object();
  bool valid = true;
  for (const auto& [block, functions] : values) {
    blocks[block] = json::array();
    for (std::size_t i = 0; i < functions.size(); ++i) {
      const double value = functions[i];
      const bool finite = std::isfinite(value) && value >= 0.0;
      valid = valid && finite;
      blocks[block].push_back({{"function_index", i},
          {"value", finite ? json(value) : json(nullptr)},
          {"availability", finite ? "measured" : "invalid"}});
    }
  }
  return {{"availability", valid ? "measured" : "invalid"},
      {"quantity", "per_function_kinetic_quadratic"},
      {"units", "hartree_times_response_amplitude_scale_squared"},
      {"reduction", "none"}, {"active_gate", false}, {"target", nullptr},
      {"satisfied", nullptr}, {"valid_for_saved_state", saved && valid},
      {"phase", phase}, {"blocks", blocks},
      {"sample", {{"attempt_directory", attempt}, {"iteration", iteration},
                  {"threshold", threshold}, {"k", k}}}};
}

// Validate actual scalars too: a corrupt component cannot hide behind a
// top-level "measured" flag. Missing records on old checkpoints are handled
// by the caller and never manufacture evidence.
inline bool valid_function_kinetic_record(const nlohmann::json& record) {
  if (record.value("availability", "") != "measured" ||
      !record.contains("blocks") || !record.at("blocks").is_object()) return false;
  for (const auto& block : record.at("blocks").items()) {
    if (!block.value().is_array()) return false;
    std::size_t index = 0;
    for (const auto& function : block.value()) {
      if (function.value("availability", "") != "measured" ||
          !function.contains("value") || !function.at("value").is_number() ||
          function.value("function_index", -1) != static_cast<int>(index++)) return false;
      const double value = function.at("value").get<double>();
      if (!std::isfinite(value) || value < 0.0) return false;
    }
  }
  return true;
}
} // namespace molresponse_v3
