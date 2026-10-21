#pragma once

#include <madness/chem/FunctionMetricSummary.h>
#include <nlohmann/json.hpp>
#include <cmath>
#include <string>

namespace molresponse_v3 {
// Root labels survive permutations, but the present eigensolver also mixes
// roots. Labels alone do not establish physical state correspondence. Report
// the positive trace of each stored root without inventing temporal changes.
inline nlohmann::json es_kinetic_evidence(const nlohmann::json& blocks,
    int iteration, int slot, int root_count, int bookkeeping_index, const std::string& model,
    const std::string& shell, double threshold, int k, bool locked) {
  using json = nlohmann::json;
  const auto trace = madness::function_metric_trace(blocks);
  const char* units = "hartree_times_response_amplitude_scale_squared";
  return {{"identity", {{"slot", slot}, {"root_count", root_count}, {"bookkeeping_index", bookkeeping_index},
              {"correspondence", "unvalidated"}, {"representation", "real_response_amplitudes"}, {"model", model}, {"shell", shell},
              {"reference_generation", nullptr}, {"reference_identity_reason", "native_bundle_without_named_reference_binding"}}},
      {"root_locked", locked},
      {"kinetic_energy", {{"quantity", "per_function_kinetic_quadratic"}, {"units", units},
          {"availability", trace.at("availability")}, {"blocks", blocks}, {"trace", trace},
          {"collection_definition", "all stored spin and X/Y blocks of this root; positive sum without occupation weights"},
          {"phase", "solver_step_output"},
          {"sampled_this_iteration", true}, {"normalization", "as_stored_no_rescaling"}, {"active_gate", false}, {"target", nullptr},
          {"sample", {{"iteration", iteration}, {"threshold", threshold}, {"k", k}, {"attempt_directory", ""}}}}},
      {"delta_kinetic_energy", {{"quantity", "absolute_kinetic_trace_change"}, {"units", units},
          {"availability", "unavailable"}, {"value", nullptr}, {"reason", "state_correspondence_unvalidated"},
          {"active_gate", false}, {"target", nullptr}, {"satisfied", nullptr}}}};
}

inline bool valid_es_kinetic_evidence(const nlohmann::json& report) {
  try {
    const auto& kinetic = report.at("kinetic_energy");
    const auto& identity = report.at("identity");
    const auto& sample = kinetic.at("sample");
    const auto& blocks = kinetic.at("blocks");
    const auto model = identity.at("model").get<std::string>();
    const auto shell = identity.at("shell").get<std::string>();
    if ((model != "full" && model != "tda") ||
        (shell != "closed_shell" && shell != "open_shell")) return false;
    const bool full = model == "full", open = shell == "open_shell";
    if (blocks.size() != (full ? 2U : 1U) * (open ? 2U : 1U) ||
        !blocks.contains("x_alpha") || blocks.at("x_alpha").empty()) return false;
    for (const auto* spin : {"alpha", "beta"}) {
      if (std::string(spin) == "beta" && !open) continue;
      const auto x = std::string("x_") + spin, y = std::string("y_") + spin;
      if (!blocks.contains(x) || (full && (!blocks.contains(y) || blocks.at(x).size() != blocks.at(y).size())))
        return false;
    }
    if (!sample.at("iteration").is_number_integer() || sample.at("iteration").get<int>() < 0 ||
        !sample.at("k").is_number_integer() || sample.at("k").get<int>() < 1 ||
        !sample.at("threshold").is_number() || !std::isfinite(sample.at("threshold").get<double>()) ||
        sample.at("threshold").get<double>() <= 0.0 || !sample.at("attempt_directory").is_string() ||
        !identity.at("slot").is_number_integer() || identity.at("slot").get<int>() < 0 ||
        !identity.at("root_count").is_number_integer() || identity.at("root_count").get<int>() < 1 ||
        identity.at("slot").get<int>() >= identity.at("root_count").get<int>() ||
        !identity.at("bookkeeping_index").is_number_integer() || identity.at("bookkeeping_index").get<int>() < -1 ||
        !report.at("root_locked").is_boolean()) return false;
    auto expected = es_kinetic_evidence(kinetic.at("blocks"), sample.at("iteration"), identity.at("slot"),
        identity.at("root_count"), identity.at("bookkeeping_index"), identity.at("model"), identity.at("shell"),
        sample.at("threshold"), sample.at("k"), report.at("root_locked"));
    expected["kinetic_energy"]["sample"]["attempt_directory"] = sample.at("attempt_directory");
    return expected.at("kinetic_energy").at("availability") == "measured" && expected == report;
  } catch (const std::exception&) { return false; }
}
} // namespace molresponse_v3
