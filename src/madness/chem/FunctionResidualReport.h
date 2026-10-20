#pragma once

#include <madness/chem/FunctionResidual.h>
#include <nlohmann/json.hpp>
#include <cmath>

namespace madness {
// The containing aggregate supplies units, sample identity and sampling phase.
// Components are diagnostics and never introduce individual stopping gates.
inline nlohmann::json function_residual_record(const FunctionResidualBlocks& values) {
    using json = nlohmann::json;
    json blocks = json::object();
    bool valid = true;
    std::size_t count = 0;
    for (const auto& [block, functions] : values) {
        blocks[block] = json::array();
        for (std::size_t i = 0; i < functions.size(); ++i) {
            const double value = functions[i];
            const bool finite = std::isfinite(value) && value >= 0.0;
            valid = valid && finite;
            ++count;
            blocks[block].push_back({{"function_index", i},
                {"value", finite ? json(value) : json(nullptr)},
                {"availability", finite ? "measured" : "invalid"}});
        }
    }
    return {{"quantity", "raw_fixed_point_function_L2"},
        {"availability", valid && count ? "measured" : "invalid"},
        {"normalization", "none"}, {"reduction", "none"},
        {"active_gate", false}, {"blocks", blocks}};
}

inline bool valid_function_residual_record(const nlohmann::json& record) {
    // Malformed optional evidence is rejected; older records may omit it entirely.
    try {
        if (record.at("availability") != "measured" ||
            record.at("quantity") != "raw_fixed_point_function_L2" ||
            record.at("normalization") != "none" || record.at("reduction") != "none" ||
            record.at("active_gate") != false || !record.at("blocks").is_object()) return false;
        std::size_t count = 0;
        for (const auto& block : record.at("blocks").items()) {
            if (!block.value().is_array()) return false;
            std::size_t index = 0;
            for (const auto& function : block.value()) {
                if (function.at("availability") != "measured" ||
                    !function.at("function_index").is_number_integer() ||
                    function.at("function_index") != index++ ||
                    !function.at("value").is_number()) return false;
                const double value = function.at("value").get<double>();
                if (!std::isfinite(value) || value < 0.0) return false;
                ++count;
            }
        }
        return count != 0;
    } catch (const nlohmann::json::exception&) {
        return false;
    }
}
} // namespace madness
