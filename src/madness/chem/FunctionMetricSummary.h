#pragma once

#include <nlohmann/json.hpp>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace madness {
// Scalar evidence only: callers supply the physical quantity, units, sample
// identity and phase. No MPI, solver policy, weighting or tolerance lives here.
// Block names and contiguous indices identify stored functions within a sample.
inline nlohmann::json function_metric_maximum(const nlohmann::json& blocks) {
    using json = nlohmann::json;
    json result = {{"availability", "unavailable"}, {"reason", "empty_collection"},
        {"value", nullptr}, {"limiting_function", nullptr}, {"reduction", "maximum"}};
    try {
        if (!blocks.is_object()) throw std::runtime_error("invalid blocks");
        for (const auto& block : blocks.items()) {
            if (!block.value().is_array()) throw std::runtime_error("invalid block");
            std::size_t index = 0;
            for (const auto& function : block.value()) {
                if (function.at("availability") != "measured" ||
                    !function.at("function_index").is_number_integer() ||
                    function.at("function_index") != index ||
                    !function.at("value").is_number())
                    throw std::runtime_error("invalid component");
                const double value = function.at("value").get<double>();
                if (!std::isfinite(value) || value < 0.0)
                    throw std::runtime_error("invalid value");
                // Ties select the first block/index in deterministic map order.
                if (result.at("value").is_null() || value > result.at("value").get<double>()) {
                    result["value"] = value;
                    result["limiting_function"] = {{"block", block.key()}, {"function_index", index}};
                }
                ++index;
            }
        }
        if (!result.at("value").is_null()) {
            result["availability"] = "measured";
            result["reason"] = "";
        }
    } catch (const std::exception&) {
        result["availability"] = "invalid";
        result["reason"] = "invalid_component";
        result["value"] = nullptr;
        result["limiting_function"] = nullptr;
    }
    return result;
}

// Only call after establishing physical correspondence (same reference,
// equation, representation and sampling convention). Equal array sizes alone
// do not establish orbital identity across rotations or different solves.
inline nlohmann::json function_metric_absolute_change(const nlohmann::json& current,
                                                       const nlohmann::json& previous) {
    using json = nlohmann::json;
    auto result = function_metric_maximum(current);
    const auto prior = function_metric_maximum(previous);
    result["blocks"] = json::object();
    if (result.at("availability") == "invalid" || prior.at("availability") == "invalid") {
        result["availability"] = "invalid";
        result["reason"] = "invalid_component";
    } else if (result.at("availability") != "measured" || prior.at("availability") != "measured") {
        result["availability"] = "unavailable";
        result["reason"] = "empty_collection";
    } else {
        bool matched = current.size() == previous.size();
        for (const auto& block : current.items())
            matched = matched && previous.contains(block.key()) &&
                      previous.at(block.key()).size() == block.value().size();
        if (matched) {
            auto changes = current;
            for (auto& block : changes.items())
                for (std::size_t i = 0; i < block.value().size(); ++i)
                    block.value()[i]["value"] = std::abs(
                        current.at(block.key())[i].at("value").get<double>() -
                        previous.at(block.key())[i].at("value").get<double>());
            result = function_metric_maximum(changes);
            result["blocks"] = std::move(changes);
            return result;
        }
        result["availability"] = "unavailable";
        result["reason"] = "incomparable_function_layout";
    }
    result["value"] = nullptr;
    result["limiting_function"] = nullptr;
    return result;
}
} // namespace madness
