#pragma once
// The results.json form of the Raman properties. Pure JSON, no runtime, so the
// shape is unit-tested without a ground state (raman_analysis.hpp computes the
// numbers and calls these).
#include <madness/chem/RunFiles.h>

#include <cmath>
#include <string>

namespace molresponse_v3 {

/// One raman_spectrum item: the spectrum at one laser frequency.
///
/// `hessian` holds the sources of the Hessian (the static nuclear states);
/// `gradient` is the polarizability_gradient property (null when none), whose
/// item at this frequency adds its sources. `verdict` is the frequency's
/// completeness record {rows_found, rows_expected, all_rows_converged}, null
/// when no Hessian was computed; `rows` its mode rows, null when the spectrum
/// was not assembled. A complete frequency of a molecule with no vibrational
/// mode has an empty spectrum.
inline nlohmann::json raman_spectrum_item(double frequency, const nlohmann::json& hessian,
                                          const nlohmann::json& gradient, const nlohmann::json& verdict,
                                          const nlohmann::json& rows) {
    namespace it = madness::iteration;
    nlohmann::json item{{"frequency", frequency}};
    const auto add = [&](const nlohmann::json& from) {
        for (const auto& source : from.value(it::results_sources_key, nlohmann::json::array()))
            it::add_results_source(item, source.at("node").get<std::string>());
    };
    add(hessian);
    if (gradient.is_object())
        for (const auto& value : gradient.value("values", nlohmann::json::array()))
            if (std::abs(value.value("frequency", -1.0) - frequency) < 1e-9) add(value);
    if (verdict.is_null()) {
        item["status"] = "unavailable";
        item["availability"] = {{"spectrum", "qualified Hessian and derivatives are required"}};
        return item;
    }
    const bool whole = verdict.at("rows_found") == verdict.at("rows_expected") &&
                       verdict.at("all_rows_converged").get<bool>();
    if (whole || !rows.is_null()) item["value"] = rows.is_null() ? nlohmann::json::array() : rows;
    item["status"] = whole ? "complete" : rows.is_null() ? "unavailable" : "partial";
    if (!whole) item["availability"] = {{"spectrum", "polarizability derivatives are incomplete"}};
    return item;
}

/// The vibrational_analysis property: `value` holds the Hessian, frequencies,
/// normal modes and reduced masses; null `value` when no Hessian was computed.
inline nlohmann::json vibrational_analysis_property(const nlohmann::json& hessian, const nlohmann::json& units,
                                                    const nlohmann::json& value) {
    auto property = hessian; // its sources
    property["value_type"] = "real";
    if (value.is_null()) {
        property["status"] = "unavailable";
        property["availability"] = {{"hessian", "nuclear response inputs have not met the final request"}};
        return property;
    }
    property["units"] = units;
    property["status"] = "complete";
    property["value"] = value;
    return property;
}

} // namespace molresponse_v3
