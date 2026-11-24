#pragma once

// Per-update history: one JSON row per decision, in one format for the ground
// state, driven response and excited states. It records; it never measures or
// judges. No MADNESS runtime or MPI: the caller chooses the writing process.
#include <madness/chem/MetricRecords.h>
#include <madness/chem/IterationMetrics.h>

#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>

namespace madness::iteration {

// A cheap property value reported for information. It is never a result.
struct Diagnostic {
    std::string name;      // for example "polarizability"
    std::string subject;   // the state it was computed from
    std::string component; // for example "zz", or empty
    double value = 0.0;
};

inline const char* action_name(Action action) {
    switch (action) {
    case Action::update: return "update";
    case Action::finish_protocol: return "finish_protocol";
    case Action::complete: return "complete";
    case Action::stop: return "stop";
    }
    throw std::invalid_argument("unknown iteration action");
}

inline const char* check_result_name(CheckResult result) {
    switch (result) {
    case CheckResult::met: return "met";
    case CheckResult::unmet: return "unmet";
    case CheckResult::unavailable: return "unavailable";
    case CheckResult::not_applicable: return "not_applicable";
    case CheckResult::invalid: return "invalid";
    }
    throw std::invalid_argument("unknown check result");
}

// One check exactly as the decision holds it. `target` is the effective target
// of the level; `limit` is what the value was compared with (target, or target
// times 1 + scale under the scaled rule). Absent fields were not available.
inline nlohmann::json check_record(const MetricCheck& check) {
    nlohmann::json entry{{"subject", check.subject}, {"metric", metric_name(check.metric)},
                         {"result", check_result_name(check.result)}};
    if (check.result != CheckResult::not_applicable) {
        entry["rule"] = check.relative > 0.0 ? "scaled" : "plain";
        entry["target"] = check.absolute;
    }
    if (check.limit) entry["limit"] = *check.limit;
    if (check.scale) entry["scale"] = *check.scale;
    if (!check.reason.empty()) entry["reason"] = check.reason;
    return entry;
}

// Build one row from what the solver already has after an update. The first
// five fields are the response history's existing fields, unchanged.
inline nlohmann::json history_row(const IterationMetrics& metrics, const Decision& decision,
                                  const std::map<Metric, std::string>& definitions,
                                  const std::vector<Diagnostic>& diagnostics = {}) {
    namespace records = madness::checkpoint::records;
    nlohmann::json row{{"task_id", metrics.info.task_id},
                       {"protocol", records::protocol_record(metrics.info.protocol)},
                       {"iteration", metrics.info.iteration},
                       {"metrics", records::measurements(metrics, definitions)},
                       {"decision", decision.reason},
                       {"action", action_name(decision.action)},
                       {"checks", nlohmann::json::array()}};
    for (const auto& check : decision.checks) row["checks"].push_back(check_record(check));
    if (!diagnostics.empty()) {
        auto& values = row["diagnostics"] = nlohmann::json::array();
        for (const auto& diagnostic : diagnostics)
            values.push_back({{"name", diagnostic.name}, {"subject", diagnostic.subject},
                              {"component", diagnostic.component}, {"value", diagnostic.value}});
    }
    return row;
}

// Append one row as one line. Creates the directory if needed.
inline void append_history(const std::filesystem::path& file, const nlohmann::json& row) {
    if (file.has_parent_path()) std::filesystem::create_directories(file.parent_path());
    std::ofstream history(file, std::ios::app);
    history.exceptions(std::ios::failbit | std::ios::badbit);
    history << row.dump() << '\n';
    history.close();
}

} // namespace madness::iteration
