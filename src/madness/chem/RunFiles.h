#pragma once

// The two files at the top of a run folder, written directly by rank zero of
// the main World: results.json (one entry per application, each written by its
// own application) and timings.json (one entry per launch). No lock: within a
// launch one process writes them; separate launches must not share a folder.
#include <madness/chem/ConvergenceDecision.h>
#include <madness/chem/DerivedResults.h>
#include <madness/chem/RootCoordination.h>
#include <madness/chem/RunSummary.h>
#include <madness/misc/info.h>
#include <madness/world/thread.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace madness::iteration {

/// When this process started, as UTC text. A launch is identified by it.
inline const std::chrono::system_clock::time_point launch_clock = std::chrono::system_clock::now();
inline std::string launch_started_at() {
    const auto time = std::chrono::system_clock::to_time_t(launch_clock);
    std::tm utc{};
    gmtime_r(&time, &utc);
    char text[32];
    std::strftime(text, sizeof text, "%Y-%m-%dT%H:%M:%S", &utc);
    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(launch_clock.time_since_epoch()).count() % 1000;
    char fraction[8];
    std::snprintf(fraction, sizeof fraction, ".%03dZ", static_cast<int>(ms));
    return std::string(text) + fraction;
}

/// Collective. A run folder written before the states moved to node folders
/// has a checkpoints/ folder; it is refused, never read.
inline void require_current_run_folder(World& world, const std::filesystem::path& root) {
    madness::checkpoint::detail::on_root(world, [&] {
        if (std::filesystem::exists(std::filesystem::symlink_status(root / "checkpoints")))
            throw std::runtime_error("the run folder " + root.string() +
                                     " was written by an older version of madqc; use a new output prefix");
        return std::string{};
    });
}

/// The archive I/O backend this launch writes saved states with, for the
/// provenance in results.json. madqc sets it when a deck turns HDF5 on.
inline std::string& archive_io_backend() {
    static std::string backend = "native";
    return backend;
}

/// Who wrote results.json and how it ran: replaced on each write, so it
/// describes the launch that wrote the file last.
inline nlohmann::json results_provenance(World& world) {
    return {{"code_version", madness::info::git_commit()}, {"launched_at", launch_started_at()},
            {"mpi_processes", world.size()}, {"threads_per_process", ThreadPool::size() + 1},
            {"archive_io", archive_io_backend()}};
}

inline constexpr int results_version = 5;

/// Rank zero. Write `results`, with its derived quantities computed afresh
/// (add_derived_quantities), to results.json and its text summary to
/// summary.txt beside it, so the two never disagree.
inline void replace_results_files(const std::filesystem::path& root, nlohmann::json results) {
    namespace cp = madness::checkpoint;
    add_derived_quantities(results);
    cp::detail::replace_json(root / "results.json", results);
    const auto path = root / "summary.txt";
    const auto temp = root / ".summary.txt.tmp";
    if (std::filesystem::is_regular_file(std::filesystem::symlink_status(temp))) std::filesystem::remove(temp);
    cp::detail::write_exclusive(temp, results_summary(results));
    cp::records::require(!std::filesystem::is_symlink(std::filesystem::symlink_status(path)),
                         "symlink record destination");
    std::filesystem::rename(temp, path);
}

/// The final-level target of each metric, as the convergence decision module
/// resolved it: {metric name: target}.
inline nlohmann::json results_targets(const ResolvedSettings& settings) {
    nlohmann::json targets = nlohmann::json::object();
    for (const auto& criterion : settings.criteria) targets[metric_name(criterion.metric)] = criterion.effective;
    return targets;
}

/// The key under which a writer lists the states a value was computed from,
/// while the entry is assembled. It never reaches the file: each item is
/// `{"node": folder name}` (every subject of that node) or
/// `{"node": folder name, "subject": id}` (one subject, e.g. an excited root).
inline constexpr const char* results_sources_key = "sources";

/// Add a source to a value's list, once.
inline void add_results_source(nlohmann::json& value, const std::string& node, const std::string& subject = {}) {
    auto& sources = value[results_sources_key];
    if (!sources.is_array()) sources = nlohmann::json::array();
    nlohmann::json source{{"node", node}};
    if (!subject.empty()) source["subject"] = subject;
    if (std::find(sources.begin(), sources.end(), source) == sources.end()) sources.push_back(std::move(source));
}

/// The metrics a value's `convergence` reports when its entry names no
/// targets.
inline const std::vector<std::string>& convergence_metrics() {
    static const std::vector<std::string> names{"bsh_residual", "density_change"};
    return names;
}

/// The metrics a `convergence` reports: every metric with a target in the
/// entry's `targets`, flat or named by kind (driven_states, excited_states).
/// A metric a value's states did not measure is still left out of its
/// convergence (a driven state has no energy_change).
inline std::vector<std::string> results_convergence_metrics(const nlohmann::json& targets) {
    std::vector<std::string> names;
    const auto add = [&](const std::string& name) {
        if (std::find(names.begin(), names.end(), name) == names.end()) names.push_back(name);
    };
    if (targets.is_object())
        for (const auto& [key, target] : targets.items()) {
            if (target.is_number()) add(key);
            else if (target.is_object())
                for (const auto& [metric, value] : target.items())
                    if (value.is_number()) add(metric);
        }
    return names.empty() ? convergence_metrics() : names;
}

/// The worst measured value of each of `metrics` over `sources`, from the
/// entry's state_metrics. A metric no source measured is left out.
inline nlohmann::json worst_convergence(const nlohmann::json& state_metrics, const nlohmann::json& sources,
                                        const std::vector<std::string>& metrics = convergence_metrics()) {
    nlohmann::json worst = nlohmann::json::object();
    for (const auto& source : sources) {
        const auto node = state_metrics.find(source.at("node").get<std::string>());
        if (node == state_metrics.end()) continue;
        const auto& subjects = node->at("subjects");
        for (const auto& [subject, measured] : subjects.items()) {
            if (source.contains("subject") && source.at("subject") != subject) continue;
            const auto& values = measured.at("values");
            for (const auto& metric : metrics) {
                const auto value = values.find(metric);
                if (value == values.end()) continue;
                if (!worst.contains(metric) || value->get<double>() > worst.at(metric).get<double>())
                    worst[metric] = *value;
            }
        }
    }
    return worst;
}

/// `x` rounded to `digits` significant figures. Convergence and targets are
/// written this way; state_metrics and property values keep full precision.
inline double round_significant(double x, int digits = 3) {
    if (x == 0.0 || !std::isfinite(x)) return x;
    char text[32];
    std::snprintf(text, sizeof text, "%.*e", digits - 1, x);
    return std::strtod(text, nullptr);
}

/// Every floating-point number in `j`, at any depth, rounded to three
/// significant figures.
inline nlohmann::json round_results_numbers(nlohmann::json j) {
    std::vector<nlohmann::json*> pending{&j};
    while (!pending.empty()) {
        auto* node = pending.back();
        pending.pop_back();
        if (node->is_number_float()) *node = round_significant(node->get<double>());
        else if (node->is_structured())
            for (auto& element : *node) pending.push_back(&element);
    }
    return j;
}

/// The written convergence of a value or entry: the worst measured values
/// over its sources, to three significant figures.
inline nlohmann::json results_convergence(const nlohmann::json& state_metrics, const nlohmann::json& value,
                                          const std::vector<std::string>& metrics) {
    return round_results_numbers(
        worst_convergence(state_metrics, value.value(results_sources_key, nlohmann::json::array()), metrics));
}

/// Turn one value (a single-answer property or an item of `values`) into its
/// written form: its sources become `convergence` (unless the entry carries
/// the convergence for all its values), and an empty availability is dropped.
inline void finish_results_value(nlohmann::json& value, const nlohmann::json& state_metrics,
                                 bool with_convergence = true,
                                 const std::vector<std::string>& metrics = convergence_metrics()) {
    if (with_convergence) value["convergence"] = results_convergence(state_metrics, value, metrics);
    value.erase(results_sources_key);
    if (value.contains("availability") && value.at("availability").empty()) value.erase("availability");
}

/// The written form of an entry's properties. With `with_convergence` false
/// the entry carries one convergence for every value, and the properties none.
/// `metrics` are those the convergence reports (results_convergence_metrics).
inline nlohmann::json finish_results_properties(nlohmann::json properties, const nlohmann::json& state_metrics,
                                                bool with_convergence = true,
                                                const std::vector<std::string>& metrics = convergence_metrics()) {
    for (auto& [name, property] : properties.items()) {
        if (property.contains("values")) {
            for (auto& item : property.at("values"))
                finish_results_value(item, state_metrics, with_convergence, metrics);
            property.erase(results_sources_key);
            if (property.contains("availability") && property.at("availability").empty()) property.erase("availability");
        } else {
            finish_results_value(property, state_metrics, with_convergence, metrics);
        }
    }
    return properties;
}

/// Collective. Replace this application's entry in results.json, keyed by its
/// calculation_id, and leave the other entries alone. A file of another
/// version is replaced, not merged. Each property's sources become its
/// convergence (finish_results_properties). An entry that lists sources itself
/// (every value from the same states, as the ground state) carries one
/// convergence instead, and its properties none. A convergence reports every
/// metric that has a target. Convergence and targets are written to three
/// significant figures.
inline void write_results_entry(World& world, const std::filesystem::path& root, nlohmann::json entry) {
    namespace cp = madness::checkpoint;
    const auto provenance = results_provenance(world);
    const auto state_metrics = entry.value("state_metrics", nlohmann::json::object());
    const auto metrics = results_convergence_metrics(entry.value("targets", nlohmann::json::object()));
    const bool entry_convergence = entry.contains(results_sources_key);
    if (entry_convergence) {
        entry["convergence"] = results_convergence(state_metrics, entry, metrics);
        entry.erase(results_sources_key);
    }
    if (entry.contains("properties"))
        entry["properties"] =
            finish_results_properties(entry.at("properties"), state_metrics, !entry_convergence, metrics);
    if (entry.contains("targets")) entry["targets"] = round_results_numbers(entry.at("targets"));
    cp::detail::on_root(world, [&] {
        const auto path = root / "results.json";
        nlohmann::json results{{"version", results_version}, {"calculations", nlohmann::json::array()}};
        if (std::filesystem::exists(std::filesystem::symlink_status(path))) {
            auto saved = cp::detail::read_json(path);
            if (saved.value("version", 0) == results_version) {
                cp::records::require(saved.at("calculations").is_array(), "unsupported results file");
                results = std::move(saved);
            }
        }
        const auto id = cp::records::id(entry.at("calculation_id"));
        nlohmann::json calculations = nlohmann::json::array();
        for (auto& other : results.at("calculations"))
            if (other.at("calculation_id") != id) calculations.push_back(std::move(other));
        calculations.push_back(entry);
        results["calculations"] = std::move(calculations);
        results["provenance"] = provenance;
        replace_results_files(root, results);
        return std::string{};
    });
}

/// Collective. Add properties to this application's entry in results.json,
/// which must have been written already. Their sources become convergence
/// from that entry's state_metrics, unless the entry carries one convergence
/// for all its values.
inline void add_results_properties(World& world, const std::filesystem::path& root, const std::string& id,
                                   const nlohmann::json& properties) {
    namespace cp = madness::checkpoint;
    cp::detail::on_root(world, [&] {
        const auto path = root / "results.json";
        auto results = cp::detail::read_json(path);
        cp::records::require(results.value("version", 0) == results_version, "unsupported results file");
        bool found = false;
        for (auto& entry : results.at("calculations"))
            if (entry.at("calculation_id") == id) {
                entry.at("properties").update(finish_results_properties(
                    properties, entry.value("state_metrics", nlohmann::json::object()), !entry.contains("convergence"),
                    results_convergence_metrics(entry.value("targets", nlohmann::json::object()))));
                found = true;
            }
        cp::records::require(found, "results.json has no entry " + id);
        replace_results_files(root, results);
        return std::string{};
    });
}

/// Collective. Add one application run's timing record to this launch's entry
/// in timings.json and return the run's part number in that entry. The first
/// record of a launch appends the entry. Each run's top-level steps (scf,
/// response) go under the entry's top-level step, tagged with the part number.
/// A run publishes a snapshot of its whole record each time: passing the part
/// back replaces that run's steps; a new run (no part) appends after the others.
/// With more than one part, the top-level wall interval is the time from the
/// launch start to this write; a failure is kept on the entry.
inline std::size_t write_timing_entry(World& world, const std::filesystem::path& root, const nlohmann::json& entry,
                                      std::optional<std::size_t> part = std::nullopt) {
    namespace cp = madness::checkpoint;
    using json = nlohmann::json;
    const auto written = cp::detail::on_root(world, [&] {
        cp::records::require(entry.at("started_at").is_string() && entry.at("execution").is_object() &&
                             entry.at("timing").is_object(), "invalid timing entry");
        const auto path = root / "timings.json";
        auto timings = std::filesystem::exists(std::filesystem::symlink_status(path))
            ? cp::detail::read_json(path) : json{{"version", 2}, {"launches", json::array()}};
        cp::records::require(timings.at("version") == 2 && timings.at("launches").is_array(),
                             "unsupported timings file");
        auto& launches = timings.at("launches");
        const bool same_launch = !launches.empty() && launches.back().at("started_at") == entry.at("started_at");
        cp::records::require(same_launch || !part, "a timing part belongs to another launch");
        std::size_t number = 0;
        json steps = entry.at("timing").value("children", json::array());
        if (!same_launch) {
            for (auto& step : steps) step["part"] = number;
            auto launch = entry;
            launch["timing"]["children"] = std::move(steps);
            launches.push_back(std::move(launch));
        } else {
            auto& launch = launches.back();
            auto& tree = launch.at("timing");
            auto& children = tree["children"];
            if (!children.is_array()) children = json::array();
            std::size_t position = children.size(), parts = 0;
            if (part) {
                number = *part;
                json kept = json::array();
                for (auto& child : children) {
                    if (child.at("part") == number) position = std::min(position, kept.size());
                    else kept.push_back(std::move(child));
                }
                children = std::move(kept);
                position = std::min(position, children.size());
            } else {
                for (const auto& child : children) number = std::max(number, child.at("part").get<std::size_t>() + 1);
            }
            for (auto& step : steps) step["part"] = number;
            children.insert(children.begin() + static_cast<std::ptrdiff_t>(position), steps.begin(), steps.end());
            for (const auto& child : children) parts = std::max(parts, child.at("part").get<std::size_t>() + 1);
            tree["seconds"] = parts > 1
                ? std::chrono::duration<double>(std::chrono::system_clock::now() - launch_clock).count()
                : entry.at("timing").at("seconds").get<double>();
            if (entry.contains("failure")) launch["failure"] = entry.at("failure");
        }
        cp::detail::replace_json(path, timings);
        return std::to_string(number);
    });
    return std::stoul(written);
}
} // namespace madness::iteration
