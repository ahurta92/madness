#pragma once

// One node's saved states in one fixed folder:
//   <folder>/description.json  what a restart needs; measurements, never verdicts
//   <folder>/history.jsonl     per-update rows; appended by the caller, never touched here
//   <folder>/archive/          the functions; written and read by the caller
// A save replaces the archive and the description together. No MADNESS runtime,
// MPI or lock: the caller chooses the one process that touches the folder.
// Known limit: no fsync; a save survives the process dying, not a power loss.
#include <madness/chem/CheckpointRecords.h>
#include <madness/chem/IterationMetrics.h>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <functional>
#include <nlohmann/json.hpp>
#include <optional>

namespace madness::iteration {

// What a restart needs to know about a saved node. No target, no converged
// flag, no launch identifier, no generation.
struct NodeDescription {
    std::string kind;                    // ground_state, driven, excited_set, second_order
    Protocol protocol;                   // the level the states were saved at
    std::size_t updates_at_protocol = 0;
    nlohmann::json settings;             // what must match: the molecule and the model, never tolerances
    nlohmann::json metrics;              // node_metrics(...) of the saved states
    nlohmann::json lineage;              // what the node was computed from; empty for the ground state
    std::string saved_at;                // UTC; filled by save when empty
    std::string code_version;
    bool operator==(const NodeDescription&) const = default;
};

inline bool node_kind(const std::string& kind) {
    return kind == "ground_state" || kind == "driven" || kind == "excited_set" || kind == "second_order";
}

inline nlohmann::json description_record(const NodeDescription& d) {
    return {{"kind", d.kind}, {"protocol", madness::checkpoint::records::protocol_record(d.protocol)},
            {"updates_at_protocol", d.updates_at_protocol}, {"settings", d.settings}, {"metrics", d.metrics},
            {"lineage", d.lineage}, {"saved_at", d.saved_at}, {"code_version", d.code_version}};
}

// Throws on a malformed record.
inline NodeDescription read_description(const nlohmann::json& r) {
    namespace records = madness::checkpoint::records;
    NodeDescription d{r.at("kind"), records::read_protocol(r.at("protocol")),
                      records::count(r.at("updates_at_protocol")), r.at("settings"), r.at("metrics"),
                      r.at("lineage"), r.at("saved_at"), r.at("code_version")};
    records::require(node_kind(d.kind), "unknown node kind: " + d.kind);
    records::require(d.metrics.is_array(), "node metrics must be an array");
    return d;
}

// The measurements in the form records::measurements writes, with each
// comparison reduced to the level and subject it compared with: the task and
// launch it came from are not stored, and stored_metrics rebinds them.
inline nlohmann::json node_metrics(const IterationMetrics& metrics,
                                   const std::map<Metric, std::string>& definitions) {
    auto entries = madness::checkpoint::records::measurements(metrics, definitions);
    for (auto& entry : entries)
        if (entry.contains("comparison")) {
            const auto& c = entry.at("comparison");
            entry["comparison"] = nlohmann::json{{"protocol", c.at("protocol")}, {"subject", c.at("subject")}};
        }
    return entries;
}

// The stored measurements as measured values of `as_task` (its task, launch and
// context), at the saved level after `updates_at_protocol` updates. A stored
// comparison at the same level is taken as the update before; with no update
// at the level it cannot be, and the change check then finds no comparison.
inline IterationMetrics stored_metrics(const NodeDescription& description, const IterationInfo& as_task) {
    namespace records = madness::checkpoint::records;
    IterationMetrics result;
    result.info = {as_task.task_id, as_task.attempt_id, description.protocol, description.updates_at_protocol,
                   as_task.context_id};
    const auto& info = result.info;
    for (const auto& entry : description.metrics) {
        const auto id = entry.at("subject").get<std::string>();
        auto state = std::find_if(result.states.begin(), result.states.end(),
                                  [&](const StateMetrics& s) { return s.id == id; });
        if (state == result.states.end()) state = result.states.insert(result.states.end(), StateMetrics{id, {}, {}});
        const auto status = entry.at("status").get<std::string>();
        MetricValue m = status == "measured"       ? MetricValue::measured(entry.at("value").get<double>())
                      : status == "invalid"        ? MetricValue::invalid(entry.at("reason"))
                      : status == "not_applicable" ? MetricValue::not_applicable(entry.at("reason"))
                                                   : MetricValue::unavailable(entry.at("reason"));
        if (entry.contains("scale")) {
            const auto& scale = entry.at("scale");
            if (scale.at("status") == "measured") m.scale = Scale::measured(scale.at("value"));
            else m.scale.reason = scale.at("reason");
        }
        if (entry.contains("comparison")) {
            const auto& c = entry.at("comparison");
            const auto level = records::read_protocol(c.at("protocol"));
            const auto before = level != info.protocol ? 0 : info.iteration > 0 ? info.iteration - 1 : info.iteration;
            m.comparison = Comparison{info.task_id, info.attempt_id, level, before, c.at("subject").get<std::string>()};
        }
        const auto e = entry.value("evidence", nlohmann::json::object());
        const auto number = [&](const char* key) {
            return e.contains(key) ? std::optional<double>(e.at(key).get<double>()) : std::nullopt;
        };
        m.evidence = {number("current"), number("previous"), number("current_norm"), number("previous_norm"), {},
                      e.contains("limiting_id") ? std::optional(e.at("limiting_id").get<std::string>()) : std::nullopt};
        for (const auto& c : e.value("components", nlohmann::json::array()))
            m.evidence.components.push_back({c.at("id"), c.at("residual_norm"),
                c.contains("function_norm") ? std::optional(c.at("function_norm").get<double>()) : std::nullopt});
        state->metrics[metric_from_name(entry.at("metric"))] = std::move(m);
    }
    return result;
}

class Node {
public:
    explicit Node(std::filesystem::path folder) : folder_(std::move(folder)) {}

    // True when a readable description and its archive folder are present.
    bool exists() const { return inspect().has_value() && std::filesystem::is_directory(archive()); }

    // Finishes or discards an interrupted save, then reads the description file
    // and nothing else. Nothing when it is missing or unreadable.
    std::optional<NodeDescription> inspect() const {
        recover();
        try {
            std::ifstream file(folder_ / "description.json");
            return read_description(nlohmann::json::parse(file));
        } catch (const std::exception&) {
            return std::nullopt;
        }
    }

    // Replace the saved states. `write_archive` is given an empty folder; if it
    // throws, the previous save is intact. The complete new description, under
    // its temporary name, is the commit point: from then on the save finishes,
    // here or in the next inspect.
    void save(NodeDescription description,
              const std::function<void(const std::filesystem::path&)>& write_archive) {
        namespace fs = std::filesystem;
        madness::checkpoint::records::require(node_kind(description.kind), "unknown node kind: " + description.kind);
        if (description.saved_at.empty()) description.saved_at = utc_now();
        fs::create_directories(folder_);
        recover();
        const auto fresh = folder_ / "archive.new", staged = folder_ / "description.json.tmp";
        try {
            fs::create_directory(fresh);
            write_archive(fresh);
            std::ofstream file(staged);
            file.exceptions(std::ios::failbit | std::ios::badbit);
            file << description_record(description).dump(2) << '\n';
            file.close();
        } catch (...) {
            std::error_code ignored;
            fs::remove_all(fresh, ignored), fs::remove(staged, ignored);
            throw;
        }
        fs::rename(staged, folder_ / "description.json.new");
        recover();
    }

    std::filesystem::path archive() const { return folder_ / "archive"; }
    std::filesystem::path history() const { return folder_ / "history.jsonl"; }

private:
    // With a committed new description: swap in the new archive, then the
    // description. Then remove every temporary name. Never touches history.
    void recover() const {
        namespace fs = std::filesystem;
        const auto fresh = folder_ / "archive.new", old = folder_ / "archive.old",
                   committed = folder_ / "description.json.new";
        if (fs::exists(committed)) {
            if (fs::is_directory(fresh)) {
                fs::remove_all(old);
                if (fs::exists(archive())) fs::rename(archive(), old);
                fs::rename(fresh, archive());
            }
            fs::rename(committed, folder_ / "description.json");
        }
        fs::remove_all(old);
        fs::remove_all(fresh);
        fs::remove(folder_ / "description.json.tmp");
    }

    static std::string utc_now() {
        const auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        std::tm utc{};
        gmtime_r(&now, &utc);
        char text[32];
        std::strftime(text, sizeof text, "%Y-%m-%dT%H:%M:%SZ", &utc);
        return text;
    }

    std::filesystem::path folder_;
};

} // namespace madness::iteration
