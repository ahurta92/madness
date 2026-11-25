#pragma once

// The two files at the top of a run folder, written directly by rank zero of
// the main World: results.json (one entry per application, each written by its
// own application) and timings.json (one entry per launch). No lock: within a
// launch one process writes them; separate launches must not share a folder.
#include <madness/chem/RootCoordination.h>
#include <madness/misc/info.h>
#include <madness/world/thread.h>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <optional>
#include <string>

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

/// Collective. Replace this application's entry in results.json, keyed by its
/// calculation_id, and leave the other entries alone. A file of another
/// version is replaced, not merged.
inline void write_results_entry(World& world, const std::filesystem::path& root, const nlohmann::json& entry) {
    namespace cp = madness::checkpoint;
    const auto provenance = results_provenance(world);
    cp::detail::on_root(world, [&] {
        const auto path = root / "results.json";
        nlohmann::json results{{"version", 3}, {"calculations", nlohmann::json::array()}};
        if (std::filesystem::exists(std::filesystem::symlink_status(path))) {
            auto saved = cp::detail::read_json(path);
            if (saved.value("version", 0) == 3) {
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
        cp::detail::replace_json(path, results);
        return std::string{};
    });
}

/// Collective. Add properties to this application's entry in results.json,
/// which must have been written already.
inline void add_results_properties(World& world, const std::filesystem::path& root, const std::string& id,
                                   const nlohmann::json& properties) {
    namespace cp = madness::checkpoint;
    cp::detail::on_root(world, [&] {
        const auto path = root / "results.json";
        auto results = cp::detail::read_json(path);
        cp::records::require(results.value("version", 0) == 3, "unsupported results file");
        bool found = false;
        for (auto& entry : results.at("calculations"))
            if (entry.at("calculation_id") == id) { entry.at("properties").update(properties); found = true; }
        cp::records::require(found, "results.json has no entry " + id);
        cp::detail::replace_json(path, results);
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
