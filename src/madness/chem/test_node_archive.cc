// Unit tests for the node archive. Plain values and a temporary folder: no
// World, no MPI, no solver. Each block follows one line of the step 4a test list.
#include <madness/chem/ConvergenceDecision.h>
#include <madness/chem/NodeArchive.h>

#include <cstdlib>
#include <iostream>
#include <iterator>
#include <set>
#include <stdexcept>
#include <string>
#include <unistd.h>

namespace {
using namespace madness::iteration;
namespace fs = std::filesystem;
using json = nlohmann::json;
int failures = 0;
void check(bool ok, const std::string& message) {
    if (!ok) { ++failures; std::cerr << "FAILED: " << message << '\n'; }
}

std::string read(const fs::path& file) {
    std::ifstream in(file, std::ios::binary);
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}
void write(const fs::path& file, const std::string& text) { std::ofstream(file, std::ios::binary) << text; }
std::set<std::string> names(const fs::path& folder) {
    std::set<std::string> result;
    for (const auto& entry : fs::directory_iterator(folder)) result.insert(entry.path().filename());
    return result;
}
// A writer that puts one file with the given text into the archive folder.
auto writer(const std::string& text) {
    return [text](const fs::path& folder) {
        check(fs::is_directory(folder) && fs::is_empty(folder), "the writer is given an empty folder");
        write(folder / "orbitals.00000", text);
    };
}

const std::map<Metric, std::string> definitions{{Metric::bsh_residual, "t.bsh.v1"},
                                                {Metric::density_change, "t.density.v1"}};

// Two states after the third update at one level, each compared with update 2.
IterationMetrics sample() {
    IterationMetrics metrics;
    metrics.info = {"task", {1e-4, 6}, 3, "context"};
    for (const std::string id : {"a", "b"}) {
        StateMetrics state{id, Comparison{"task", {1e-4, 6}, 2, id}, {}};
        state.metrics[Metric::bsh_residual] = MetricValue::measured(id == "a" ? 1e-5 : 3e-5, 1.5);
        state.metrics[Metric::density_change] = MetricValue::measured(2e-5, 0.5);
        metrics.states.push_back(state);
    }
    return metrics;
}

NodeDescription description(const std::string& version) {
    return {"driven", {1e-4, 6}, 3, json{{"molecule", "h2"}, {"model", "hf"}}, node_metrics(sample(), definitions),
            json{{"ground_state", "gs"}}, "2026-11-23T00:00:00Z", version};
}
} // namespace

int main() {
    const auto root = fs::temp_directory_path() / ("test_node_archive." + std::to_string(::getpid()));
    fs::remove_all(root);
    const auto first = description("v1"), second = description("v2");

    // A folder that does not exist: no node, and inspect returns nothing.
    {
        Node node(root / "missing");
        check(!node.exists() && !node.inspect(), "no node in a missing folder");
        check(!fs::exists(root / "missing"), "inspect creates nothing");
    }

    const auto folder = root / "node";
    Node node(folder);

    // After a save: exists, the same description, the writer's files.
    node.save(first, writer("first"));
    check(node.exists(), "exists after a save");
    check(node.inspect() == first, "inspect returns the saved description");
    check(read(node.archive() / "orbitals.00000") == "first", "the archive holds what the writer wrote");
    check(node.archive() == folder / "archive" && node.history() == folder / "history.jsonl", "fixed paths");

    // A save leaves an existing history file unchanged, byte for byte.
    const std::string rows = "{\"iteration\":1}\n{\"iteration\":2}\n";
    write(node.history(), rows);

    // A second save replaces the first: one archive, one description.
    node.save(second, writer("second"));
    check(node.inspect() == second, "the second description replaces the first");
    check(names(node.archive()) == std::set<std::string>{"orbitals.00000"}
              && read(node.archive() / "orbitals.00000") == "second", "the second archive replaces the first");
    check(names(folder) == std::set<std::string>{"archive", "description.json", "history.jsonl"},
          "only the archive, the description and the history remain");
    check(read(node.history()) == rows, "the history is unchanged after a save");

    // The writer throws during a save: the previous save is intact.
    try {
        node.save(first, [](const fs::path& f) { write(f / "partial", "x"); throw std::runtime_error("disk full"); });
        check(false, "the writer's exception reaches the caller");
    } catch (const std::runtime_error&) {}
    check(node.inspect() == second && read(node.archive() / "orbitals.00000") == "second"
              && names(node.archive()).size() == 1, "a failed save leaves the previous save intact");
    check(names(folder) == std::set<std::string>{"archive", "description.json", "history.jsonl"},
          "a failed save leaves no temporary names");
    check(read(node.history()) == rows, "the history is unchanged after a failed save");

    // Leftovers from an interrupted save are removed and do not change what is found.
    const auto new_record = description_record(first).dump();
    const auto interrupted = [&](bool committed, int step) {
        fs::create_directory(folder / "archive.new");
        write(folder / "archive.new" / "orbitals.00000", "new");
        write(folder / (committed ? "description.json.new" : "description.json.tmp"), committed ? new_record : "{tru");
        if (step >= 1) fs::rename(folder / "archive", folder / "archive.old");
        if (step >= 2) fs::rename(folder / "archive.new", folder / "archive");
    };
    // Killed before the commit point: the old save, untouched.
    interrupted(false, 0);
    check(node.inspect() == second && read(node.archive() / "orbitals.00000") == "second",
          "before the commit point inspect finds the old save");
    check(names(folder) == std::set<std::string>{"archive", "description.json", "history.jsonl"},
          "inspect removes the uncommitted leftovers");
    interrupted(false, 0);
    node.save(second, writer("second"));
    check(node.inspect() == second && names(folder).size() == 3, "save removes the uncommitted leftovers");
    // Killed after the commit point, at each step of the swap: the new save, complete.
    for (int step = 0; step <= 2; ++step) {
        node.save(second, writer("second"));
        interrupted(true, step);
        check(node.inspect() == first && read(node.archive() / "orbitals.00000") == "new"
                  && names(node.archive()).size() == 1, "after the commit point inspect finds the new save, step "
                  + std::to_string(step));
        check(names(folder) == std::set<std::string>{"archive", "description.json", "history.jsonl"},
              "no temporary names remain, step " + std::to_string(step));
    }
    check(read(node.history()) == rows, "the history is unchanged by recovery");

    // A missing or unreadable description: no node, even with an archive folder.
    write(folder / "description.json", "{\"kind\": ");
    check(!node.exists() && !node.inspect(), "an unreadable description means no node");
    write(folder / "description.json", "{\"kind\": \"driven\"}");
    check(!node.exists() && !node.inspect(), "an incomplete description means no node");
    fs::remove(folder / "description.json");
    check(fs::is_directory(node.archive()) && !node.exists() && !node.inspect(),
          "a missing description means no node");

    // stored_metrics meets saved_completion: accepted when met, rejected when tighter.
    {
        const ModelMetric criterion{.support = Support::supported, .default_criterion = true, .default_target = 1e-4};
        const ModelMetrics model{{Metric::bsh_residual, criterion}, {Metric::density_change, criterion}};
        const ExpectedTask task{"task", "context", {"a", "b"}};
        const IterationInfo now{"task", {}, 0, "context"};
        const auto stored = stored_metrics(first, now);
        check(stored.info.task_id == "task" && stored.info.protocol == first.protocol
                  && stored.info.iteration == 3 && stored.states.size() == 2, "stored metrics take the current task");
        const auto met = resolve({.targets = {{Metric::bsh_residual, 1e-4}, {Metric::density_change, 1e-4}}},
                                 model, {1e-4, 6}, true);
        const auto tight = resolve({.targets = {{Metric::bsh_residual, 1e-6}, {Metric::density_change, 1e-4}}},
                                   model, {1e-4, 6}, true);
        check(bool(saved_completion(met, task, stored, true)), "met targets: saved_completion accepts");
        check(!saved_completion(tight, task, stored, true), "tighter targets: saved_completion rejects");
        auto zero = first;
        zero.updates_at_protocol = 0;
        check(!saved_completion(met, task, stored_metrics(zero, now), true),
              "a change at the level with no update before it has no comparison");
    }

    // The description has no target, converged flag, launch identifier or generation.
    node.save(first, writer("first"));
    const auto saved = json::parse(read(folder / "description.json"));
    std::set<std::string> keys;
    for (const auto& [key, value] : saved.items()) keys.insert(key);
    check(keys == std::set<std::string>{"kind", "protocol", "updates_at_protocol", "settings", "metrics",
                                        "lineage", "saved_at", "code_version", "level_finished", "details"},
          "the description's fields");
    const auto text = saved.dump();
    for (const auto* word : {"target", "converged", "attempt", "generation", "launch", "limit", "satisfied"})
        check(text.find(word) == std::string::npos, std::string("no ") + word + " in the description");
    auto unstamped = first;
    unstamped.saved_at.clear();
    node.save(unstamped, writer("first"));
    check(node.inspect()->saved_at.size() == 20, "save stamps the time when none is given");

    fs::remove_all(root);
    if (failures) { std::cerr << failures << " check(s) failed\n"; return 1; }
    std::cout << "test_node_archive passed\n";
    return 0;
}
