// Unit tests for the per-update history row. Plain values only: no World, no
// MPI, no solver. Each block follows one line of the step 2 test list.
#include <madness/chem/ConvergenceDecision.h>
#include <madness/chem/IterationHistory.h>

#include <cstdlib>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

namespace {
using namespace madness::iteration;
using json = nlohmann::json;
int failures = 0;
void check(bool ok, const std::string& message) {
    if (!ok) { ++failures; std::cerr << "FAILED: " << message << '\n'; }
}

const std::map<Metric, std::string> definitions{{Metric::bsh_residual, "t.bsh.v1"},
    {Metric::density_change, "t.density.v1"}, {Metric::energy_change, "t.energy.v1"}};

// Two states after the third update at one level, the second compared with update 2.
IterationMetrics sample() {
    IterationMetrics metrics;
    metrics.info = {"task", {1e-4, 6}, 3, "context"};
    const Comparison previous{"task", {1e-4, 6}, 2, ""};
    for (const std::string id : {"a", "b"}) {
        StateMetrics state{id, previous, {}};
        state.comparison->subject = id;
        state.metrics[Metric::bsh_residual] = MetricValue::measured(id == "a" ? 1e-5 : 3e-4, 1.5);
        state.metrics[Metric::density_change] = MetricValue::measured(2e-4, 0.5);
        state.metrics[Metric::energy_change] = MetricValue::not_applicable("no_energy");
        metrics.states.push_back(state);
    }
    return metrics;
}

// The decision step 1 returns for `sample`: scaled BSH and density, energy written but not applicable.
Decision sample_decision(const IterationMetrics& metrics, std::size_t max_updates = 10, bool final = true) {
    const ModelMetrics model{
        {Metric::bsh_residual, {Support::supported, true, 1e-4, LimitRule::scaled, {}}},
        {Metric::density_change, {Support::supported, true, 1e-4, LimitRule::scaled, {}}},
        {Metric::energy_change, {Support::not_applicable, false, {}, LimitRule::plain, {}}}};
    Request request;
    request.targets = {{Metric::energy_change, 1e-5}};
    request.max_updates = max_updates;
    const auto settings = resolve(request, model, metrics.info.protocol, final);
    return decide(settings, {"task", "context", {"a", "b"}}, {metrics.info.iteration, metrics.info.iteration, true},
                  metrics, final);
}

const json* find_check(const json& row, const std::string& subject, const std::string& metric) {
    for (const auto& entry : row.at("checks"))
        if (entry.at("subject") == subject && entry.at("metric") == metric) return &entry;
    return nullptr;
}

std::filesystem::path scratch_directory() {
    auto base = std::filesystem::temp_directory_path() / ("test_iteration_history." + std::to_string(::getpid()));
    std::filesystem::remove_all(base);
    return base;
}
} // namespace

int main() {
    const auto metrics = sample();
    const auto decision = sample_decision(metrics);
    const auto row = history_row(metrics, decision, definitions);

    // A row holds the existing fields with their existing names.
    check(row.at("task_id") == "task", "task_id");
    check(row.at("protocol") == json{{"threshold", 1e-4}, {"k", 6}}, "protocol");
    check(row.at("iteration") == 3, "iteration");
    check(row.at("metrics") == madness::checkpoint::records::measurements(metrics, definitions),
          "metrics in the form measurements writes");
    // A change names the earlier state it was measured against: its task,
    // level, update number and subject, and no launch.
    {
        const auto& density = row.at("metrics").at(1);
        check(density.at("metric") == "density_change" &&
                  density.at("comparison") == json{{"task_id", "task"}, {"protocol", {{"threshold", 1e-4}, {"k", 6}}},
                                                   {"iteration", 2}, {"subject", "a"}},
              "comparison record holds the task, level, update and subject");
        bool residual_compared = false;
        for (const auto& entry : row.at("metrics"))
            residual_compared = residual_compared || (entry.at("metric") == "bsh_residual" && entry.contains("comparison"));
        check(!residual_compared, "a residual has no comparison");
    }
    check(row.at("decision") == decision.reason && decision.reason == "criteria_unmet", "decision is the reason");

    // action matches the decision for each of the four actions.
    {
        auto copy = decision;
        for (const auto& [action, name] : std::vector<std::pair<Action, std::string>>{
                 {Action::update, "update"}, {Action::finish_protocol, "finish_protocol"},
                 {Action::complete, "complete"}, {Action::stop, "stop"}}) {
            copy.action = action;
            check(history_row(metrics, copy, definitions).at("action") == name, "action " + name);
        }
        check(row.at("action") == "update", "action of the decided sample");
        check(history_row(metrics, sample_decision(metrics, 3), definitions).at("action") == "stop",
              "action at the final iteration limit");
        check(history_row(metrics, sample_decision(metrics, 3, false), definitions).at("action") == "finish_protocol",
              "action at an intermediate iteration limit");
    }

    // checks holds one entry per state and metric, with the limit, the scale
    // and the result taken from the decision.
    {
        check(row.at("checks").size() == decision.checks.size() && decision.checks.size() == 6,
              "one check per state and metric");
        std::set<std::pair<std::string, std::string>> pairs;
        for (const auto& entry : row.at("checks")) pairs.emplace(entry.at("subject"), entry.at("metric"));
        check(pairs.size() == 6, "no duplicate check");
        for (const auto& expected : decision.checks) {
            const auto* entry = find_check(row, expected.subject, metric_name(expected.metric));
            check(entry != nullptr, "check present: " + expected.subject);
            if (!entry) continue;
            check(entry->at("result") == check_result_name(expected.result), "result from the decision");
            check(entry->contains("limit") == expected.limit.has_value() &&
                  (!expected.limit || entry->at("limit") == *expected.limit), "limit from the decision");
            check(entry->contains("scale") == expected.scale.has_value() &&
                  (!expected.scale || entry->at("scale") == *expected.scale), "scale from the decision");
        }
        const auto* a = find_check(row, "a", "bsh_residual");
        const auto* b = find_check(row, "b", "bsh_residual");
        check(a && a->at("result") == "met" && a->at("rule") == "scaled" && a->at("target") == 1e-4 &&
              a->at("limit") == 1e-4 * 2.5 && a->at("scale") == 1.5, "met scaled check values");
        check(b && b->at("result") == "unmet", "unmet check");
    }

    // A metric that is not applicable appears as not applicable, with no limit.
    {
        const auto* entry = find_check(row, "a", "energy_change");
        check(entry && entry->at("result") == "not_applicable" && !entry->contains("limit") &&
              !entry->contains("target") && !entry->contains("scale"), "not applicable check");
    }

    // With no diagnostics, the row has no diagnostics field.
    check(!row.contains("diagnostics"), "no diagnostics field");

    // With diagnostics, each one appears with its name, subject, component and value.
    {
        const std::vector<Diagnostic> values{{"polarizability", "a", "zz", 6.4353},
                                             {"excitation_energy", "b", "", 0.4658}};
        const auto with = history_row(metrics, decision, definitions, values);
        check(with.at("diagnostics").size() == 2, "two diagnostics");
        for (std::size_t i = 0; i < values.size(); ++i) {
            const auto& entry = with.at("diagnostics").at(i);
            check(entry == json{{"name", values[i].name}, {"subject", values[i].subject},
                                {"component", values[i].component}, {"value", values[i].value}},
                  "diagnostic " + values[i].name);
        }
        auto without = with;
        without.erase("diagnostics");
        check(without == row, "diagnostics add one field and change nothing else");
    }

    // Appending two rows gives two lines, each parsing back to the row written;
    // a directory that does not exist is created.
    {
        const auto base = scratch_directory();
        const auto file = base / "states" / "ground_state" / "history.jsonl";
        check(!std::filesystem::exists(file.parent_path()), "directory absent before append");
        const auto second = history_row(metrics, sample_decision(metrics, 3), definitions,
                                        {{"polarizability", "a", "zz", 1.0 / 3.0}});
        append_history(file, row);
        check(std::filesystem::is_directory(file.parent_path()), "append creates the directory");
        append_history(file, second);
        std::ifstream input(file);
        std::vector<std::string> lines;
        for (std::string line; std::getline(input, line);) lines.push_back(line);
        check(lines.size() == 2, "two lines");
        if (lines.size() == 2) {
            check(json::parse(lines[0]) == row, "first line parses back");
            check(json::parse(lines[1]) == second, "second line parses back, values exact");
        }
        std::filesystem::remove_all(base);
    }

    if (failures) { std::cerr << failures << " failure(s)\n"; return 1; }
    std::cout << "test_iteration_history passed\n";
    return 0;
}
