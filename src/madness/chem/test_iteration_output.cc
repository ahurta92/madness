// Unit tests for step 3: the timing record of one launch, the verbosity type,
// and the printing rules that had no test. Plain values only: no World, no MPI,
// no solver.
#include <madness/chem/ConvergenceDecision.h>
#include <madness/chem/IterationHistory.h>
#include <madness/chem/TimingRecords.h>

#include <algorithm>
#include <cctype>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
using namespace madness::iteration;
using json = nlohmann::json;
int failures = 0;
void check(bool ok, const std::string& message) {
    if (!ok) { ++failures; std::cerr << "FAILED: " << message << '\n'; }
}
std::string lower(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c) { return std::tolower(c); });
    return text;
}

const std::map<Metric, std::string> definitions{{Metric::bsh_residual, "t.bsh.v1"},
    {Metric::density_change, "t.density.v1"}, {Metric::energy_change, "t.energy.v1"}};

// Two states after the second update; state b has no BSH measurement.
IterationMetrics sample() {
    IterationMetrics metrics;
    metrics.info = {"task", {1e-4, 6}, 2, "context"};
    const Comparison previous{"task", {1e-4, 6}, 1, ""};
    for (const std::string id : {"a", "b"}) {
        StateMetrics state{id, previous, {}};
        state.comparison->subject = id;
        state.metrics[Metric::bsh_residual] = id == "a" ? MetricValue::measured(3e-2, 1.5)
                                                        : MetricValue::unavailable("not_measured");
        state.metrics[Metric::density_change] = MetricValue::measured(2e-2, 0.5);
        state.metrics[Metric::energy_change] = MetricValue::not_applicable("no_energy");
        metrics.states.push_back(state);
    }
    return metrics;
}

// The decision for `sample` at its iteration limit, final level.
Decision limit_decision(const IterationMetrics& metrics) {
    const ModelMetrics model{
        {Metric::bsh_residual, {Support::supported, true, 1e-4, LimitRule::plain, {}}},
        {Metric::density_change, {Support::supported, true, 1e-4, LimitRule::scaled, {}}},
        {Metric::energy_change, {Support::not_applicable, false, {}, LimitRule::plain, {}}}};
    Request request;
    request.max_updates = 2;
    const auto settings = resolve(request, model, metrics.info.protocol, true);
    return decide(settings, {"task", "context", {"a", "b"}}, {2, 2, true}, metrics, true);
}

// Everything a print function could change, in a comparable form.
json snapshot(const IterationMetrics& metrics, const Decision& decision) {
    auto row = history_row(metrics, decision, definitions);
    row["checkpoint"] = decision.checkpoint;
    row["state_count"] = metrics.states.size();
    return row;
}
} // namespace

int main() {
    // timing_entry gives the record the ground state and response built by hand.
    {
        const json tree{{"name", "run"}, {"seconds", 1.5}, {"children", json::array()}};
        const int ranks = 2, threads = 9; // World::size() and ThreadPool::size() + 1 are int
        json by_hand{{"started_at", "2026-11-23T20:01:05.123Z"},
            {"execution", {{"mpi_ranks", ranks}, {"threads_per_rank", threads}}},
            {"timing", tree}};
        const auto record = timing_entry("2026-11-23T20:01:05.123Z", {2, 9}, tree);
        check(record == by_hand, "timing_entry equals the hand-built record");
        check(record.dump() == by_hand.dump(), "timing_entry prints as the hand-built record");
        std::vector<std::string> keys;
        for (const auto& [key, value] : record.items()) keys.push_back(key);
        check(keys == std::vector<std::string>{"execution", "started_at", "timing"}, "field names");
        check(record.at("execution").size() == 2, "execution has two fields");

        // A failure appears when given, and is absent when not.
        by_hand["failure"] = "solver threw";
        const auto failed = timing_entry("2026-11-23T20:01:05.123Z", {2, 9}, tree, "solver threw");
        check(failed == by_hand, "failure recorded under the hand-built name");
        check(!record.contains("failure"), "no failure field for a normal end");
        check(!timing_entry("2026-11-23T20:01:05.123Z", {}, tree, "").contains("failure"), "empty failure is absent");
        check(timing_entry("x", {}, tree).at("execution") == json{{"mpi_ranks", 1}, {"threads_per_rank", 1}},
              "default execution is one rank, one thread");
    }

    // prints: true at and above the needed level, false below.
    {
        const std::vector<Verbosity> levels{Verbosity::silent, Verbosity::progress,
                                            Verbosity::breakdown, Verbosity::detail};
        for (std::size_t c = 0; c < levels.size(); ++c)
            for (std::size_t n = 0; n < levels.size(); ++n)
                check(prints(levels[c], levels[n]) == (c >= n),
                      "prints(" + std::to_string(c) + ", " + std::to_string(n) + ")");
    }

    // The one mapping from print_level, for 0, 1, 2, 3 and 10 (and outside them).
    {
        check(verbosity(0) == Verbosity::silent, "print_level 0");
        check(verbosity(1) == Verbosity::progress, "print_level 1");
        check(verbosity(2) == Verbosity::breakdown, "print_level 2");
        check(verbosity(3) == Verbosity::detail, "print_level 3");
        check(verbosity(10) == Verbosity::detail, "print_level 10");
        check(verbosity(-1) == Verbosity::silent, "print_level -1");
        check(verbosity(4) == Verbosity::detail, "print_level 4");
        // The table reproduces the guards each application used at the step 2 tag.
        for (int p = -1; p <= 12; ++p) {
            const auto v = verbosity(p);
            const std::string at = " at print_level " + std::to_string(p);
            // Ground state: print_level > 0, > 1, > 2. Response: PrintLevel >= Normal, i.e. print_level >= 1.
            for (const auto level : {needed::iteration_metrics, needed::iteration_timing, needed::checkpoint_saved,
                                     needed::restart_summary, needed::checkpoint_timing})
                check(prints(v, level) == (p > 0), "progress output" + at);
            check(prints(v, needed::energy_breakdown) == (p > 1), "energy breakdown" + at);
            check(prints(v, needed::final_orbital_timing) == (p > 2), "final orbital timing" + at);
            // Warnings print at every level (SCF at print_level 0 is the one intended change).
            check(prints(v, needed::targets_below_floor), "below-floor warning" + at);
            check(prints(v, needed::unconverged_diagnostics), "NOT CONVERGED block" + at);
        }
    }

    const auto metrics = sample();
    auto decision = limit_decision(metrics);
    check(decision.action == Action::stop && decision.reason == "iteration_limit", "sample stops at the limit");
    const std::vector<Diagnostic> diagnostics{{"polarizability", "b", "zz", 3.84}};

    // A value that is unavailable is never printed as zero, in the NOT CONVERGED block,
    // and the block gives its reason.
    {
        std::ostringstream block;
        write_unconverged_diagnostics(block, metrics, decision, diagnostics);
        const auto text = block.str();
        std::istringstream lines(text);
        std::string line;
        bool seen = false;
        while (std::getline(lines, line))
            if (line.find("BSH residual") != std::string::npos) {
                seen = true;
                check(line.find(" - ") != std::string::npos, "unavailable BSH shows a dash");
                check(line.find("0.0e+00") == std::string::npos, "unavailable BSH is not zero");
                check(line.find("unavailable") != std::string::npos, "unavailable BSH says unavailable");
                check(line.find("unavailable (not measured)") != std::string::npos,
                      "unavailable BSH gives its reason");
            }
        check(seen, "the block lists the BSH check of state b");
    }

    // "Checkpoint saved" comes only from write_checkpoint_saved, never from the
    // table or the restart summary, whatever the decision asks for.
    {
        for (const bool checkpoint : {false, true})
            for (const auto action : {Action::update, Action::finish_protocol, Action::complete, Action::stop}) {
                auto asked = decision;
                asked.checkpoint = checkpoint;
                asked.action = action;
                std::ostringstream table;
                write_iteration_metrics(table, "SCF", metrics, asked);
                write_unconverged_diagnostics(table, metrics, asked, diagnostics);
                Assessment assessment{true, false, "", asked.checks};
                write_restart_summary(table, "SCF", metrics, metrics.info.protocol, assessment, false);
                write_restart_summary(table, "SCF", metrics, metrics.info.protocol, assessment, true);
                check(lower(table.str()).find("checkpoint saved") == std::string::npos,
                      "no checkpoint notice outside write_checkpoint_saved");
            }
        std::ostringstream notice;
        write_checkpoint_saved(notice, decision);
        check(notice.str().rfind("Checkpoint saved.", 0) == 0, "the notice says saved");
    }

    // The printer computes nothing: no print function changes a value, a limit
    // or a decision, and the iteration table shows the decision's own limits.
    {
        const auto before = snapshot(metrics, decision);
        const Assessment assessment{true, false, "", decision.checks};
        std::ostringstream sink;
        write_iteration_metrics(sink, "response", metrics, decision);
        write_metric_table(sink, metrics, decision.checks);
        write_unconverged_diagnostics(sink, metrics, decision, diagnostics);
        write_checkpoint_saved(sink, decision);
        write_restart_summary(sink, "response", metrics, {1e-5, 8}, assessment, false);
        write_iteration_timing(sink, "response", metrics.info, {{{"evaluate", 0.25}}, 0.5});
        write_checkpoint_timing(sink, "response", "task", 0.125);
        write_value_table(sink, "values", {{"x", 1.0}});
        check(snapshot(metrics, decision) == before, "printing leaves metrics and decision unchanged");
        for (const auto& c : decision.checks)
            if (c.limit) {
                std::ostringstream text;
                text << std::scientific << std::setprecision(4) << *c.limit;
                check(sink.str().find(text.str()) != std::string::npos, "the decision's limit is printed");
            }
    }

    if (failures) { std::cerr << failures << " check(s) failed\n"; return 1; }
    std::cout << "test_iteration_output: all checks passed\n";
    return 0;
}
