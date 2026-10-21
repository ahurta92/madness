// Pure scalar diagnostics: no World, numerical solver or external fixture.
#include <madness/chem/FunctionMetricSummary.h>
#include <madness/chem/ExplicitCriteria.h>
#include "../solvers/function_metric_diagnostics.hpp"
#include "../solvers/iterate.hpp"
#include "../solvers/convergence_policy.hpp"
#include <cstdio>
#include <limits>
using nlohmann::json;
namespace {
int failed = 0;
void expect(bool ok, const char* name) {
    std::printf("[%s] %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++failed;
}
json blocks(std::initializer_list<double> values) {
    json result = {{"x_alpha", json::array()}};
    int i = 0;
    for (double value : values)
        result["x_alpha"].push_back({{"function_index", i++}, {"value", value},
                                    {"availability", "measured"}});
    return result;
}
}
struct DiagnosticSolver {
    struct State { int iter = 0; bool diverged = false; };
    int diverge_at = 0;
    bool reported_converged = false;
    State step(State s) const { ++s.iter; s.diverged = s.iter == diverge_at; return s; }
    [[nodiscard]] static bool converged(const State& s) { return s.iter >= 2; }
    void print_final(const State&, bool converged) { reported_converged = converged; }
};
int main() {
    expect(molresponse_v3::fd_exploded(std::numeric_limits<double>::quiet_NaN(), 1., 2.),
           "nonfinite BSH evidence trips invalid-state protection");
    DiagnosticSolver solver;
    molresponse_v3::solvers::IteratePolicy forced;
    forced.max_iters = 8;
    forced.forced_iterations = 4;
    expect(molresponse_v3::solvers::iterate(solver, DiagnosticSolver::State{}, forced).iter == 4,
           "forced diagnostic executes exactly supplied iteration count");
    solver.diverge_at = 3;
    expect(molresponse_v3::solvers::iterate(solver, DiagnosticSolver::State{}, forced).iter == 3,
           "forced diagnostic cannot bypass divergent-state protection");
    expect(!solver.reported_converged, "divergent diagnostic cannot retain an earlier convergence verdict");
    const json targets = {{"trace_change", {{"value", 0.1}, {"units", "hartree"}}}};
    json measurements = {{"trace_change", {{"value", 0.05}, {"units", "hartree"}, {"availability", "measured"}}}};
    expect(madness::evaluate_explicit_criteria(measurements, targets).at("status") == "satisfied",
           "explicit criterion evaluates supplied target");
    measurements["trace_change"]["value"] = 0.2;
    expect(madness::evaluate_explicit_criteria(measurements, targets).at("status") == "unmet",
           "unmet explicit criterion is not success");
    expect(madness::evaluate_explicit_criteria(json::object(), targets).at("status") == "unavailable",
           "missing required criterion is not successful accuracy");
    measurements["trace_change"]["units"] = "dimensionless";
    expect(madness::evaluate_explicit_criteria(measurements, targets).at("status") == "invalid",
           "incompatible criterion units rejected");
    expect(madness::evaluate_explicit_criteria(measurements, json::object()).at("status") == "not_requested",
           "no explicit targets introduces no numerical default");
    using madness::function_metric_maximum;
    using madness::function_metric_absolute_change;
    using madness::function_metric_trace;
    using madness::function_metric_trace_change;
    auto trace = function_metric_trace(blocks({1., 3.}));
    expect(trace.value("value", -1.) == 4., "unweighted trace sums diagonal contributions");
    expect(function_metric_trace_change(blocks({2., 2.}), blocks({1., 3.})).value("value", -1.) == 0.,
           "unitary redistribution preserves trace change");
    expect(function_metric_trace_change(blocks({2., 3.}), blocks({1., 3.})).value("value", -1.) == 1.,
           "changed subspace trace is detected");
    expect(function_metric_trace_change(blocks({4.}), blocks({1., 3.})).value("availability", "") == "unavailable",
           "equal trace cannot compare different dimensions");
    expect(function_metric_trace(blocks({std::numeric_limits<double>::infinity()})).value("availability", "") == "invalid",
           "nonfinite trace contribution rejected");
    auto r = function_metric_maximum(blocks({3., 4.}));
    expect(r.value("value", -1.) == 4., "maximum differs from quadrature (5) and RMS");
    expect(r.value("limiting_function", json()) == json({{"block", "x_alpha"}, {"function_index", 1}}),
           "identify the limiting stored function");
    r = function_metric_absolute_change(blocks({2., 2.}), blocks({1., 3.}));
    expect(r.value("value", -1.) == 1., "difference each function before reducing; equal totals do not imply zero");
    r = function_metric_absolute_change(blocks({2., 3.}), blocks({1., 3.}));
    expect(r.value("value", -1.) == 1., "unchanged maximum kinetic value does not imply zero change");
    r = function_metric_absolute_change(blocks({1., 3.}), blocks({1., 3.}));
    expect(r.value("value", -1.) == 0., "genuine comparable equal values have measured zero change");
    expect(function_metric_maximum(blocks({})).value("availability", "") == "unavailable",
           "empty required collection is unavailable");
    expect(function_metric_absolute_change(blocks({1.}), blocks({1., 2.})).value("availability", "") == "unavailable",
           "different component counts are incomparable");
    auto different = blocks({1.}); different["y_alpha"] = different["x_alpha"]; different.erase("x_alpha");
    expect(function_metric_absolute_change(blocks({1.}), different).value("availability", "") == "unavailable",
           "different stored blocks are incomparable");
    for (double bad : {-1., std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
        expect(function_metric_maximum(blocks({bad})).value("availability", "") == "invalid",
               "negative or nonfinite component cannot produce a maximum");
        expect(function_metric_absolute_change(blocks({1.}), blocks({bad})).value("availability", "") == "invalid",
               "invalid prior component cannot produce a kinetic change");
    }
    auto bad = blocks({1., 2.}); bad["x_alpha"][1]["function_index"] = 0;
    expect(function_metric_maximum(bad).value("availability", "") == "invalid", "duplicate function index rejected");
    bad = blocks({1.}); bad["x_alpha"][0]["value"] = nullptr;
    expect(function_metric_maximum(bad).value("availability", "") == "invalid", "null is not zero");
    expect(function_metric_maximum(json::array()).value("availability", "") == "invalid", "malformed block collection rejected");

    auto density = molresponse_v3::density_history_record(0.125, true, 0.01, 0.1, 2.0,
        {{"iteration", 1}, {"attempt_directory", "a"}});
    expect(density.at("value") == 0.125 && std::abs(density.at("effective_limit").get<double>() - 0.21) < 1e-15,
           "density history keeps full precision and contemporaneous normalization");
    expect(molresponse_v3::density_history_record(0., false, 0.01, 0.1, 2., json::object()).at("value").is_null(),
           "first density comparison remains unavailable");
    using molresponse_v3::kinetic_change_record;
    using molresponse_v3::bsh_maximum_comparison;
    json previous = {{"blocks", blocks({1., 3.})}, {"availability", "measured"},
        {"quantity", "per_function_kinetic_quadratic"}, {"units", "response_kinetic"},
        {"phase", "initial_solver_iterate"}, {"valid_for_saved_state", false},
        {"sample", {{"attempt_directory", "attempt-a"}, {"iteration", 0}, {"threshold", 1e-4}, {"k", 6}}}};
    auto current = previous;
    current["blocks"] = blocks({2., 2.});
    current["sample"]["iteration"] = 1;
    current["phase"] = "solver_iterate_after_acceleration";
    auto change = kinetic_change_record(current, previous);
    expect(change.value("value", -1.) == 0., "response adapter compares adjacent traces, invariant to rotation");
    expect(change.value("active_gate", true) == false && change.contains("target") && change["target"].is_null(),
           "kinetic comparison has no invented target or active gate");
    expect(kinetic_change_record(previous).value("availability", "") == "unavailable", "first sample has no temporal difference");
    expect(kinetic_change_record(current, current).value("availability", "") == "unavailable", "duplicate reevaluation is not iteration convergence");
    for (const char* field : {"attempt_directory", "threshold", "k"}) {
        auto other = previous;
        other["sample"][field] = field == std::string("attempt_directory") ? json("attempt-b") : json(12);
        expect(kinetic_change_record(current, other).value("availability", "") == "unavailable",
               "changed attempt or representation invalidates comparison");
    }
    auto wrong_phase = current; wrong_phase["phase"] = "unrelated_evaluation";
    expect(kinetic_change_record(wrong_phase, previous).value("availability", "") == "unavailable",
           "incompatible current sampling phase cannot compare traces");
    auto other = previous; other["units"] = "different";
    expect(kinetic_change_record(current, other).value("availability", "") == "unavailable", "different units cannot be compared");
    auto metric = json{{"per_function", {{"blocks", blocks({3., 4.})}}},
        {"units", "response_L2"}, {"phase", "solver_step_input_before_acceleration"},
        {"sample", current["sample"]}, {"valid_for_saved_state", false},
        {"absolute_limit", 0.5}, {"relative_factor", 2.}, {"normalization_scale", 2.},
        {"normalization_phase", "solver_step_output_after_acceleration"}, {"effective_limit", 4.5}};
    auto candidate = bsh_maximum_comparison(metric);
    expect(candidate.value("value", -1.) == 4. && candidate.value("satisfied", false), "max passes 4.5 while quadrature 5 fails");
    metric["absolute_limit"] = 0.; metric["effective_limit"] = 4.;
    expect(bsh_maximum_comparison(metric).value("satisfied", false), "response inclusive equality preserved");
    metric.erase("normalization_scale");
    candidate = bsh_maximum_comparison(metric);
    expect(candidate.contains("satisfied") && candidate["satisfied"].is_null(), "missing historical norm never borrows final limit");
    metric["normalization_scale"] = 2.; metric["effective_limit"] = 99.;
    expect(bsh_maximum_comparison(metric).value("availability", "") == "invalid", "inconsistent limit arithmetic rejected");

    metric["effective_limit"] = 4.;
    metric["normalization_phase"] = "solver_step_input_before_acceleration";
    expect(bsh_maximum_comparison(metric).value("availability", "") == "invalid", "wrong-phase normalization cannot support a verdict");
    metric["normalization_phase"] = "solver_step_output_after_acceleration";
    metric["effective_limit"] = std::nextafter(4., 5.);
    expect(bsh_maximum_comparison(metric).value("availability", "") == "measured",
           "one-ulp multiply-add rounding difference is valid evidence");
    using molresponse_v3::valid_metric_diagnostics;
    metric["effective_limit"] = 4.;
    metric["maximum_comparison"] = bsh_maximum_comparison(metric);
    json density_entry = {{"density_history", json::array({density})}};
    expect(valid_metric_diagnostics(density_entry), "density history round-trips");
    density_entry["density_history"][0]["effective_limit"] = 99.;
    expect(!valid_metric_diagnostics(density_entry), "corrupted density limit is rejected");
    json entry = {{"function_bsh_history", json::array({metric})}};
    expect(valid_metric_diagnostics(entry), "valid counterfactual round-trips");
    entry["function_bsh_history"][0]["maximum_comparison"]["value"] = 0.;
    expect(!valid_metric_diagnostics(entry), "tampered maximum cannot hide behind valid raw components");
    previous["change_from_previous"] = kinetic_change_record(previous);
    current["change_from_previous"] = kinetic_change_record(current, previous);
    entry = {{"function_kinetic_history", json::array({previous, current})}};
    expect(valid_metric_diagnostics(entry), "valid kinetic history round-trips");
    entry["function_kinetic_history"][1]["change_from_previous"]["current_trace"] = 99.;
    expect(!valid_metric_diagnostics(entry), "tampered trace change rejected");
    entry["function_kinetic_history"][1]["change_from_previous"] = kinetic_change_record(current, previous);
    entry["function_kinetic_history"][1]["sample"]["k"] = 8;
    expect(!valid_metric_diagnostics(entry), "changed segment cannot retain measured differences");

    auto input = previous;
    input["phase"] = "solver_step_input_before_acceleration";
    entry = {{"function_kinetic_history", json::array({previous, current})},
        {"convergence_measurements", {{"kinetic_energy", input},
            {"delta_kinetic_energy", previous["change_from_previous"]}}}};
    expect(valid_metric_diagnostics(entry), "historical input agrees with penultimate iterate");
    entry["convergence_measurements"]["kinetic_energy"]["blocks"]["x_alpha"][0]["value"] = 77.;
    expect(!valid_metric_diagnostics(entry), "historical raw kinetic value must agree with delta source");
    entry["convergence_measurements"]["kinetic_energy"] = input;
    entry["convergence_measurements"]["kinetic_energy"]["sample"]["iteration"] = 7;
    expect(!valid_metric_diagnostics(entry), "historical kinetic sample must agree with delta source");
    expect(valid_metric_diagnostics(json::object()), "old records without new diagnostics remain usable");
    return failed ? 1 : 0;
}
