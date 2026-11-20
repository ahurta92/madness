// Shared-controller contract. Scalar fixtures deliberately need no World or solver.
#include <madness/chem/ConvergenceSettings.h>
#include <madness/chem/IterationMetrics.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
using namespace madness::iteration;
int failures = 0;
void check(bool ok, const std::string& message) {
    if (!ok) { ++failures; std::cerr << "FAILED: " << message << '\n'; }
}
template <class F> void rejects(F&& call, const std::string& message) {
    try { call(); check(false, message); }
    catch (const std::invalid_argument&) {}
}
Request request() {
    Request r;
    r.criteria = {Metric::bsh_residual, Metric::energy_change};
    r.tolerances = {{Metric::bsh_residual, {1e-5, 0.0}},
                    {Metric::energy_change, {1e-6, 0.0}}};
    r.max_updates = 2;
    return r;
}
ModelMetrics model() {
    return {{Metric::bsh_residual, {Support::supported, 1e-5}},
            {Metric::energy_change, {Support::supported, 1e-6}},
            {Metric::density_change, {Support::supported, 1e-4}},
            {Metric::kinetic_change, {Support::unsupported, std::nullopt}}};
}
ExpectedTask expected() { return {"scf", "model1", {"ground"}}; }
IterationMetrics iteration_metrics() {
    IterationMetrics e;
    e.info = {"scf", "run1", {1e-6, 8}, 0, "model1"};
    StateMetrics s;
    s.id = "ground";
    s.comparison = Comparison{"scf", "run0", {1e-6, 8}, 7, "ground"};
    s.metrics[Metric::bsh_residual] = MetricValue::measured(5e-6, 1.0);
    s.metrics[Metric::energy_change] = MetricValue::measured(5e-7, 75.0);
    e.states.push_back(s);
    return e;
}
void resolution() {
    auto r = request();
    auto m = model();
    m.at(Metric::bsh_residual).intermediate = LegacyFloor{5e-4, "5 * threshold", {{"threshold", 1e-4}, {"factor", 5.0}}};
    auto coarse = resolve(r, m, {1e-4, 6}, false);
    auto fine = resolve(r, m, {1e-6, 8}, true);
    check(coarse.criteria[0].absolute == 5e-4 && fine.criteria[0].absolute == 1e-5,
          "legacy floor relaxes intermediate absolute limits only");
    check(coarse.criteria[0].adjustment.has_value() && !fine.criteria[0].adjustment,
          "resolved report retains only applied intermediate heuristic");
    check(coarse.criteria[1].absolute == 1e-6 && !coarse.criteria[1].adjustment,
          "a metric without a legacy rule keeps its requested limit");
    r.tolerances.at(Metric::bsh_residual) = {1e-3, 0.25};
    coarse = resolve(r, m, {1e-4, 6}, false);
    check(coarse.criteria[0].absolute == 1e-3 && coarse.criteria[0].relative == 0.25,
          "intermediate heuristic cannot tighten a requested limit or change relative tolerance");
    r = request(); r.tolerances.erase(Metric::energy_change);
    check(resolve(r, m, {1e-6, 8}, true).criteria[1].absolute == 1e-6,
          "omitted absolute tolerance uses the declared model default");
    r.tolerances[Metric::energy_change] = {};
    r.tolerances[Metric::energy_change].relative = 0.05;
    const auto relative_only = resolve(r, m, {1e-6, 8}, true);
    check(relative_only.criteria[1].absolute == 1e-6 && relative_only.criteria[1].relative == 0.05,
          "relative-only override preserves the model absolute default");
    r.criteria = {Metric::bsh_residual}; r.tolerances[Metric::density_change] = {1e-9, 0.0};
    check(resolve(r, m, {1e-6, 8}, true).criteria.size() == 1,
          "a tolerance override does not activate another criterion");
    r.criteria.clear(); rejects([&] { resolve(r, m, {1e-6, 8}, true); }, "empty selection rejected");
    r = request(); r.criteria.push_back(Metric::bsh_residual);
    rejects([&] { resolve(r, m, {1e-6, 8}, true); }, "duplicate selection rejected");
    r.criteria = {Metric::kinetic_change};
    rejects([&] { resolve(r, m, {1e-6, 8}, true); }, "unsupported applicable metric rejected");
    r = request(); m.at(Metric::energy_change).default_absolute.reset();
    r.tolerances.erase(Metric::energy_change);
    rejects([&] { resolve(r, m, {1e-6, 8}, true); }, "missing absolute default rejected");
    r = request(); r.tolerances[Metric::energy_change].relative = -1;
    rejects([&] { resolve(r, m, {1e-6, 8}, true); }, "negative relative tolerance rejected");
    r = request(); r.tolerances[Metric::energy_change].absolute = std::numeric_limits<double>::infinity();
    rejects([&] { resolve(r, m, {1e-6, 8}, true); }, "nonfinite tolerance rejected");
    rejects([] { metric_from_name("dconv"); }, "legacy names are not aliases");
    check(metric_from_name("energy_change") == Metric::energy_change, "explicit metric name accepted");
}
void assessments() {
    auto r = request(); auto m = model(); auto e = iteration_metrics();
    auto settings = resolve(r, m, e.info.protocol, true);
    check(assess(settings, expected(), e).all_met, "same prepared sample satisfies both criteria");
    e.states[0].comparison.reset();
    check(!assess(settings, expected(), e).all_met, "zero or small change without predecessor cannot pass");
    e.states[0].metrics[Metric::energy_change].comparison = Comparison{"scf", "run0", {1e-6, 8}, 7, "ground"};
    check(assess(settings, expected(), e).all_met, "metric-specific predecessor can supply the comparison");
    e.states[0].metrics[Metric::energy_change].comparison = Comparison{"scf", "run1", {1e-6, 8}, 0, "ground"};
    check(!assess(settings, expected(), e).all_met, "loaded state cannot be its own predecessor");
    e = iteration_metrics(); e.states[0].metrics[Metric::bsh_residual] = MetricValue::unavailable("not_measured");
    check(!assess(settings, expected(), e).all_met, "missing metric cannot be a measured zero");
    e = iteration_metrics(); e.states[0].metrics[Metric::bsh_residual].value = std::numeric_limits<double>::quiet_NaN();
    check(!assess(settings, expected(), e).all_met, "nonfinite metric cannot converge");
    e.states[0].metrics[Metric::bsh_residual].value = -1e-6;
    check(!assess(settings, expected(), e).all_met, "negative residual cannot converge");
    e = iteration_metrics(); e.states[0].metrics[Metric::kinetic_change] = MetricValue::invalid("diagnostic_failure");
    check(assess(settings, expected(), e).all_met, "unselected diagnostic does not become a stopping gate");
    m.at(Metric::energy_change).support = Support::not_applicable;
    e.states[0].metrics.erase(Metric::energy_change);
    check(assess(resolve(r, m, e.info.protocol, true), expected(), e).all_met,
          "model-declared inapplicability excludes the criterion");
    m = model(); e = iteration_metrics();
    e.states[0].metrics[Metric::energy_change] = MetricValue::not_applicable("claimed_by_sample");
    check(!assess(settings, expected(), e).all_met, "sample cannot waive an applicable selected metric");
    e = iteration_metrics(); e.states.clear();
    check(!assess(settings, expected(), e).valid, "empty bundle cannot pass vacuously");
    e = iteration_metrics(); e.states.push_back(e.states.front());
    check(!assess(settings, expected(), e).valid, "duplicate subject invalidates bundle");
    e = iteration_metrics(); auto task = expected(); task.states.emplace_back("root2");
    check(!assess(settings, task, e).valid, "missing required root invalidates bundle");
    e = iteration_metrics(); e.info.context_id = "different_model";
    check(!assess(settings, expected(), e).valid, "different scientific context cannot pass");
    e = iteration_metrics(); e.info.protocol.k = 6;
    check(!assess(settings, expected(), e).valid, "same threshold with wrong k cannot pass");
}
void scales() {
    auto r = request(); auto e = iteration_metrics();
    r.criteria = {Metric::energy_change}; r.tolerances[Metric::energy_change] = {0.1, 0.01};
    auto settings = resolve(r, model(), e.info.protocol, true);
    e.states[0].metrics[Metric::energy_change] = MetricValue::measured(0.85, 75.0);
    auto a = assess(settings, expected(), e);
    check(a.all_met && a.checks[0].limit && std::abs(*a.checks[0].limit - 0.85) < 1e-14,
          "combined tolerance includes equality and reports its effective limit");
    auto& metric = e.states[0].metrics[Metric::energy_change];
    metric = MetricValue::measured(0.11, 0.0);
    check(!assess(settings, expected(), e).all_met, "zero scale has no artificial floor");
    metric = MetricValue::measured(0.01);
    check(!assess(settings, expected(), e).all_met, "relative request cannot fall back to absolute-only");
    metric.evidence.current = -75.0; metric.evidence.previous = -74.99;
    check(assess(settings, expected(), e).all_met, "negative energy evidence reconstructs magnitude scale");
    metric.scale = Scale::measured(std::numeric_limits<double>::infinity());
    check(!assess(settings, expected(), e).all_met, "invalid scale is not silently replaced by evidence");
    r.tolerances[Metric::energy_change].relative = 0.0;
    check(assess(resolve(r, model(), e.info.protocol, true), expected(), e).all_met,
          "absolute-only test does not require a usable scale");
}
void decisions() {
    auto r = request(); auto e = iteration_metrics(); auto settings = resolve(r, model(), e.info.protocol, true);
    Progress p;
    auto d = decide(settings, expected(), p, e);
    check(d.action == Action::complete && d.checkpoint && d.info.iteration == 0,
          "converged live sample requests endpoint publication before any update");
    p.main_updates = 2; p.updates_at_protocol = 2; e.info.iteration = 2;
    check(decide(settings, expected(), p, e).action == Action::complete,
          "criteria passing at cap still completes");
    e.states[0].metrics[Metric::bsh_residual].value = 1.0;
    d = decide(settings, expected(), p, e);
    check(d.action == Action::stop && d.reason == "iteration_limit" && d.checkpoint,
          "unmet final endpoint stops with valid checkpoint");
    p.main_updates = 0; e.info.iteration = 0;
    d = decide(settings, expected(), p, e);
    check(d.action == Action::update && !d.checkpoint, "ordinary unmet sample updates without default archive write");
    settings.final_protocol = false;
    p.main_updates = 2; p.updates_at_protocol = 2; e.info.iteration = 2;
    check(decide(settings, expected(), p, e).action == Action::finish_protocol,
          "intermediate cap permits advancement without completion");
    p.main_updates = 0; e.info.iteration = 0; e = iteration_metrics();
    check(decide(settings, expected(), p, e).action == Action::finish_protocol,
          "intermediate targets passing early finish only that protocol");
    e.states[0].metrics[Metric::bsh_residual].value = 1.0;
    settings.checkpoint_interval = 2; p.updates_at_protocol = 2;
    check(decide(settings, expected(), p, e).checkpoint, "optional interval requests publication");
    p.updates_at_protocol = 0;
    check(!decide(settings, expected(), p, e).checkpoint, "zero updates does not trigger interval save");
    p.state_valid = false;
    d = decide(settings, expected(), p, e);
    check(d.action == Action::stop && !d.checkpoint, "invalid state never replaces usable checkpoint");
    p = {}; settings.final_protocol = true; e = iteration_metrics();
    const auto saved = saved_completion(settings, expected(), e);
    check(saved && saved->action == Action::complete && !saved->checkpoint && saved->info.attempt_id == "run1",
          "saved qualifying evidence completes without numerical work or relabeling");
    r.tolerances[Metric::energy_change].absolute = 1e-8;
    check(!saved_completion(resolve(r, model(), e.info.protocol, true), expected(), e),
          "unmet tighter request needs fresh numerical evaluation");
}
void coupled_and_terminal() {
    auto r = request(); r.criteria = {Metric::bsh_residual};
    auto e = iteration_metrics(); e.states[0].id = "root1";
    auto second = e.states[0]; second.id = "root2";
    second.metrics[Metric::bsh_residual].value = 1.0;
    e.states.push_back(second);
    auto task = expected(); task.states = {"root1", "root2"};
    auto settings = resolve(r, model(), e.info.protocol, true);
    const auto choice = decide(settings, task, {}, e);
    check(choice.action == Action::update && choice.checks.size() == 2
          && choice.checks[0].result == CheckResult::met && choice.checks[1].result == CheckResult::unmet,
          "one passing root cannot complete a coupled bundle with an unmet root");
    check(!saved_completion(settings, task, e), "coupled reuse requires every root to pass");
    std::swap(e.states[0], e.states[1]);
    check(decide(settings, task, {}, e).action == Action::update, "root order does not alter membership or policy");

    e = iteration_metrics(); e.states[0].metrics[Metric::bsh_residual].value = 1.0;
    Progress progress;
    int evaluations = 0, updates = 0;
    for (;;) {
        ++evaluations;
        e.info.iteration = progress.main_updates;
        auto decision = decide(settings, expected(), progress, e);
        if (decision.action != Action::update) break;
        ++updates; ++progress.main_updates; ++progress.updates_at_protocol;
    }
    check(evaluations == 3 && updates == 2,
          "two allowed updates have a third terminal evaluation and no third update");
    auto invalid = settings;
    invalid.criteria.clear();
    rejects([&] { assess(invalid, expected(), e); }, "empty resolved criteria cannot pass");
    invalid = settings; invalid.criteria[0].applicable = false;
    rejects([&] { assess(invalid, expected(), e); }, "all-inapplicable resolved criteria cannot pass vacuously");
}
void iteration_identity() {
    auto r = request(); auto e = iteration_metrics();
    auto settings = resolve(r, model(), e.info.protocol, true);
    auto& comparison = e.states[0].comparison;
    comparison = Comparison{"scf", "run1", {1e-6, 8}, 0, "ground"};
    check(!assess(settings, expected(), e).all_met, "same task/attempt/protocol/iteration cannot compare with itself");
    comparison->iteration = 1;
    check(!assess(settings, expected(), e).all_met, "a future iteration cannot be a predecessor in the same solve");
    comparison->iteration = 0;
    comparison->protocol = {1e-4, 6};
    check(assess(settings, expected(), e).all_met, "same iteration at a different protocol has distinct identity");
    comparison->protocol = e.info.protocol; comparison->attempt_id = "run0";
    check(assess(settings, expected(), e).all_met, "same iteration in a prior attempt has distinct identity");
    comparison->task_id = "another_task";
    check(!assess(settings, expected(), e).all_met, "comparison cannot silently come from another task");
    e = iteration_metrics();
    const auto stale = decide(settings, expected(), {1, 1, true}, e);
    check(stale.action == Action::stop && !stale.checkpoint && stale.reason == "iteration_context_mismatch",
          "metrics for iteration zero cannot drive iteration one");
    const auto reused = saved_completion(settings, expected(), e);
    check(reused && reused->info == e.info, "saved evidence retains its full producing iteration identity");
}
void application_owned_loop() {
    // Catch extra KAIN updates, lost terminal evaluation, and proposal rebuilding.
    auto r = request(); r.criteria = {Metric::bsh_residual};
    auto e = iteration_metrics(); auto settings = resolve(r, model(), e.info.protocol, true);
    int evaluations = 0, proposals = 0, advances = 0, publications = 0;
    double state = 1.0;
    Progress progress;
    for (;;) {
        ++evaluations; ++proposals; // app.evaluate caches a proposal for THIS state.
        e.info.iteration = progress.main_updates;
        e.states[0].metrics[Metric::bsh_residual] = MetricValue::measured(state);
        const auto d = decide(settings, expected(), progress, e);
        if (d.checkpoint) ++publications;
        if (d.action != Action::update) break;
        ++advances; state = 0.0; // app retains its own history and uses cached raw work.
        ++progress.main_updates; ++progress.updates_at_protocol;
    }
    check(evaluations == 2 && proposals == 2 && advances == 1 && publications == 1,
          "functional loop evaluates current state then stops before another update");
    e.info.iteration = 0; r.max_updates = 0; settings = resolve(r, model(), e.info.protocol, true);
    e.states[0].metrics[Metric::bsh_residual] = MetricValue::measured(1.0);
    check(decide(settings, expected(), {}, e).action == Action::stop,
          "zero allowance still evaluates its endpoint but never updates");
}
}
int main() {
    resolution(); assessments(); scales(); decisions(); coupled_and_terminal(); iteration_identity(); application_owned_loop();
    if (failures) std::cerr << failures << " shared-controller checks failed\n";
    return failures ? 1 : 0;
}
