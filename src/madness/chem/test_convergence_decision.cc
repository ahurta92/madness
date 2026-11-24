// Unit tests for the convergence decision module. Plain values only: no World,
// no MPI, no solver. Each block follows one line of the step 1 test list.
#include <madness/chem/ConvergenceDecision.h>
#include <madness/chem/IterationOutput.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {
using namespace madness::iteration;
int failures = 0;
void check(bool ok, const std::string& message) {
    if (!ok) { ++failures; std::cerr << "FAILED: " << message << '\n'; }
}
// Requires an invalid_argument whose message contains every given fragment.
template <class F> void rejects(F&& call, const std::vector<std::string>& fragments, const std::string& message) {
    try { call(); check(false, message + " (no error)"); }
    catch (const std::invalid_argument& error) {
        for (const auto& fragment : fragments)
            check(std::string(error.what()).find(fragment) != std::string::npos,
                  message + " (message lacks '" + fragment + "': " + error.what() + ")");
    }
}
bool near(double a, double b) { return std::abs(a - b) <= 1e-15 * std::max(1.0, std::abs(b)); }

const Protocol coarse{1e-4, 6}, fine{1e-6, 8};

// A ground-state-like model: three defaults, kinetic on request.
ModelMetrics ground() {
    return {{Metric::bsh_residual, {Support::supported, true, 1e-4, LimitRule::plain, {}}},
            {Metric::density_change, {Support::supported, true, 1e-4, LimitRule::plain, 7.0}},
            {Metric::energy_change, {Support::supported, true, 1e-5, LimitRule::plain, {}}},
            {Metric::kinetic_change, {Support::supported, false, 1e-5, LimitRule::plain, {}}}};
}
// A driven-response-like model: scaled BSH and density, energy not applicable.
ModelMetrics driven() {
    return {{Metric::bsh_residual, {Support::supported, true, 5e-4, LimitRule::scaled, {}}},
            {Metric::density_change, {Support::supported, true, 5e-4, LimitRule::scaled, {}}},
            {Metric::energy_change, {Support::not_applicable, false, {}, LimitRule::plain, {}}},
            {Metric::kinetic_change, {Support::supported, false, {}, LimitRule::plain, {}}}};
}
const ResolvedCriterion* find(const ResolvedSettings& s, Metric metric) {
    for (const auto& c : s.criteria) if (c.metric == metric) return &c;
    return nullptr;
}

const ExpectedTask task{"scf", "ctx", {"ground"}};
IterationMetrics bsh_only(double value, std::size_t iteration = 1, Protocol protocol = fine) {
    return {{"scf", protocol, iteration, "ctx"},
            {{"ground", {}, {{Metric::bsh_residual, MetricValue::measured(value)}}}}};
}
ResolvedSettings bsh_settings(double target, Protocol protocol = fine, std::size_t cap = 5) {
    Request r; r.targets[Metric::bsh_residual] = target; r.max_updates = cap;
    return resolve(r, {{Metric::bsh_residual, {Support::supported, false, {}, LimitRule::plain, 0.0}}}, protocol, true);
}
MetricCheck one(const ResolvedCriterion& c, const MetricValue& value) {
    return check_metric(c, value, std::nullopt, {"scf", fine, 1, "ctx"});
}
} // namespace

int main() {
    // ---- Resolving criteria ----
    {
        const auto s = resolve({}, ground(), fine, true);
        check(s.criteria.size() == 3 && find(s, Metric::bsh_residual) && find(s, Metric::density_change)
                  && find(s, Metric::energy_change) && !find(s, Metric::kinetic_change),
              "no written targets: criteria are exactly the defaults");
        check(near(find(s, Metric::bsh_residual)->requested, 1e-4) && near(find(s, Metric::energy_change)->requested, 1e-5),
              "default targets are used when nothing is written");
        const auto d = resolve({}, driven(), coarse, true);
        check(d.criteria.size() == 2 && !find(d, Metric::energy_change), "driven defaults: BSH and density only");
    }
    {
        Request r; r.targets[Metric::kinetic_change] = 3e-5;
        const auto s = resolve(r, ground(), fine, true);
        check(s.criteria.size() == 4 && find(s, Metric::kinetic_change) && near(find(s, Metric::kinetic_change)->requested, 3e-5),
              "written target for a supported non-default metric adds it");
        r.targets[Metric::bsh_residual] = 2e-4;
        check(near(find(resolve(r, ground(), fine, true), Metric::bsh_residual)->requested, 2e-4), "written target replaces the default");
    }
    {
        auto m = ground(); m.erase(Metric::kinetic_change);
        Request r; r.targets[Metric::kinetic_change] = 1e-5;
        rejects([&] { resolve(r, m, fine, true, "ground state"); }, {"kinetic_change", "ground state"},
                "target for a metric absent from the model is unsupported");
        m[Metric::kinetic_change].support = Support::unsupported;
        rejects([&] { resolve(r, m, fine, true, "ground state"); }, {"kinetic_change", "not supported", "ground state"},
                "target for an unsupported metric names metric and application");
    }
    {
        // A written target for a not-applicable metric is reported, not judged.
        Request r; r.targets[Metric::energy_change] = 1e-6; r.max_updates = 5;
        const auto s = resolve(r, driven(), fine, true, "driven response");
        check(s.criteria.size() == 2 && !find(s, Metric::energy_change)
                  && s.not_applicable == std::vector<Metric>{Metric::energy_change},
              "not-applicable target is not a criterion and is listed as not applicable");
        auto unscaled = s; for (auto& c : unscaled.criteria) c.rule = LimitRule::plain;
        IterationMetrics m{{"scf", fine, 1, "ctx"}, {{"ground", {}, {}}}};
        const auto d = decide(unscaled, task, {1, 1, true}, m, true);
        bool reported = false;
        for (const auto& c : d.checks)
            reported = reported || (c.metric == Metric::energy_change && c.result == CheckResult::not_applicable);
        check(reported && d.action == Action::update, "not-applicable row is reported and never blocks or passes");
        check(resolve({}, driven(), fine, true).not_applicable.empty(), "nothing written: nothing reported");
    }
    {
        auto m = ground(); m.at(Metric::energy_change).default_target.reset();
        rejects([&] { resolve({}, m, fine, true); }, {"energy_change"}, "default criterion with no target is an error");
        Request r; r.targets[Metric::energy_change] = 1e-6;
        check(near(find(resolve(r, m, fine, true), Metric::energy_change)->requested, 1e-6),
              "a written target supplies a default criterion that has no default");
    }
    for (double bad : {-1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
        Request r; r.targets[Metric::bsh_residual] = bad;
        rejects([&] { resolve(r, ground(), fine, true); }, {"bsh_residual"}, "invalid target rejected");
        Request f; f.floor_factor = bad;
        rejects([&] { resolve(f, ground(), fine, true); }, {"floor_factor"}, "invalid user floor factor rejected");
        auto m = ground(); m.at(Metric::bsh_residual).floor_factor = bad;
        rejects([&] { resolve({}, m, fine, true); }, {"bsh_residual"}, "invalid model floor factor rejected");
    }
    rejects([&] { resolve({}, ground(), {0.0, 6}, true); }, {"protocol"}, "invalid protocol rejected");
    rejects([&] { resolve({}, {}, fine, true); }, {"no convergence criteria"}, "a model with no criteria is an error");

    // ---- The floor ----
    {
        // Intermediate level: the floor raises BSH and density targets.
        const auto c = *find(resolve({}, ground(), fine, false), Metric::bsh_residual); // 1e-4 vs 5e-6
        check(near(c.floor, 5e-6) && near(c.effective, 1e-4), "target above the floor is unchanged");
        Request r; r.targets[Metric::bsh_residual] = 1e-6;
        const auto raised = *find(resolve(r, ground(), fine, false), Metric::bsh_residual);
        check(near(raised.requested, 1e-6) && near(raised.floor, 5e-6) && near(raised.effective, 5e-6),
              "intermediate: a target below the floor becomes the floor");
        Request equal; equal.targets[Metric::bsh_residual] = 0.625; // 5 * 0.125, exact in binary
        const auto at_floor = *find(resolve(equal, ground(), {0.125, 6}, false), Metric::bsh_residual);
        check(at_floor.floor == 0.625 && at_floor.effective == 0.625, "a target equal to the floor is kept");
        Request below; below.targets[Metric::bsh_residual] = 0.5;
        check(find(resolve(below, ground(), {0.125, 6}, false), Metric::bsh_residual)->effective == 0.625,
              "a target just below the floor is raised");
    }
    {
        // Final level: the request is judged as written; the floor is only reported.
        Request r; r.targets[Metric::bsh_residual] = 1e-6; r.max_updates = 5;
        const ModelMetrics bsh{{Metric::bsh_residual, {Support::supported, false, {}, LimitRule::plain, {}}}};
        const auto s = resolve(r, bsh, fine, true);
        const auto& c = *find(s, Metric::bsh_residual);
        check(near(c.effective, 1e-6) && near(c.floor, 5e-6), "final level: written target kept, floor reported");
        const auto fin = decide(s, task, {1, 1, true}, bsh_only(4e-6), true);
        check(fin.action == Action::update && fin.checks[0].limit && near(*fin.checks[0].limit, 1e-6),
              "final level: a value between target and floor is unmet");
        const auto mid = decide(resolve(r, bsh, fine, false), task, {1, 1, true}, bsh_only(4e-6), false);
        check(mid.action == Action::finish_protocol && near(*mid.checks[0].limit, 5e-6),
              "intermediate level: the same value meets the floor");
    }
    {
        // Energy and kinetic change have no floor at any level.
        Request r; r.targets[Metric::energy_change] = 1e-9; r.targets[Metric::kinetic_change] = 1e-9;
        for (bool final_level : {false, true}) {
            const auto s = resolve(r, ground(), coarse, final_level);
            for (auto metric : {Metric::energy_change, Metric::kinetic_change})
                check(find(s, metric)->floor == 0.0 && find(s, metric)->effective == 1e-9,
                      "energy and kinetic targets are never floored");
        }
        check(has_floor(Metric::bsh_residual) && has_floor(Metric::density_change)
                  && !has_floor(Metric::energy_change) && !has_floor(Metric::kinetic_change),
              "only BSH residual and density change have a floor");
    }
    {
        auto m = ground();
        const auto app = resolve({}, m, fine, false);
        check(near(find(app, Metric::density_change)->floor, 7e-6), "application factor overrides the shared default");
        check(near(find(app, Metric::bsh_residual)->floor, shared_floor_factor * 1e-6) && shared_floor_factor == 5.0,
              "shared default factor is 5");
        Request r; r.floor_factor = 2.0;
        const auto user = resolve(r, m, fine, false);
        check(near(find(user, Metric::density_change)->floor, 2e-6) && near(find(user, Metric::bsh_residual)->floor, 2e-6),
              "user factor overrides application and shared factors");
        r.floor_factor = 0.0;
        check(near(find(resolve(r, m, fine, false), Metric::bsh_residual)->floor, 0.0), "a zero user factor is allowed");
    }
    {
        Request r; r.targets[Metric::bsh_residual] = 1e-5;
        const auto at_coarse = *find(resolve(r, ground(), coarse, false), Metric::bsh_residual);
        const auto at_fine = *find(resolve(r, ground(), fine, false), Metric::bsh_residual);
        check(near(at_coarse.floor, 5e-4) && near(at_fine.floor, 5e-6)
                  && near(at_coarse.effective, 5e-4) && near(at_fine.effective, 1e-5),
              "the same request at two thresholds gives two floors");
    }
    {
        // The start-of-task warning lists final targets below the floor, unchanged.
        Request r; r.targets[Metric::bsh_residual] = 1e-6;
        std::ostringstream warned, silent;
        write_targets_below_floor(warned, "SCF", "ground_state", resolve(r, ground(), fine, true));
        write_targets_below_floor(silent, "SCF", "ground_state", resolve({}, ground(), fine, true));
        check(warned.str().find("Warning: SCF ground_state") != std::string::npos
                  && warned.str().find("requested 1.0000e-06, floor 5.0000e-06") != std::string::npos
                  && silent.str().empty(), "warning names each target below the final floor");
    }

    // ---- Limits ----
    {
        ResolvedCriterion plain{Metric::bsh_residual, 1e-3, 0.0, 1e-3, LimitRule::plain};
        const auto equal = one(plain, MetricValue::measured(1e-3));
        check(equal.result == CheckResult::met && equal.limit && near(*equal.limit, 1e-3)
                  && near(equal.absolute, 1e-3) && equal.relative == 0.0,
              "plain: value equal to T is met; relative recorded as zero");
        check(one(plain, MetricValue::measured(std::nextafter(1e-3, 1.0))).result == CheckResult::unmet,
              "plain: value just above T is unmet");
        check(one(plain, MetricValue::measured(5e-4, 100.0)).limit == 1e-3, "plain ignores any measured scale");
    }
    {
        ResolvedCriterion scaled{Metric::bsh_residual, 1e-3, 0.0, 1e-3, LimitRule::scaled};
        const auto zero = one(scaled, MetricValue::measured(1e-3, 0.0));
        check(zero.result == CheckResult::met && zero.limit && near(*zero.limit, 1e-3)
                  && near(zero.absolute, 1e-3) && near(zero.relative, 1e-3) && zero.scale && *zero.scale == 0.0,
              "scaled, scale 0: limit T; relative recorded as T");
        check(one(scaled, MetricValue::measured(std::nextafter(1e-3, 1.0), 0.0)).result == CheckResult::unmet,
              "scaled, scale 0: just above T is unmet");
        const auto large = one(scaled, MetricValue::measured(0.05, 99.0));
        check(large.result == CheckResult::met && large.limit && near(*large.limit, 0.1), "scaled, large scale: T*(1+scale)");
        check(one(scaled, MetricValue::measured(0.1 + 1e-12, 99.0)).result == CheckResult::unmet, "scaled: above T*(1+scale) is unmet");
        auto derived = MetricValue::measured(0.002);
        derived.evidence.components = {{"x", 0.002, 1.0}, {"y", 0.001, 3.0}};
        const auto from_evidence = one(scaled, derived);
        check(from_evidence.scale && *from_evidence.scale == 3.0 && near(*from_evidence.limit, 4e-3),
              "scaled: the scale may come from the application's evidence");
    }
    {
        ResolvedCriterion scaled{Metric::bsh_residual, 1e-3, 0.0, 1e-3, LimitRule::scaled};
        const auto missing = one(scaled, MetricValue::measured(1e-9));
        check(missing.result == CheckResult::unmet && missing.reason == "scale_unavailable" && !missing.limit,
              "scaled with unavailable scale is unmet with a reason, never plain");
        auto bad = MetricValue::measured(1e-9); bad.scale = Scale::measured(-1.0);
        check(one(scaled, bad).result == CheckResult::invalid, "scaled with an invalid scale is invalid");
    }
    {
        ResolvedCriterion plain{Metric::bsh_residual, 1.0, 0.0, 1.0, LimitRule::plain};
        for (const auto& value : {MetricValue::unavailable("nope"), MetricValue::invalid("bad"),
                                  MetricValue::not_applicable("na"), MetricValue::measured(-1.0),
                                  MetricValue::measured(std::numeric_limits<double>::quiet_NaN())})
            check(one(plain, value).result != CheckResult::met, "unavailable or invalid measurement is never met");
        ResolvedCriterion density{Metric::density_change, 1.0, 0.0, 1.0, LimitRule::plain};
        check(one(density, MetricValue::measured(0.0)).reason == "comparison_unavailable",
              "a change metric without a valid comparison is not met");
    }

    // ---- The decision ----
    {
        const auto s = bsh_settings(1e-3, fine, 5);
        const auto met = bsh_only(1e-4);
        auto d = decide(s, task, {1, 1, true}, met, false);
        check(d.action == Action::finish_protocol && d.reason == "intermediate_criteria_met" && d.checkpoint,
              "all met at an intermediate level finishes the level");
        d = decide(s, task, {1, 1, true}, met, true);
        check(d.action == Action::complete && d.reason == "final_criteria_met", "all met at the final level completes");
    }
    {
        const auto s = bsh_settings(1e-3, fine, 3);
        const auto unmet = bsh_only(1.0, 3);
        const auto mid = decide(s, task, {3, 3, true}, unmet, false);
        const auto fin = decide(s, task, {3, 3, true}, unmet, true);
        check(mid.action == Action::finish_protocol && mid.reason == "iteration_limit", "cap at intermediate finishes the level");
        check(fin.action == Action::stop && fin.reason == "iteration_limit", "cap at final stops");
        const auto below = decide(s, task, {2, 2, true}, bsh_only(1.0, 2), true);
        check(below.action == Action::update && below.reason == "criteria_unmet" && !below.checkpoint, "unmet below cap updates");
        auto cadence = s; cadence.checkpoint_interval = 2;
        check(decide(cadence, task, {2, 2, true}, bsh_only(1.0, 2), true).checkpoint, "update on the cadence checkpoints");
    }
    {
        const auto s = bsh_settings(1e-3, fine, 5);
        for (bool final_level : {true, false}) {
            const auto d = decide(s, task, {1, 1, false}, bsh_only(0.0), final_level);
            check(d.action == Action::stop && d.reason == "invalid_state", "invalid state stops, whatever the measurements");
        }
        auto other_task = bsh_only(0.0); other_task.info.task_id = "other";
        auto other_level = bsh_only(0.0); other_level.info.protocol = coarse;
        auto other_context = bsh_only(0.0); other_context.info.context_id = "other";
        for (const auto& m : {other_task, other_level, other_context}) {
            const auto d = decide(s, task, {1, 1, true}, m, true);
            check(d.action == Action::stop && d.reason == "iteration_context_mismatch", "foreign task/level stops with context reason");
        }
        const auto d = decide(s, task, {2, 2, true}, bsh_only(0.0, 1), true);
        check(d.action == Action::stop && d.reason == "iteration_context_mismatch", "foreign iteration stops with context reason");
        auto stranger = bsh_only(0.0); stranger.states[0].id = "someone";
        check(decide(s, task, {1, 1, true}, stranger, true).reason == "subject_membership_mismatch", "foreign subject stops");
    }
    {
        // Sweep inputs: no outcome is ever a stall, and every action is one of four.
        const auto s = bsh_settings(1e-3, fine, 4);
        for (std::size_t n = 0; n <= 5; ++n)
            for (double value : {0.0, 1e-3, 1.0, 1e3})
                for (bool final_level : {true, false})
                    for (bool valid : {true, false}) {
                        const auto d = decide(s, task, {n, n, valid}, bsh_only(value, n), final_level);
                        check(d.reason.find("stall") == std::string::npos
                                  && (d.action == Action::update || d.action == Action::finish_protocol
                                      || d.action == Action::complete || d.action == Action::stop),
                              "no input produces a stall outcome");
                    }
    }

    // ---- Saved states ----
    {
        const auto saved = bsh_only(1e-4, 7);
        const auto loose = bsh_settings(1e-3, fine);
        const auto completion = saved_completion(loose, task, saved, true);
        check(completion && completion->action == Action::complete && completion->reason == "saved_metrics_satisfy_request",
              "saved measurements meeting the current limits complete at the final level");
        check(!saved_completion(loose, task, saved, false), "saved measurements at an intermediate level return nothing");
        const auto tight = bsh_settings(1e-5, fine);
        check(!saved_completion(tight, task, saved, true), "met an older, looser request but not the current one: nothing");
        auto foreign = saved; foreign.info.task_id = "other";
        check(!saved_completion(loose, task, foreign, true), "invalid saved measurements return nothing");
    }

    if (failures == 0) std::cout << "test_convergence_decision passed\n";
    return failures == 0 ? 0 : 1;
}
