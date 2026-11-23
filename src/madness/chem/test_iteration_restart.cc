// Restart uses the same measurement assessment as a live solve, per task.
#include <madness/chem/LadderPlanner.h>
#include <iostream>
#include <stdexcept>

namespace {
using namespace madness::iteration;
int failures = 0;
void check(bool ok, const char* why) { if (!ok) { ++failures; std::cerr << "FAILED: " << why << '\n'; } }
template<class F> void rejects(F&& f, const char* why) {
    try { f(); check(false, why); } catch (const std::invalid_argument&) {}
}
}
int main() {
    const std::vector<Protocol> ladder{{1e-4, 6}, {1e-6, 8}};
    // Floor factor zero: these cases are about restart planning, not the floor.
    Request request{{{Metric::bsh_residual, 1e-5}}, 0.0, 3, 0};
    ModelMetrics model{{Metric::bsh_residual, {Support::supported, true, 1e-5, LimitRule::plain, {}}}};
    ExpectedTask task{"state_z_0.05", "ground_generation_1", {"z"}};
    SavedTask saved;
    saved.metrics.info = {task.task_id, "old_attempt", ladder.back(), 20, task.context_id};
    StateMetrics subject; subject.id = "z";
    subject.metrics[Metric::bsh_residual] = MetricValue::measured(2e-6, 0.5);
    saved.metrics.states = {subject};
    saved.protocol_complete = true;
    saved.updates_at_protocol = 20;
    auto plan = plan_restart(ladder, request, model, task, saved);
    check(plan.action == RestartAction::reuse && plan.start_protocol == 1 && plan.completion
          && !plan.completion->checkpoint && plan.completion->info.attempt_id == "old_attempt",
          "saved final evidence reused with original identity and no archive request");
    request.targets[Metric::bsh_residual] = 3e-6;
    check(plan_restart(ladder, request, model, task, saved).action == RestartAction::reuse,
          "tighter already-satisfied request needs no evaluation");
    request.targets[Metric::bsh_residual] = 1e-7;
    auto floored = request; floored.floor_factor.reset();
    check(plan_restart(ladder, floored, model, task, saved).action == RestartAction::resume,
          "no floor at the final level: a 1e-7 request is not met by the saved 2e-6");
    plan = plan_restart(ladder, request, model, task, saved);
    check(plan.action == RestartAction::resume && plan.start_protocol == 1
          && plan.main_updates == 0 && plan.max_updates == 3 && plan.updates_at_protocol == 20,
          "unmet final criteria resume final protocol with fresh allowance and separate history");
    saved.metrics.info.protocol = ladder.front(); saved.protocol_complete = false;
    plan = plan_restart(ladder, request, model, task, saved);
    check(plan.action == RestartAction::resume && plan.start_protocol == 0,
          "interrupted intermediate stage resumes there");
    saved.protocol_complete = true;
    plan = plan_restart(ladder, request, model, task, saved);
    check(plan.action == RestartAction::resume && plan.start_protocol == 1 && plan.updates_at_protocol == 0,
          "completed intermediate stage advances even when its old metrics miss final tolerances");
    saved.metrics.info.protocol = {1e-5, 7}; saved.conversion_supported = true;
    plan = plan_restart(ladder, request, model, task, saved);
    check(plan.action == RestartAction::guess && plan.start_protocol == 0 && plan.reset_comparisons,
          "removed protocol uses supported conversion as a fresh first-stage guess");
    saved.conversion_supported = false;
    rejects([&] { plan_restart(ladder, request, model, task, saved); }, "unsupported conversion fails strictly");
    saved.metrics.info.protocol = {1e-6, 6};
    rejects([&] { plan_restart(ladder, request, model, task, saved); }, "wrong k is not a matching final checkpoint");
    saved.metrics.info.protocol = {1e-7, 8}; saved.conversion_supported = true;
    check(plan_restart(ladder, request, model, task, saved).action == RestartAction::guess,
          "finer removed threshold still enters the first requested stage as a guess");
    saved.metrics.info.protocol = ladder.back();
    saved.use = SavedEvidenceUse::guess_only;
    plan = plan_restart(ladder, request, model, task, saved);
    check(plan.action == RestartAction::resume && plan.start_protocol == 1 && plan.reset_comparisons,
          "changed operator settings keep matching protocol but require fresh evidence");
    saved.use = SavedEvidenceUse::invalid;
    rejects([&] { plan_restart(ladder, request, model, task, saved); }, "invalid selected checkpoint never falls back");
    saved.use = SavedEvidenceUse::compatible;
    saved.metrics.info.context_id = "stale_ground";
    rejects([&] { plan_restart(ladder, request, model, task, saved); }, "compatible claim with wrong context rejected");
    saved.metrics.info.context_id = task.context_id;
    saved.bundle_expanded = true; task.states.emplace_back("new_root");
    saved.metrics.info.protocol = ladder.front(); saved.protocol_complete = true;
    plan = plan_restart(ladder, request, model, task, saved);
    check(plan.action == RestartAction::resume && plan.start_protocol == 0 && plan.reset_comparisons,
          "added coupled roots renew the solve at its saved stage, without advancing old completed bundle");
    saved.bundle_expanded = false;
    rejects([&] { plan_restart(ladder, request, model, task, saved); }, "unexplained missing bundle member rejected");
    const auto fresh = plan_restart(ladder, request, model, task, std::nullopt);
    check(fresh.action == RestartAction::fresh && fresh.start_protocol == 0 && fresh.main_updates == 0,
          "new independent state starts first regardless of other saved states");
    auto appended = ladder; appended.push_back({1e-7, 10});
    task.states = {"z"}; saved.metrics.info.protocol = ladder.back();
    plan = plan_restart(appended, request, model, task, saved);
    check(plan.start_protocol == 2 && plan.action == RestartAction::resume,
          "appended final stage is mandatory after completed former final stage");
    rejects([&] { plan_restart({}, request, model, task, saved); }, "empty ladder rejected");
    auto duplicate = ladder; duplicate.push_back(ladder.back());
    rejects([&] { plan_restart(duplicate, request, model, task, saved); }, "ambiguous duplicate protocols rejected");
    return failures ? 1 : 0;
}
