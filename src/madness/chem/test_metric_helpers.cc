#include <madness/chem/MetricHelpers.h>
#include <madness/chem/ConvergenceSettings.h>
#include <madness/chem/SCFEvaluationPolicy.h>
#include <cmath>
#include <iostream>
#include <limits>

using namespace madness::iteration;
int main() {
    int failures = 0;
    auto check = [&](bool ok, const char* why) {
        if (!ok) { ++failures; std::cerr << why << '\n'; }
    };
    const Comparison previous{"scf", {1e-6, 8}, 2, "ground"};
    auto residual = maximum_residual({{"alpha/0", 0.1, 1.0}, {"beta/0", 0.2, 2.0}});
    check(residual.value == 0.2 && residual.scale.value == 2.0 &&
          residual.evidence.limiting_id == "beta/0" && residual.evidence.components.size() == 2,
          "residual is the unweighted maximum, retaining every component and scale");
    auto energy = scalar_change(-3.0, -2.0, previous);
    check(energy.value == 1.0 && energy.scale.value == 3.0 && energy.comparison &&
          energy.evidence.current == -3.0 && energy.evidence.previous == -2.0,
          "scalar change retains physical values and max absolute scale");
    auto kinetic = scalar_change(2.0, std::nullopt, std::nullopt);
    check(kinetic.status == Availability::unavailable && !kinetic.value &&
          !kinetic.comparison && kinetic.evidence.current == 2.0,
          "no predecessor is unavailable, never a false zero");
    auto density = maximum_field_change({{2.0, 3.0, 0.4}, {5.0, 4.0, 0.3}}, previous);
    check(density.value == 0.4 && density.scale.value == 5.0 &&
          density.evidence.current_norm == 5.0 && density.evidence.previous_norm == 4.0,
          "density uses the largest spin-field change and both samples for scale");
    check(maximum_field_change({{2.0, {}, {}}}, {}).status == Availability::unavailable,
          "missing density history cannot pass");
    check(maximum_residual({}).status == Availability::invalid, "empty residual is invalid");
    check(maximum_residual({{"a", -1.0, 1.0}}).status == Availability::invalid,
          "negative norm cannot become measured evidence");
    check(maximum_residual({{"a", 1.0, 1.0}, {"a", 0.5, 1.0}}).status == Availability::invalid,
          "component identities are unique");
    check(maximum_residual({{"a", 1.0, {}}}).scale.status == Availability::unavailable,
          "missing function norm has no invented relative scale");
    check(scalar_change(std::numeric_limits<double>::infinity(), -2.0, previous).status == Availability::invalid,
          "nonfinite physical quantity is invalid");
    check(maximum_field_change({{2.0, 3.0, -0.1}}, previous).status == Availability::invalid,
          "invalid density evidence is not hidden by maximum reduction");
    const ResolvedCriterion scaled{Metric::density_change, 0.1, 0.0, 0.1, LimitRule::scaled};
    const auto outcome = check_metric(scaled, density, {}, {"scf", {1e-6, 8}, 3, "ctx"});
    check(outcome.result == CheckResult::met && outcome.limit && std::abs(*outcome.limit - 0.6) < 1e-15,
          "constructed evidence works with the shared scaled rule");
    for (std::size_t updates : {0u, 1u}) {
        IterationMetrics state{{"scf", {1e-6, 8}, updates, "ctx"},
            {{"ground", {}, {{Metric::bsh_residual, MetricValue::measured(0.1, 1.0)}}}}};
        ResolvedSettings settings{{1e-6, 8}, 1, 0,
            {{Metric::bsh_residual, 0.5, 0.0, 0.5, LimitRule::plain}}, {}};
        int prepared = 0;
        auto endpoint = madness::scf_continuation::assess_prepared_endpoint(
            [&] { return decide(settings, {"scf", "ctx", {"ground"}}, {updates, updates, true}, state, true); },
            [&] { ++prepared; state.states[0].metrics[Metric::bsh_residual] = MetricValue::measured(0.8, 1.0); },
            [](const Decision& decision) { return decision.action != Action::update; });
        check(prepared == 1 && endpoint.action == (updates == 0 ? Action::update : Action::stop),
              "shared decision revokes a preliminary pass after final preparation");
    }
    return failures ? 1 : 0;
}
