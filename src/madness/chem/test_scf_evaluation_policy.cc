#include <madness/chem/SCFEvaluationPolicy.h>
#include <iostream>

int main() {
    using madness::convergence::StopReason;
    int failures = 0;
    auto check = [&](bool ok, const char* message) {
        if (!ok) { ++failures; std::cerr << message << '\n'; }
    };
    // The real stop policy sees a preliminary pass and a failed prepared state.
    // Test both remaining allowance and exhaustion, with one preparation only.
    for (int updates : {0, 1}) {
        double cached_residual = 0.1;
        int preparations = 0, assessments = 0;
        struct Assessment { StopReason rung_stop; double residual; };
        auto assess = [&] {
            ++assessments;
            return Assessment{madness::convergence::rung_exit({
                madness::convergence::Decider::engine_check, cached_residual < 0.5,
                false, false, updates, 1, 0, true}), cached_residual};
        };
        const auto result = madness::scf_continuation::assess_prepared_endpoint(assess, [&] {
            ++preparations;
            cached_residual = 0.8; // fresh equation work after final preparation
        });
        check(preparations == 1 && assessments == 2, "endpoint evaluated once after preparation");
        check(result.residual == 0.8, "continuation retains fresh endpoint residual");
        check(result.rung_stop == (updates == 0 ? StopReason::in_progress : StopReason::iteration_limit),
              "failed endpoint revokes preliminary convergence and respects remaining updates");
    }
    int preparations = 0;
    struct Assessment { StopReason rung_stop; };
    const auto result = madness::scf_continuation::assess_prepared_endpoint(
        [] { return Assessment{StopReason::in_progress}; }, [&] { ++preparations; });
    check(result.rung_stop == StopReason::in_progress && preparations == 0,
          "ordinary continuation does not run endpoint preparation");
    return failures ? 1 : 0;
}
