#pragma once

#include <madness/chem/ConvergenceSettings.h>
#include <utility>

namespace madness::scf_continuation {

/// A possible stop is provisional until endpoint preparation and reevaluation.
/// Both callbacks capture application-owned state and comparison history. The
/// preparation callback replaces cached equation work; assessment never advances
/// the solver. A failed final check may therefore return in_progress.
template <class Assess, class PrepareAndEvaluate, class IsEndpoint>
auto assess_prepared_endpoint(Assess&& assess, PrepareAndEvaluate&& prepare_and_evaluate,
                              IsEndpoint&& is_endpoint) {
    auto decision = assess();
    if (is_endpoint(decision)) {
        prepare_and_evaluate();
        decision = assess();
    }
    return decision;
}

template <class Assess, class PrepareAndEvaluate>
auto assess_prepared_endpoint(Assess&& assess, PrepareAndEvaluate&& prepare_and_evaluate) {
    return assess_prepared_endpoint(std::forward<Assess>(assess),
        std::forward<PrepareAndEvaluate>(prepare_and_evaluate), [](const auto& decision) {
            return decision.rung_stop != convergence::StopReason::in_progress;
        });
}

} // namespace madness::scf_continuation
