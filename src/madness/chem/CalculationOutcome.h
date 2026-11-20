#pragma once

#include <madness/chem/IterationMetrics.h>
#include <stdexcept>
#include <string>

namespace madness {

/// Application completion used by workflow handoffs. This is not restart
/// evidence: the producing application must qualify its state and finish any
/// required publication before returning it.
struct CalculationOutcome {
    bool complete = false;
    std::string reason = "not_run";

    /// For a single iterative task. Applications with several tasks also own
    /// the prerequisite and requested-property completion checks.
    static CalculationOutcome from_iteration(const iteration::Decision& decision) {
        if (decision.action != iteration::Action::complete &&
            decision.action != iteration::Action::stop)
            throw std::invalid_argument("application outcome requires a terminal decision");
        return {decision.action == iteration::Action::complete, decision.reason};
    }
};

} // namespace madness
