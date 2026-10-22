// Table of cases for plan_scf_continuation (DESIGN B5, B6): where a named
// moldft stage continues from. Pure: no MPI, no files.

#include <madness/chem/ScfContinuation.h>

#include <iostream>
#include <string>
#include <vector>

namespace {
using namespace madness::scf_continuation;

struct Case {
  std::string name;
  SavedRung saved;
  std::vector<RequestedRung> request;
  bool saved_meets_request;
  Action action;
  std::size_t start_rung;
  int iterations_done;
  std::string unmet;
};

const Rung coarse{1e-4, 6};
const Rung fine{1e-6, 8};
const Rung middle{1e-5, 7};

} // namespace

int main() {
  const std::vector<RequestedRung> ladder = {{coarse, 30}, {fine, 30}};
  const std::vector<Case> cases = {
      {"converged final rung that meets the request is loaded",
       {fine, 5, RungStop::converged},
       ladder,
       true,
       Action::load,
       1,
       0,
       ""},
      {"interrupted inside rung 2 resumes rung 2 with its count",
       {fine, 3, RungStop::in_progress},
       ladder,
       false,
       Action::resume,
       1,
       3,
       ""},
      {"interrupted inside rung 1 resumes rung 1 with its count",
       {coarse, 4, RungStop::in_progress},
       ladder,
       false,
       Action::resume,
       0,
       4,
       ""},
      {"finished rung 1 is never redone: continue at rung 2",
       {coarse, 7, RungStop::converged},
       ladder,
       false,
       Action::resume,
       1,
       0,
       ""},
      {"rung 1 ended at its cap: advance to rung 2",
       {coarse, 30, RungStop::iteration_limit},
       ladder,
       false,
       Action::resume,
       1,
       0,
       ""},
      {"rung 1 stalled: advance to rung 2",
       {coarse, 9, RungStop::stalled},
       ladder,
       false,
       Action::resume,
       1,
       0,
       ""},
      {"final rung used its cap, targets unmet: load and report",
       {fine, 30, RungStop::iteration_limit},
       ladder,
       false,
       Action::load,
       1,
       0,
       "unmet: iteration limit reached at rung 2"},
      {"raised cap gives the final rung more iterations",
       {fine, 30, RungStop::iteration_limit},
       {{coarse, 30}, {fine, 40}},
       false,
       Action::resume,
       1,
       30,
       ""},
      {"killed exactly at the cap counts as finished by the cap",
       {fine, 30, RungStop::in_progress},
       ladder,
       false,
       Action::load,
       1,
       0,
       "unmet: iteration limit reached at rung 2"},
      {"converged final rung, tighter request, cap left: resume there",
       {fine, 5, RungStop::converged},
       ladder,
       false,
       Action::resume,
       1,
       5,
       ""},
      {"rungs matched by value: a new middle rung after a finished coarse one",
       {coarse, 7, RungStop::converged},
       {{coarse, 30}, {middle, 30}, {fine, 30}},
       false,
       Action::resume,
       1,
       0,
       ""},
      {"rungs matched by value: coarser requested rungs are skipped",
       {middle, 2, RungStop::in_progress},
       {{coarse, 30}, {middle, 30}, {fine, 30}},
       false,
       Action::resume,
       1,
       2,
       ""},
      {"saved rung absent from the request: first tighter requested rung",
       {middle, 2, RungStop::in_progress},
       ladder,
       false,
       Action::resume,
       1,
       0,
       ""},
      {"saved state tighter than every requested rung: load if it meets",
       {Rung{1e-7, 9}, 4, RungStop::converged},
       ladder,
       true,
       Action::load,
       1,
       0,
       ""},
  };
  int failures = 0;
  for (const auto& c : cases) {
    const auto plan
        = plan_scf_continuation(c.saved, c.request, c.saved_meets_request);
    const bool ok = plan.action == c.action && plan.start_rung == c.start_rung
                    && plan.iterations_done == c.iterations_done
                    && plan.unmet == c.unmet;
    if (!ok) {
      ++failures;
      std::cerr << "FAILED: " << c.name << " -> "
                << (plan.action == Action::load ? "load" : "resume") << " rung "
                << plan.start_rung << " done " << plan.iterations_done
                << " unmet '" << plan.unmet << "' why: " << plan.why << '\n';
    }
  }
  const auto resumed
      = plan_scf_continuation({fine, 3, RungStop::in_progress}, ladder, false);
  if (describe(resumed, ladder)
      != "resumed at rung 2 (thresh 1e-06, k 8), iteration 4 of that rung") {
    ++failures;
    std::cerr << "FAILED: describe -> " << describe(resumed, ladder) << '\n';
  }
  std::cout << "scf continuation failures: " << failures << '\n';
  return failures ? 1 : 0;
}
