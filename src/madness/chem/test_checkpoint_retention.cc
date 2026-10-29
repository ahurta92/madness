// Checkpoint cadence and retention (DESIGN A3): the pure decisions in
// CheckpointRetention.h, as tables of cases. No MPI, no files.

#include <madness/chem/CheckpointRetention.h>

#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace {
using namespace madness::checkpoint;

int failures = 0;
void check(bool ok, const std::string& what) {
  if (!ok) {
    ++failures;
    std::cerr << "FAILED: " << what << '\n';
  }
}
} // namespace

int main() {
  // publish_now: decision point, total at the rung, this run, maxiter, K.
  struct CadenceCase {
    const char* name;
    bool decision_point;
    int at_rung, this_run, cap, every;
    bool expected;
  };
  const std::vector<CadenceCase> cadence = {
      {"default: inside a rung, nothing", false, 3, 3, 30, 0, false},
      {"default: a rung end or stop", true, 3, 3, 30, 0, true},
      {"K = 3: the third iteration", false, 3, 3, 30, 3, true},
      {"K = 3: the fourth iteration", false, 4, 4, 30, 3, false},
      {"K = 3: the sixth iteration", false, 6, 6, 30, 3, true},
      {"K = 3 across --continue: total 6, second of this run", false, 6, 2, 30,
       3, true},
      {"K = 3 across --continue: total 5, first of this run", false, 5, 1, 30,
       3, false},
      {"K = 1: every iteration", false, 1, 1, 30, 1, true},
      {"maxiter used up without a stop sample", false, 1, 1, 1, 0, true},
      {"below maxiter, default", false, 29, 29, 30, 0, false},
      {"a stop with K set", true, 4, 4, 30, 3, true},
  };
  for (const auto& c : cadence)
    check(publish_now(c.decision_point, c.at_rung, c.this_run, c.cap, c.every)
              == c.expected,
          std::string("publish_now: ") + c.name);

  // generations_to_delete (D-070): g1..g5 of one stage.
  const std::vector<std::string> store = {"g-1", "g-2", "g-3", "g-4", "g-5"};
  struct RetentionCase {
    const char* name;
    std::string latest;
    std::set<std::string> referenced;
    std::vector<std::string> expected;
  };
  const std::vector<RetentionCase> retention = {
      {"every generation but latest goes",
       "g-5",
       {},
       {"g-1", "g-2", "g-3", "g-4"}},
      {"a latest in the middle keeps only itself",
       "g-3",
       {},
       {"g-1", "g-2", "g-4", "g-5"}},
      {"a referenced generation is kept",
       "g-5",
       {"g-3"},
       {"g-1", "g-2", "g-4"}},
      {"everything referenced or latest: nothing goes",
       "g-5",
       {"g-1", "g-2", "g-3", "g-4"},
       {}},
  };
  for (const auto& c : retention)
    check(generations_to_delete(store, c.latest, c.referenced) == c.expected,
          std::string("generations_to_delete: ") + c.name);
  check(generations_to_delete({"g-4"}, "g-4").empty(),
        "generations_to_delete: only latest: no candidates, no scan needed");
  check(generations_to_delete({}, "g-1").empty(),
        "generations_to_delete: an empty store");

  std::cout << "checkpoint retention failures: " << failures << '\n';
  return failures ? 1 : 0;
}
