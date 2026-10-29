#pragma once

// When a stage publishes a checkpoint and which older generations it keeps
// (DESIGN A3). The adapter's observer asks these two questions at every save
// boundary; no engine counts iterations for them.
//
// Pure: no World, MPI, archives or files. The adapter reads the manifests and
// the storage layer deletes; these functions only decide.

#include <set>
#include <string>
#include <vector>

namespace madness::checkpoint {

/// Publish a generation at this save boundary? (A3, Adrian 2026-10-23)
///   - always at a decision point: a rung end or a stop (`decision_point`);
///   - within a rung only when the deck opted in with
///     checkpoint.every_iterations K: every K-th iteration at the rung,
///     counted across --continue runs (`iterations_at_rung`), so the cadence
///     holds over a restart;
///   - when this sample used the rung's whole maxiter in this run. An engine
///     can end a rung on such a sample without a separate stop sample (SCF
///     tests no exit on its first iteration), and the rung's state must not
///     be lost.
[[nodiscard]] inline bool publish_now(bool decision_point,
                                      int iterations_at_rung,
                                      int iterations_this_run, int cap,
                                      int every_iterations) {
  if (decision_point)
    return true;
  if (cap > 0 && iterations_this_run >= cap)
    return true;
  return every_iterations > 0 && iterations_at_rung > 0
         && iterations_at_rung % every_iterations == 0;
}

/// The generations retention deletes after a successful publish (A3, D-070:
/// one copy of each saved function): every one of this stage's generations
/// except the one `latest` names and any another stage's manifest references
/// (A1). With no `referenced` set this gives the candidates, so a caller can
/// skip reading other stages when there are none.
[[nodiscard]] inline std::vector<std::string>
generations_to_delete(const std::vector<std::string>& generations,
                      const std::string& latest,
                      const std::set<std::string>& referenced = {}) {
  std::vector<std::string> doomed;
  for (const auto& generation : generations)
    if (generation != latest && !referenced.contains(generation))
      doomed.push_back(generation);
  return doomed;
}

} // namespace madness::checkpoint
