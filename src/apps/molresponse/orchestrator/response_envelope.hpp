#pragma once
// Task-record envelope for the response task (unified task record, Phase 0 spec
// §4.1), derived from the response_metadata document that the workflow already
// writes. The SCF task has carried `precision` / `convergence` since Phase 0; the
// response task did not, so a regression test had no single key that says "this
// run converged". This header gives it one.
//
//   precision   = {k, thresh, protocol: [thresh per rung, coarse -> fine],
//                  protocol_key} of the FINEST rung
//   convergence = {status, iterations, n_states, n_unconverged, stop_reason}
//
// FINEST means physically finest, not "highest `protocols[*].index`": every
// production writer registers its rung with index = -1 (fd_save_load.hpp,
// es_save_load.hpp, vbc_save_load.hpp, dalton_import.hpp), so the index ties
// across the whole registry and cannot order anything. Rungs are therefore
// ordered by accuracy — thresh descending, then k ascending on a tie — and the
// last one is the finest.
//
// Only canonical ladder keys take part. `protocol_key(thresh, k)`
// (ResponseProtocol.hpp) formats "%.0e_k%d", e.g. "1e-06_k8"; a dalton.dir seed
// additionally registers a SYNTHETIC sibling "<key>_dseed" at the same
// (thresh, k) (dalton_import.hpp) whose FD entries are placeholders
// (converged = false, iter = 0). That is a seed bundle, not a rung, so keys that
// do not match the canonical pattern are excluded from both the finest pick and
// the `protocol` array — otherwise a fully converged seeded run would report
// protocol_key "…_dseed", status "unconverged", iterations 0, and a duplicate
// threshold in `protocol`.
//
// status is "converged" iff at least one state exists at the finest rung, every
// FD / excited-state / VBC entry there reports converged == true, and
// run_summary.stop_reason == "complete"; "unconverged" otherwise; "unknown" when
// the metadata names no ladder protocol (nothing ran). iterations is the maximum
// over those states (`iter`, or `metrics.iters` for VBC entries), -1 when none.

#include <nlohmann/json.hpp>

#include <algorithm>
#include <regex>
#include <string>
#include <vector>

namespace molresponse_v3 {

struct ResponseEnvelope {
  nlohmann::json precision;    // null when no protocol ran
  nlohmann::json convergence;  // always an object
};

/// True for a canonical protocol ladder key ("%.0e_k%d", e.g. "1e-06_k8").
/// False for the dalton seed's synthetic "<key>_dseed" sibling, and for
/// anything else that is not a rung of the ladder.
inline bool is_ladder_protocol_key(const std::string &key) {
  static const std::regex canonical(R"(^[0-9]e[-+][0-9]{2}_k[0-9]+$)");
  return std::regex_match(key, canonical);
}

inline ResponseEnvelope response_task_envelope(const nlohmann::json &md) {
  using json = nlohmann::json;
  ResponseEnvelope env;
  const std::string stop_reason =
      (md.contains("run_summary") && md["run_summary"].is_object())
          ? md["run_summary"].value("stop_reason", std::string("unknown"))
          : std::string("unknown");
  env.convergence = {{"status", "unknown"}, {"iterations", -1}, {"n_states", 0},
                     {"n_unconverged", 0}, {"stop_reason", stop_reason}};
  if (!md.contains("protocols") || !md["protocols"].is_object() || md["protocols"].empty())
    return env;

  struct Rung { double thresh; int k; std::string key; };
  std::vector<Rung> rungs;
  for (const auto &[key, p] : md["protocols"].items()) {
    if (!p.is_object()) continue;                // a malformed entry is not a rung
    if (!is_ladder_protocol_key(key)) continue;  // seed bundles are not rungs
    rungs.push_back({p.value("thresh", 0.0), p.value("k", 0), key});
  }
  if (rungs.empty()) return env;
  // Coarse -> fine: the looser threshold first; on a tie the larger k is finer.
  std::sort(rungs.begin(), rungs.end(), [](const Rung &a, const Rung &b) {
    if (a.thresh != b.thresh) return a.thresh > b.thresh;
    return a.k < b.k;
  });
  const Rung &fine = rungs.back();
  json protocol = json::array();
  for (const auto &r : rungs) protocol.push_back(r.thresh);
  env.precision = {{"k", fine.k}, {"thresh", fine.thresh},
                   {"protocol", protocol}, {"protocol_key", fine.key}};

  int n_states = 0, n_bad = 0, iterations = 0;
  auto visit = [&](const json &st) {
    if (!st.is_object()) return;
    ++n_states;
    if (!st.value("converged", false)) ++n_bad;
    int it = st.value("iter", 0);
    if (st.contains("metrics") && st["metrics"].is_object())
      it = std::max(it, st["metrics"].value("iters", 0));
    iterations = std::max(iterations, it);
  };
  if (md.contains("fd_states") && md["fd_states"].is_object())
    for (const auto &[name, by_key] : md["fd_states"].items())
      if (by_key.is_object() && by_key.contains(fine.key) && by_key[fine.key].is_object())
        for (const auto &[fkey, st] : by_key[fine.key].items()) visit(st);
  if (md.contains("excited_states") && md["excited_states"].is_object() &&
      md["excited_states"].contains(fine.key))
    visit(md["excited_states"][fine.key]);
  if (md.contains("vbc_states") && md["vbc_states"].is_object())
    for (const auto &[name, by_key] : md["vbc_states"].items())
      if (by_key.is_object() && by_key.contains(fine.key)) visit(by_key[fine.key]);

  env.convergence["iterations"] = n_states ? iterations : -1;
  env.convergence["n_states"] = n_states;
  env.convergence["n_unconverged"] = n_bad;
  env.convergence["status"] =
      (n_states > 0 && n_bad == 0 && stop_reason == "complete") ? "converged" : "unconverged";
  return env;
}

} // namespace molresponse_v3
