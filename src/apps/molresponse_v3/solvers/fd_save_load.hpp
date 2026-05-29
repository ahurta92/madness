#ifndef MOLRESPONSE_V3_SOLVERS_FD_SAVE_LOAD_HPP
#define MOLRESPONSE_V3_SOLVERS_FD_SAVE_LOAD_HPP

// =========================================================================
// FD save/load (Inc 13c-i: save only; load + restart precedence land in
// 13c-ii / 13c-iii).
//
// FD persistence is per-point — one Storage archive per (perturbation,
// protocol, freq) triple. The unit of persistence is therefore simpler than
// the ES bundle: a single ResponseStateX<Shell> or ResponseStateXY<Shell>
// archive plus a metadata entry in the unified response_metadata.json.
//
// Disk layout (per doc 13):
//
//   <calc dir>/
//   ├── response_metadata.json                      ← 13b aggregator
//   ├── dipole_x__1e-06_k8__f0.05700.00000          ← FD point (parallel
//   ├── dipole_x__1e-06_k8__f0.05700.00001            archive, .000N per rank)
//   └── es_bundle__1e-06_k8/                        ← ES bundle
//       ├── roots.json
//       └── root_0.00000 ...
//
// response_filename(pert, protocol_key, freq) is the join key for both the
// archive name and the metadata path (fd_states/<pert>/<key>/<fkey>) so a
// property layer can find FD and ES inputs by string compare on protocol_key.
//
// Collective discipline: state.responses[0].save() is collective on `world`.
// The metadata upsert is rank-0 only (filesystem op) and the function fences
// after both finish. Caller does NOT need to add its own fence.
// =========================================================================

#include "../Perturbations.hpp"          // Perturbation
#include "../ResponseProtocol.hpp"        // protocol_key()
#include "../kernels/tags.hpp"            // Static/Full/TDA, ClosedShell/OpenShell
#include "fd_solver.hpp"
#include "response_metadata.hpp"
#include "response_state.hpp"
#include "state_metrics.hpp"              // measure_state (per-state mem/iters)

#include <madness/external/nlohmann_json/json.hpp>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <algorithm>
#include <filesystem>
#include <optional>
#include <string>
#include <type_traits>
#include <vector>

namespace molresponse_v3 {

namespace detail_fd_save_load {

template <typename Type>
inline const char *type_tag() {
  if constexpr (std::is_same_v<Type, Static>) return "static";
  else if constexpr (std::is_same_v<Type, Full>) return "full";
  else if constexpr (std::is_same_v<Type, TDA>)  return "tda";
  else return "unknown";
}

template <typename Shell>
inline const char *shell_tag() {
  if constexpr (std::is_same_v<Shell, ClosedShell>) return "closed_shell";
  else if constexpr (std::is_same_v<Shell, OpenShell>) return "open_shell";
  else return "unknown";
}

} // namespace detail_fd_save_load

/// Canonical archive name for an FD point. Per doc 13:
///   <pert>__<protocol_key>__<freq_key>      dipole_x__1e-06_k8__f0.05700
inline std::string response_filename(const std::string &pert,
                                     const std::string &protocol_key_str,
                                     double freq) {
  return pert + "__" + protocol_key_str + "__"
       + ResponseMetadata::freq_key(freq);
}

/// Save one FD point. Writes the Storage archive (collective) and upserts
/// the entry into response_metadata.json (rank-0). Assumes the active
/// FunctionDefaults<3> (k, thresh) is the protocol the state was solved at —
/// the caller's protocol step has already set this.
///
/// Invariant: state.responses.size() == 1. The skeleton (and present v3 FD
/// driver) solve one perturbation/freq per State; multi-channel save is a
/// later concern.
template <typename Type, typename Shell>
void save_fd_state(madness::World &world,
                   const typename FDSolver<Type, Shell>::State &state,
                   const std::string &dir,
                   const Perturbation &pert,
                   double freq,
                   bool converged) {
  MADNESS_CHECK(state.responses.size() == 1);

  if (world.rank() == 0) {
    std::filesystem::create_directories(dir);
  }
  world.gop.fence();

  const double      thresh = madness::FunctionDefaults<3>::get_thresh();
  const int         k_now  = madness::FunctionDefaults<3>::get_k();
  const std::string key    = protocol_key(thresh, k_now);
  const std::string fkey   = ResponseMetadata::freq_key(freq);
  const std::string pdesc  = pert.description();
  const std::string archive_basename = response_filename(pdesc, key, freq);
  const std::string archive_path     = dir + "/" + archive_basename;

  // (1) Collective binary save — same per-state primitive ES uses per-root.
  state.responses[0].save(world, archive_path);

  // (1b) Collective per-state metrics (coeffs/bytes/rss/iters) — every rank.
  const StateMetrics metrics =
      measure_state(world, state.responses[0], state.iter);

  // (2) Rank-0 metadata upsert.
  if (world.rank() == 0) {
    auto meta = ResponseMetadata::load_or_create(
        dir + "/response_metadata.json");

    // Register the protocol if this is the first artifact at this key.
    // `index` is unknown to this writer (it's a property of the caller's
    // ramp); a later orchestrator can overwrite it with the real ordering.
    if (!meta.json()["protocols"].contains(key)) {
      meta.set_protocol(key, thresh, k_now, /*index=*/-1);
    }

    const double bsh_res =
        state.last_bsh_residual.empty()
            ? 0.0
            : *std::max_element(state.last_bsh_residual.begin(),
                                state.last_bsh_residual.end());

    nlohmann::json entry = {
        {"freq",         freq},
        {"type",         detail_fd_save_load::type_tag<Type>()},
        {"shell",        detail_fd_save_load::shell_tag<Shell>()},
        {"converged",    converged},
        {"diverged",     state.diverged},
        {"iter",         state.iter},
        {"bsh_residual", bsh_res},
        {"archive",      archive_basename},
        {"metrics",      metrics.to_json()},
    };
    meta.set_fd_state(pdesc, key, fkey, entry);
    meta.save();

    madness::print("[SAVE] fd_state: pert=", pdesc,
                   "  protocol_key=", key,
                   "  freq=", freq,
                   "  archive=", archive_basename,
                   "  bsh_res=", bsh_res,
                   "  converged=", converged);
  }
  world.gop.fence();
}

/// Exact-match load (Inc 13c-ii). Looks up the archive at
///   dir/response_filename(pert.description(), protocol_key(), freq)
/// using the ACTIVE FunctionDefaults<3> (k, thresh) — caller is expected to
/// have called set_response_protocol() first so the key matches.
///
/// Returns a fresh State with one response slot: responses[0] = loaded
/// Storage; iter=0; diverged=false; last_bsh_residual seeded from the
/// metadata entry's bsh_residual (for diagnostic continuity — the solver
/// recomputes on the first step). Density-residual history is not
/// persisted; the first-iter Δρ guard inside FDSolver::step() handles a
/// missing rho_alpha_prev (same path as ESSolver after es_load).
///
/// Throws if either the archive or the metadata entry is missing.
/// Cross-protocol fallback (load at a lower protocol than the active one)
/// is deferred to try_load_fd_state in 13c-iii.
template <typename Type, typename Shell>
typename FDSolver<Type, Shell>::State
load_fd_state(madness::World &world,
              const std::string &dir,
              const Perturbation &pert,
              double freq) {
  using State   = typename FDSolver<Type, Shell>::State;
  using Storage = typename FDSolver<Type, Shell>::Storage;

  const std::string key   = protocol_key();
  const std::string fkey  = ResponseMetadata::freq_key(freq);
  const std::string pdesc = pert.description();
  const std::string archive_basename = response_filename(pdesc, key, freq);
  const std::string archive_path     = dir + "/" + archive_basename;

  // Rank-0 reads the metadata entry (required); broadcasts the small
  // numeric diagnostics. The archive load below is collective.
  double bsh_res        = 0.0;
  int    iter_at_save   = 0;
  int    converged_int  = 0;   // bool via int for gop.broadcast
  if (world.rank() == 0) {
    const std::string meta_path = dir + "/response_metadata.json";
    if (!std::filesystem::exists(meta_path)) {
      throw std::runtime_error(
          "load_fd_state: missing " + meta_path);
    }
    auto meta = ResponseMetadata::load_or_create(meta_path);
    const auto &j = meta.json();
    if (!j["fd_states"].contains(pdesc) ||
        !j["fd_states"][pdesc].contains(key) ||
        !j["fd_states"][pdesc][key].contains(fkey)) {
      throw std::runtime_error(
          "load_fd_state: no fd_states/" + pdesc + "/" + key + "/" + fkey +
          " in " + meta_path);
    }
    const auto &e = j["fd_states"][pdesc][key][fkey];
    bsh_res       = e.value("bsh_residual", 0.0);
    iter_at_save  = e.value("iter", 0);
    converged_int = e.value("converged", false) ? 1 : 0;
  }
  world.gop.broadcast(bsh_res,       0);
  world.gop.broadcast(iter_at_save,  0);
  world.gop.broadcast(converged_int, 0);

  // Collective binary load — same primitive ES uses per-root.
  State s;
  s.responses.resize(1);
  s.responses[0] = Storage::load(world, archive_path);

  s.last_bsh_residual     = {bsh_res};
  s.last_density_residual = {};      // not persisted; first-iter recomputes
  s.rho_alpha_prev        .clear();  // first-iter Δρ guard handles empty
  s.iter                  = 0;
  s.diverged              = false;

  if (world.rank() == 0) {
    madness::print("[LOAD] fd_state: pert=", pdesc,
                   "  protocol_key=", key,
                   "  freq=", freq,
                   "  archive=", archive_basename,
                   "  bsh_res_at_save=", bsh_res,
                   "  iter_at_save=", iter_at_save,
                   "  converged_at_save=", converged_int != 0);
  }
  return s;
}

/// Result of try_load_fd_state. When `exact == false`, the caller is
/// responsible for re-projecting `state` to the active (k, thresh) — the
/// natural place is the first `prepare(...)` call inside iterate_protocol,
/// which already does this for every step.
template <typename Type, typename Shell>
struct FDRestartResult {
  typename FDSolver<Type, Shell>::State state;
  std::string source_protocol_key;
  bool        exact = false;
};

/// Restart precedence (Inc 13c-iii). Reads response_metadata.json and picks
/// the best on-disk bundle for (pert, freq):
///
///   1. exact match at the active protocol_key   → exact = true
///   2. else any saved (thresh, k) that's COARSER-OR-EQUAL to active
///      (saved_thresh >= active_thresh AND saved_k <= active_k):
///      pick max k, then min thresh — closest to active = best initial guess
///   3. else nullopt — caller seeds fresh
///
/// Returns nullopt if no metadata file or no compatible record. The
/// archive load itself is collective, so this function must be called by
/// every rank.
template <typename Type, typename Shell>
std::optional<FDRestartResult<Type, Shell>>
try_load_fd_state(madness::World &world,
                  const std::string &dir,
                  const Perturbation &pert,
                  double freq) {
  using State   = typename FDSolver<Type, Shell>::State;
  using Storage = typename FDSolver<Type, Shell>::Storage;

  const std::string active_key = protocol_key();
  const double active_thresh   = madness::FunctionDefaults<3>::get_thresh();
  const int    active_k        = madness::FunctionDefaults<3>::get_k();
  const std::string fkey       = ResponseMetadata::freq_key(freq);
  const std::string pdesc      = pert.description();

  // Rank-0 picks the source_key. Empty string ↔ no compatible match.
  std::string source_key;
  if (world.rank() == 0) {
    const std::string meta_path = dir + "/response_metadata.json";
    if (std::filesystem::exists(meta_path)) {
      auto meta = ResponseMetadata::load_or_create(meta_path);
      const auto &j = meta.json();
      if (j["fd_states"].contains(pdesc) &&
          j.contains("protocols") && j["protocols"].is_object()) {
        struct Cand { std::string key; double thresh; int k; };
        std::vector<Cand> cands;
        for (const auto &[key, ent] : j["fd_states"][pdesc].items()) {
          if (!ent.contains(fkey))             continue;
          if (!j["protocols"].contains(key))   continue;
          const double t  = j["protocols"][key].value("thresh", 0.0);
          const int    kk = j["protocols"][key].value("k", 0);
          if (t >= active_thresh && kk <= active_k) {
            cands.push_back({key, t, kk});
          }
        }
        // Exact match wins.
        for (const auto &c : cands) {
          if (c.key == active_key) { source_key = c.key; break; }
        }
        // Else closest-to-active: max k, then min thresh.
        if (source_key.empty() && !cands.empty()) {
          auto best = std::max_element(
              cands.begin(), cands.end(),
              [](const Cand &a, const Cand &b) {
                if (a.k != b.k)      return a.k < b.k;       // higher k wins
                return a.thresh > b.thresh;                  // smaller thresh wins
              });
          source_key = best->key;
        }
      }
    }
  }
  world.gop.broadcast_serializable(source_key, 0);

  if (source_key.empty()) {
    if (world.rank() == 0) {
      madness::print("[LOAD] try_load_fd_state: no compatible bundle for pert=",
                     pdesc, " freq=", freq, " in ", dir,
                     " (active_key=", active_key, ")");
    }
    return std::nullopt;
  }

  // Collective load.
  const std::string archive_basename = response_filename(pdesc, source_key, freq);
  const std::string archive_path     = dir + "/" + archive_basename;

  FDRestartResult<Type, Shell> r;
  r.state.responses.resize(1);
  r.state.responses[0] = Storage::load(world, archive_path);
  r.state.last_bsh_residual    .clear();
  r.state.last_density_residual.clear();
  r.state.rho_alpha_prev       .clear();
  r.state.iter     = 0;
  r.state.diverged = false;
  r.source_protocol_key = source_key;
  r.exact = (source_key == active_key);

  if (world.rank() == 0) {
    madness::print("[LOAD] try_load_fd_state: pert=", pdesc,
                   "  freq=", freq,
                   "  source_protocol_key=", source_key,
                   "  active=", active_key,
                   "  exact=", r.exact,
                   "  archive=", archive_basename);
  }
  return r;
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_FD_SAVE_LOAD_HPP
