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

#include <madness/external/nlohmann_json/json.hpp>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <algorithm>
#include <filesystem>
#include <string>
#include <type_traits>

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
        {"iter",         state.iter},
        {"bsh_residual", bsh_res},
        {"archive",      archive_basename},
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

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_FD_SAVE_LOAD_HPP
