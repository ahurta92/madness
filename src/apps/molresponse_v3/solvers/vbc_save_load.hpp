#ifndef MOLRESPONSE_V3_SOLVERS_VBC_SAVE_LOAD_HPP
#define MOLRESPONSE_V3_SOLVERS_VBC_SAVE_LOAD_HPP

// =========================================================================
// VBC persistence (beta-ii-b). The VBC quadratic source is a
// ResponseStateXY<ClosedShell> (same Storage as an FD Full state); we save it
// as a parallel archive under the calc dir plus a vbc_states/<vbc_id>/<key>
// entry in response_metadata.json. VBC is rebuilt from its converged FD inputs
// at each protocol, so there is no try_load_vbc — reconcile's Skip (converged at
// this protocol_key) is what avoids recomputation; restart is handled by the
// FD inputs already being on disk.
// =========================================================================

#include "../ResponseProtocol.hpp"   // protocol_key()
#include "response_metadata.hpp"
#include "response_state.hpp"
#include "state_metrics.hpp"

#include <madness/external/nlohmann_json/json.hpp>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <filesystem>
#include <string>

namespace molresponse_v3 {

/// Filesystem-safe archive basename for a VBC node at a protocol:
///   "vbc:dipole_x__dipole_y@fB_fC"  ->  "vbc_dipole_x__dipole_y_fB_fC__<key>"
inline std::string vbc_archive_basename(const std::string &vbc_id,
                                        const std::string &key) {
  std::string base = vbc_id;
  for (char &ch : base)
    if (ch == ':' || ch == '@') ch = '_';
  return base + "__" + key;
}

/// Save one VBC source: collective archive write + rank-0 metadata upsert.
/// Assumes the active FunctionDefaults<3> (k, thresh) is the build protocol.
inline void save_vbc_state(madness::World &world,
                           const ResponseStateXY<ClosedShell> &vbc,
                           const std::string &dir,
                           const std::string &vbc_id,
                           bool converged) {
  if (world.rank() == 0) std::filesystem::create_directories(dir);
  world.gop.fence();

  const double      thresh = madness::FunctionDefaults<3>::get_thresh();
  const int         k_now  = madness::FunctionDefaults<3>::get_k();
  const std::string key    = protocol_key(thresh, k_now);
  const std::string base   = vbc_archive_basename(vbc_id, key);
  const std::string archive_path = dir + "/" + base;

  vbc.save(world, archive_path);  // collective
  const StateMetrics metrics = measure_state(world, vbc, /*iter=*/0);

  if (world.rank() == 0) {
    auto meta = ResponseMetadata::load_or_create(dir + "/response_metadata.json");
    if (!meta.json()["protocols"].contains(key))
      meta.set_protocol(key, thresh, k_now, /*index=*/-1);
    nlohmann::json entry = {
        {"converged", converged},
        {"diverged",  false},
        {"archive",   base},
        {"metrics",   metrics.to_json()},
    };
    meta.set_vbc_state(vbc_id, key, entry);
    meta.save();
    madness::print("[SAVE] vbc_state: id=", vbc_id, "  protocol_key=", key,
                   "  archive=", base, "  converged=", converged);
  }
  world.gop.fence();
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_VBC_SAVE_LOAD_HPP
