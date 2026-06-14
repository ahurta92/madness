#ifndef MOLRESPONSE_V3_SOLVERS_NODE_SUBWORLDS_HPP
#define MOLRESPONSE_V3_SOLVERS_NODE_SUBWORLDS_HPP

// =========================================================================
// Node-aligned subworld creation — Inc 1 of the ES node-φ state-parallel
// design (docs/21).
//
// One subworld per PHYSICAL NODE, via MPI_Comm_split_type(SHARED): every rank
// on the same host lands in the same subworld. This is the prerequisite for
// the Stage-2 design — with node-aligned subworlds the ground state can be
// Cloud-copied into each node-subworld ONCE and left DISTRIBUTED over that
// node's ranks, giving exactly one φ copy per node (per-rank φ identical to
// the current single-World cost). That sidesteps the v2 per-RANK replication
// OOM without needing shared-memory Function storage (which MADNESS lacks).
//
// MacroTaskQ's default create_worlds uses color = rank % nsubworld, which
// INTERLEAVES subworlds across nodes (safempi Split). We want the opposite:
// node-CONTIGUOUS subworlds, which Split_type(SHARED) gives directly.
//
// This header is placement/diagnostics only — no solver behavior. The build
// fan-out that uses these subworlds is a later increment (docs/21 §5).
// =========================================================================

#include <madness/world/world.h>
#include <madness/world/worldgop.h>
#include <madness/world/safempi.h>
#include <madness/world/ranks_and_hosts.h>   // ranks_per_host, get_hostname
#include <madness/world/print.h>

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace molresponse_v3 {

/// Diagnostics describing how the universe split into node-aligned subworlds.
struct NodeSubworldInfo {
  std::string hostname;        ///< this rank's host
  int         universe_rank = 0;
  int         universe_size = 0;
  int         subworld_rank = 0;   ///< rank within this node's subworld
  int         subworld_size = 0;   ///< ranks on this node (= ranks in subworld)
  int         n_nodes       = 0;   ///< distinct hosts in the universe
};

/// Create one subworld per physical node (MPI_COMM_TYPE_SHARED). All ranks on
/// the same host land in the same subworld; the `Key = universe.rank()` keeps
/// the within-node ordering deterministic. Collective on `universe`. Fills
/// `info` if non-null. The returned World must be destroyed (reset) before
/// `finalize()` and before the universe.
inline std::shared_ptr<madness::World>
make_node_aligned_subworld(madness::World &universe,
                           NodeSubworldInfo *info = nullptr) {
  SafeMPI::Intracomm node_comm = universe.mpi.comm().Split_type(
      SafeMPI::Intracomm::SHARED_SPLIT_TYPE, /*Key=*/universe.rank());
  auto subworld = std::make_shared<madness::World>(node_comm);
  universe.gop.fence();

  if (info) {
    auto rph = madness::ranks_per_host(universe);   // collective: host -> [ranks]
    info->hostname      = madness::get_hostname();
    info->universe_rank = universe.rank();
    info->universe_size = universe.size();
    info->subworld_rank = subworld->rank();
    info->subworld_size = subworld->size();
    info->n_nodes       = static_cast<int>(rph.size());
  }
  return subworld;
}

/// Print (rank 0) the host -> ranks map and check the "one subworld per node"
/// invariant: on every rank, the subworld size must equal the number of
/// universe ranks sharing that rank's host. Returns true iff the invariant
/// holds on ALL ranks. Collective on `universe`.
inline bool
verify_one_subworld_per_node(madness::World &universe,
                             madness::World &subworld) {
  auto rph = madness::ranks_per_host(universe);   // full map on every rank
  const std::string host = madness::get_hostname();
  const int expected = static_cast<int>(rph.count(host) ? rph[host].size() : 0);
  const int actual   = subworld.size();
  int ok = (expected == actual && expected > 0) ? 1 : 0;
  universe.gop.fence();
  universe.gop.sum(&ok, 1);                        // total #ranks that passed
  const bool all_ok = (ok == universe.size());

  if (universe.rank() == 0) {
    madness::print("=== node-aligned subworld map ===");
    madness::print("  universe ranks =", universe.size(),
                   "   distinct nodes =", static_cast<int>(rph.size()));
    for (const auto &kv : rph)
      madness::print("  node", kv.first, ":", static_cast<int>(kv.second.size()),
                     "ranks ->", kv.second);
    madness::print("  one-subworld-per-node:", all_ok ? "PASS" : "FAIL",
                   " (", ok, "/", universe.size(), "ranks consistent)");
  }
  return all_ok;
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_NODE_SUBWORLDS_HPP
