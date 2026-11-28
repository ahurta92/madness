#ifndef MOLRESPONSE_V3_CALC_NODE_SNAPSHOT_HPP
#define MOLRESPONSE_V3_CALC_NODE_SNAPSHOT_HPP

// ===========================================================================
// The scheduler's snapshot: what is saved, built from the node folders at the
// start of every scheduler pass (CalcManager::run), in the shape the pure
// scheduler functions read (calc_manager.hpp):
//
//   protocols/<key>                     {thresh, k}
//   fd_states/<pert>/<key>/<freq key>   a driven node
//   excited_states/<key>                an excited-state set
//   vbc_states/<vbc id>/<key>           a VBC quadratic source
//
// <key> is the protocol key of the level the node is saved at; a node holds
// one level. `converged` is not read from anywhere: the caller judges each
// node against the current request (Session::node_meets_request) and passes
// the verdict in.
//
// Pure: no World, no filesystem.
// ===========================================================================

#include "calc_manager.hpp"

#include <madness/chem/NodeArchive.h>

#include <nlohmann/json.hpp>

#include <string>
#include <vector>

namespace molresponse_v3 {

/// One saved node as the snapshot sees it: its folder name under
/// states/response/, its description, and whether it meets the current
/// request (judged by the caller).
struct SnapshotNode {
  std::string                        name;
  madness::iteration::NodeDescription description;
  bool                               converged = false;
};

namespace detail_snapshot {

/// The maximum BSH residual stored for a driven node, 0 when it has none.
inline double stored_bsh_residual(const madness::iteration::NodeDescription &d) {
  namespace it = madness::iteration;
  try {
    const auto frame = it::stored_metrics(d, {"snapshot", d.protocol, 0, ""});
    if (frame.states.empty()) return 0.0;
    const auto found = frame.states.front().metrics.find(it::Metric::bsh_residual);
    if (found == frame.states.front().metrics.end()) return 0.0;
    return found->second.value.value_or(0.0);
  } catch (const std::exception &) {
    return 0.0;
  }
}

/// A driven node's name split into its perturbation and frequency key
/// ("dipole_z-f0.10000"): a perturbation description never holds a hyphen,
/// so the first one separates them (driven_node_name).
inline std::pair<std::string, std::string> split_driven_name(const std::string &name) {
  const auto dash = name.find('-');
  if (dash == std::string::npos) return {};
  return {name.substr(0, dash), name.substr(dash + 1)};
}

} // namespace detail_snapshot

/// The snapshot of `nodes`. A driven node needs settings.perturbation and
/// settings.frequency; an excited-state set settings.excited_model and
/// details.roots; a VBC source settings.kind "VBC" and the lineage entries of
/// roles B and C. A node without them, and every other kind, is left out.
/// Two excited-state sets saved at the same level share one entry; the later
/// in `nodes` wins (the shape keys sets by level alone).
inline nlohmann::json node_snapshot(const std::vector<SnapshotNode> &nodes) {
  using nlohmann::json;
  json s{{"protocols", json::object()},
         {"fd_states", json::object()},
         {"excited_states", json::object()},
         {"vbc_states", json::object()}};
  for (const auto &node : nodes) {
    const auto &d = node.description;
    const std::string key = protocol_key(d.protocol.thresh, d.protocol.k);
    const std::string folder = "states/response/" + node.name;
    try {
      if (d.kind == "driven") {
        const auto pert = d.settings.at("perturbation").get<std::string>();
        const double freq = d.settings.at("frequency").get<double>();
        s["fd_states"][pert][key][freq_key(freq)] = {
            {"converged", node.converged},
            {"freq", freq},
            {"frequency", freq},
            {"iter", d.updates_at_protocol},
            {"bsh_residual", detail_snapshot::stored_bsh_residual(d)},
            {"diverged", false},
            {"node", folder}};
      } else if (d.kind == "excited_set") {
        const auto &roots = d.details.at("roots");
        s["excited_states"][key] = {
            {"converged", node.converged},
            {"type", d.settings.at("excited_model").get<std::string>()},
            {"shell", "closed_shell"},
            {"n_roots", roots.size()},
            {"roots", roots},
            {"iter", d.updates_at_protocol},
            {"diverged", false},
            {"node", folder}};
      } else if (d.kind == "second_order" && d.settings.is_object() &&
                 d.settings.value("kind", std::string()) == "VBC") {
        std::string b, c;
        for (const auto &input : d.lineage) {
          const auto role = input.at("role").get<std::string>();
          if (role == "B") b = input.at("name").get<std::string>();
          if (role == "C") c = input.at("name").get<std::string>();
        }
        const auto [pb, fb] = detail_snapshot::split_driven_name(b);
        const auto [pc, fc] = detail_snapshot::split_driven_name(c);
        if (pb.empty() || pc.empty()) continue;
        // vbc_node_id from the two inputs' names.
        const std::string id = "vbc:" + pb + "__" + pc + "@" + fb + "_" + fc;
        s["vbc_states"][id][key] = {{"converged", node.converged},
                                    {"diverged", false},
                                    {"node", folder}};
      } else {
        continue;
      }
    } catch (const std::exception &) {
      continue;  // a description without the fields above is not scheduled from
    }
    s["protocols"][key] = {{"thresh", d.protocol.thresh}, {"k", d.protocol.k}};
  }
  return s;
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_CALC_NODE_SNAPSHOT_HPP
