#pragma once

// Which saved ground-state Fock matrix a prepare at level (thresh, k) may read.
//
// A ground-state node's archive/scf.fock.json holds one level: the one its
// description records, written in the same save as the orbitals. Its entry is
// named "thresh: <thresh> k: <k>" with six decimals, so levels below 1e-6 can
// share a name (1e-7 and 1e-8 both print as 0.000000). With the node's level
// known, the decision is therefore made on the numbers, never on the name.
// A standalone moldft file (no node) may hold several levels and is matched by
// name, as before.
#include <madness/chem/IterationMetrics.h>
#include <nlohmann/json.hpp>
#include <cstdio>
#include <optional>
#include <string>

namespace molresponse_v3 {

// The entry name SCF::fock_record writes for a level.
inline std::string fock_level_name(double thresh, int k) {
  return "thresh: " + std::to_string(thresh) + " k: " + std::to_string(k);
}

struct FockChoice {
  std::string entry;   // the entry to read; empty: rebuild from the orbitals
  std::string reason;  // why nothing is read, when `entry` is empty
};

// `file` is the saved Fock JSON; `node_level` the level in the node
// description, or nothing for a standalone file.
inline FockChoice choose_saved_fock(const nlohmann::json& file,
                                    const std::optional<madness::iteration::Protocol>& node_level,
                                    double thresh, int k) {
  if (node_level) {
    if (node_level->thresh != thresh || node_level->k != k) {
      char saved[64];
      std::snprintf(saved, sizeof saved, "thresh %.3g k %d", node_level->thresh, node_level->k);
      return {"", std::string("the node was saved at ") + saved};
    }
    if (!file.is_object() || file.size() != 1)
      return {"", "the node's Fock file does not hold exactly one level"};
    return {file.begin().key(), ""};
  }
  const auto name = fock_level_name(thresh, k);
  if (file.is_object() && file.contains(name)) return {name, ""};
  std::string levels;
  if (file.is_object())
    for (const auto& item : file.items()) levels += (levels.empty() ? "" : ", ") + item.key();
  return {"", "the file holds " + (levels.empty() ? std::string("no level") : levels)};
}

}  // namespace molresponse_v3
