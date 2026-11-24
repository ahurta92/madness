#pragma once
// Saved response states of a madqc run folder, for the inspection tools and
// the kernel tests. A tool lists the node folders under <run>/states/response/,
// reads each description through the node module, picks by kind and settings,
// and loads through the node loaders. It never rebuilds a node name and never
// writes: a node is read without finishing or discarding a save in progress,
// so a tool can look at a run that is still going.
#include "../GroundState.hpp"
#include "../ResponseProtocol.hpp"
#include "shared_es_checkpoint.hpp"

#include <madness/chem/SCFCheckpoint.h>  // ground_state_node
#include <madness/chem/molecule.h>

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace molresponse_v3::saved_states {
namespace fs = std::filesystem;
namespace it = madness::iteration;
namespace cp = madness::checkpoint;
using json = nlohmann::json;

// One saved node (or seed) and its description.
struct SavedNode {
  fs::path folder;
  it::NodeDescription description;
};

// The node at `folder`. Nothing when the folder holds neither a description
// nor an archive (a launch that has not saved yet); a folder that holds one of
// them but does not read as a node stops. Collective.
inline std::optional<SavedNode> read_saved_node(madness::World& world, const fs::path& folder) {
  const auto present = cp::detail::on_root(world, [&] {
    return std::string(fs::exists(folder / "description.json") || fs::exists(it::Node(folder).archive()) ? "1" : "");
  });
  if (present.empty()) return std::nullopt;
  auto d = response_shared::read_other_node(world, folder);
  cp::records::require(d.has_value(), "half-saved or unreadable node (needs a readable description.json and "
                                      "archive/): " + folder.string());
  return SavedNode{folder, std::move(*d)};
}

// Every saved node under <run>/states/response/, in folder-name order. Stops
// when the run folder has no states/response/. Collective.
inline std::vector<SavedNode> response_nodes(madness::World& world, const fs::path& run) {
  const auto root = fs::absolute(run).lexically_normal() / "states" / "response";
  const auto names = cp::detail::on_root(world, [&] {
    cp::records::require(fs::is_directory(root), "no saved response states: " + root.string() + " is not a folder");
    std::vector<std::string> found;
    for (const auto& entry : fs::directory_iterator(root))
      if (entry.is_directory()) found.push_back(entry.path().filename().string());
    std::sort(found.begin(), found.end());
    return json(found).dump();
  });
  std::vector<SavedNode> nodes;
  for (const auto& name : json::parse(names))
    if (auto node = read_saved_node(world, root / name.get<std::string>())) nodes.push_back(std::move(*node));
  return nodes;
}

// The completely saved nodes under <run>/states/response/, in folder-name
// order, each read without changing it. A folder with a save under way or
// left unfinished, or with no readable description, is passed over, and a run
// folder with no states/response/ has none. Not collective: for the root
// process (the response scheduler's snapshot).
inline std::vector<SavedNode> completed_response_nodes(const fs::path& run) {
  const auto root = fs::absolute(run).lexically_normal() / "states" / "response";
  std::vector<SavedNode> nodes;
  if (!fs::is_directory(root)) return nodes;
  std::vector<fs::path> folders;
  for (const auto& entry : fs::directory_iterator(root))
    if (entry.is_directory()) folders.push_back(entry.path());
  std::sort(folders.begin(), folders.end());
  for (const auto& folder : folders) {
    if (!response_shared::completely_saved(folder)) continue;
    if (auto d = response_shared::read_node_description(folder)) nodes.push_back({folder, std::move(*d)});
  }
  return nodes;
}

// The node or seed at `folder`, which must hold a saved state. Collective.
inline SavedNode require_node(madness::World& world, const fs::path& folder) {
  auto node = read_saved_node(world, fs::absolute(folder));
  cp::records::require(node.has_value(), "no saved node in " + folder.string());
  return std::move(*node);
}

// A seed's settings keep the state's own settings under "state".
inline const json& state_settings(const SavedNode& node) {
  const auto& s = node.description.settings;
  return node.description.kind == "seed" ? s.at("state") : s;
}
inline std::string layout(const SavedNode& node) { return state_settings(node).at("layout").get<std::string>(); }
inline bool closed_shell(const SavedNode& node) { return layout(node).rfind("closed_shell_", 0) == 0; }
inline bool full_layout(const SavedNode& node) { return layout(node).ends_with("_full"); }

inline bool is_driven(const SavedNode& node) { return node.description.kind == "driven"; }
// An excited-state set, or a seed for one.
inline bool is_excited(const SavedNode& node) {
  const auto& d = node.description;
  return d.kind == "excited_set" ||
         (d.kind == "seed" && d.details.value("seed_for", std::string()) == "excited_set");
}

// Half the cubic cell edge the ground state was solved in, from the ground
// model the node's settings hold.
inline double cell_half_width(const SavedNode& node) {
  const auto& s = node.description.settings;
  const auto& ground = s.contains("ground_settings") ? s.at("ground_settings") : s.at("ground_model");
  return ground.at("representation_settings").at("cell").at("values").at(1).get<double>();
}

// One line per loaded node: folder, kind, protocol, save time.
inline void print_loaded(madness::World& world, const SavedNode& node) {
  if (world.rank() != 0) return;
  const auto& d = node.description;
  madness::print("  loaded", node.folder.string(), " kind=", d.kind,
                 " protocol=", protocol_key(d.protocol.thresh, d.protocol.k), " saved_at=", d.saved_at);
}

inline it::Protocol active_protocol() {
  return {madness::FunctionDefaults<3>::get_thresh(), madness::FunctionDefaults<3>::get_k()};
}

// A driven state's functions at the active level. Storage is
// ResponseStateX<ClosedShell> for a static layout, ResponseStateXY for a full
// one. Collective.
template <class Storage>
Storage load_driven(madness::World& world, const SavedNode& node) {
  cp::records::require(is_driven(node), "not a driven-state node: " + node.folder.string());
  return response_shared::load_driven_functions<Storage>(world, it::Node(node.folder).archive(),
                                                         node.description, active_protocol());
}

// A driven state as (X, Y) at the active level. A static state (X only) gives
// Y = X, as a dynamic state started from a static one does. Collective.
inline ResponseStateXY<ClosedShell> load_driven_xy(madness::World& world, const SavedNode& node) {
  cp::records::require(closed_shell(node), "not a closed-shell driven state: " + node.folder.string());
  if (full_layout(node)) return load_driven<ResponseStateXY<ClosedShell>>(world, node);
  auto x = load_driven<ResponseStateX<ClosedShell>>(world, node);
  ResponseStateXY<ClosedShell> xy;
  xy.x_alpha = madness::copy(world, x.x_alpha);
  xy.y_alpha = std::move(x.x_alpha);
  return xy;
}

// The saved closed-shell driven state of `perturbation` (e.g. "dipole_z") at
// `frequency` (to the 1e-5 the node names hold), or nothing. Collective.
inline std::optional<SavedNode> find_driven(madness::World& world, const fs::path& run,
                                            const std::string& perturbation, double frequency) {
  for (auto& node : response_nodes(world, run)) {
    if (!is_driven(node) || !closed_shell(node)) continue;
    const auto& s = node.description.settings;
    if (s.at("perturbation").get<std::string>() == perturbation &&
        std::abs(s.at("frequency").get<double>() - frequency) < 5e-6)
      return std::move(node);
  }
  return std::nullopt;
}

// The saved closed-shell excited-state set of the run with a full (RPA) or TDA
// layout, or nothing. Collective.
inline std::optional<SavedNode> find_excited(madness::World& world, const fs::path& run, bool full) {
  for (auto& node : response_nodes(world, run))
    if (node.description.kind == "excited_set" && closed_shell(node) && full_layout(node) == full)
      return std::move(node);
  return std::nullopt;
}

// The ground state of a run folder, read from its node without changing it:
//   <run>/states/ground_state/            the node: description (level, box) and archive,
//                                         with the Fock matrix of its level when the
//                                         description's details name it (details.fock)
//   <run>/analysis/ground_state/          the SCF stage's analysis files, where the
//                                         node layout keeps moldft.calc_info.json
//                                         (the molecule)
// The analysis folder is the workflow's analysis/<calculation id> folder
// (Drivers.hpp, Applications.hpp); no library function names it.
struct GroundNode {
  SavedNode node;
  fs::path archive;        // the name the native loader takes (no .00000 suffix)
  madness::Molecule molecule;
  std::string fock_json;   // empty when there is none: the Fock matrix is built from the orbitals
  double L = 0.0;          // half the cubic cell edge
  [[nodiscard]] it::Protocol protocol() const { return node.description.protocol; }
};

// Collective. Stops when the run has no saved ground state, the node is
// half-saved or unreadable, or the molecule cannot be found.
inline GroundNode read_ground(madness::World& world, const fs::path& run) {
  const auto folder = madness::ground_state_node(fs::absolute(run).lexically_normal());
  auto node = read_saved_node(world, folder);
  cp::records::require(node.has_value(), "no saved ground state in " + folder.string());
  cp::records::require(node->description.kind == "ground_state",
                       "not a ground-state node: " + folder.string());
  GroundNode ground;
  ground.archive = it::Node(folder).archive() / "scf.restartdata";
  ground.L = node->description.details.at("orbitals").at("L").get<double>();
  const auto analysis = fs::absolute(run).lexically_normal() / "analysis" / "ground_state";
  const auto found = json::parse(cp::detail::on_root(world, [&] {
    const auto info = analysis / "moldft.calc_info.json";
    cp::records::require(fs::is_regular_file(info), "no ground-state molecule: " + info.string() + " is missing");
    json calc_info;
    std::ifstream(info) >> calc_info;
    json molecule;
    if (calc_info.contains("tasks") && calc_info["tasks"].is_array() && !calc_info["tasks"].empty())
      molecule = calc_info["tasks"][0].value("molecule", json());
    else
      molecule = calc_info.value("molecule", json());
    cp::records::require(!molecule.is_null(), "no molecule in " + info.string());
    // A node saved before the Fock file existed names none: the matrix is then
    // built from the orbitals.
    const auto name = node->description.details.value("fock", std::string());
    const auto fock = it::Node(folder).archive() / name;
    return json{{"molecule", molecule},
                {"fock", !name.empty() && fs::is_regular_file(fock) ? fock.string() : std::string()}}.dump();
  }));
  ground.molecule.from_json(found.at("molecule"));
  ground.fock_json = found.at("fock").get<std::string>();
  ground.node = std::move(*node);
  return ground;
}

// The ground-state orbitals at level (thresh, k): FunctionDefaults are set to
// that level first, so the node's own level loads them without projection.
// Prints the loaded-node line. Collective.
inline GroundState load_ground(madness::World& world, const GroundNode& ground, double thresh, int k) {
  set_response_protocol(world, ground.L, thresh, k);
  print_loaded(world, ground.node);
  auto gs = GroundState::from_archive(world, ground.archive.string(), ground.molecule);
  gs.set_fock_level(ground.protocol());  // its Fock file is read only at the node's level
  return gs;
}

// An excited-state set's (or seed's) roots at the active level. Collective.
template <class Type>
typename ESSolver<Type, ClosedShell>::State load_excited(madness::World& world, const SavedNode& node) {
  cp::records::require(is_excited(node), "not an excited-state set or seed: " + node.folder.string());
  return response_shared::load_excited_functions<Type, ClosedShell>(world, it::Node(node.folder).archive(),
                                                                   node.description, active_protocol());
}
}  // namespace molresponse_v3::saved_states
