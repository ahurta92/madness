#ifndef MOLRESPONSE_V3_ORCHESTRATOR_RESPONSE_WORKFLOW_HPP
#define MOLRESPONSE_V3_ORCHESTRATOR_RESPONSE_WORKFLOW_HPP

// -----------------------------------------------------------------------------
// run_response — the single public entry point (doc 16, L1).
//
// One self-contained call: given a checkpoint + a response plan + settings, load
// the ground state, drive CalcManager, assemble Tier-A properties, and return a
// structured Output. It owns no scheduling, ownership, or file-format logic —
// those live in CalcManager and the persistence layer. main.cpp, the test
// runner, and (later, R3) madqc all build a ResponseWorkflowInput and call this,
// so the engine is testable, scriptable, and Python-bindable through one seam.
//
// R0a scope: wraps EXACTLY today's flow (no behavior change). The timing /
// diagnostics / exports Output slots are reserved here and filled by R1 / R2.
// -----------------------------------------------------------------------------

#include "../GroundState.hpp"
#include "../ResponseProtocol.hpp"          // set_response_protocol
#include "../ResponsePropertyPlanner.hpp"   // ResponsePlan
#include "../calc/calc_executor.hpp"        // ExecutorSettings/Context, FdResponseExecutor, assemble_*
#include "../calc/calc_manager.hpp"         // CalcManager
#include "../solvers/response_metadata.hpp"

#include <madness/external/nlohmann_json/json.hpp>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <filesystem>
#include <fstream>
#include <optional>
#include <string>
#include <vector>

namespace molresponse_v3 {

// ---------------------------------------------------------------------------
// Public contract (doc 16). Input is World-free and self-contained; Output is
// pure JSON so madqc/tests/Python consume one object.
// ---------------------------------------------------------------------------

struct ResponseWorkflowInput {
  std::string         archive_file;   // SCF/moldft checkpoint (ground state)
  std::vector<double> protocols;      // coarse->fine truncation-threshold ladder
  // R0a: the Tier-A plan is pre-built by the caller (plan_one/merge_plans, plus
  // any es-full override). FUTURE: carry std::vector<ResponsePropertyRequest>
  // and lower to a plan inside run_response once the TDA/Full choice is a
  // request field (needed by madqc, R3).
  ResponsePlan        plan;
  ExecutorSettings    settings;       // convergence policy + ES/seed/accept knobs (World-free)
  std::optional<madness::Molecule> molecule;  // else built from the archive dir
};

struct ResponseWorkflowOutput {
  nlohmann::json properties;   // Tier-A tensors (mirror of metadata "properties")
  nlohmann::json metadata;     // full response_metadata.json
  nlohmann::json timing;       // R1 — reserved
  nlohmann::json diagnostics;  // R1 — reserved
  nlohmann::json exports;      // R2 — reserved
  nlohmann::json debug_log;    // optional iteration-level
  int            rc = 0;
};

// ---------------------------------------------------------------------------
// Small helpers factored out of the old driver (the archive-adjacent lookups).
// ---------------------------------------------------------------------------

/// Build a Molecule from moldft/mad.calc_info.json next to the archive. Returns
/// an empty Molecule if none is found (harmless for pure-alpha; nuclear/Raman
/// needs it for the per-atom expansion). Mirrors the old test-runner logic.
inline madness::Molecule molecule_from_archive_dir(const std::string &archive_file) {
  madness::Molecule molecule;
  const auto dir = std::filesystem::path(archive_file).parent_path();
  for (const auto &name : {"moldft.calc_info.json", "mad.calc_info.json"}) {
    const auto candidate = dir / name;
    if (!std::filesystem::exists(candidate)) continue;
    std::ifstream ifs(candidate);
    nlohmann::json j;
    ifs >> j;
    nlohmann::json mol_json;
    if (j.contains("tasks") && j["tasks"].is_array() && !j["tasks"].empty())
      mol_json = j["tasks"][0]["molecule"];
    else if (j.contains("molecule"))
      mol_json = j["molecule"];
    if (!mol_json.is_null()) molecule.from_json(mol_json);
    break;
  }
  return molecule;
}

/// Resolve the moldft/mad.fock.json path next to the archive ("" if none).
inline std::string fock_json_from_archive_dir(const std::string &archive_file) {
  const auto dir = std::filesystem::path(archive_file).parent_path();
  for (const auto &name : {"moldft.fock.json", "mad.fock.json"}) {
    const auto candidate = dir / name;
    if (std::filesystem::exists(candidate)) return candidate.string();
  }
  return {};
}

// ---------------------------------------------------------------------------
// run_response — load ground -> plan -> CalcManager -> assemble -> Output.
// ---------------------------------------------------------------------------

inline ResponseWorkflowOutput
run_response(madness::World &world, const ResponseWorkflowInput &in) {
  ResponseWorkflowOutput out;
  MADNESS_CHECK(!in.protocols.empty());

  // 1. Protocol + ground state. set_response_protocol to the coarsest rung up
  //    front (the executor re-prepares per protocol during the solve).
  const auto header = GroundState::read_archive_header(world, in.archive_file);
  set_response_protocol(world, header.L, in.protocols.front());

  const madness::Molecule molecule =
      in.molecule ? *in.molecule : molecule_from_archive_dir(in.archive_file);
  auto gs = GroundState::from_archive(world, in.archive_file, molecule);
  if (world.rank() == 0) gs.print_info();

  const std::string fock_json = fock_json_from_archive_dir(in.archive_file);
  const double cur_thresh = madness::FunctionDefaults<3>::get_thresh();
  auto coulop = madness::poperatorT(
      madness::CoulombOperatorPtr(world, gs.params().lo(), 0.001 * cur_thresh));
  gs.prepare(world, 0.001 * cur_thresh, coulop, fock_json);

  // 2. Drive the calc manager.
  CalcManager::Policy mgr_policy;
  mgr_policy.max_iters_per_step = in.settings.max_iters;
  CalcManager mgr(in.plan, in.settings.calc_dir, mgr_policy);
  mgr.build(molecule.natom());

  ExecutorContext ctx(world, gs, header.L, fock_json, in.settings);
  FdResponseExecutor exec(ctx);
  mgr.run(world, exec);

  // 3. Tier-A property assembly (off the solve path). beta OR raman fill
  //    plan.vbc and share assemble_beta; alpha only for the plain-FD path.
  //    es-only runs (plan.es non-empty, no vbc) assemble no scalar here.
  if (!in.plan.vbc.empty())
    assemble_beta(ctx, in.plan, in.protocols.back());
  else if (in.plan.es.empty())
    assemble_alpha(ctx, in.plan, in.protocols.back());

  // 4. Collect the Output from the aggregate metadata (rank 0 authoritative).
  if (world.rank() == 0) {
    auto meta = ResponseMetadata::load_or_create(
        in.settings.calc_dir + "/response_metadata.json");
    out.metadata = meta.json();
    if (out.metadata.contains("properties"))
      out.properties = out.metadata["properties"];
  }
  return out;
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_ORCHESTRATOR_RESPONSE_WORKFLOW_HPP
