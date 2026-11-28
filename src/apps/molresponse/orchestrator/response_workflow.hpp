#ifndef MOLRESPONSE_V3_ORCHESTRATOR_RESPONSE_WORKFLOW_HPP
#define MOLRESPONSE_V3_ORCHESTRATOR_RESPONSE_WORKFLOW_HPP

// -----------------------------------------------------------------------------
// run_response_with_ground — the response workflow (doc 16, L1).
//
// Given a prepared ground state, a response plan and settings that carry the
// response session, drive CalcManager, assemble Tier-A properties, and return
// a structured Output. It owns no scheduling, ownership, or file-format logic —
// those live in CalcManager and the session. The madqc adapter
// (madqc_adapter.hpp) builds the ResponseWorkflowInput and calls it.
// -----------------------------------------------------------------------------

#include "../GroundState.hpp"
#include "../ResponseProtocol.hpp"          // set_response_protocol
#include "../ResponsePropertyPlanner.hpp"   // ResponsePlan
#include "../calc/calc_executor.hpp"        // ExecutorSettings/Context, FdResponseExecutor, assemble_*
#include "../solvers/raman_analysis.hpp"    // full Raman: Hessian + normal-mode spectrum
#include "../calc/calc_manager.hpp"         // CalcManager

#include <nlohmann/json.hpp>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <cstdlib>
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
  // The Tier-A plan is pre-built by the caller (plan_one/merge_plans, plus
  // any es-full override).
  ResponsePlan        plan;
  ExecutorSettings    settings;       // the session, convergence policy and ES knobs
};

struct ResponseWorkflowOutput {
  nlohmann::json timing;       // R1 — stage wall/cpu times (plan_build, solve, assemble)
  nlohmann::json diagnostics;  // the scheduler's stop reason and dropped work
  nlohmann::json exports;      // R2 — reserved
  nlohmann::json debug_log;    // optional iteration-level
  int            rc = 0;
};

// ---------------------------------------------------------------------------
// Stage timing (R1a, doc 16 L3). Rank-0 wall/cpu deltas — the "shape of the
// run" (is time in solving or assembly? which stage dominates?). wall_time() /
// cpu_time() are LOCAL (not collective); stage boundaries are bracketed by the
// collectives inside each stage, so rank-0's wall is representative. Per-rank-max
// reduction is a later refinement.
// ---------------------------------------------------------------------------
namespace detail_workflow {
struct StageTimer {
  double w0 = madness::wall_time();
  double c0 = madness::cpu_time();
  nlohmann::json lap() const {
    return nlohmann::json{{"wall_s", madness::wall_time() - w0},
                          {"cpu_s", madness::cpu_time() - c0}};
  }
};
} // namespace detail_workflow

// ---------------------------------------------------------------------------
// run_response_with_ground — plan -> CalcManager -> assemble -> Output.
// ---------------------------------------------------------------------------

/// Core (stages 2–4): DAG build → CalcManager solve → Tier-A assembly → collect.
/// `gs` must already be loaded AND prepared at the coarsest protocol; `L` is the
/// cubic-cell half-edge and `fock_json` the ground-state node's Fock file ("" = compute).
/// Fills `Output.timing` for plan_build/solve/assemble (the caller owns load +
/// total). Called by the madqc adapter.
inline ResponseWorkflowOutput
run_response_with_ground_impl(madness::World &world, GroundState &gs, double L,
                         const std::string &fock_json,
                         const ResponseWorkflowInput &in) {
  ResponseWorkflowOutput out;
  MADNESS_CHECK(!in.protocols.empty());
  // DFT + quadratic response gate. The quadratic sources (V^BC for beta and
  // Raman, the derived-FD legs + residue for 2PA) need the SECOND xc kernel
  // g''_xc[rho^B rho^C] in addition to f_xc. g''_xc is exactly zero for HF
  // and is not implemented for any functional — running would silently drop
  // it and return plausible-looking wrong numbers. Refuse instead. Linear
  // response (alpha) and excited states are fully supported with f_xc.
  if (gs.scf().xc.is_dft() &&
      (!in.plan.vbc.empty() || !in.plan.derived_fd.empty()))
    throw std::runtime_error(
        "quadratic response (beta / Raman / 2PA) on a DFT ground state needs "
        "the second xc kernel g''_xc, which is not implemented. Use an HF "
        "ground state for quadratic properties, or request linear response / "
        "excited states (f_xc supported) with this functional.");
  nlohmann::json timing;

  // 2a. Build the DAG.
  detail_workflow::StageTimer t_build;
  MADNESS_CHECK_THROW(in.settings.shared_response, "response workflow requires a response session");
  in.settings.shared_response->timing().begin("plan_build");
  CalcManager::Policy mgr_policy;
  mgr_policy.fd_subworlds       = in.settings.fd_subworlds;   // F2
  mgr_policy.fd_subworld_ranks  = in.settings.fd_subworld_ranks;
  mgr_policy.state_plans        = in.settings.state_plans;     // D-016
  mgr_policy.verbosity          = in.settings.verbosity;
  {
    auto session=in.settings.shared_response;
    mgr_policy.saved_states=[session]{return detail_exec::session_snapshot(*session);};
    mgr_policy.collect_subworlds=[session]{session->collect_subworlds();};
    mgr_policy.dependencies_ready=[session](const CalcNode& node) {
      return std::all_of(node.prerequisites.begin(),node.prerequisites.end(),[&](const auto& id) {
        return session->complete(id);
      });
    };
    mgr_policy.qualified_bundle=[session](const std::string& id){return session->excited_at(id).complete();};
    mgr_policy.plan_derived=[session,&in](const CalcNode& node) {
      const auto& bundle=session->excited_at(node.prerequisites.front());
      auto& binding=session->bind_derived(node.id,{
          {"frequency",node.freq},{"perturbation",node.pert.description()},
          {"layout",detail_exec::is_static_freq(node.freq)?"closed_shell_static":"closed_shell_full"},
          {"exchange_tensor",in.settings.policy.exchange_tensor},
          {"derived_from",{{"bundle",bundle.info().task_id},{"root",node.es_root_id},{"factor",node.es_freq_factor}}}});
      StatePlan entry;
      entry.kind=binding.reusable()?StatePlan::Kind::load : binding.has_saved()?StatePlan::Kind::resume:StatePlan::Kind::fresh;
      entry.start_thresh=binding.controls().protocol.thresh;
      if(binding.reusable()) {entry.met=true;entry.loaded_key=protocol_key(binding.controls().protocol.thresh,binding.controls().protocol.k);}
      return entry;
    };
  }

  CalcManager mgr(in.plan, in.settings.calc_dir, mgr_policy);
  mgr.build(gs.molecule().natom());
  timing["plan_build"] = t_build.lap();
  in.settings.shared_response->timing().end();

  // 2b. Drive the calc manager (the solve).
  detail_workflow::StageTimer t_solve;
  ExecutorContext ctx(world, gs, L, fock_json, in.settings);
  FdResponseExecutor exec(ctx);

  // F2 (doc 32 §5): the per-subworld solve closure. Only built when an archive is
  // available — GroundState's only construction path is from_archive (doc 32 §3),
  // so the madqc in-memory-GS path can't fan out (it leaves fan_out null → the
  // single-World path). Each node-subworld loads its own GS and solves its owned
  // FD items into their node folders. This is exactly the F1-proven block,
  // lifted into the live run().
  CalcManager::SubworldSolve fan_out;
  const bool fan_requested =
      in.settings.fd_subworlds > 0 || in.settings.fd_subworld_ranks > 0;
  if (fan_requested && in.archive_file.empty() && world.rank() == 0)
    madness::print("F2: fd_subworlds =", in.settings.fd_subworlds,
                   " fd_subworld_ranks =", in.settings.fd_subworld_ranks,
                   "requested but no archive_file — staying single-World "
                   "(subworlds reload the ground state from the archive)");
  if (fan_requested && !in.archive_file.empty()) {
    const std::string  archive = in.archive_file;
    const std::string  fock = fock_json;
    const auto         fock_level = gs.fock_level();
    const double       L_   = L;
    const ExecutorSettings base = in.settings;
    fan_out = [archive, fock, fock_level, L_, base](
                  madness::World &sub, const std::vector<WorkItem> &items,
                  double thresh, int gid, const std::string &log_prefix) {
      std::set<std::string> states;
      for(const auto& item:items) {
        if(item.node->kind==CalcKind::VBC) {
          states.insert(fd_node_id(item.node->pert,item.node->freq));
          states.insert(fd_node_id(item.node->pert_c,item.node->freq_c));
        } else states.insert(item.node->id);
      }
      auto shared=base.shared_response->for_world(sub,states,gid);
      try {
        shared->timing().begin("reference_load");
        set_response_protocol(sub,L_,thresh,-1,false);
        // The molecule and occupation come from the ground-state node.
        const auto& ground_settings=base.shared_response->ground().identity.at("settings");
        auto gs_s=GroundState::from_archive(sub,archive,madness::saved_molecule(ground_settings),nullptr,nullptr,
                                            madness::saved_occupation(ground_settings));
        if(fock_level) gs_s.set_fock_level(*fock_level);
        shared->timing().end();
        ExecutorSettings settings=base;
        settings.shared_response=shared;
        settings.log_prefix=log_prefix;
        ExecutorContext context(sub,gs_s,L_,fock,settings);
        FdResponseExecutor executor(context);
        for(const auto& item:items) executor.run_protocol(item);
        shared->timing().end();shared->timing().end();shared->publish_timings();
      } catch(const std::exception& error) {
        try {shared->publish_timings(error.what());} catch(...) {}
        throw;
      }
    };
  }
  const auto sched_diag=mgr.run(world,exec,fan_out);
  timing["solve"] = t_solve.lap();

  // 3. Tier-A property assembly (off the solve path). A mixed request may want
  //    several properties at once (e.g. alpha + beta), so assemble each that the
  //    plan supports — not one-XOR-the-other. assemble_alpha self-guards (returns
  //    if no dipole FD; it ignores nuclear FD) so it is safe to always attempt;
  //    beta/raman come from VBC. (ES-derived scalar properties are a later step.)
  detail_workflow::StageTimer t_assemble;
  in.settings.shared_response->timing().begin("properties");
  assemble_alpha(ctx, in.plan, in.protocols.back());
  for(const auto& bundle:in.plan.es) {
    if(bundle.tda) assemble_shared_excited<TDA>(ctx,bundle.n_roots,in.protocols.back());
    else assemble_shared_excited<Full>(ctx,bundle.n_roots,in.protocols.back());
  }
  if (!in.plan.vbc.empty())
    assemble_beta(ctx, in.plan, in.protocols.back());
  // Full Raman: the plan carries the all-atoms nuclear sentinel only when the
  // request asked for every coordinate, so its presence IS the gate for the
  // normal-mode analysis (Hessian from the nuclear response states, then the
  // per-mode spectrum from the raman rows assemble_beta just wrote).
  if (plan_requests_full_raman(in.plan)) {
    in.settings.shared_response->timing().begin("vibrational_analysis");
    assemble_raman_spectrum(ctx,in.plan,in.protocols.back());
    in.settings.shared_response->timing().end();
  }
  // 2PA: derived FD legs are only planned when the request asked for the
  // two-photon contraction (request.tpa), so their presence IS the gate.
  // assemble_tpa self-guards (returns if no ES bundle; ClosedShell only).
  if (!in.plan.derived_fd.empty())
    assemble_tpa(ctx, in.plan, in.protocols.back());
  // The two-photon states of earlier excitation energies go once the current
  // roots' states are saved (step 4d); every process group has finished.
  if (!in.plan.derived_fd.empty()) {
    const auto removed = in.settings.shared_response->prune_two_photon();
    if (world.rank() == 0 && !removed.empty() &&
        madness::iteration::prints(in.settings.verbosity, madness::iteration::needed::checkpoint_saved))
      for (const auto& name : removed)
        std::cout << "Removed: " << name << ", a two-photon state no current root needs\n";
  }
  timing["assemble"] = t_assemble.lap();
  {
    auto& session=*in.settings.shared_response;
    session.timing().end();
    session.timing().begin("publish_results");
    session.set_run({{"stop_reason",sched_diag.at("stop_reason")},
                     {"dropped_work",sched_diag.value("dropped_beta",nlohmann::json::array())}});
    session.report();
    session.timing().end();
    session.timing().end(); // calculation
    session.timing().end(); // run
    session.publish_timings();
  }
  if (!in.settings.shared_response->complete())
    throw std::runtime_error("response convergence criteria unmet; native endpoint and unavailable property recorded");

  // 4. The Output: the stage times and the scheduler's verdict. The values
  // are in results.json; the archive I/O backend is its provenance.archive_io.
  out.timing = std::move(timing);
  out.diagnostics = sched_diag;
  return out;
}

inline ResponseWorkflowOutput
run_response_with_ground(madness::World& world,GroundState& gs,double L,
                         const std::string& fock_json,const ResponseWorkflowInput& in) {
  MADNESS_CHECK_THROW(in.settings.shared_response, "response workflow requires a response session");
  const auto work=[&] {return run_response_with_ground_impl(world,gs,L,fock_json,in);};
  return in.settings.shared_response->with_failure_timing(work);
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_ORCHESTRATOR_RESPONSE_WORKFLOW_HPP
