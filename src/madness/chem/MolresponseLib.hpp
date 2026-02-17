/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680


  $Id$
 */
#pragma once

#include "ResponseParameters.hpp"
#include "madness_exception.h"
#include <apps/molresponse_v2/DerivedStatePlanner.hpp>
#include <apps/molresponse_v2/FrequencyLoop.hpp>
#include <apps/molresponse_v2/GroundStateData.hpp>
#include <apps/molresponse_v2/PropertyManager.hpp>
#include <apps/molresponse_v2/ResponseDebugLogger.hpp>
#include <apps/molresponse_v2/ResponseManager.hpp>
#include <apps/molresponse_v2/ResponseRecord.hpp>
#include <apps/molresponse_v2/StateGenerator.hpp>
#include <apps/molresponse_v2/StateParallelPlanner.hpp>
#include <apps/molresponse_v2/VBCMacrotask.hpp>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <madness/chem/InputWriter.hpp>
#include <madness/chem/ParameterManager.hpp>
#include <madness/chem/Results.h>
#include <madness/mra/macrotaskq.h>

/*
Developer Overview
- Orchestration file for molresponse, structured as three stages:
  1) bootstrap + planning, 2) state solves, 3) property assembly.
- Stage 1: read ground/checkpoint context, generate linear states, build derived
  request plan, and build state-parallel ownership/execution plan.
- Stage 2: solve linear states over protocol thresholds; fresh runs start
  state-oriented, then later thresholds can fan out by frequency-point.
  Restart runs with complete protocol-0 saved data may promote fanout to ti==0.
  When subgroup mode is enabled, work runs in macrotask subworlds with shard
  metadata/log merge and synchronization barriers between protocol steps.
- Stage 2c: evaluate derived-state dependency gate and execute ready requests
  before property assembly (serial or subgroup lanes).
- Stage 3: compute requested properties (alpha/beta/raman) with optional
  property-group execution and broadcast of assembled outputs.
*/

struct molresponse_lib {
  struct Results {
    nlohmann::json metadata;
    nlohmann::json properties;
    nlohmann::json vibrational_analysis;
    nlohmann::json raman_spectra;
    nlohmann::json debug_log;
  };

  static constexpr const char *label() { return "molresponse"; }

private:
  struct GroundContext {
    // Molecule read from checkpoint; reused across all response stages.
    Molecule molecule;
    // Ground-state object used to build and solve response states.
    GroundStateData ground;
    // Runtime manager configured from response input.
    ResponseManager response_manager;
    // Archive and Fock filenames resolved for this run directory.
    std::string archive_file;
    std::string fock_json_file;
  };

  struct PlannedStates {
    // Linear states generated directly from requested perturbations.
    GeneratedStateData generated_states;
    // Derived-state requests (currently VBC-driven scaffolding/execution).
    DerivedStatePlan derived_state_plan;
    // Ownership and subgroup execution plan for Stage 2.
    StateParallelPlan state_parallel_plan;
  };

  struct SolvedStates {
    // Carries full planning context forward into property stage.
    PlannedStates planned_states;
    // Merged response metadata after linear + derived execution.
    nlohmann::json metadata;
    // Combined debug log (including subgroup shards when used).
    nlohmann::json debug_log;
  };

  struct PropertyStageOutput {
    // Aggregated property JSON (alpha/beta/raman blocks).
    nlohmann::json properties;
    // Vibrational artifacts needed for Raman post-processing/output.
    VibrationalResults vibrational_analysis;
    RamanResults raman_spectra;
  };

  class JsonStateSolvePersistence final : public StateSolvePersistence {
  public:
    JsonStateSolvePersistence(World &world, const std::string &meta_file,
                              const std::string &debug_file)
        : response_record_(world, meta_file), debug_logger_(debug_file) {}

    void
    initialize_states(const std::vector<LinearResponseDescriptor> &states) {
      response_record_.initialize_states(states);
    }

    void print_summary() const { response_record_.print_summary(); }

    [[nodiscard]] bool is_saved(const LinearResponsePoint &pt) const override {
      return response_record_.is_saved(pt.perturbationDescription(),
                                       pt.threshold(), pt.frequency());
    }

    [[nodiscard]] bool
    is_converged(const LinearResponsePoint &pt) const override {
      return response_record_.is_converged(pt.perturbationDescription(),
                                           pt.threshold(), pt.frequency());
    }

    void record_status(const LinearResponsePoint &pt, bool c) override {
      response_record_.record_status(pt, c);
    }

    void record_timing(const LinearResponsePoint &pt, double wall_seconds,
                       double cpu_seconds) override {
      response_record_.record_timing(pt, wall_seconds, cpu_seconds);
    }

    ResponseDebugLogger &logger() override { return debug_logger_; }

    void flush_debug_log(World &world) override {
      if (debug_logger_.enabled() && world.rank() == 0) {
        debug_logger_.write_to_disk();
      }
    }

    [[nodiscard]] nlohmann::json metadata_json() const {
      return response_record_.to_json();
    }

    [[nodiscard]] nlohmann::json debug_log_json() const {
      return debug_logger_.to_json();
    }

  private:
    ResponseRecord2 response_record_;
    ResponseDebugLogger debug_logger_;
  };

  static std::string group_shard_file(const std::string &filename, size_t gid) {
    std::filesystem::path in(filename);
    const auto parent = in.parent_path();
    const auto stem = in.stem().string();
    const auto ext = in.extension().string();
    const auto grouped = ext.empty()
                             ? stem + ".group" + std::to_string(gid)
                             : stem + ".group" + std::to_string(gid) + ext;
    if (parent.empty()) {
      return grouped;
    }
    return (parent / grouped).string();
  }

  static std::string group_console_file(size_t gid) {
    return group_shard_file("response_console.log", gid);
  }

  static std::string group_derived_timing_file(size_t gid) {
    return group_shard_file("derived_request_timings.json", gid);
  }

  class ScopedRankLogRedirect final {
  public:
    ScopedRankLogRedirect(bool enabled, const std::string &filename,
                          const std::string &header = "")
        : enabled_(enabled) {
      if (!enabled_) {
        return;
      }
      sink_.open(filename, std::ios::out | std::ios::app);
      if (!sink_) {
        return;
      }

      old_cout_ = std::cout.rdbuf(sink_.rdbuf());
      old_cerr_ = std::cerr.rdbuf(sink_.rdbuf());
      active_ = true;
      if (!header.empty()) {
        std::cout << "\n=== " << header << " ===\n";
      }
    }

    ~ScopedRankLogRedirect() { restore(); }
    ScopedRankLogRedirect(const ScopedRankLogRedirect &) = delete;
    ScopedRankLogRedirect &operator=(const ScopedRankLogRedirect &) = delete;

  private:
    void restore() {
      if (!active_) {
        return;
      }
      std::cout.flush();
      std::cerr.flush();
      std::cout.rdbuf(old_cout_);
      std::cerr.rdbuf(old_cerr_);
      active_ = false;
    }

    bool enabled_ = false;
    bool active_ = false;
    std::ofstream sink_;
    std::streambuf *old_cout_ = nullptr;
    std::streambuf *old_cerr_ = nullptr;
  };

  static void write_json_file(const std::string &filename,
                              const nlohmann::json &json_data) {
    std::ofstream out(filename);
    if (!out) {
      throw std::runtime_error("Cannot open " + filename + " for writing");
    }
    out << std::setw(2) << json_data << "\n";
  }

  static nlohmann::json read_json_file_or_object(const std::string &filename) {
    if (!std::filesystem::exists(filename)) {
      return nlohmann::json::object();
    }
    std::ifstream in(filename);
    if (!in) {
      return nlohmann::json::object();
    }
    nlohmann::json parsed = nlohmann::json::object();
    in >> parsed;
    return parsed;
  }

  static void merge_state_metadata_json(nlohmann::json &merged,
                                        const nlohmann::json &shard) {
    if (!merged.is_object()) {
      merged = nlohmann::json::object();
    }
    if (!merged.contains("states") || !merged["states"].is_object()) {
      merged["states"] = nlohmann::json::object();
    }
    if (!shard.is_object() || !shard.contains("states") ||
        !shard["states"].is_object()) {
      return;
    }

    for (const auto &[state_id, shard_state] : shard["states"].items()) {
      auto &dst_state = merged["states"][state_id];
      if (!dst_state.is_object()) {
        dst_state = nlohmann::json::object();
      }
      if (!dst_state.contains("protocols") ||
          !dst_state["protocols"].is_object()) {
        dst_state["protocols"] = nlohmann::json::object();
      }

      if (shard_state.contains("final_saved") &&
          shard_state["final_saved"].is_boolean()) {
        const bool lhs = dst_state.contains("final_saved") &&
                         dst_state["final_saved"].is_boolean() &&
                         dst_state["final_saved"].get<bool>();
        const bool rhs = shard_state["final_saved"].get<bool>();
        dst_state["final_saved"] = (lhs || rhs);
      }

      if (!shard_state.contains("protocols") ||
          !shard_state["protocols"].is_object()) {
        continue;
      }

      for (const auto &[protocol_key, shard_proto] :
           shard_state["protocols"].items()) {
        auto &dst_proto = dst_state["protocols"][protocol_key];
        if (!dst_proto.is_object()) {
          dst_proto = nlohmann::json::object();
        }
        if (!dst_proto.contains("saved") || !dst_proto["saved"].is_object()) {
          dst_proto["saved"] = nlohmann::json::object();
        }
        if (!dst_proto.contains("converged") ||
            !dst_proto["converged"].is_object()) {
          dst_proto["converged"] = nlohmann::json::object();
        }
        if (!dst_proto.contains("timings") ||
            !dst_proto["timings"].is_object()) {
          dst_proto["timings"] = nlohmann::json::object();
        }

        auto merge_flag_map = [&](const char *name) {
          if (!shard_proto.contains(name) || !shard_proto[name].is_object()) {
            return;
          }
          for (const auto &[freq_key, shard_value] :
               shard_proto[name].items()) {
            const bool rhs =
                shard_value.is_boolean() && shard_value.get<bool>();
            const bool lhs = dst_proto[name].contains(freq_key) &&
                             dst_proto[name][freq_key].is_boolean() &&
                             dst_proto[name][freq_key].get<bool>();
            dst_proto[name][freq_key] = (lhs || rhs);
          }
        };

        merge_flag_map("saved");
        merge_flag_map("converged");
        if (shard_proto.contains("timings") &&
            shard_proto["timings"].is_object()) {
          for (const auto &[freq_key, timing_value] :
               shard_proto["timings"].items()) {
            if (!dst_proto["timings"].contains(freq_key)) {
              dst_proto["timings"][freq_key] = timing_value;
            }
          }
        }
      }
    }
  }

  static void merge_debug_log_recursive(nlohmann::json &dst,
                                        const nlohmann::json &src) {
    if (!src.is_object()) {
      dst = src;
      return;
    }

    if (!dst.is_object()) {
      dst = nlohmann::json::object();
    }

    for (const auto &[key, src_value] : src.items()) {
      if (!dst.contains(key)) {
        dst[key] = src_value;
        continue;
      }

      auto &dst_value = dst[key];
      if (dst_value.is_object() && src_value.is_object()) {
        merge_debug_log_recursive(dst_value, src_value);
      } else if (dst_value.is_array() && src_value.is_array()) {
        for (const auto &entry : src_value) {
          dst_value.push_back(entry);
        }
      } else {
        dst_value = src_value;
      }
    }
  }

  static void merge_debug_log_json(nlohmann::json &merged,
                                   const nlohmann::json &shard) {
    if (!shard.is_object()) {
      return;
    }
    if (!merged.is_object()) {
      merged = nlohmann::json::object();
    }
    merge_debug_log_recursive(merged, shard);
  }

  static bool point_ready_in_metadata(const nlohmann::json &metadata,
                                      const LinearResponsePoint &pt,
                                      bool require_saved,
                                      bool require_converged) {
    if (!metadata.is_object() || !metadata.contains("states") ||
        !metadata["states"].is_object()) {
      return false;
    }

    const auto state_key = pt.perturbationDescription();
    const auto protocol_key = ResponseRecord2::protocol_key(pt.threshold());
    const auto freq_key = ResponseRecord2::freq_key(pt.frequency());

    const auto states_it = metadata["states"].find(state_key);
    if (states_it == metadata["states"].end()) {
      return false;
    }
    if (!states_it->contains("protocols") ||
        !(*states_it)["protocols"].is_object()) {
      return false;
    }
    const auto protos_it = (*states_it)["protocols"].find(protocol_key);
    if (protos_it == (*states_it)["protocols"].end()) {
      return false;
    }

    auto check_flag = [&](const char *flag_name, bool required) {
      if (!required) {
        return true;
      }
      if (!protos_it->contains(flag_name) ||
          !(*protos_it)[flag_name].is_object()) {
        return false;
      }
      const auto values_it = (*protos_it)[flag_name].find(freq_key);
      if (values_it == (*protos_it)[flag_name].end()) {
        return false;
      }
      return values_it->is_boolean() && values_it->get<bool>();
    };

    return check_flag("saved", require_saved) &&
           check_flag("converged", require_converged);
  }

  static GroundContext
  make_ground_context(World &world, const CalculationParameters &calc_params,
                      const std::shared_ptr<SCF> &scf_calc,
                      const std::filesystem::path &outdir) {
    auto indir = scf_calc->work_dir;
    auto rel = std::filesystem::relative(indir, outdir);
    auto prox = std::filesystem::proximate(indir, outdir);

    if (world.rank() == 0) {
      print("Running MolresponseLib::run_response() in directory: ", outdir);
      print("Ground state archive: ", indir);
      print("Relative path: ", rel);
      print("Proximate path: ", prox);
    }

    const auto &prefix = calc_params.prefix();
    std::string archive_name = prefix + ".restartdata";
    std::string fock_json_file = (prox / (prefix + ".fock.json")).string();
    std::string moldft_checkpt = prox / "moldft.calc_info.json";
    auto relative_archive = prox / archive_name;

    if (!std::filesystem::exists(moldft_checkpt)) {
      if (world.rank() == 0) {
        print("Error: Missing ground-state checkpoint file: ", moldft_checkpt);
      }
      throw std::runtime_error("Missing ground-state checkpoint file");
    }

    auto read_molecule = [](const std::string &scf_ckpt) -> Molecule {
      std::ifstream ifs(scf_ckpt);
      nlohmann::json j;
      ifs >> j;
      ifs.close();
      Molecule mol;
      if (j.contains("molecule")) {
        mol.from_json(j["molecule"]);
      } else {
        throw std::runtime_error(
            "Molecule information missing from checkpoint JSON.");
      }
      return mol;
    };

    Molecule molecule = read_molecule(moldft_checkpt);
    print("Read molecule with ", molecule.natom(), " atoms.");

    GroundStateData ground(world, relative_archive.string(), molecule);
    ResponseManager response_manager(world, calc_params);

    return GroundContext{std::move(molecule), std::move(ground),
                         std::move(response_manager), relative_archive.string(),
                         std::move(fock_json_file)};
  }

  static PlannedStates
  plan_required_states(World &world, const CalculationParameters &calc_params,
                       const GroundContext &ctx,
                       const ResponseParameters &response_params) {
    StateGenerator state_generator(ctx.molecule, calc_params.protocol(),
                                   ctx.ground.isSpinRestricted(),
                                   response_params);
    auto generated_states = state_generator.generateStates();

    if (world.rank() == 0) {
      GeneratedStateData::print_generated_state_map(generated_states.state_map);
    }
    DerivedStatePlan derived_state_plan =
        DerivedStatePlanner::build_vbc_driven_quadratic_plan(
            response_params, ctx.molecule, ctx.ground.isSpinRestricted(),
            calc_params.protocol());
    StateParallelPlan state_parallel_plan = StateParallelPlanner::build(
        response_params, world.size(), generated_states.states);
    if (world.rank() == 0 && !derived_state_plan.requests.empty()) {
      print("🧩 Planned ", derived_state_plan.requests.size(),
            " VBC-derived quadratic-state requests.");
    }
    if (world.rank() == 0 && response_params.state_parallel() != "off") {
      print(
          "State-parallel plan: mode=", state_parallel_plan.effective_mode,
          " requested_groups=", state_parallel_plan.requested_groups,
          " mapping_groups=", state_parallel_plan.mapping_groups,
          " state_owner_groups=", state_parallel_plan.state_owner_groups,
          " effective_point_groups=",
          state_parallel_plan.effective_point_groups,
          " point_parallel_start_protocol_index=",
          state_parallel_plan.point_parallel_start_protocol_index,
          " frequency_policy=", state_parallel_plan.frequency_partition_policy,
          " reason=", state_parallel_plan.reason);
    }
    world.gop.fence();
    return PlannedStates{std::move(generated_states),
                         std::move(derived_state_plan),
                         std::move(state_parallel_plan)};
  }

  struct RuntimePointOwnershipPolicy {
    size_t point_parallel_start_protocol_index = 0;
    bool restart_protocol0_saved_complete = false;
    bool restart_point_parallel_promoted = false;
  };

  static std::vector<size_t>
  build_owner_by_state_index(const std::vector<LinearResponseDescriptor> &states,
                             const StateParallelPlan &state_parallel_plan) {
    std::vector<size_t> owner_by_state_index(states.size(), 0);
    for (const auto &assignment : state_parallel_plan.assignments) {
      if (assignment.state_index < owner_by_state_index.size()) {
        owner_by_state_index[assignment.state_index] = assignment.owner_group;
      }
    }
    return owner_by_state_index;
  }

  static RuntimePointOwnershipPolicy compute_runtime_point_ownership_policy(
      World &world, const CalculationParameters &calc_params,
      const StateParallelPlan &state_parallel_plan,
      const std::vector<LinearResponseDescriptor> &linear_states,
      bool owner_group_schedule) {
    RuntimePointOwnershipPolicy runtime_policy;
    runtime_policy.point_parallel_start_protocol_index =
        state_parallel_plan.point_parallel_start_protocol_index;

    const bool has_protocol_thresholds = !calc_params.protocol().empty();
    const bool restart_promotion_candidate =
        owner_group_schedule && has_protocol_thresholds &&
        state_parallel_plan.point_parallel_start_protocol_index == 1;
    if (restart_promotion_candidate) {
      if (world.rank() == 0) {
        const nlohmann::json existing_metadata =
            read_json_file_or_object("response_metadata.json");
        runtime_policy.restart_protocol0_saved_complete = true;
        for (const auto &state : linear_states) {
          for (size_t freq_idx = 0; freq_idx < state.num_frequencies();
               ++freq_idx) {
            LinearResponsePoint pt{state, 0, freq_idx};
            if (!point_ready_in_metadata(existing_metadata, pt,
                                         /*require_saved=*/true,
                                         /*require_converged=*/false)) {
              runtime_policy.restart_protocol0_saved_complete = false;
              break;
            }
          }
          if (!runtime_policy.restart_protocol0_saved_complete) {
            break;
          }
        }
      }
      world.gop.broadcast_serializable(
          runtime_policy.restart_protocol0_saved_complete, 0);
    }

    runtime_policy.point_parallel_start_protocol_index =
        state_parallel_plan.effective_point_parallel_start_protocol_index(
            owner_group_schedule, has_protocol_thresholds,
            runtime_policy.restart_protocol0_saved_complete);
    runtime_policy.restart_point_parallel_promoted =
        runtime_policy.point_parallel_start_protocol_index !=
        state_parallel_plan.point_parallel_start_protocol_index;
    return runtime_policy;
  }

  static void build_local_state_workset(
      const std::vector<LinearResponseDescriptor> &linear_states,
      const std::vector<size_t> &owner_by_state_index,
      const PointOwnershipScheduler &point_scheduler, size_t subgroup_id,
      size_t mapping_groups, std::vector<size_t> &local_state_indices,
      std::vector<LinearResponseDescriptor> &local_states) {
    local_state_indices.clear();
    local_states.clear();
    local_state_indices.reserve(
        linear_states.size() / std::max<size_t>(1, mapping_groups));
    local_states.reserve(linear_states.size());

    const size_t point_owner_groups = point_scheduler.owner_groups();
    for (size_t state_index = 0; state_index < linear_states.size();
         ++state_index) {
      const bool owns_state_at_first_protocol =
          owner_by_state_index[state_index] == subgroup_id;
      bool owns_any_point_after_first_protocol = false;
      if (point_owner_groups > 1) {
        const auto &state = linear_states[state_index];
        for (size_t freq_index = 0; freq_index < state.num_frequencies();
             ++freq_index) {
          if (point_scheduler.owner_group(state_index, freq_index) ==
              subgroup_id) {
            owns_any_point_after_first_protocol = true;
            break;
          }
        }
      }
      if (owns_state_at_first_protocol) {
        local_state_indices.push_back(state_index);
      }
      if (owns_state_at_first_protocol || owns_any_point_after_first_protocol) {
        local_states.push_back(linear_states[state_index]);
      }
    }
  }

  struct StateSolveScheduleContext {
    // True when subgroup solve path is enabled and requested.
    bool subgroup_parallel_requested = false;
    // True when deterministic owner-lane scheduling is active.
    bool owner_group_schedule = false;
    // Planner output (group counts, mode, protocol switch threshold).
    const StateParallelPlan &state_parallel_plan;
    // Full set of linear states generated in Stage 1.
    const std::vector<LinearResponseDescriptor> &linear_states;
    // state_index -> owner lane for protocol ranges using state ownership.
    std::vector<size_t> owner_by_state_index;
    // Deterministic (state,freq) -> owner lane mapping for point ownership.
    PointOwnershipScheduler point_scheduler;
    // Runtime protocol index where state-ownership switches to point-ownership.
    size_t runtime_point_parallel_start_protocol_index = 0;
    // Restart diagnostics propagated to stage metadata.
    bool restart_point_parallel_promoted = false;
    bool restart_protocol0_saved_complete = false;

    StateSolveScheduleContext(
        bool subgroup_parallel_requested_, bool owner_group_schedule_,
        const StateParallelPlan &state_parallel_plan_,
        const std::vector<LinearResponseDescriptor> &linear_states_,
        std::vector<size_t> owner_by_state_index_, size_t effective_point_groups,
        size_t runtime_point_parallel_start_protocol_index_,
        bool restart_point_parallel_promoted_,
        bool restart_protocol0_saved_complete_)
        : subgroup_parallel_requested(subgroup_parallel_requested_),
          owner_group_schedule(owner_group_schedule_),
          state_parallel_plan(state_parallel_plan_),
          linear_states(linear_states_),
          owner_by_state_index(std::move(owner_by_state_index_)),
          point_scheduler(linear_states_, effective_point_groups),
          runtime_point_parallel_start_protocol_index(
              runtime_point_parallel_start_protocol_index_),
          restart_point_parallel_promoted(restart_point_parallel_promoted_),
          restart_protocol0_saved_complete(
              restart_protocol0_saved_complete_) {}

    [[nodiscard]] size_t point_owner_groups() const {
      return point_scheduler.owner_groups();
    }
  };

  // Runtime policy gate shared by serial and subgroup paths.
  // For protocol indices below the runtime switch threshold we keep all
  // frequencies of a state together; after the switch threshold we fan out by
  // independent state-frequency points.
  static bool
  use_state_ownership_for_protocol_runtime(const StateSolveScheduleContext &ctx,
                                           size_t protocol_index) {
    if (ctx.state_parallel_plan.mapping_groups <= 1) {
      return true;
    }
    return protocol_index < ctx.runtime_point_parallel_start_protocol_index;
  }

  static void
  print_state_solve_execution_mode(World &world,
                                   const StateParallelPlan &state_parallel_plan,
                                   bool subgroup_parallel_requested,
                                   bool owner_group_schedule) {
    if (world.rank() == 0 && subgroup_parallel_requested) {
      print("State ownership mapping is active across ",
            state_parallel_plan.mapping_groups,
            " groups (protocol-0 owner groups=",
            state_parallel_plan.state_owner_groups,
            "); executing owner-group solves in parallel subworlds.");
    } else if (world.rank() == 0 && owner_group_schedule) {
      print("State ownership mapping is active across ",
            state_parallel_plan.mapping_groups,
            " groups (protocol-0 owner groups=",
            state_parallel_plan.state_owner_groups,
            "); executing deterministic owner-group solve passes "
            "serially on the universe communicator.");
    } else if (world.rank() == 0 && state_parallel_plan.mapping_groups > 1) {
      print("State ownership mapping is active across ",
            state_parallel_plan.mapping_groups,
            " groups; falling back to plain serial state loop.");
    }
  }

  // Build runtime scheduling context for stage-2 linear solves:
  // - picks subgroup vs serial-lane mode flags,
  // - constructs owner maps/schedulers,
  // - applies restart-aware protocol switch policy.
  static StateSolveScheduleContext
  build_state_solve_schedule_context(World &world,
                                     const CalculationParameters &calc_params,
                                     const PlannedStates &planned_states) {
    const auto &state_parallel_plan = planned_states.state_parallel_plan;
    const bool subgroup_parallel_requested =
        state_parallel_plan.execution_enabled &&
        state_parallel_plan.subgroup_parallel_enabled &&
        state_parallel_plan.execution_groups > 1;
    const bool owner_group_schedule =
        state_parallel_plan.mapping_groups > 1 &&
        state_parallel_plan.effective_mode != "serial";
    print_state_solve_execution_mode(world, state_parallel_plan,
                                     subgroup_parallel_requested,
                                     owner_group_schedule);

    const auto &linear_states = planned_states.generated_states.states;
    auto owner_by_state_index =
        build_owner_by_state_index(linear_states, state_parallel_plan);
    const RuntimePointOwnershipPolicy runtime_policy =
        compute_runtime_point_ownership_policy(
            world, calc_params, state_parallel_plan, linear_states,
            owner_group_schedule);
    if (runtime_policy.restart_point_parallel_promoted && world.rank() == 0) {
      print("State-parallel restart detected: protocol-0 points are saved; "
            "enabling point ownership from protocol index 0.");
    }

    return StateSolveScheduleContext(
        subgroup_parallel_requested, owner_group_schedule, state_parallel_plan,
        linear_states, std::move(owner_by_state_index),
        state_parallel_plan.effective_point_groups,
        runtime_policy.point_parallel_start_protocol_index,
        runtime_policy.restart_point_parallel_promoted,
        runtime_policy.restart_protocol0_saved_complete);
  }

  static void execute_serial_state_solve(
      World &world, const CalculationParameters &calc_params, GroundContext &ctx,
      const StateSolveScheduleContext &schedule_ctx,
      nlohmann::json &state_metadata_json, nlohmann::json &debug_log_json) {
    // Serial execution path:
    // 1) initialize metadata/log persistence,
    // 2) iterate protocol thresholds and solve pending work,
    // 3) emit in-memory metadata/debug JSON.
    JsonStateSolvePersistence persistence(world, "response_metadata.json",
                                          "response_log.json");
    persistence.initialize_states(schedule_ctx.linear_states);

    if (world.rank() == 0) {
      persistence.print_summary();
    }
    world.gop.fence();

    auto needs_solving_at_protocol = [&](double protocol_thresh,
                                         size_t thresh_index) {
      // A protocol threshold is considered "active" when at least one point is
      // missing on disk, or not converged at the final threshold.
      const bool at_final_protocol = protocol_thresh == calc_params.protocol().back();

      for (const auto &state : schedule_ctx.linear_states) {
        for (size_t freq_idx = 0; freq_idx < state.num_frequencies();
             ++freq_idx) {
          LinearResponsePoint pt{state, thresh_index, freq_idx};
          const bool is_saved = persistence.is_saved(pt);
          const bool should_solve =
              !is_saved || (at_final_protocol && !persistence.is_converged(pt));

          if (world.rank() == 0) {
            print("Checking state ", pt.perturbationDescription(),
                  " at thresh ", protocol_thresh, " freq ", pt.frequency(),
                  " is_saved=", is_saved, " at_final_protocol=",
                  at_final_protocol, " should_solve=", should_solve);
          }

          if (should_solve) {
            return true;
          }
        }
      }

      return false;
    };

    const auto &protocol = calc_params.protocol();
    for (size_t ti = 0; ti < protocol.size(); ++ti) {
      const double thresh = protocol[ti];

      if (!needs_solving_at_protocol(thresh, ti)) {
        if (world.rank() == 0) {
          madness::print("✓ All states converged at thresh", thresh,
                         "skipping to next protocol.");
        }
        continue;
      }

      ctx.response_manager.setProtocol(world, ctx.ground.getL(), thresh);
      ctx.ground.prepareOrbitals(world, FunctionDefaults<3>::get_k(), thresh);
      ctx.ground.computePreliminaries(world, *ctx.response_manager.getCoulombOp(),
                                      ctx.response_manager.getVtol(),
                                      ctx.fock_json_file);

      const bool at_final_protocol = (ti + 1 == protocol.size());
      auto solve_state = [&](size_t state_index) {
        auto &state = schedule_ctx.linear_states[state_index];
        computeFrequencyLoop(world, ctx.response_manager, state, ti, ctx.ground,
                             persistence, at_final_protocol);
        persistence.flush_debug_log(world);
      };
      auto solve_state_frequency = [&](size_t state_index, size_t freq_index) {
        const auto &state = schedule_ctx.linear_states[state_index];
        LinearResponseDescriptor single_frequency_state(
            state.perturbation, {state.frequency(freq_index)}, state.thresholds,
            state.spin_restricted);
        computeFrequencyLoop(world, ctx.response_manager, single_frequency_state,
                             ti, ctx.ground, persistence, at_final_protocol);
        persistence.flush_debug_log(world);
      };

      if (schedule_ctx.owner_group_schedule &&
          !use_state_ownership_for_protocol_runtime(schedule_ctx, ti)) {
        // Point mode: route independent (state,freq) points to deterministic
        // owner lanes.
        const size_t point_owner_groups = schedule_ctx.point_owner_groups();
        for (size_t gid = 0; gid < point_owner_groups; ++gid) {
          if (world.rank() == 0) {
            print("State-frequency solve lane ", gid, "/",
                  point_owner_groups - 1, " at protocol thresh ", thresh);
          }
          for (size_t state_index = 0;
               state_index < schedule_ctx.linear_states.size(); ++state_index) {
            const auto &state = schedule_ctx.linear_states[state_index];
            for (size_t freq_index = 0; freq_index < state.num_frequencies();
                 ++freq_index) {
              if (schedule_ctx.point_scheduler.owner_group(state_index,
                                                           freq_index) == gid) {
                solve_state_frequency(state_index, freq_index);
              }
            }
          }
        }
      } else if (schedule_ctx.owner_group_schedule) {
        // State mode: keep all frequencies for an owned state on one lane.
        for (size_t gid = 0;
             gid < schedule_ctx.state_parallel_plan.state_owner_groups;
             ++gid) {
          if (world.rank() == 0) {
            print("State solve lane ", gid, "/",
                  schedule_ctx.state_parallel_plan.state_owner_groups - 1,
                  " at protocol thresh ", thresh);
          }
          for (size_t state_index = 0;
               state_index < schedule_ctx.linear_states.size(); ++state_index) {
            if (schedule_ctx.owner_by_state_index[state_index] == gid) {
              solve_state(state_index);
            }
          }
        }
      } else {
        for (size_t state_index = 0; state_index < schedule_ctx.linear_states.size();
             ++state_index) {
          solve_state(state_index);
        }
      }
    }

    state_metadata_json = persistence.metadata_json();
    debug_log_json = persistence.debug_log_json();
  }

  static bool execute_subgroup_state_solve(
      World &world, const CalculationParameters &calc_params, GroundContext &ctx,
      const StateSolveScheduleContext &schedule_ctx,
      nlohmann::json &state_metadata_json, nlohmann::json &debug_log_json) {
    // Subgroup execution path:
    // 1) create macrotask subworlds,
    // 2) solve local ownership shards in each subgroup,
    // 3) merge subgroup metadata/debug shards on world rank 0 and broadcast.
    const auto &state_parallel_plan = schedule_ctx.state_parallel_plan;
    try {
      auto subworld_ptr = MacroTaskQ::create_worlds(
          world, state_parallel_plan.execution_groups);
      if (!subworld_ptr) {
        throw std::runtime_error("subworld creation returned null");
      }
      World &subworld = *subworld_ptr;
      const size_t subgroup_id = static_cast<size_t>(
          world.rank() % static_cast<int>(state_parallel_plan.execution_groups));
      ScopedRankLogRedirect subgroup_console_redirect(
          subworld.rank() == 0, group_console_file(subgroup_id),
          "stage2-linear subgroup=" + std::to_string(subgroup_id));

      auto old_pmap3 = FunctionDefaults<3>::get_pmap();
      auto restore_pmap = [&]() { FunctionDefaults<3>::set_pmap(old_pmap3); };

      // Molresponse state solves are strictly 3D; only swap the 3D default
      // pmap.
      FunctionDefaults<3>::set_default_pmap(subworld);
      try {
        std::vector<size_t> local_state_indices;
        std::vector<LinearResponseDescriptor> local_states;
        build_local_state_workset(
            schedule_ctx.linear_states, schedule_ctx.owner_by_state_index,
            schedule_ctx.point_scheduler, subgroup_id,
            state_parallel_plan.mapping_groups, local_state_indices, local_states);

        const std::string metadata_shard_file =
            group_shard_file("response_metadata.json", subgroup_id);
        const std::string debug_shard_file =
            group_shard_file("response_log.json", subgroup_id);
        const std::string fock_shard_file =
            group_shard_file(ctx.fock_json_file, subgroup_id);

        JsonStateSolvePersistence local_persistence(
            subworld, metadata_shard_file, debug_shard_file);
        local_persistence.initialize_states(local_states);
        if (subworld.rank() == 0) {
          print("State-parallel subgroup ", subgroup_id, " owns ",
                local_state_indices.size(), " protocol-0 states and ",
                local_states.size(), " active states across all protocols.");
          local_persistence.print_summary();
        }
        subworld.gop.fence();

        if (!local_states.empty()) {
          GroundStateData local_ground(subworld, ctx.archive_file, ctx.molecule);
          ResponseManager local_response_manager(subworld, calc_params);

          auto local_needs_solving_at_protocol = [&](double protocol_thresh,
                                                     size_t thresh_index) {
            // Local "needs solve" probe mirrors the serial path, but filtered
            // by this subgroup's ownership mode (state lane or point lane).
            const bool at_final_protocol =
                protocol_thresh == calc_params.protocol().back();

            auto point_needs_solving =
                [&](const LinearResponseDescriptor &state, size_t freq_idx) {
                  LinearResponsePoint pt{state, thresh_index, freq_idx};
                  const bool is_saved = local_persistence.is_saved(pt);
                  const bool should_solve =
                      !is_saved || (at_final_protocol &&
                                    !local_persistence.is_converged(pt));
                  return should_solve;
                };

            const bool use_state_ownership_for_thresh =
                use_state_ownership_for_protocol_runtime(schedule_ctx,
                                                        thresh_index);
            if (use_state_ownership_for_thresh) {
              for (const auto state_index : local_state_indices) {
                const auto &state = schedule_ctx.linear_states[state_index];
                for (size_t freq_idx = 0; freq_idx < state.num_frequencies();
                     ++freq_idx) {
                  if (point_needs_solving(state, freq_idx)) {
                    return true;
                  }
                }
              }
            } else {
              for (size_t state_index = 0;
                   state_index < schedule_ctx.linear_states.size();
                   ++state_index) {
                const auto &state = schedule_ctx.linear_states[state_index];
                for (size_t freq_idx = 0; freq_idx < state.num_frequencies();
                     ++freq_idx) {
                  if (schedule_ctx.point_scheduler.owner_group(state_index,
                                                               freq_idx) !=
                      subgroup_id) {
                    continue;
                  }
                  if (point_needs_solving(state, freq_idx)) {
                    return true;
                  }
                }
              }
            }
            return false;
          };

          const auto &protocol = calc_params.protocol();
          for (size_t ti = 0; ti < protocol.size(); ++ti) {
            const double thresh = protocol[ti];
            if (!local_needs_solving_at_protocol(thresh, ti)) {
              if (subworld.rank() == 0) {
                print("Subgroup ", subgroup_id, " has no pending states at thresh ",
                      thresh, "; skipping.");
              }
              // Keep protocol progression globally lock-step across all
              // subgroups because ti>0 points may depend on ti-1 files
              // produced by other owners.
              world.gop.fence();
              continue;
            }

            local_response_manager.setProtocol(subworld, local_ground.getL(),
                                               thresh);
            local_ground.prepareOrbitals(subworld, FunctionDefaults<3>::get_k(),
                                         thresh);
            local_ground.computePreliminaries(
                subworld, *local_response_manager.getCoulombOp(),
                local_response_manager.getVtol(), fock_shard_file);

            const bool at_final_protocol = (ti + 1 == protocol.size());
            const bool use_state_ownership_for_ti =
                use_state_ownership_for_protocol_runtime(schedule_ctx, ti);
            if (use_state_ownership_for_ti) {
              for (const auto state_index : local_state_indices) {
                auto &state = schedule_ctx.linear_states[state_index];
                computeFrequencyLoop(subworld, local_response_manager, state, ti,
                                     local_ground, local_persistence,
                                     at_final_protocol);
                local_persistence.flush_debug_log(subworld);
              }
            } else {
              if (subworld.rank() == 0) {
                print("Subgroup ", subgroup_id,
                      " solving independent state-frequency points at thresh ",
                      thresh, ".");
              }
              for (size_t state_index = 0;
                   state_index < schedule_ctx.linear_states.size();
                   ++state_index) {
                const auto &state = schedule_ctx.linear_states[state_index];
                for (size_t freq_idx = 0; freq_idx < state.num_frequencies();
                     ++freq_idx) {
                  if (schedule_ctx.point_scheduler.owner_group(state_index,
                                                               freq_idx) !=
                      subgroup_id) {
                    continue;
                  }
                  LinearResponseDescriptor single_frequency_state(
                      state.perturbation, {state.frequency(freq_idx)},
                      state.thresholds, state.spin_restricted);
                  computeFrequencyLoop(subworld, local_response_manager,
                                       single_frequency_state, ti, local_ground,
                                       local_persistence, at_final_protocol);
                  local_persistence.flush_debug_log(subworld);
                }
              }
            }
            // ti+1 must not begin until every subgroup has finished ti.
            world.gop.fence();
          }
        } else if (subworld.rank() == 0) {
          print("Subgroup ", subgroup_id,
                " has no owned states; skipping subgroup solve loop.");
        }

        local_persistence.flush_debug_log(subworld);
        subworld.gop.fence();
        restore_pmap();
      } catch (...) {
        restore_pmap();
        throw;
      }

      world.gop.fence();
      if (world.rank() == 0) {
        nlohmann::json merged_metadata = nlohmann::json::object();
        nlohmann::json merged_debug_log = nlohmann::json::object();
        for (size_t gid = 0; gid < state_parallel_plan.execution_groups; ++gid) {
          const std::string metadata_file =
              group_shard_file("response_metadata.json", gid);
          const std::string debug_file =
              group_shard_file("response_log.json", gid);
          merge_state_metadata_json(merged_metadata,
                                    read_json_file_or_object(metadata_file));
          merge_debug_log_json(merged_debug_log,
                               read_json_file_or_object(debug_file));
        }

        write_json_file("response_metadata.json", merged_metadata);
        write_json_file("response_log.json", merged_debug_log);
        state_metadata_json = std::move(merged_metadata);
        debug_log_json = std::move(merged_debug_log);
      }

      std::string metadata_dump;
      std::string debug_dump;
      if (world.rank() == 0) {
        metadata_dump = state_metadata_json.dump();
        debug_dump = debug_log_json.dump();
      }
      world.gop.broadcast_serializable(metadata_dump, 0);
      world.gop.broadcast_serializable(debug_dump, 0);
      if (world.rank() != 0) {
        state_metadata_json = metadata_dump.empty()
                                  ? nlohmann::json::object()
                                  : nlohmann::json::parse(metadata_dump);
        debug_log_json = debug_dump.empty()
                             ? nlohmann::json::object()
                             : nlohmann::json::parse(debug_dump);
      }
      world.gop.fence();

      return true;
    } catch (const std::exception &ex) {
      if (world.rank() == 0) {
        print("State-parallel subgroup execution failed: ", ex.what(),
              " Falling back to serial state loop.");
      }
      return false;
    }
  }

  struct DerivedRequestTiming {
    bool success = false;
    double wall_seconds = 0.0;
    double cpu_seconds = 0.0;
  };

  static DerivedRequestTiming
  run_derived_request(World &exec_world, const GroundStateData &ground,
                      const DerivedStateRequest &req,
                      SimpleVBCComputer &vbc_computer, long &completed,
                      long &failed) {
    // Single derived-request runner used by both subgroup and serial fallback
    // execution paths.
    const double start_wall = madness::wall_time();
    const double start_cpu = madness::cpu_time();
    DerivedRequestTiming timing;
    try {
      VBCResponseState vbc_state =
          DerivedStatePlanner::make_vbc_state(req, ground.isSpinRestricted());
      vbc_computer.compute_and_save(vbc_state);
      timing.success = true;
      if (exec_world.rank() == 0) {
        ++completed;
      }
    } catch (const std::exception &ex) {
      if (exec_world.rank() == 0) {
        ++failed;
        print("Derived-state request ", req.derived_state_id,
              " failed during VBC solve: ", ex.what());
      }
    }
    timing.wall_seconds = madness::wall_time() - start_wall;
    timing.cpu_seconds = madness::cpu_time() - start_cpu;
    return timing;
  }

  struct DerivedExecutionResult {
    // Dependency readiness report for all derived requests.
    DerivedStateGateReport dependency_gate;
    // Execution summary JSON written into response metadata.
    nlohmann::json execution;
  };

  struct FinalProtocolState {
    double threshold = 0.0;
    size_t threshold_index = 0;
  };

  static DerivedExecutionResult execute_derived_state_requests(
      World &world, const CalculationParameters &calc_params, GroundContext &ctx,
      const PlannedStates &planned_states, const nlohmann::json &state_metadata,
      bool subgroup_parallel_requested, bool owner_group_schedule,
      size_t final_ti, double final_thresh) {
    // Derived stage:
    // 1) evaluate dependency gate against final linear-state readiness,
    // 2) execute ready requests in subgroup mode when available,
    // 3) fall back to deterministic serial lanes when needed,
    // 4) return both gate diagnostics and execution summary.
    const auto &state_parallel_plan = planned_states.state_parallel_plan;

    DerivedStateGateReport derived_gate =
        DerivedStatePlanner::evaluate_dependency_gate(
            planned_states.derived_state_plan, planned_states.generated_states,
            final_ti, [&](const LinearResponsePoint &pt) {
              return point_ready_in_metadata(state_metadata, pt,
                                             /*require_saved=*/true,
                                             /*require_converged=*/true);
            });
    if (world.rank() == 0 && derived_gate.total_requests > 0) {
      print("Derived-state dependency gate: ready ",
            derived_gate.ready_requests, "/", derived_gate.total_requests,
            " blocked=", derived_gate.blocked_requests);
    }

    // Execution summary persisted under metadata["derived_state_planner"].
    nlohmann::json derived_execution = {
        {"attempted", false},
        {"mode", "none_ready"},
        {"execution_groups", 1},
        {"ready_requests", derived_gate.ready_requests},
        {"blocked_requests", derived_gate.blocked_requests},
        {"completed_requests", 0},
        {"failed_requests", 0},
        {"total_wall_seconds", 0.0},
        {"total_cpu_seconds", 0.0},
        {"request_timings", nlohmann::json::array()}};

    // Plan indices that passed the dependency gate at the final threshold.
    std::vector<size_t> ready_request_indices;
    ready_request_indices.reserve(planned_states.derived_state_plan.requests.size());
    for (size_t i = 0; i < planned_states.derived_state_plan.requests.size();
         ++i) {
      if (i < derived_gate.entries.size() && derived_gate.entries[i].ready) {
        ready_request_indices.push_back(i);
      }
    }

    if (!ready_request_indices.empty()) {
      const size_t derived_owner_groups =
          subgroup_parallel_requested
              ? state_parallel_plan.execution_groups
              : std::max<size_t>(1, state_parallel_plan.mapping_groups);
      // Deterministic lane assignment for ready derived requests.
      std::vector<size_t> owner_by_ready_index(ready_request_indices.size(), 0);
      for (size_t i = 0; i < ready_request_indices.size(); ++i) {
        owner_by_ready_index[i] = i % derived_owner_groups;
      }

      derived_execution["attempted"] = true;
      derived_execution["execution_groups"] = derived_owner_groups;

      bool ran_subgroup_derived = false;
      if (subgroup_parallel_requested && derived_owner_groups > 1) {
        derived_execution["mode"] = "owner_group_subworld";
        try {
          auto subworld_ptr = MacroTaskQ::create_worlds(
              world, state_parallel_plan.execution_groups);
          if (!subworld_ptr) {
            throw std::runtime_error("derived subworld creation returned null");
          }

          World &subworld = *subworld_ptr;
          const size_t subgroup_id = static_cast<size_t>(
              world.rank() %
              static_cast<int>(state_parallel_plan.execution_groups));
          ScopedRankLogRedirect subgroup_console_redirect(
              subworld.rank() == 0, group_console_file(subgroup_id),
              "stage2-derived subgroup=" + std::to_string(subgroup_id));
          const std::string timing_shard_file =
              group_derived_timing_file(subgroup_id);

          auto old_pmap3 = FunctionDefaults<3>::get_pmap();
          auto restore_pmap = [&]() { FunctionDefaults<3>::set_pmap(old_pmap3); };

          long local_completed = 0;
          long local_failed = 0;
          double local_wall_seconds = 0.0;
          double local_cpu_seconds = 0.0;
          nlohmann::json local_request_timings = nlohmann::json::array();
          FunctionDefaults<3>::set_default_pmap(subworld);
          try {
            std::vector<size_t> local_ready_positions;
            local_ready_positions.reserve(
                ready_request_indices.size() /
                std::max<size_t>(1, derived_owner_groups));
            for (size_t pos = 0; pos < ready_request_indices.size(); ++pos) {
              if (owner_by_ready_index[pos] == subgroup_id) {
                local_ready_positions.push_back(pos);
              }
            }

            if (!local_ready_positions.empty()) {
              GroundStateData local_ground(subworld, ctx.archive_file,
                                           ctx.molecule);
              ResponseManager local_response_manager(subworld, calc_params);
              local_response_manager.setProtocol(subworld, local_ground.getL(),
                                                 final_thresh);
              local_ground.prepareOrbitals(
                  subworld, FunctionDefaults<3>::get_k(), final_thresh);

              SimpleVBCComputer local_vbc_computer(subworld, local_ground);
              for (const auto pos : local_ready_positions) {
                const auto &req =
                    planned_states.derived_state_plan
                        .requests[ready_request_indices[pos]];
                const auto timing = run_derived_request(
                    subworld, local_ground, req, local_vbc_computer,
                    local_completed, local_failed);
                if (subworld.rank() == 0) {
                  local_wall_seconds += timing.wall_seconds;
                  local_cpu_seconds += timing.cpu_seconds;
                  local_request_timings.push_back(
                      {{"derived_state_id", req.derived_state_id},
                       {"owner_group", subgroup_id},
                       {"success", timing.success},
                       {"wall_seconds", timing.wall_seconds},
                       {"cpu_seconds", timing.cpu_seconds}});
                }
              }
            } else if (subworld.rank() == 0) {
              print("Derived-state subgroup ", subgroup_id,
                    " has no owned ready requests.");
            }

            if (subworld.rank() == 0) {
              write_json_file(timing_shard_file, local_request_timings);
            }
            subworld.gop.fence();
            restore_pmap();
          } catch (...) {
            restore_pmap();
            throw;
          }

          if (subworld.rank() != 0) {
            local_completed = 0;
            local_failed = 0;
            local_wall_seconds = 0.0;
            local_cpu_seconds = 0.0;
          }
          world.gop.sum(local_completed);
          world.gop.sum(local_failed);
          world.gop.sum(local_wall_seconds);
          world.gop.sum(local_cpu_seconds);
          derived_execution["completed_requests"] = local_completed;
          derived_execution["failed_requests"] = local_failed;
          derived_execution["total_wall_seconds"] = local_wall_seconds;
          derived_execution["total_cpu_seconds"] = local_cpu_seconds;
          if (world.rank() == 0) {
            print("Derived-state timing summary: completed=", local_completed,
                  " failed=", local_failed, " wall=", local_wall_seconds,
                  "s cpu=", local_cpu_seconds, "s");
          }

          world.gop.fence();
          nlohmann::json merged_request_timings = nlohmann::json::array();
          if (world.rank() == 0) {
            for (size_t gid = 0; gid < derived_owner_groups; ++gid) {
              const auto shard =
                  read_json_file_or_object(group_derived_timing_file(gid));
              if (!shard.is_array()) {
                continue;
              }
              for (const auto &entry : shard) {
                merged_request_timings.push_back(entry);
              }
            }
          }
          std::string request_timing_dump;
          if (world.rank() == 0) {
            request_timing_dump = merged_request_timings.dump();
          }
          world.gop.broadcast_serializable(request_timing_dump, 0);
          if (world.rank() != 0) {
            merged_request_timings = request_timing_dump.empty()
                                         ? nlohmann::json::array()
                                         : nlohmann::json::parse(
                                               request_timing_dump);
          }
          derived_execution["request_timings"] =
              std::move(merged_request_timings);
          world.gop.fence();
          ran_subgroup_derived = true;
        } catch (const std::exception &ex) {
          if (world.rank() == 0) {
            print("Derived-state subgroup execution failed: ", ex.what(),
                  " Falling back to serial derived-state loop.");
          }
        }
      }

      if (!ran_subgroup_derived) {
        // Deterministic serial fallback keeps a predictable owner-lane order
        // for easier debugging and reproducible logs.
        derived_execution["mode"] =
            (owner_group_schedule && derived_owner_groups > 1)
                ? "owner_group_serial_lanes"
                : "serial";
        if (world.rank() == 0 &&
            derived_execution["mode"].get<std::string>() ==
                "owner_group_serial_lanes") {
          print("Derived-state ownership mapping is active across ",
                derived_owner_groups,
                " groups; executing deterministic derived-state passes "
                "serially on the universe communicator.");
        }

        long local_completed = 0;
        long local_failed = 0;
        double local_wall_seconds = 0.0;
        double local_cpu_seconds = 0.0;
        nlohmann::json request_timings = nlohmann::json::array();
        SimpleVBCComputer serial_vbc_computer(world, ctx.ground);
        if (owner_group_schedule && derived_owner_groups > 1) {
          for (size_t gid = 0; gid < derived_owner_groups; ++gid) {
            if (world.rank() == 0) {
              print("Derived-state solve lane ", gid, "/",
                    derived_owner_groups - 1);
            }
            for (size_t pos = 0; pos < ready_request_indices.size(); ++pos) {
              if (owner_by_ready_index[pos] != gid) {
                continue;
              }
              const auto &req = planned_states.derived_state_plan
                                    .requests[ready_request_indices[pos]];
              const auto timing = run_derived_request(
                  world, ctx.ground, req, serial_vbc_computer, local_completed,
                  local_failed);
              if (world.rank() == 0) {
                local_wall_seconds += timing.wall_seconds;
                local_cpu_seconds += timing.cpu_seconds;
                request_timings.push_back({{"derived_state_id", req.derived_state_id},
                                           {"owner_group", gid},
                                           {"success", timing.success},
                                           {"wall_seconds", timing.wall_seconds},
                                           {"cpu_seconds", timing.cpu_seconds}});
              }
            }
          }
        } else {
          for (const auto idx : ready_request_indices) {
            const auto &req = planned_states.derived_state_plan.requests[idx];
            const auto timing =
                run_derived_request(world, ctx.ground, req, serial_vbc_computer,
                                    local_completed, local_failed);
            if (world.rank() == 0) {
              local_wall_seconds += timing.wall_seconds;
              local_cpu_seconds += timing.cpu_seconds;
              request_timings.push_back({{"derived_state_id", req.derived_state_id},
                                         {"owner_group", 0},
                                         {"success", timing.success},
                                         {"wall_seconds", timing.wall_seconds},
                                         {"cpu_seconds", timing.cpu_seconds}});
            }
          }
        }

        if (world.rank() != 0) {
          local_completed = 0;
          local_failed = 0;
          local_wall_seconds = 0.0;
          local_cpu_seconds = 0.0;
        }
        world.gop.sum(local_completed);
        world.gop.sum(local_failed);
        world.gop.sum(local_wall_seconds);
        world.gop.sum(local_cpu_seconds);
        derived_execution["completed_requests"] = local_completed;
        derived_execution["failed_requests"] = local_failed;
        derived_execution["total_wall_seconds"] = local_wall_seconds;
        derived_execution["total_cpu_seconds"] = local_cpu_seconds;
        if (world.rank() == 0) {
          print("Derived-state timing summary: completed=", local_completed,
                " failed=", local_failed, " wall=", local_wall_seconds,
                "s cpu=", local_cpu_seconds, "s");
        }
        derived_execution["request_timings"] = std::move(request_timings);
      }
    }

    std::string derived_execution_dump;
    if (world.rank() == 0) {
      derived_execution_dump = derived_execution.dump();
    }
    world.gop.broadcast_serializable(derived_execution_dump, 0);
    if (world.rank() != 0) {
      derived_execution =
          derived_execution_dump.empty()
              ? nlohmann::json::object()
              : nlohmann::json::parse(derived_execution_dump);
    }

    return DerivedExecutionResult{std::move(derived_gate),
                                  std::move(derived_execution)};
  }

  // Stage 2 post-linear finalize:
  // 1) configure ground/response context at the final protocol threshold,
  // 2) verify every linear point is converged in merged state metadata.
  static FinalProtocolState prepare_and_validate_final_protocol_state(
      World &world, const CalculationParameters &calc_params, GroundContext &ctx,
      const PlannedStates &planned_states,
      const nlohmann::json &state_metadata_json) {
    const auto &protocol = calc_params.protocol();
    FinalProtocolState final_state{protocol.back(), protocol.size() - 1};

    ctx.response_manager.setProtocol(world, ctx.ground.getL(),
                                     final_state.threshold);
    ctx.ground.prepareOrbitals(world, FunctionDefaults<3>::get_k(),
                               final_state.threshold);
    ctx.ground.computePreliminaries(world, *ctx.response_manager.getCoulombOp(),
                                    ctx.response_manager.getVtol(),
                                    ctx.fock_json_file);

    bool all_are_converged = true;
    for (const auto &state : planned_states.generated_states.states) {
      for (size_t fi = 0; fi < state.num_frequencies(); ++fi) {
        LinearResponsePoint pt{state, final_state.threshold_index, fi};
        if (!point_ready_in_metadata(state_metadata_json, pt,
                                     /*require_saved=*/false,
                                     /*require_converged=*/true)) {
          all_are_converged = false;
          break;
        }
      }
      if (!all_are_converged) {
        break;
      }
    }
    MADNESS_ASSERT(all_are_converged);
    return final_state;
  }

  // Assemble stage-2 metadata consumed by downstream property execution and
  // restart diagnostics.
  static nlohmann::json build_state_stage_metadata(
      const PlannedStates &planned_states, const nlohmann::json &state_metadata,
      const StateSolveScheduleContext &schedule_ctx,
      const DerivedExecutionResult &derived_result) {
    auto metadata = state_metadata;
    metadata["state_parallel_planner"] =
        planned_states.state_parallel_plan.to_json();
    metadata["state_parallel_runtime"] = {
        {"effective_point_groups", schedule_ctx.point_owner_groups()},
        {"effective_point_parallel_start_protocol_index",
         schedule_ctx.runtime_point_parallel_start_protocol_index},
        {"restart_protocol0_saved_complete",
         schedule_ctx.restart_protocol0_saved_complete},
        {"restart_point_parallel_promoted",
         schedule_ctx.restart_point_parallel_promoted}};
    if (schedule_ctx.state_parallel_plan.execution_groups > 1) {
      nlohmann::json group_console_logs = nlohmann::json::array();
      for (size_t gid = 0;
           gid < schedule_ctx.state_parallel_plan.execution_groups; ++gid) {
        group_console_logs.push_back(group_console_file(gid));
      }
      metadata["state_parallel_runtime"]["group_console_logs"] =
          std::move(group_console_logs);
    }
    metadata["derived_state_planner"] = {
        {"note",
         "Stage 2c: ready VBC-derived requests are executed before property "
         "assembly; blocked requests remain gated."},
        {"plan", planned_states.derived_state_plan.to_json()},
        {"dependency_gate", derived_result.dependency_gate.to_json()},
        {"execution", derived_result.execution}};
    return metadata;
  }

  // Stage-2 orchestrator:
  // 2a) build runtime ownership policy, 2b) solve linear states,
  // 2c) execute dependency-gated derived requests, 2d) publish stage metadata.
  static SolvedStates
  solve_all_states(World &world, const CalculationParameters &calc_params,
                   GroundContext &ctx,
                   const ResponseParameters &response_params,
                   PlannedStates planned_states) {
    (void)response_params;

    // Stage 2a setup: build runtime ownership schedule and restart-aware policy.
    const StateSolveScheduleContext schedule_ctx =
        build_state_solve_schedule_context(world, calc_params, planned_states);

    nlohmann::json state_metadata_json = nlohmann::json::object();
    nlohmann::json debug_log_json = nlohmann::json::object();

    // Stage 2b linear solve: subgroup path first, serial fallback second.
    bool ran_subgroup_path = false;
    if (schedule_ctx.subgroup_parallel_requested) {
      ran_subgroup_path = execute_subgroup_state_solve(
          world, calc_params, ctx, schedule_ctx, state_metadata_json,
          debug_log_json);
    }
    if (!ran_subgroup_path) {
      execute_serial_state_solve(world, calc_params, ctx, schedule_ctx,
                                 state_metadata_json, debug_log_json);
    }

    const FinalProtocolState final_state = prepare_and_validate_final_protocol_state(
        world, calc_params, ctx, planned_states, state_metadata_json);

    // Stage 2c derived solve: dependency-gated execution + summary metadata.
    const DerivedExecutionResult derived_result = execute_derived_state_requests(
        world, calc_params, ctx, planned_states, state_metadata_json,
        schedule_ctx.subgroup_parallel_requested, schedule_ctx.owner_group_schedule,
        final_state.threshold_index, final_state.threshold);

    // Stage 2 metadata assembly for downstream property stage.
    auto metadata = build_state_stage_metadata(
        planned_states, state_metadata_json, schedule_ctx, derived_result);

    return SolvedStates{std::move(planned_states), std::move(metadata),
                        std::move(debug_log_json)};
  }

  struct PropertyContext {
    // Universe world used for property assembly.
    World &world;
    // Parsed response inputs (frequencies, directions, requested properties).
    const ResponseParameters &response_params;
    // Ground-state objects and managers prepared in Stage 1.
    const GroundContext &ground_ctx;
    // Solved-state metadata/map produced by Stage 2.
    const SolvedStates &solved_states;
    // Property accumulator/writer.
    PropertyManager &properties;
    // Optional SCF handle needed by Raman/hessian routines.
    std::shared_ptr<SCF> scf_calc;
  };

  enum class PropertyType { Alpha, Beta, Raman };

  inline static PropertyType parse_property_name(const std::string &raw) {
    auto key = raw;
    if (key.size() >= 2 && ((key.front() == '"' && key.back() == '"') ||
                            (key.front() == '\'' && key.back() == '\''))) {
      key = key.substr(1, key.size() - 2);
    }
    if (key == "polarizability")
      return PropertyType::Alpha;
    if (key == "hyperpolarizability")
      return PropertyType::Beta;
    if (key == "raman")
      return PropertyType::Raman;
    MADNESS_EXCEPTION(std::string("Unknown property: " + key).c_str(), 0);
  }

  inline static void compute_polarizability(PropertyContext &ctx) {
    if (ctx.world.rank() == 0)
      madness::print("▶️ Computing polarizability α...");

    compute_alpha(
        ctx.world, ctx.solved_states.planned_states.generated_states.state_map,
        ctx.ground_ctx.ground, ctx.response_params.dipole_frequencies(),
        ctx.response_params.dipole_directions(), ctx.properties);

    ctx.properties.save();
  }

  inline static void compute_hyperpolarizability(PropertyContext &ctx) {
    if (ctx.world.rank() == 0)
      madness::print("▶️ Computing hyperpolarizability β...");

    auto dip_dirs = ctx.response_params.dipole_directions();
    ::compute_hyperpolarizability(ctx.world, ctx.ground_ctx.ground,
                                  ctx.response_params.dipole_frequencies(),
                                  dip_dirs, ctx.properties);

    ctx.properties.save();
  }

  inline static void compute_raman(PropertyContext &ctx,
                                   VibrationalResults &vib,
                                   RamanResults &raman) {
    if (ctx.world.rank() == 0) {
      madness::print("------------------------- Raman Computation "
                     "-------------------------");
    }

    vib = compute_hessian(
        ctx.world, ctx.solved_states.planned_states.generated_states.state_map,
        ctx.ground_ctx.ground, ctx.response_params.dipole_directions(),
        ctx.scf_calc);

    const double csg_factor = 142.9435756;
    Tensor<double> normal_modes = *vib.normalmodes_atomic;
    if (ctx.world.rank() == 0) {
      print("normal modes in atomic coordinates (au): \n", normal_modes);
    }
    auto vib_freq = *vib.frequencies * constants::au2invcm;

    std::vector<int> mode;
    for (int i = 0; i < vib_freq.size(); ++i) {

      if (abs(vib_freq(i)) > 1e-2) {
        mode.push_back(i);
        raman.vibrational_frequencies.push_back(vib_freq(i));
      }
    }
    auto nnmodes = normal_modes(_, Slice(mode[0], -1, 1));

    raman.normal_modes = nnmodes;

    if (ctx.world.rank() == 0) {
      print(mode);
      print(mode[0], mode[mode.size() - 1]);
    }
    auto alpha_derivatives =
        compute_Raman(ctx.world, ctx.ground_ctx.ground,
                      ctx.response_params.dipole_frequencies(),
                      ctx.response_params.dipole_directions(),
                      ctx.response_params.nuclear_directions(), ctx.properties);
    raman.polarization_frequencies = ctx.response_params.dipole_frequencies();
    ctx.properties.save();
    ctx.world.gop.fence();
    auto compute_alpha2 = [](const Tensor<double> &alpha) {
      auto alpha_mean = 0.0;
      for (int i = 0; i < 3; ++i) {
        alpha_mean += alpha(i, i);
      }
      alpha_mean *= (1.0 / 3.0);
      return alpha_mean * alpha_mean;
    };
    auto compute_beta2 = [](const Tensor<double> &alpha) {
      auto beta2 = 0.0;
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          beta2 +=
              0.5 * (3 * alpha(i, j) * alpha(i, j) - alpha(i, i) * alpha(j, j));
        }
      }
      return beta2;
    };

    // Raman Intensity linearly polarized light
    auto RamanIntensityL = [](const double alpha2, const double beta2) {
      return 45 * alpha2 + 4 * beta2; // in a.u.
    };
    auto DepolarizationRatio = [](const double alpha2, const double beta2) {
      return 3 * beta2 / (45 * alpha2 + 4 * beta2);
    };
    bool debug = ctx.world.rank() == 0 && ctx.response_params.print_level() > 1;
    if (debug) {
      for (size_t freq_idx = 0;
           freq_idx < raman.polarization_frequencies.size(); ++freq_idx) {
        print("Alpha derivatives for frequency ",
              raman.polarization_frequencies[freq_idx], " (a.u.): \n",
              alpha_derivatives[freq_idx]);
      };
    }
    // Compute Normal Mode
    for (size_t freq_idx = 0; freq_idx < raman.polarization_frequencies.size();
         ++freq_idx) {

      const auto &alpha_dxyz = alpha_derivatives[freq_idx];
      double pol_freq = raman.polarization_frequencies[freq_idx];
      // Save polarizability derivatives to Results
      raman.polarizability_derivatives.push_back(alpha_dxyz);
      auto alpha_qi = inner(alpha_dxyz, nnmodes);
      raman.polarizability_derivatives_normal_modes.push_back(alpha_qi);
      if (debug) {
        print("Alpha derivative projected onto normal modes for frequency ",
              pol_freq, " (a.u.): \n", alpha_qi);
      }
    }
    // 1 a.u. of photon energy → wavelength (nm)
    auto wavelength_nm_from_au = [](double omega_au) {
      // lambda(nm) = 1239.841973 eV*nm / (27.211386 eV * omega_au)
      //            ≈ 45.56335 / omega_au
      const double AU_TO_NM_FACTOR = 45.5633525316;
      return (omega_au > 0.0) ? (AU_TO_NM_FACTOR / omega_au) : 0.0;
    };
    struct ColumnSpec {
      int width;
      int precision;
    };

    static const ColumnSpec COL_FREQ{10, 2};  // "Freq." like Dalton (xx.xx)
    static const ColumnSpec COL_FLOAT{11, 6}; // six decimals for the rest
                                              //
    auto print_header = [&](std::ostream &os, double omega_au) {
      const double lambda_nm = wavelength_nm_from_au(omega_au);
      os << "     Raman related properties for freq.  " << std::fixed
         << std::setprecision(6) << std::setw(9) << omega_au
         << " au  = " << std::setw(9) << std::setprecision(2) << lambda_nm
         << " nm\n";
      os << "     "
            "----------------------------------------------------------"
            "----"
            "-\n\n";
      os << " Mode    Freq.     Alpha**2   Beta(a)**2   Pol.Int.   "
            "Depol.Int.  Dep. Ratio \n\n";
    };

    auto print_row = [&](std::ostream &os,
                         const RamanResults::RamanModeRow &r) {
      // Compute missing fields if not provided
      double pol = r.pol_int ? *r.pol_int : (45.0 * r.alpha2 + 4.0 * r.beta2);
      double depi = r.depol_int ? *r.depol_int : (3.0 * r.beta2);
      double rho;
      if (r.dep_ratio) {
        rho = *r.dep_ratio;
      } else {
        rho = (pol != 0.0) ? (depi / pol) : 0.0; // guard div-by-zero
      }

      // Mode index
      os << std::setw(5) << r.mode
         << " "
         // Freq.
         << std::fixed << std::setw(COL_FREQ.width)
         << std::setprecision(COL_FREQ.precision) << r.freq_cm1
         << "  "
         // Alpha**2
         << std::setw(COL_FLOAT.width) << std::setprecision(COL_FLOAT.precision)
         << r.alpha2
         // Beta(a)**2
         << std::setw(COL_FLOAT.width) << std::setprecision(COL_FLOAT.precision)
         << r.beta2
         // Pol.Int.
         << std::setw(COL_FLOAT.width) << std::setprecision(COL_FLOAT.precision)
         << pol
         // Depol.Int.
         << std::setw(COL_FLOAT.width) << std::setprecision(COL_FLOAT.precision)
         << depi
         // Dep. Ratio
         << std::setw(COL_FLOAT.width) << std::setprecision(COL_FLOAT.precision)
         << rho << "\n";
    };

    auto print_block =
        [&](std::ostream &os, double omega_au,
            const std::vector<RamanResults::RamanModeRow> &rows) {
          print_header(os, omega_au);
          // print in reverse order to match Dalton output
          // for (auto const& x : range | std::views::reverse)
          for (auto r_it = rows.rbegin(); r_it != rows.rend(); ++r_it)
            print_row(os, *r_it);
          os << "\n";
        };
    using RamanModeRow = RamanResults::RamanModeRow;

    for (size_t freq_idx = 0; freq_idx < raman.polarization_frequencies.size();
         ++freq_idx) {
      auto alpha_qi = raman.polarizability_derivatives_normal_modes[freq_idx];
      vector<RamanModeRow> raman_rows;
      for (size_t i = 0; i < mode.size(); ++i) {
        RamanModeRow row;
        double vib_freq_i = vib_freq(mode[i]);
        row.mode = static_cast<int>(i) + 1;
        row.freq_cm1 = vib_freq_i;
        auto alpha_i = copy(alpha_qi(_, i));
        auto alpha = alpha_i.reshape(3, 3);
        auto alpha2_au = compute_alpha2(alpha);
        row.alpha2 = alpha2_au * csg_factor;
        auto beta2_au = compute_beta2(alpha);
        row.beta2 = beta2_au * csg_factor;
        auto depol_int_au = RamanIntensityL(alpha2_au, beta2_au);
        row.pol_int = depol_int_au * csg_factor;
        row.dep_ratio = DepolarizationRatio(alpha2_au, beta2_au);
        row.depol_int = (*row.dep_ratio) * (*row.pol_int);
        ctx.world.gop.fence();
        raman_rows.push_back(row);
      }
      raman.raman_spectra[raman.polarization_frequencies[freq_idx]] =
          raman_rows;
      ctx.world.gop.fence();
      if (ctx.world.rank() == 0) {
        print_block(std::cout, raman.polarization_frequencies[freq_idx],
                    raman_rows);
      }

      ctx.world.gop.fence();
    }
  }

  static PropertyStageOutput compute_requested_properties(
      World &world, const ResponseParameters &response_params,
      const GroundContext &ground_ctx, const SolvedStates &solved_states,
      const std::shared_ptr<SCF> &scf_calc) {
    VibrationalResults vib;
    RamanResults raman;
    PropertyManager properties(world, "properties.json");

    PropertyContext prop_ctx{world,         response_params, ground_ctx,
                             solved_states, properties,      scf_calc};

    for (const std::string &prop : response_params.requested_properties()) {
      auto prop_type = parse_property_name(prop);
      switch (prop_type) {
      case PropertyType::Alpha:
        compute_polarizability(prop_ctx);
        break;
      case PropertyType::Beta:
        compute_hyperpolarizability(prop_ctx);
        break;
      case PropertyType::Raman:
        compute_raman(prop_ctx, vib, raman);
        break;
      }
    }

    return PropertyStageOutput{properties.to_json(), std::move(vib),
                               std::move(raman)};
  }

  static PropertyStageOutput compute_requested_properties_with_property_group(
      World &world, const CalculationParameters &calc_params,
      const ResponseParameters &response_params,
      const GroundContext &ground_ctx, const SolvedStates &solved_states,
      const std::shared_ptr<SCF> &scf_calc) {
    const auto &state_parallel_plan =
        solved_states.planned_states.state_parallel_plan;
    const bool subgroup_property_mode =
        state_parallel_plan.execution_enabled &&
        state_parallel_plan.subgroup_parallel_enabled &&
        state_parallel_plan.execution_groups > 1;

    if (!subgroup_property_mode) {
      return compute_requested_properties(world, response_params, ground_ctx,
                                          solved_states, scf_calc);
    }

    const size_t property_group =
        response_params.state_parallel_property_group();
    if (property_group >= state_parallel_plan.execution_groups) {
      if (world.rank() == 0) {
        print("State-parallel property subgroup id ", property_group,
              " is out of range for execution_groups=",
              state_parallel_plan.execution_groups,
              ". Falling back to world property stage.");
      }
      return compute_requested_properties(world, response_params, ground_ctx,
                                          solved_states, scf_calc);
    }

    if (world.rank() == 0) {
      print("Property stage: executing on subgroup ", property_group, "/",
            state_parallel_plan.execution_groups - 1);
    }

    try {
      auto subworld_ptr = MacroTaskQ::create_worlds(
          world, state_parallel_plan.execution_groups);
      if (!subworld_ptr) {
        throw std::runtime_error("property subworld creation returned null");
      }

      World &subworld = *subworld_ptr;
      const size_t subgroup_id = static_cast<size_t>(
          world.rank() %
          static_cast<int>(state_parallel_plan.execution_groups));
      const bool in_property_group = (subgroup_id == property_group);

      auto old_pmap3 = FunctionDefaults<3>::get_pmap();
      auto restore_pmap = [&]() { FunctionDefaults<3>::set_pmap(old_pmap3); };

      std::string properties_dump;
      std::string vib_dump;
      std::string raman_dump;

      // Property assembly uses 3D response objects; only swap the 3D pmap.
      FunctionDefaults<3>::set_default_pmap(subworld);
      try {
        if (in_property_group) {
          const auto &protocol = calc_params.protocol();
          const double final_thresh = protocol.empty()
                                          ? FunctionDefaults<3>::get_thresh()
                                          : protocol.back();
          const std::string property_fock_file =
              group_shard_file(ground_ctx.fock_json_file, property_group);

          GroundStateData local_ground(subworld, ground_ctx.archive_file,
                                       ground_ctx.molecule);
          ResponseManager local_response_manager(subworld, calc_params);
          local_response_manager.setProtocol(subworld, local_ground.getL(),
                                             final_thresh);
          local_ground.prepareOrbitals(subworld, FunctionDefaults<3>::get_k(),
                                       final_thresh);
          local_ground.computePreliminaries(
              subworld, *local_response_manager.getCoulombOp(),
              local_response_manager.getVtol(), property_fock_file);

          GroundContext local_ground_ctx{
              ground_ctx.molecule, std::move(local_ground),
              std::move(local_response_manager), ground_ctx.archive_file,
              property_fock_file};

          PropertyStageOutput subgroup_output = compute_requested_properties(
              subworld, response_params, local_ground_ctx, solved_states,
              scf_calc);
          if (subworld.rank() == 0) {
            properties_dump = subgroup_output.properties.dump();
            vib_dump = subgroup_output.vibrational_analysis.to_json().dump();
            raman_dump = subgroup_output.raman_spectra.to_json().dump();
          }
        }
        subworld.gop.fence();
        restore_pmap();
      } catch (...) {
        restore_pmap();
        throw;
      }

      const int property_world_root = static_cast<int>(property_group);
      world.gop.broadcast_serializable(properties_dump, property_world_root);
      world.gop.broadcast_serializable(vib_dump, property_world_root);
      world.gop.broadcast_serializable(raman_dump, property_world_root);

      PropertyStageOutput stage_output;
      stage_output.properties = properties_dump.empty()
                                    ? nlohmann::json::object()
                                    : nlohmann::json::parse(properties_dump);
      if (!vib_dump.empty()) {
        stage_output.vibrational_analysis.from_json(
            nlohmann::json::parse(vib_dump));
      }
      if (!raman_dump.empty()) {
        stage_output.raman_spectra.from_json(nlohmann::json::parse(raman_dump));
      }
      world.gop.fence();
      return stage_output;
    } catch (const std::exception &ex) {
      if (world.rank() == 0) {
        print("State-parallel property subgroup execution failed: ", ex.what(),
              " Falling back to world property stage.");
      }
      return compute_requested_properties(world, response_params, ground_ctx,
                                          solved_states, scf_calc);
    }
  }

public:
  /**
   * @brief Run the full molecular response & property workflow.
   *
   * @param world      The MADNESS world communicator
   * @param params     Unified parameters containing response and molecule info
   * @param outdir     Directory where all outputs will be written
   * @return Results   Structured JSON fragments: metadata + properties
   */
  inline static Results run_response(World &world, const Params &params,
                                     const std::shared_ptr<SCF> &scf_calc,
                                     const std::filesystem::path &outdir) {
    const auto &calc_params = params.get<CalculationParameters>();
    const auto &rp_copy = params.get<ResponseParameters>();

    if (world.rank() == 0) {
      json response_input_json = {};
      response_input_json["response"] =
          rp_copy.to_json_if_precedence("defined");
      print("response_input_json: ", response_input_json.dump(4));
      std::ofstream ofs("response.in");
      write_json_to_input_file(response_input_json, {"response"}, ofs);
      ofs.close();
    }
    world.gop.fence();
    commandlineparser parser;
    parser.set_keyval("input", "response.in");
    if (world.rank() == 0) {
      ::print("input filename: ", parser.value("input"));
    }

    auto response_params = ResponseParameters(world, parser);
    GroundContext ground_ctx =
        make_ground_context(world, calc_params, scf_calc, outdir);

    PlannedStates planned_states =
        plan_required_states(world, calc_params, ground_ctx, response_params);
    SolvedStates solved_states =
        solve_all_states(world, calc_params, ground_ctx, response_params,
                         std::move(planned_states));
    PropertyStageOutput property_stage =
        compute_requested_properties_with_property_group(
            world, calc_params, response_params, ground_ctx, solved_states,
            scf_calc);

    // finalize & stats
    if (world.rank() == 0) {
      madness::print(
          "\n✅ Molecular response & property calculation complete.");
    }
    world.gop.fence();
    world.gop.fence();
    print_stats(world);

    // aggregate JSON results
    Results results;
    results.metadata = std::move(solved_states.metadata);
    results.properties = std::move(property_stage.properties);
    results.debug_log = std::move(solved_states.debug_log);
    results.vibrational_analysis =
        property_stage.vibrational_analysis.to_json();
    results.raman_spectra = property_stage.raman_spectra.to_json();
    return results;
  }
}; // namespace molresponse_lib
