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
#include <apps/molresponse_v2/ExcitedStateBundleSolver.hpp>
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
#include <numeric>
#include <optional>
#include <memory>
#include <utility>
#include <madness/chem/InputWriter.hpp>
#include <madness/chem/ParameterManager.hpp>
#include <madness/chem/Results.h>
#include <madness/mra/macrotaskq.h>

/*
Developer Overview
- Orchestration file for molresponse, structured as three stages:
  1) bootstrap + planning, 2) state solves, 3) property assembly.
- Terminology used in this file:
  - perturbation channel: one perturbation descriptor (e.g. Dipole_x)
  - frequency series: ordered frequencies for one perturbation channel
  - channel point: one (channel, frequency, protocol) solve target
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
                              const std::string &debug_file,
                              nlohmann::json baseline_metadata =
                                  nlohmann::json::object())
        : response_record_(world, meta_file), debug_logger_(debug_file),
          baseline_metadata_(std::move(baseline_metadata)) {}

    void
    initialize_states(const std::vector<LinearResponseDescriptor> &states) {
      response_record_.initialize_states(states);
    }

    void print_summary() const { response_record_.print_summary(); }

    [[nodiscard]] bool is_saved(const LinearResponsePoint &pt) const override {
      return response_record_.is_saved(pt.perturbationDescription(),
                                       pt.threshold(), pt.frequency()) ||
             point_ready_in_baseline(pt, /*require_converged=*/false);
    }

    [[nodiscard]] bool
    is_converged(const LinearResponsePoint &pt) const override {
      return response_record_.is_converged(pt.perturbationDescription(),
                                           pt.threshold(), pt.frequency()) ||
             point_ready_in_baseline(pt, /*require_converged=*/true);
    }

    void record_status(const LinearResponsePoint &pt, bool c) override {
      response_record_.record_status(pt, c);
    }

    void record_timing(const LinearResponsePoint &pt, double wall_seconds,
                       double cpu_seconds) override {
      response_record_.record_timing(pt, wall_seconds, cpu_seconds);
    }

    void record_restart_provenance(
        const LinearResponsePoint &pt, const std::string &source_kind,
        bool loaded_from_disk, bool promoted_from_static,
        const std::optional<double> &source_protocol,
        const std::optional<double> &source_frequency) override {
      response_record_.record_restart_provenance(
          pt, source_kind, loaded_from_disk, promoted_from_static,
          source_protocol, source_frequency);
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
    [[nodiscard]] bool
    point_ready_in_baseline(const LinearResponsePoint &pt,
                            bool require_converged) const {
      if (!baseline_metadata_.is_object() ||
          !baseline_metadata_.contains("states") ||
          !baseline_metadata_["states"].is_object()) {
        return false;
      }
      const std::string state_key = pt.perturbationDescription();
      const std::string protocol_key = ResponseRecord2::protocol_key(pt.threshold());
      const std::string freq_key = ResponseRecord2::freq_key(pt.frequency());

      const auto states_it = baseline_metadata_["states"].find(state_key);
      if (states_it == baseline_metadata_["states"].end() ||
          !states_it->contains("protocols") ||
          !(*states_it)["protocols"].is_object()) {
        return false;
      }
      const auto protos_it = (*states_it)["protocols"].find(protocol_key);
      if (protos_it == (*states_it)["protocols"].end()) {
        return false;
      }
      if (!protos_it->contains("saved") || !(*protos_it)["saved"].is_object()) {
        return false;
      }
      const auto saved_it = (*protos_it)["saved"].find(freq_key);
      if (saved_it == (*protos_it)["saved"].end() ||
          !saved_it->is_boolean() || !saved_it->get<bool>()) {
        return false;
      }
      if (!require_converged) {
        return true;
      }
      if (!protos_it->contains("converged") ||
          !(*protos_it)["converged"].is_object()) {
        return false;
      }
      const auto converged_it = (*protos_it)["converged"].find(freq_key);
      return converged_it != (*protos_it)["converged"].end() &&
             converged_it->is_boolean() && converged_it->get<bool>();
    }

    ResponseRecord2 response_record_;
    ResponseDebugLogger debug_logger_;
    nlohmann::json baseline_metadata_;
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

  static nlohmann::json
  broadcast_json_object(World &world, nlohmann::json payload,
                        int root_rank = 0) {
    std::string payload_dump;
    if (world.rank() == root_rank) {
      payload_dump = payload.dump();
    }
    world.gop.broadcast_serializable(payload_dump, root_rank);
    if (payload_dump.empty()) {
      return nlohmann::json::object();
    }
    return nlohmann::json::parse(payload_dump, nullptr,
                                 /*allow_exceptions=*/true);
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
        if (!dst_proto.contains("restart_provenance") ||
            !dst_proto["restart_provenance"].is_object()) {
          dst_proto["restart_provenance"] = nlohmann::json::object();
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
        if (shard_proto.contains("restart_provenance") &&
            shard_proto["restart_provenance"].is_object()) {
          for (const auto &[freq_key, provenance_value] :
               shard_proto["restart_provenance"].items()) {
            if (!dst_proto["restart_provenance"].contains(freq_key)) {
              dst_proto["restart_provenance"][freq_key] = provenance_value;
            }
          }
          if (src_protocol.contains("iterations") &&
              src_protocol["iterations"].is_number_unsigned()) {
            const auto src_iterations =
                src_protocol["iterations"].get<size_t>();
            const auto dst_iterations =
                (dst_protocol.contains("iterations") &&
                 dst_protocol["iterations"].is_number_unsigned())
                    ? dst_protocol["iterations"].get<size_t>()
                    : static_cast<size_t>(0);
            dst_protocol["iterations"] =
                std::max(dst_iterations, src_iterations);
          }
          if (src_protocol.contains("residual_norms") &&
              src_protocol["residual_norms"].is_array()) {
            if (!dst_protocol.contains("residual_norms") ||
                !dst_protocol["residual_norms"].is_array() ||
                dst_protocol["residual_norms"].empty()) {
              dst_protocol["residual_norms"] = src_protocol["residual_norms"];
            }
          }
          if (src_protocol.contains("iteration_max_residuals") &&
              src_protocol["iteration_max_residuals"].is_array()) {
            if (!dst_protocol.contains("iteration_max_residuals") ||
                !dst_protocol["iteration_max_residuals"].is_array() ||
                dst_protocol["iteration_max_residuals"].size() <
                    src_protocol["iteration_max_residuals"].size()) {
              dst_protocol["iteration_max_residuals"] =
                  src_protocol["iteration_max_residuals"];
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
    //print("Read molecule with ", molecule.natom(), " atoms.");

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
          " channel_owner_groups=", state_parallel_plan.channel_owner_groups,
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

  struct ProtocolExecutionPolicy {
    bool use_channel_series_mode = true;
    size_t active_groups = 1;
    size_t pending_channels = 0;
    size_t pending_points = 0;

    [[nodiscard]] nlohmann::json to_json() const {
      return {{"use_channel_series_mode", use_channel_series_mode},
              // Legacy key kept for downstream tooling compatibility.
              {"use_state_mode", use_channel_series_mode},
              {"active_groups", active_groups},
              {"pending_channels", pending_channels},
              // Legacy key kept for downstream tooling compatibility.
              {"pending_states", pending_channels},
              {"pending_points", pending_points}};
    }
  };

  struct RuntimePointOwnershipPolicy {
    size_t point_parallel_start_protocol_index = 0;
    bool restart_protocol0_saved_complete = false;
    bool restart_point_parallel_promoted = false;
    std::vector<ProtocolExecutionPolicy> protocol_policies;
  };

  static std::vector<size_t>
  build_owner_by_channel_index(const std::vector<LinearResponseDescriptor> &states,
                             const StateParallelPlan &state_parallel_plan) {
    std::vector<size_t> owner_by_channel_index(states.size(), 0);
    for (const auto &assignment : state_parallel_plan.channel_assignments) {
      if (assignment.channel_index < owner_by_channel_index.size()) {
        owner_by_channel_index[assignment.channel_index] =
            assignment.owner_group;
      }
    }
    return owner_by_channel_index;
  }

  static std::vector<ProtocolExecutionPolicy>
  build_protocol_execution_policy(World &world,
                                  const CalculationParameters &calc_params,
                                  const StateParallelPlan &state_parallel_plan,
                                  const std::vector<LinearResponseDescriptor>
                                      &linear_states,
                                  bool owner_group_schedule,
                                  size_t point_parallel_start_protocol_index) {
    const auto &protocol = calc_params.protocol();
    std::vector<ProtocolExecutionPolicy> protocol_policies(protocol.size());
    if (protocol.empty()) {
      return protocol_policies;
    }

    if (world.rank() == 0) {
      const nlohmann::json existing_metadata =
          read_json_file_or_object("response_metadata.json");
      const size_t point_owner_groups =
          std::max<size_t>(1, state_parallel_plan.effective_point_groups);
      const size_t channel_owner_groups =
          std::max<size_t>(1, state_parallel_plan.channel_owner_groups);
      for (size_t ti = 0; ti < protocol.size(); ++ti) {
        ProtocolExecutionPolicy policy;
        policy.use_channel_series_mode =
            !owner_group_schedule || ti < point_parallel_start_protocol_index;
        const bool at_final_protocol = (ti + 1 == protocol.size());

        for (const auto &state : linear_states) {
          bool state_has_pending = false;
          for (size_t freq_idx = 0; freq_idx < state.num_frequencies();
               ++freq_idx) {
            LinearResponsePoint pt{state, ti, freq_idx};
            const bool point_ready = point_ready_in_metadata(
                existing_metadata, pt, /*require_saved=*/true,
                /*require_converged=*/at_final_protocol);
            if (!point_ready) {
              ++policy.pending_points;
              state_has_pending = true;
            }
          }
          if (state_has_pending) {
            ++policy.pending_channels;
          }
        }

        // Phase-0 policy:
        // - if protocol-0 still has pending points, keep strict state ownership
        //   regardless of requested point-start protocol.
        // - if protocol-0 is complete (restart), runtime may use point mode.
        if (owner_group_schedule && ti == 0 && policy.pending_points > 0) {
          policy.use_channel_series_mode = true;
        }

        if (!owner_group_schedule || policy.pending_points == 0) {
          policy.active_groups = 1;
        } else if (policy.use_channel_series_mode) {
          // Keep full state-owner lanes active so fixed owner mapping remains
          // stable across restarts/partial completion.
          policy.active_groups = channel_owner_groups;
        } else {
          policy.active_groups = std::max<size_t>(
              1, std::min(point_owner_groups, policy.pending_points));
        }
        protocol_policies[ti] = policy;
      }
    }

    nlohmann::json payload = nlohmann::json::array();
    if (world.rank() == 0) {
      for (const auto &policy : protocol_policies) {
        payload.push_back(policy.to_json());
      }
    }
    payload = broadcast_json_object(world, std::move(payload), 0);
    if (!payload.is_array()) {
      return protocol_policies;
    }
    protocol_policies.clear();
    protocol_policies.reserve(payload.size());
    for (const auto &entry : payload) {
      ProtocolExecutionPolicy policy;
      if (entry.is_object()) {
        policy.use_channel_series_mode =
            entry.value("use_channel_series_mode",
                        entry.value("use_state_mode", true));
        policy.active_groups = std::max<size_t>(1, entry.value("active_groups", 1));
        policy.pending_channels =
            entry.value("pending_channels",
                        entry.value("pending_states", size_t(0)));
        policy.pending_points = entry.value("pending_points", size_t(0));
      }
      protocol_policies.push_back(policy);
    }
    return protocol_policies;
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
    runtime_policy.protocol_policies = build_protocol_execution_policy(
        world, calc_params, state_parallel_plan, linear_states,
        owner_group_schedule, runtime_policy.point_parallel_start_protocol_index);
    return runtime_policy;
  }

  static void build_local_channel_workset(
      const std::vector<LinearResponseDescriptor> &linear_states,
      const std::vector<size_t> &owner_by_channel_index,
      const PointOwnershipScheduler &point_scheduler, size_t subgroup_id,
      size_t mapping_groups, std::vector<size_t> &local_channel_indices,
      std::vector<LinearResponseDescriptor> &local_channels) {
    local_channel_indices.clear();
    local_channels.clear();
    local_channel_indices.reserve(
        linear_states.size() / std::max<size_t>(1, mapping_groups));
    local_channels.reserve(linear_states.size());

    const size_t point_owner_groups = point_scheduler.owner_groups();
    for (size_t channel_index = 0; channel_index < linear_states.size();
         ++channel_index) {
      const bool owns_channel_at_first_protocol =
          owner_by_channel_index[channel_index] == subgroup_id;
      bool owns_any_point_after_first_protocol = false;
      if (point_owner_groups > 1) {
        const auto &channel = linear_states[channel_index];
        for (size_t freq_index = 0; freq_index < channel.num_frequencies();
             ++freq_index) {
          if (point_scheduler.owner_group(channel_index, freq_index) ==
              subgroup_id) {
            owns_any_point_after_first_protocol = true;
            break;
          }
        }
      }
      if (owns_channel_at_first_protocol) {
        local_channel_indices.push_back(channel_index);
      }
      if (owns_channel_at_first_protocol || owns_any_point_after_first_protocol) {
        local_channels.push_back(linear_states[channel_index]);
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
    // channel_index -> owner lane for protocol ranges using channel-series
    // ownership.
    std::vector<size_t> owner_by_channel_index;
    // Deterministic (channel,freq) -> owner lane mapping for channel-point
    // ownership.
    PointOwnershipScheduler point_scheduler;
    // Runtime protocol index where channel-series ownership switches to
    // channel-point ownership.
    size_t runtime_point_parallel_start_protocol_index = 0;
    // Restart diagnostics propagated to stage metadata.
    bool restart_point_parallel_promoted = false;
    bool restart_protocol0_saved_complete = false;
    // Runtime per-protocol mode + active owner groups after restart analysis.
    std::vector<ProtocolExecutionPolicy> runtime_protocol_policies;

    StateSolveScheduleContext(
        bool subgroup_parallel_requested_, bool owner_group_schedule_,
        const StateParallelPlan &state_parallel_plan_,
        const std::vector<LinearResponseDescriptor> &linear_states_,
        std::vector<size_t> owner_by_channel_index_, size_t effective_point_groups,
        size_t runtime_point_parallel_start_protocol_index_,
        bool restart_point_parallel_promoted_,
        bool restart_protocol0_saved_complete_,
        std::vector<ProtocolExecutionPolicy> runtime_protocol_policies_)
        : subgroup_parallel_requested(subgroup_parallel_requested_),
          owner_group_schedule(owner_group_schedule_),
          state_parallel_plan(state_parallel_plan_),
          linear_states(linear_states_),
          owner_by_channel_index(std::move(owner_by_channel_index_)),
          point_scheduler(linear_states_, effective_point_groups),
          runtime_point_parallel_start_protocol_index(
              runtime_point_parallel_start_protocol_index_),
          restart_point_parallel_promoted(restart_point_parallel_promoted_),
          restart_protocol0_saved_complete(
              restart_protocol0_saved_complete_),
          runtime_protocol_policies(std::move(runtime_protocol_policies_)) {}

    [[nodiscard]] size_t point_owner_groups() const {
      return point_scheduler.owner_groups();
    }

    [[nodiscard]] ProtocolExecutionPolicy
    protocol_policy(size_t protocol_index) const {
      if (protocol_index < runtime_protocol_policies.size()) {
        return runtime_protocol_policies[protocol_index];
      }
      ProtocolExecutionPolicy fallback;
      fallback.use_channel_series_mode =
          protocol_index < runtime_point_parallel_start_protocol_index;
      fallback.active_groups = fallback.use_channel_series_mode
                                   ? std::max<size_t>(
                                         1, state_parallel_plan.channel_owner_groups)
                                   : std::max<size_t>(1, point_owner_groups());
      return fallback;
    }
  };

  // Runtime policy gate shared by serial and subgroup paths.
  // For protocol indices below the runtime switch threshold we keep all
  // frequencies of a channel together; after the switch threshold we fan out by
  // independent channel-frequency points.
  static bool
  use_channel_series_ownership_for_protocol_runtime(const StateSolveScheduleContext &ctx,
                                           size_t protocol_index) {
    if (ctx.state_parallel_plan.mapping_groups <= 1) {
      return true;
    }
    return ctx.protocol_policy(protocol_index).use_channel_series_mode;
  }

  static size_t active_owner_groups_for_protocol_runtime(
      const StateSolveScheduleContext &ctx, size_t protocol_index) {
    if (ctx.state_parallel_plan.mapping_groups <= 1) {
      return 1;
    }
    return std::max<size_t>(1, ctx.protocol_policy(protocol_index).active_groups);
  }

  template <typename PointFn>
  static void
  for_each_owned_point_in_point_mode(const StateSolveScheduleContext &schedule_ctx,
                                     size_t lane_id, PointFn &&fn) {
    for (size_t state_index = 0; state_index < schedule_ctx.linear_states.size();
         ++state_index) {
      const auto &state = schedule_ctx.linear_states[state_index];
      for (size_t freq_index = 0; freq_index < state.num_frequencies();
           ++freq_index) {
        if (schedule_ctx.point_scheduler.owner_group(state_index, freq_index) ==
            lane_id) {
          fn(state_index, freq_index);
        }
      }
    }
  }

  template <typename StateFn>
  static void
  for_each_owned_channel_in_channel_series_mode(const StateSolveScheduleContext &schedule_ctx,
                                     size_t lane_id, StateFn &&fn) {
    for (size_t state_index = 0; state_index < schedule_ctx.linear_states.size();
         ++state_index) {
      if (schedule_ctx.owner_by_channel_index[state_index] == lane_id) {
        fn(state_index);
      }
    }
  }

  template <typename PersistenceT>
  static bool
  point_needs_solving(PersistenceT &persistence, const LinearResponsePoint &pt,
                      bool at_final_protocol) {
    const bool is_saved = persistence.is_saved(pt);
    return !is_saved || (at_final_protocol && !persistence.is_converged(pt));
  }

  static bool point_needs_solving_from_metadata(
      const nlohmann::json &metadata, const LinearResponsePoint &pt,
      bool at_final_protocol) {
    const bool is_saved = point_ready_in_metadata(metadata, pt,
                                                  /*require_saved=*/true,
                                                  /*require_converged=*/false);
    if (!is_saved) {
      return true;
    }
    if (!at_final_protocol) {
      return false;
    }
    return !point_ready_in_metadata(metadata, pt,
                                    /*require_saved=*/true,
                                    /*require_converged=*/true);
  }

  template <typename PointNeedsFn>
  static bool any_state_point_needs_solving(const std::vector<size_t> &state_indices,
                                            const StateSolveScheduleContext &schedule_ctx,
                                            size_t thresh_index,
                                            PointNeedsFn &&point_needs_solving_fn) {
    for (const auto state_index : state_indices) {
      const auto &state = schedule_ctx.linear_states[state_index];
      for (size_t freq_idx = 0; freq_idx < state.num_frequencies(); ++freq_idx) {
        if (point_needs_solving_fn(state, freq_idx)) {
          return true;
        }
      }
    }
    return false;
  }

  static LinearResponseDescriptor
  make_single_frequency_state(const LinearResponseDescriptor &state,
                              size_t freq_index) {
    return LinearResponseDescriptor(state.perturbation,
                                    {state.frequency(freq_index)},
                                    state.thresholds, state.spin_restricted);
  }

  template <typename NeedsFn, typename SkipFn, typename PrepareFn,
            typename ExecuteFn, typename FinalizeFn>
  static void run_protocol_threshold_loop(const std::vector<double> &protocol,
                                          NeedsFn &&needs_solving_at_protocol,
                                          SkipFn &&on_skip_protocol,
                                          PrepareFn &&prepare_protocol,
                                          ExecuteFn &&execute_protocol_work,
                                          FinalizeFn &&finalize_protocol) {
    for (size_t ti = 0; ti < protocol.size(); ++ti) {
      const double thresh = protocol[ti];
      if (!needs_solving_at_protocol(thresh, ti)) {
        on_skip_protocol(thresh, ti);
        continue;
      }

      prepare_protocol(thresh, ti);
      const bool at_final_protocol = (ti + 1 == protocol.size());
      execute_protocol_work(thresh, ti, at_final_protocol);
      finalize_protocol(thresh, ti);
    }
  }

  template <typename PersistenceT, typename StateT>
  static void run_frequency_loop_with_flush(
      World &exec_world, ResponseManager &response_manager,
      const StateT &state_for_compute, size_t thresh_index,
      GroundStateData &ground_state, PersistenceT &persistence,
      bool at_final_protocol) {
    computeFrequencyLoop(exec_world, response_manager, state_for_compute,
                         thresh_index, ground_state, persistence,
                         at_final_protocol);
    persistence.flush_debug_log(exec_world);
  }

  // Shared work dispatcher for ownership-aware linear state scheduling.
  // It handles both modes:
  // - state mode: one lane owns all frequencies of each state
  // - point mode: lanes own independent (state,frequency) points
  // Optional `single_lane_state_indices` allows subgroup execution to avoid
  // rescanning the full owner map.
  template <typename LaneBeginFn, typename SolveStateFn, typename SolvePointFn>
  static void dispatch_owned_work_for_protocol(
      const StateSolveScheduleContext &schedule_ctx, size_t thresh_index,
      size_t lane_begin, size_t lane_end,
      const std::vector<size_t> *single_lane_state_indices,
      LaneBeginFn &&on_lane_begin, SolveStateFn &&solve_state,
      SolvePointFn &&solve_state_frequency) {
    const bool use_channel_series_mode =
        use_channel_series_ownership_for_protocol_runtime(schedule_ctx, thresh_index);
    const bool single_lane = (lane_end == lane_begin + 1);
    for (size_t lane_id = lane_begin; lane_id < lane_end; ++lane_id) {
      on_lane_begin(lane_id, use_channel_series_mode);
      if (use_channel_series_mode) {
        if (single_lane && single_lane_state_indices != nullptr) {
          for (const auto state_index : *single_lane_state_indices) {
            solve_state(state_index);
          }
        } else {
          for_each_owned_channel_in_channel_series_mode(
              schedule_ctx, lane_id, [&](size_t state_index) {
                solve_state(state_index);
              });
        }
      } else {
        for_each_owned_point_in_point_mode(
            schedule_ctx, lane_id, [&](size_t state_index, size_t freq_index) {
              solve_state_frequency(state_index, freq_index);
            });
      }
    }
  }

  struct PendingPointWorkItem {
    size_t state_index = 0;
    size_t freq_index = 0;
  };

  struct PendingProtocolManifest {
    // Ownership mode used for this protocol index.
    bool use_channel_series_mode = true;
    // State-level work list for state ownership mode.
    std::vector<size_t> pending_channel_indices;
    // Point-level work list for point ownership mode.
    std::vector<PendingPointWorkItem> pending_points;

    [[nodiscard]] bool has_work() const {
      return !pending_channel_indices.empty() || !pending_points.empty();
    }
  };

  template <typename PointNeedsFn>
  static PendingProtocolManifest build_pending_work_manifest(
      const StateSolveScheduleContext &schedule_ctx, size_t thresh_index,
      size_t lane_begin, size_t lane_end,
      const std::vector<size_t> *single_lane_state_indices,
      PointNeedsFn &&point_needs_solving_fn) {
    (void)single_lane_state_indices;
    PendingProtocolManifest manifest;
    manifest.use_channel_series_mode =
        use_channel_series_ownership_for_protocol_runtime(schedule_ctx, thresh_index);
    const size_t active_groups =
        active_owner_groups_for_protocol_runtime(schedule_ctx, thresh_index);
    const size_t lane_end_clamped = std::min(lane_end, active_groups);
    if (lane_begin >= lane_end_clamped) {
      return manifest;
    }

    std::vector<PendingPointWorkItem> pending_points_all;
    std::vector<size_t> pending_channel_indices_all;
    pending_points_all.reserve(schedule_ctx.linear_states.size());
    pending_channel_indices_all.reserve(schedule_ctx.linear_states.size());
    for (size_t state_index = 0; state_index < schedule_ctx.linear_states.size();
         ++state_index) {
      const auto &state = schedule_ctx.linear_states[state_index];
      bool state_has_pending = false;
      for (size_t freq_idx = 0; freq_idx < state.num_frequencies(); ++freq_idx) {
        if (point_needs_solving_fn(state, freq_idx)) {
          state_has_pending = true;
          pending_points_all.push_back({state_index, freq_idx});
        }
      }
      if (state_has_pending) {
        pending_channel_indices_all.push_back(state_index);
      }
    }

    if (manifest.use_channel_series_mode) {
      std::vector<char> owned_states(schedule_ctx.linear_states.size(), 0);
      for (const auto state_index : pending_channel_indices_all) {
        const size_t owner_lane =
            schedule_ctx.owner_by_channel_index[state_index] %
            std::max<size_t>(1, active_groups);
        if (owner_lane >= lane_begin && owner_lane < lane_end_clamped) {
          manifest.pending_channel_indices.push_back(state_index);
          owned_states[state_index] = 1;
        }
      }
      for (const auto &pending_point : pending_points_all) {
        if (owned_states[pending_point.state_index]) {
          manifest.pending_points.push_back(pending_point);
        }
      }
    } else {
      for (size_t idx = 0; idx < pending_points_all.size(); ++idx) {
        const size_t owner_lane = idx % active_groups;
        if (owner_lane >= lane_begin && owner_lane < lane_end_clamped) {
          manifest.pending_points.push_back(pending_points_all[idx]);
        }
      }
      for (const auto &pending_point : manifest.pending_points) {
        if (manifest.pending_channel_indices.empty() ||
            manifest.pending_channel_indices.back() != pending_point.state_index) {
          manifest.pending_channel_indices.push_back(pending_point.state_index);
        }
      }
    }
    return manifest;
  }

  static PendingProtocolManifest build_pending_manifest_from_metadata(
      const StateSolveScheduleContext &schedule_ctx, size_t thresh_index,
      size_t lane_begin, size_t lane_end, const nlohmann::json &metadata,
      bool at_final_protocol) {
    auto needs_solving = [&](const LinearResponseDescriptor &state,
                             size_t freq_idx) {
      LinearResponsePoint pt{state, thresh_index, freq_idx};
      return point_needs_solving_from_metadata(metadata, pt, at_final_protocol);
    };
    return build_pending_work_manifest(schedule_ctx, thresh_index, lane_begin,
                                       lane_end, nullptr, needs_solving);
  }

  static void
  print_state_solve_execution_mode(World &world,
                                   const StateParallelPlan &state_parallel_plan,
                                   bool subgroup_parallel_requested,
                                   bool owner_group_schedule) {
    if (world.rank() == 0 && subgroup_parallel_requested) {
      print("Channel-series ownership mapping is active across ",
            state_parallel_plan.mapping_groups,
            " groups (protocol-0 owner groups=",
            state_parallel_plan.channel_owner_groups,
            "); executing owner-group solves in parallel subworlds.");
    } else if (world.rank() == 0 && owner_group_schedule) {
      print("Channel-series ownership mapping is active across ",
            state_parallel_plan.mapping_groups,
            " groups (protocol-0 owner groups=",
            state_parallel_plan.channel_owner_groups,
            "); executing deterministic owner-group solve passes "
            "serially on the universe communicator.");
    } else if (world.rank() == 0 && state_parallel_plan.mapping_groups > 1) {
      print("Channel-series ownership mapping is active across ",
            state_parallel_plan.mapping_groups,
            " groups; falling back to plain serial channel loop.");
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
    auto owner_by_channel_index =
        build_owner_by_channel_index(linear_states, state_parallel_plan);
    const RuntimePointOwnershipPolicy runtime_policy =
        compute_runtime_point_ownership_policy(
            world, calc_params, state_parallel_plan, linear_states,
            owner_group_schedule);
    if (runtime_policy.restart_point_parallel_promoted && world.rank() == 0) {
      print("State-parallel restart detected: protocol-0 points are saved; "
            "enabling channel-point ownership from protocol index 0.");
    }
    if (world.rank() == 0 && !runtime_policy.protocol_policies.empty()) {
      for (size_t ti = 0; ti < runtime_policy.protocol_policies.size(); ++ti) {
        const auto &policy = runtime_policy.protocol_policies[ti];
        print("Protocol policy ti=", ti, " mode=",
              policy.use_channel_series_mode ? "channel_series"
                                             : "channel_point",
              " active_groups=", policy.active_groups,
              " pending_channels=", policy.pending_channels,
              " pending_points=", policy.pending_points);
      }
      if (runtime_policy.protocol_policies[0].use_channel_series_mode &&
          runtime_policy.protocol_policies[0].pending_points > 0) {
        print("Protocol-0 channel-series ownership (channel -> owner_group):");
        const size_t max_rows_to_print = 64;
        const size_t rows_to_print =
            std::min(max_rows_to_print, state_parallel_plan.channel_assignments.size());
        for (size_t i = 0; i < rows_to_print; ++i) {
          const auto &assignment = state_parallel_plan.channel_assignments[i];
          print("  ", assignment.channel_label, " -> ", assignment.owner_group);
        }
        if (state_parallel_plan.channel_assignments.size() > rows_to_print) {
          print("  ... ",
                state_parallel_plan.channel_assignments.size() - rows_to_print,
                " additional states omitted");
        }
      }
    }

    return StateSolveScheduleContext(
        subgroup_parallel_requested, owner_group_schedule, state_parallel_plan,
        linear_states, std::move(owner_by_channel_index),
        state_parallel_plan.effective_point_groups,
        runtime_policy.point_parallel_start_protocol_index,
        runtime_policy.restart_point_parallel_promoted,
        runtime_policy.restart_protocol0_saved_complete,
        std::move(runtime_policy.protocol_policies));
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

    const auto &protocol = calc_params.protocol();
    std::vector<std::optional<PendingProtocolManifest>> pending_manifest_by_ti(
        protocol.size());

    auto needs_solving_at_protocol = [&](double protocol_thresh,
                                         size_t thresh_index) {
      // A protocol threshold is considered "active" when at least one point is
      // missing on disk, or not converged at the final threshold.
      const bool at_final_protocol = protocol_thresh == calc_params.protocol().back();
      auto serial_point_needs_solving = [&](const LinearResponseDescriptor &state,
                                            size_t freq_idx) {
        LinearResponsePoint pt{state, thresh_index, freq_idx};
        const bool is_saved = persistence.is_saved(pt);
        const bool should_solve =
            point_needs_solving(persistence, pt, at_final_protocol);
        if (world.rank() == 0) {
          print("Checking state ", pt.perturbationDescription(), " at thresh ",
                protocol_thresh, " freq ", pt.frequency(), " is_saved=",
                is_saved, " at_final_protocol=", at_final_protocol,
                " should_solve=", should_solve);
        }
        return should_solve;
      };

      if (!schedule_ctx.owner_group_schedule) {
        std::vector<size_t> all_state_indices(schedule_ctx.linear_states.size());
        std::iota(all_state_indices.begin(), all_state_indices.end(), 0);
        return any_state_point_needs_solving(all_state_indices, schedule_ctx,
                                             thresh_index,
                                             serial_point_needs_solving);
      }

      const size_t lane_count =
          active_owner_groups_for_protocol_runtime(schedule_ctx, thresh_index);
      const auto manifest =
          build_pending_work_manifest(schedule_ctx, thresh_index, 0, lane_count,
                                      nullptr, serial_point_needs_solving);
      pending_manifest_by_ti[thresh_index] = manifest;
      return manifest.has_work();
    };

    auto solve_state = [&](size_t state_index, size_t thresh_index,
                           bool at_final_protocol) {
      auto &state = schedule_ctx.linear_states[state_index];
      run_frequency_loop_with_flush(world, ctx.response_manager, state,
                                    thresh_index, ctx.ground, persistence,
                                    at_final_protocol);
    };
    auto solve_state_frequency = [&](size_t state_index, size_t freq_index,
                                     size_t thresh_index,
                                     bool at_final_protocol) {
      const auto &state = schedule_ctx.linear_states[state_index];
      const auto single_frequency_state =
          make_single_frequency_state(state, freq_index);
      run_frequency_loop_with_flush(world, ctx.response_manager,
                                    single_frequency_state, thresh_index,
                                    ctx.ground, persistence, at_final_protocol);
    };
    run_protocol_threshold_loop(
        protocol, needs_solving_at_protocol,
        [&](double thresh, size_t /*thresh_index*/) {
          if (world.rank() == 0) {
            madness::print("✓ All states converged at thresh", thresh,
                           "skipping to next protocol.");
          }
        },
        [&](double thresh, size_t /*thresh_index*/) {
          ctx.response_manager.setProtocol(world, ctx.ground.getL(), thresh);
          ctx.ground.prepareOrbitals(world, FunctionDefaults<3>::get_k(), thresh);
          ctx.ground.computePreliminaries(world,
                                          *ctx.response_manager.getCoulombOp(),
                                          ctx.response_manager.getVtol(),
                                          ctx.fock_json_file);
        },
        [&](double thresh, size_t ti, bool at_final_protocol) {
          (void)thresh;
          if (!schedule_ctx.owner_group_schedule) {
            for (size_t state_index = 0;
                 state_index < schedule_ctx.linear_states.size(); ++state_index) {
              solve_state(state_index, ti, at_final_protocol);
            }
            return;
          }

          const size_t lane_count =
              active_owner_groups_for_protocol_runtime(schedule_ctx, ti);
          PendingProtocolManifest manifest;
          if (pending_manifest_by_ti[ti].has_value()) {
            manifest = *pending_manifest_by_ti[ti];
          } else {
            auto serial_point_needs_solving =
                [&](const LinearResponseDescriptor &state, size_t freq_idx) {
                  LinearResponsePoint pt{state, ti, freq_idx};
                  return point_needs_solving(persistence, pt, at_final_protocol);
                };
            manifest = build_pending_work_manifest(schedule_ctx, ti, 0,
                                                   lane_count, nullptr,
                                                   serial_point_needs_solving);
          }
          if (world.rank() == 0) {
            print("Protocol ", ti, " pending owned channels=",
                  manifest.pending_channel_indices.size(),
                  " pending owned points=", manifest.pending_points.size(),
                  " mode=", manifest.use_channel_series_mode
                                 ? "channel_series"
                                 : "channel_point");
          }
          if (manifest.use_channel_series_mode) {
            for (const auto state_index : manifest.pending_channel_indices) {
              solve_state(state_index, ti, at_final_protocol);
            }
          } else {
            for (const auto &work_item : manifest.pending_points) {
              solve_state_frequency(work_item.state_index, work_item.freq_index,
                                    ti, at_final_protocol);
            }
          }
        },
        [&](double /*thresh*/, size_t /*thresh_index*/) {});

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
    nlohmann::json baseline_state_metadata = nlohmann::json::object();
    if (world.rank() == 0) {
      baseline_state_metadata = read_json_file_or_object("response_metadata.json");
    }
    baseline_state_metadata =
        broadcast_json_object(world, std::move(baseline_state_metadata), 0);
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
        std::vector<size_t> local_channel_indices;
        std::vector<LinearResponseDescriptor> local_channels;
        build_local_channel_workset(
            schedule_ctx.linear_states, schedule_ctx.owner_by_channel_index,
            schedule_ctx.point_scheduler, subgroup_id, state_parallel_plan.mapping_groups,
            local_channel_indices, local_channels);

        const std::string metadata_shard_file =
            group_shard_file("response_metadata.json", subgroup_id);
        const std::string debug_shard_file =
            group_shard_file("response_log.json", subgroup_id);
        const std::string fock_shard_file =
            group_shard_file(ctx.fock_json_file, subgroup_id);

        JsonStateSolvePersistence local_persistence(
            subworld, metadata_shard_file, debug_shard_file,
            baseline_state_metadata);
        local_persistence.initialize_states(local_channels);
        if (subworld.rank() == 0) {
          print("State-parallel subgroup ", subgroup_id, " owns ",
                local_channel_indices.size(), " protocol-0 channels and ",
                local_channels.size(), " active channels across all protocols.");
          local_persistence.print_summary();
        }
        subworld.gop.fence();

        if (!local_channels.empty()) {
          GroundStateData local_ground(subworld, ctx.archive_file, ctx.molecule);
          ResponseManager local_response_manager(subworld, calc_params);

          const auto &protocol = calc_params.protocol();
          std::vector<std::optional<PendingProtocolManifest>>
              local_pending_manifest_by_ti(protocol.size());
          std::vector<char> local_manifest_ready(protocol.size(), 0);
          auto local_needs_solving_at_protocol = [&](double protocol_thresh,
                                                     size_t thresh_index) {
            if (local_manifest_ready[thresh_index] != 0) {
              return local_pending_manifest_by_ti[thresh_index]->has_work();
            }
            const bool at_final_protocol =
                protocol_thresh == calc_params.protocol().back();
            const auto manifest = build_pending_manifest_from_metadata(
                schedule_ctx, thresh_index, subgroup_id, subgroup_id + 1,
                baseline_state_metadata, at_final_protocol);
            local_pending_manifest_by_ti[thresh_index] = manifest;
            local_manifest_ready[thresh_index] = 1;
            return manifest.has_work();
          };

          run_protocol_threshold_loop(
              protocol, local_needs_solving_at_protocol,
              [&](double thresh, size_t /*thresh_index*/) {
                if (subworld.rank() == 0) {
                  print("Subgroup ", subgroup_id,
                        " has no pending channels at thresh ", thresh,
                        "; skipping.");
                }
                // Keep protocol progression globally lock-step across all
                // subgroups because ti>0 points may depend on ti-1 files
                // produced by other owners.
                world.gop.fence();
              },
              [&](double thresh, size_t /*thresh_index*/) {
                local_response_manager.setProtocol(subworld, local_ground.getL(),
                                                   thresh);
                local_ground.prepareOrbitals(subworld, FunctionDefaults<3>::get_k(),
                                             thresh);
                local_ground.computePreliminaries(
                    subworld, *local_response_manager.getCoulombOp(),
                    local_response_manager.getVtol(), fock_shard_file);
              },
              [&](double thresh, size_t ti, bool at_final_protocol) {
                PendingProtocolManifest manifest;
                if (local_pending_manifest_by_ti[ti].has_value()) {
                  manifest = *local_pending_manifest_by_ti[ti];
                } else {
                  manifest = build_pending_manifest_from_metadata(
                      schedule_ctx, ti, subgroup_id, subgroup_id + 1,
                      baseline_state_metadata, at_final_protocol);
                  local_pending_manifest_by_ti[ti] = manifest;
                  local_manifest_ready[ti] = 1;
                }
                if (subworld.rank() == 0) {
                  print("Subgroup ", subgroup_id, " protocol ", ti,
                        " pending channels=",
                        manifest.pending_channel_indices.size(),
                        " pending points=", manifest.pending_points.size(),
                        " mode=", manifest.use_channel_series_mode
                                      ? "channel_series"
                                      : "channel_point");
                }
                if (!manifest.use_channel_series_mode && subworld.rank() == 0) {
                  print("Subgroup ", subgroup_id,
                        " solving independent channel-frequency points at thresh ",
                        thresh, ".");
                }
                if (manifest.use_channel_series_mode) {
                  for (const auto state_index : manifest.pending_channel_indices) {
                    auto &state = schedule_ctx.linear_states[state_index];
                    run_frequency_loop_with_flush(
                        subworld, local_response_manager, state, ti,
                        local_ground, local_persistence, at_final_protocol);
                  }
                } else {
                  for (const auto &work_item : manifest.pending_points) {
                    const auto &state =
                        schedule_ctx.linear_states[work_item.state_index];
                    const auto single_frequency_state =
                        make_single_frequency_state(state, work_item.freq_index);
                    run_frequency_loop_with_flush(
                        subworld, local_response_manager, single_frequency_state,
                        ti, local_ground, local_persistence, at_final_protocol);
                  }
                }
              },
              [&](double /*thresh*/, size_t /*thresh_index*/) {
                // ti+1 must not begin until every subgroup has finished ti.
                world.gop.fence();
              });
        } else {
          if (subworld.rank() == 0) {
            print("Subgroup ", subgroup_id,
                  " has no owned channels; participating in protocol "
                  "synchronization only.");
          }
          // Keep protocol fences balanced with active subgroups. When a subgroup
          // owns no channels (e.g., more groups than points), it still must
          // enter
          // one global fence per protocol threshold to avoid deadlock.
          const auto &protocol = calc_params.protocol();
          for (size_t ti = 0; ti < protocol.size(); ++ti) {
            world.gop.fence();
          }
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

  struct ExcitedExecutionResult {
    // Execution summary JSON written into response metadata.
    nlohmann::json execution;
  };

  struct FinalProtocolState {
    double threshold = 0.0;
    size_t threshold_index = 0;
  };

  static void ensure_excited_protocol_placeholder_node(
      nlohmann::json &protocol_node) {
    if (!protocol_node.is_object()) {
      protocol_node = nlohmann::json::object();
    }
    if (!protocol_node.contains("saved") || !protocol_node["saved"].is_boolean()) {
      protocol_node["saved"] = false;
    }
    if (!protocol_node.contains("converged") ||
        !protocol_node["converged"].is_boolean()) {
      protocol_node["converged"] = false;
    }
    if (!protocol_node.contains("timings") ||
        !protocol_node["timings"].is_object()) {
      protocol_node["timings"] =
          nlohmann::json{{"wall_seconds", 0.0}, {"cpu_seconds", 0.0}};
    } else {
      if (!protocol_node["timings"].contains("wall_seconds") ||
          !protocol_node["timings"]["wall_seconds"].is_number()) {
        protocol_node["timings"]["wall_seconds"] = 0.0;
      }
      if (!protocol_node["timings"].contains("cpu_seconds") ||
          !protocol_node["timings"]["cpu_seconds"].is_number()) {
        protocol_node["timings"]["cpu_seconds"] = 0.0;
      }
    }
    if (!protocol_node.contains("energies") ||
        !protocol_node["energies"].is_array()) {
      protocol_node["energies"] = nlohmann::json::array();
    }
    if (!protocol_node.contains("iterations") ||
        !protocol_node["iterations"].is_number_unsigned()) {
      protocol_node["iterations"] = 0;
    }
    if (!protocol_node.contains("residual_norms") ||
        !protocol_node["residual_norms"].is_array()) {
      protocol_node["residual_norms"] = nlohmann::json::array();
    }
    if (!protocol_node.contains("iteration_max_residuals") ||
        !protocol_node["iteration_max_residuals"].is_array()) {
      protocol_node["iteration_max_residuals"] = nlohmann::json::array();
    }
  }

  static ExcitedExecutionResult execute_excited_state_bundle_stage(
      World &world, const PlannedStates &planned_states,
      const GroundContext &ground_ctx, const ResponseParameters &response_params,
      nlohmann::json &state_metadata_json) {
    const auto &plan = planned_states.excited_state_bundle_plan;
    ExcitedBundleSolverConfig solver_config;
    solver_config.archive_file = ground_ctx.archive_file;
    solver_config.output_prefix = response_params.prefix();
    solver_config.protocols = plan.protocols;
    solver_config.print_level = response_params.print_level();
    auto solver_adapter = make_excited_state_bundle_solver_adapter(solver_config);
    const std::string solver_adapter_name =
        solver_adapter ? solver_adapter->name() : "null_adapter";

    nlohmann::json execution = {
        {"attempted", false},
        {"mode", plan.enabled ? "legacy_adapter" : "disabled"},
        {"solver_adapter", solver_adapter_name},
        {"enabled", plan.enabled},
        {"owner_group", plan.owner_group},
        {"protocol_count", plan.protocols.size()},
        {"ready_protocol_placeholders", 0},
        {"restart_ready_protocols", 0},
        {"pending_protocols", 0},
        {"skipped_protocols", 0},
        {"completed_protocols", 0},
        {"failed_protocols", 0},
        {"total_wall_seconds", 0.0},
        {"total_cpu_seconds", 0.0},
        {"protocol_events", nlohmann::json::array()}};

    auto protocol_result_to_json = [](const ExcitedBundleProtocolResult &result) {
      return nlohmann::json{
          {"attempted", result.attempted},
          {"saved", result.saved},
          {"converged", result.converged},
          {"failed", result.failed},
          {"skipped", result.skipped},
          {"restart_reused", result.restart_reused},
          {"stage_status", result.stage_status},
          {"iterations", result.iterations},
          {"energies", result.energies},
          {"residual_norms", result.residual_norms},
          {"iteration_max_residuals", result.iteration_max_residuals}};
    };
    auto protocol_result_from_json = [](const nlohmann::json &node) {
      ExcitedBundleProtocolResult result;
      if (!node.is_object()) {
        return result;
      }
      result.attempted = node.value("attempted", false);
      result.saved = node.value("saved", false);
      result.converged = node.value("converged", false);
      result.failed = node.value("failed", false);
      result.skipped = node.value("skipped", false);
      result.restart_reused = node.value("restart_reused", false);
      result.stage_status = node.value("stage_status",
                                       std::string("placeholder_pending_solver"));
      result.iterations = node.value("iterations", static_cast<size_t>(0));
      if (node.contains("energies") && node["energies"].is_array()) {
        result.energies = node["energies"].get<std::vector<double>>();
      }
      if (node.contains("residual_norms") && node["residual_norms"].is_array()) {
        result.residual_norms = node["residual_norms"].get<std::vector<double>>();
      }
      if (node.contains("iteration_max_residuals") &&
          node["iteration_max_residuals"].is_array()) {
        result.iteration_max_residuals =
            node["iteration_max_residuals"].get<std::vector<double>>();
      }
      return result;
    };

    size_t ready_protocol_placeholders = 0;
    size_t restart_ready_protocols = 0;
    size_t pending_protocols = 0;
    size_t skipped_protocols = 0;
    size_t completed_protocols = 0;
    size_t failed_protocols = 0;
    size_t attempted_protocols = 0;
    double total_wall_seconds = 0.0;
    double total_cpu_seconds = 0.0;
    nlohmann::json protocol_events = nlohmann::json::array();

    if (!state_metadata_json.is_object()) {
      state_metadata_json = nlohmann::json::object();
    }
    if (world.rank() == 0) {
      if (!state_metadata_json.contains("excited_states") ||
          !state_metadata_json["excited_states"].is_object()) {
        state_metadata_json["excited_states"] = nlohmann::json::object();
      }
      auto &excited = state_metadata_json["excited_states"];
      excited["plan"] = plan.to_json();
      if (!excited.contains("protocols") || !excited["protocols"].is_object()) {
        excited["protocols"] = nlohmann::json::object();
      }
    }

    for (size_t protocol_index = 0; protocol_index < plan.protocols.size();
         ++protocol_index) {
      const double threshold = plan.protocols[protocol_index];
      const double protocol_wall_start = madness::wall_time();
      const double protocol_cpu_start = madness::cpu_time();

      nlohmann::json dispatch = nlohmann::json::object();
      if (world.rank() == 0) {
        const std::string protocol_key = ResponseRecord2::protocol_key(threshold);
        auto &excited = state_metadata_json["excited_states"];
        auto &node = excited["protocols"][protocol_key];
        ensure_excited_protocol_placeholder_node(node);
        ++ready_protocol_placeholders;
        dispatch["protocol_key"] = protocol_key;
        dispatch["saved"] = node["saved"].get<bool>();
        dispatch["converged"] = node["converged"].get<bool>();
        dispatch["energies"] = node["energies"];
        dispatch["iterations"] = node["iterations"];
        dispatch["residual_norms"] = node["residual_norms"];
        dispatch["iteration_max_residuals"] = node["iteration_max_residuals"];
        if (!plan.enabled) {
          dispatch["solver_needed"] = false;
          dispatch["stage_status"] = "disabled_skip";
        } else if (dispatch["saved"].get<bool>() &&
                   dispatch["converged"].get<bool>()) {
          dispatch["solver_needed"] = false;
          dispatch["stage_status"] = "restart_ready_skip";
        } else {
          dispatch["solver_needed"] = true;
          dispatch["stage_status"] = "placeholder_pending_solver";
        }
      }
      dispatch = broadcast_json_object(world, std::move(dispatch), 0);

      ExcitedBundleProtocolResult protocol_result;
      const bool solver_needed = dispatch.value("solver_needed", false);
      const bool node_saved = dispatch.value("saved", false);
      const bool node_converged = dispatch.value("converged", false);
      if (!solver_needed) {
        protocol_result.attempted = false;
        protocol_result.saved = node_saved;
        protocol_result.converged = node_converged;
        protocol_result.failed = false;
        protocol_result.skipped = true;
        protocol_result.restart_reused = node_saved && node_converged;
        protocol_result.stage_status = dispatch.value(
            "stage_status", std::string("placeholder_pending_solver"));
        protocol_result.iterations =
            dispatch.value("iterations", static_cast<size_t>(0));
        if (dispatch.contains("energies") && dispatch["energies"].is_array()) {
          protocol_result.energies =
              dispatch["energies"].get<std::vector<double>>();
        }
        if (dispatch.contains("residual_norms") &&
            dispatch["residual_norms"].is_array()) {
          protocol_result.residual_norms =
              dispatch["residual_norms"].get<std::vector<double>>();
        }
        if (dispatch.contains("iteration_max_residuals") &&
            dispatch["iteration_max_residuals"].is_array()) {
          protocol_result.iteration_max_residuals =
              dispatch["iteration_max_residuals"].get<std::vector<double>>();
        }
      } else {
        ExcitedBundleProtocolInput protocol_input;
        protocol_input.threshold = threshold;
        protocol_input.protocol_index = protocol_index;
        protocol_input.restart_saved = node_saved;
        protocol_input.restart_converged = node_converged;
        protocol_input.owner_group = plan.owner_group;
        protocol_input.tda = plan.tda;
        protocol_input.num_states = plan.num_states;
        protocol_input.guess_max_iter = plan.guess_max_iter;
        protocol_input.maxiter = plan.maxiter;
        protocol_input.maxsub = plan.maxsub;
        if (solver_adapter) {
          protocol_result = solver_adapter->solve_protocol(world, protocol_input);
        } else {
          protocol_result.attempted = true;
          protocol_result.saved = false;
          protocol_result.converged = false;
          protocol_result.failed = true;
          protocol_result.skipped = false;
          protocol_result.restart_reused = false;
          protocol_result.stage_status = "solver_adapter_missing";
        }
        if (protocol_result.stage_status.empty()) {
          protocol_result.stage_status = "placeholder_pending_solver";
        }
      }

      nlohmann::json protocol_result_json = nlohmann::json::object();
      if (world.rank() == 0) {
        protocol_result_json = protocol_result_to_json(protocol_result);
      }
      protocol_result_json =
          broadcast_json_object(world, std::move(protocol_result_json), 0);
      protocol_result = protocol_result_from_json(protocol_result_json);

      const double protocol_wall_seconds =
          madness::wall_time() - protocol_wall_start;
      const double protocol_cpu_seconds =
          madness::cpu_time() - protocol_cpu_start;

      if (world.rank() == 0) {
        const std::string protocol_key =
            dispatch.value("protocol_key", ResponseRecord2::protocol_key(threshold));
        auto &node = state_metadata_json["excited_states"]["protocols"][protocol_key];
        node["saved"] = protocol_result.saved;
        node["converged"] = protocol_result.converged;
        node["timings"]["wall_seconds"] = protocol_wall_seconds;
        node["timings"]["cpu_seconds"] = protocol_cpu_seconds;
        node["stage_status"] = protocol_result.stage_status;
        node["owner_group"] = plan.owner_group;
        node["iterations"] = protocol_result.iterations;
        if (!protocol_result.energies.empty()) {
          node["energies"] = protocol_result.energies;
        }
        if (!protocol_result.residual_norms.empty()) {
          node["residual_norms"] = protocol_result.residual_norms;
        }
        if (!protocol_result.iteration_max_residuals.empty()) {
          node["iteration_max_residuals"] =
              protocol_result.iteration_max_residuals;
        }

        total_wall_seconds += protocol_wall_seconds;
        total_cpu_seconds += protocol_cpu_seconds;
        if (protocol_result.attempted) {
          ++attempted_protocols;
        }
        if (protocol_result.restart_reused) {
          ++restart_ready_protocols;
        }
        if (protocol_result.skipped) {
          ++skipped_protocols;
        }
        if (!protocol_result.skipped && !protocol_result.failed &&
            !protocol_result.converged) {
          ++pending_protocols;
        }
        if (protocol_result.failed) {
          ++failed_protocols;
        } else {
          ++completed_protocols;
        }

        protocol_events.push_back(
            {{"protocol_key", protocol_key},
             {"threshold", threshold},
             {"saved", protocol_result.saved},
             {"converged", protocol_result.converged},
             {"attempted", protocol_result.attempted},
             {"failed", protocol_result.failed},
             {"skipped", protocol_result.skipped},
             {"restart_reused", protocol_result.restart_reused},
             {"stage_status", protocol_result.stage_status},
             {"iterations", protocol_result.iterations},
             {"energies", protocol_result.energies},
             {"residual_norms", protocol_result.residual_norms},
             {"iteration_max_residuals",
              protocol_result.iteration_max_residuals},
             {"wall_seconds", protocol_wall_seconds},
             {"cpu_seconds", protocol_cpu_seconds}});
      }
    }

    if (world.rank() == 0) {
      execution["attempted"] = attempted_protocols > 0;
      execution["ready_protocol_placeholders"] = ready_protocol_placeholders;
      execution["restart_ready_protocols"] = restart_ready_protocols;
      execution["pending_protocols"] = pending_protocols;
      execution["skipped_protocols"] = skipped_protocols;
      execution["completed_protocols"] = completed_protocols;
      execution["failed_protocols"] = failed_protocols;
      execution["total_wall_seconds"] = total_wall_seconds;
      execution["total_cpu_seconds"] = total_cpu_seconds;
      execution["protocol_events"] = std::move(protocol_events);

      if (plan.enabled || ready_protocol_placeholders > 0) {
        print("Excited-state bundle stage: protocols=",
              ready_protocol_placeholders,
              " pending=", pending_protocols,
              " restart_ready=", restart_ready_protocols,
              " skipped=", skipped_protocols, ".");
      }
    }

    state_metadata_json =
        broadcast_json_object(world, std::move(state_metadata_json), 0);
    execution = broadcast_json_object(world, std::move(execution), 0);
    return ExcitedExecutionResult{std::move(execution)};
  }

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
      const ExcitedExecutionResult &excited_result,
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
    nlohmann::json protocol_policy = nlohmann::json::array();
    for (const auto &policy : schedule_ctx.runtime_protocol_policies) {
      protocol_policy.push_back(policy.to_json());
    }
    metadata["state_parallel_runtime"]["protocol_execution_policy"] =
        std::move(protocol_policy);
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
    metadata["excited_state_planner"] = {
        {"note",
         "Stage scaffolding only: protocol-aware excited-state metadata is "
         "tracked with status/timing transitions, but no excited-state solve "
         "path is executed yet."},
        {"plan", planned_states.excited_state_bundle_plan.to_json()},
        {"execution", excited_result.execution}};
    return metadata;
  }

  // Stage-2 orchestrator:
  // 2a) build runtime ownership policy, 2b) solve linear states,
  // 2c) update excited-stage scaffold metadata,
  // 2d) execute dependency-gated derived requests,
  // 2e) publish stage metadata.
  static SolvedStates
  solve_all_states(World &world, const CalculationParameters &calc_params,
                   GroundContext &ctx,
                   const ResponseParameters &response_params,
                   PlannedStates planned_states) {
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

    // Stage 2c excited-state bundle scaffold: metadata placeholders only.
    const ExcitedExecutionResult excited_result =
        execute_excited_state_bundle_stage(world, planned_states, ctx,
                                           response_params, state_metadata_json);

    // Stage 2d derived solve: dependency-gated execution + summary metadata.
    const DerivedExecutionResult derived_result = execute_derived_state_requests(
        world, calc_params, ctx, planned_states, state_metadata_json,
        schedule_ctx.subgroup_parallel_requested, schedule_ctx.owner_group_schedule,
        final_state.threshold_index, final_state.threshold);

    // Stage 2 metadata assembly for downstream property stage.
    auto metadata = build_state_stage_metadata(
        planned_states, state_metadata_json, schedule_ctx, excited_result,
        derived_result);

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

    const auto requested_has = [&](PropertyType needle) {
      for (const auto &prop : response_params.requested_properties()) {
        if (parse_property_name(prop) == needle) {
          return true;
        }
      }
      return false;
    };

    const bool request_beta_components = requested_has(PropertyType::Beta);
    const bool request_raman_components = requested_has(PropertyType::Raman);

    if (request_beta_components || request_raman_components) {
      if (world.rank() == 0) {
        print("Property stage: distributed component precompute enabled "
              "(beta=",
              request_beta_components ? "on" : "off",
              ", raman=", request_raman_components ? "on" : "off",
              ") across ", state_parallel_plan.execution_groups, " subgroups.");
      }

      std::string component_claim_prefix;
      if (world.rank() == 0) {
        component_claim_prefix =
            (std::filesystem::path("property_component_claims") /
             ("run_" + iso_timestamp()))
                .string();
      }
      world.gop.broadcast_serializable(component_claim_prefix, 0);

      bool local_component_stage_failed = false;
      std::string local_component_stage_error;
      try {
        auto component_subworld_ptr = MacroTaskQ::create_worlds(
            world, state_parallel_plan.execution_groups);
        if (!component_subworld_ptr) {
          throw std::runtime_error(
              "property component subworld creation returned null");
        }

        World &component_subworld = *component_subworld_ptr;
        const size_t component_subgroup_id = static_cast<size_t>(
            world.rank() %
            static_cast<int>(state_parallel_plan.execution_groups));

        auto old_component_pmap3 = FunctionDefaults<3>::get_pmap();
        auto restore_component_pmap = [&]() {
          FunctionDefaults<3>::set_pmap(old_component_pmap3);
        };

        bool local_subgroup_component_failed = false;
        std::string local_subgroup_component_error;

        FunctionDefaults<3>::set_default_pmap(component_subworld);
        try {
          const auto &protocol = calc_params.protocol();
          const double final_thresh = protocol.empty()
                                          ? FunctionDefaults<3>::get_thresh()
                                          : protocol.back();
          const std::string component_fock_file =
              group_shard_file(ground_ctx.fock_json_file, component_subgroup_id);
          const std::string component_shard_file =
              group_shard_file("properties_components.json", component_subgroup_id);

          if (component_subworld.rank() == 0) {
            std::error_code ec;
            std::filesystem::remove(component_shard_file, ec);
          }
          component_subworld.gop.fence();

          GroundStateData component_ground(component_subworld,
                                           ground_ctx.archive_file,
                                           ground_ctx.molecule);
          ResponseManager component_response_manager(component_subworld,
                                                     calc_params);
          component_response_manager.setProtocol(component_subworld,
                                                 component_ground.getL(),
                                                 final_thresh);
          component_ground.prepareOrbitals(component_subworld,
                                          FunctionDefaults<3>::get_k(),
                                          final_thresh);
          component_ground.computePreliminaries(
              component_subworld, *component_response_manager.getCoulombOp(),
              component_response_manager.getVtol(), component_fock_file);

          PropertyManager component_properties(component_subworld,
                                              component_shard_file);

          if (request_beta_components) {
            ::compute_hyperpolarizability(
                component_subworld, component_ground,
                response_params.dipole_frequencies(),
                response_params.dipole_directions(), component_properties,
                component_claim_prefix + ".beta", component_subgroup_id);
          }
          if (request_raman_components) {
            ::compute_Raman_components(
                component_subworld, component_ground,
                response_params.dipole_frequencies(),
                response_params.dipole_directions(),
                response_params.nuclear_directions(), component_properties,
                component_claim_prefix + ".raman", component_subgroup_id);
          }
          component_properties.save();
        } catch (const std::exception &ex) {
          local_subgroup_component_failed = true;
          local_subgroup_component_error = ex.what();
        } catch (...) {
          local_subgroup_component_failed = true;
          local_subgroup_component_error =
              "unknown exception during component precompute";
        }

        long subgroup_component_failed_flag =
            local_subgroup_component_failed ? 1 : 0;
        component_subworld.gop.max(subgroup_component_failed_flag);
        restore_component_pmap();

        if (subgroup_component_failed_flag != 0) {
          local_component_stage_failed = true;
          if (local_subgroup_component_failed &&
              !local_subgroup_component_error.empty()) {
            local_component_stage_error =
                std::move(local_subgroup_component_error);
          }
        }
      } catch (const std::exception &ex) {
        local_component_stage_failed = true;
        local_component_stage_error = ex.what();
      } catch (...) {
        local_component_stage_failed = true;
        local_component_stage_error =
            "unknown exception during component stage setup";
      }

      long any_component_failure = local_component_stage_failed ? 1 : 0;
      world.gop.max(any_component_failure);
      if (any_component_failure != 0) {
        if (world.rank() == 0) {
          if (!local_component_stage_error.empty()) {
            print("State-parallel component stage failed: ",
                  local_component_stage_error,
                  " Falling back to world property stage.");
          } else {
            print("State-parallel component stage failed on at least one "
                  "rank. Falling back to world property stage.");
          }
        }
        return compute_requested_properties(world, response_params, ground_ctx,
                                            solved_states, scf_calc);
      }

      if (world.rank() == 0) {
        nlohmann::json merged_components = nlohmann::json::array();
        const auto baseline = read_json_file_or_object("properties.json");
        if (baseline.is_array()) {
          for (const auto &row : baseline) {
            merged_components.push_back(row);
          }
        }
        for (size_t gid = 0; gid < state_parallel_plan.execution_groups; ++gid) {
          const auto shard = read_json_file_or_object(
              group_shard_file("properties_components.json", gid));
          if (!shard.is_array()) {
            continue;
          }
          for (const auto &row : shard) {
            merged_components.push_back(row);
          }
        }
        write_json_file("properties.json", merged_components);
      }
      world.gop.fence();
    }

    if (world.rank() == 0) {
      print("Property stage: executing downstream assembly on subgroup ",
            property_group, "/", state_parallel_plan.execution_groups - 1);
    }

    std::string properties_dump;
    std::string vib_dump;
    std::string raman_dump;
    bool local_property_stage_failed = false;
    std::string local_property_stage_error;

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

      bool local_subgroup_failed = false;
      std::string local_subgroup_error;

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
      } catch (const std::exception &ex) {
        local_subgroup_failed = true;
        local_subgroup_error = ex.what();
      } catch (...) {
        local_subgroup_failed = true;
        local_subgroup_error =
            "unknown exception during subgroup property execution";
      }

      long subgroup_failed_flag = local_subgroup_failed ? 1 : 0;
      subworld.gop.max(subgroup_failed_flag);
      restore_pmap();

      if (subgroup_failed_flag != 0) {
        local_property_stage_failed = true;
        if (local_subgroup_failed && !local_subgroup_error.empty()) {
          local_property_stage_error = std::move(local_subgroup_error);
        }
      }
    } catch (const std::exception &ex) {
      local_property_stage_failed = true;
      local_property_stage_error = ex.what();
    } catch (...) {
      local_property_stage_failed = true;
      local_property_stage_error =
          "unknown exception during property subgroup setup";
    }

    long any_property_failure = local_property_stage_failed ? 1 : 0;
    world.gop.max(any_property_failure);
    if (any_property_failure != 0) {
      if (world.rank() == 0) {
        if (!local_property_stage_error.empty()) {
          print("State-parallel property subgroup execution failed: ",
                local_property_stage_error,
                " Falling back to world property stage.");
        } else {
          print("State-parallel property subgroup execution failed on at least "
                "one rank. Falling back to world property stage.");
        }
      }
      return compute_requested_properties(world, response_params, ground_ctx,
                                          solved_states, scf_calc);
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
  }

public:
  /**
   * @brief Run the full molecular response & property workflow.
   *
   * 1. Builds contexts and plans states
   * 2. Solves all linear and derived states
   * 3. Computes requested properties
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
