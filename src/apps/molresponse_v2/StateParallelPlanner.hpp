#pragma once

#include "ResponseState.hpp"
#include "../../madness/chem/ResponseParameters.hpp"

#include <algorithm>
#include <string>
#include <vector>

#include <madness/external/nlohmann_json/json.hpp>

/*
Developer Overview
- Purpose: Convert user knobs (off/auto/on + group count) into a deterministic
  perturbation-channel scheduling plan.
- Terminology used in this file:
  - perturbation channel: one perturbation descriptor
  - frequency series: ordered frequencies attached to one channel
  - channel point: one (channel, frequency, protocol) work item
- Strategy: Validate requested group configuration, decide effective mode, then
  assign each generated linear response channel to an owner group via
  round-robin.
- Protocol policy: In parallel mode, begin with state ownership and optionally
  switch to independent channel-point ownership at
  response.state_parallel_point_start_protocol.
- Output contract: Return a pure metadata plan (mapping + execution settings +
  per-channel ownership) consumed by MolresponseLib solve stages.
- Design note: Mapping and execution group counts are kept explicit so we can
  preserve deterministic ownership even when execution falls back to serial.
*/

using json = nlohmann::json;

struct PerturbationChannelAssignment {
  // Index into GeneratedStateData::states (one entry per perturbation channel).
  size_t channel_index = 0;
  // Human-readable channel key (for logs/JSON only).
  std::string channel_label;
  // Deterministic owner lane for this channel.
  size_t owner_group = 0;

  [[nodiscard]] json to_json() const {
    return {{"channel_index", channel_index},
            {"channel_label", channel_label},
            // Legacy key kept for downstream tooling compatibility.
            {"perturbation", channel_label},
            {"owner_group", owner_group}};
  }
};

struct StateParallelPlan {
  // Raw user knob from response.state_parallel.
  std::string requested_mode = "off";
  // Planner-selected mode used by the solver.
  std::string effective_mode = "serial";
  // Raw user knob from response.state_parallel_groups.
  size_t requested_groups = 1;
  // Number of ownership lanes used for deterministic work partitioning.
  size_t mapping_groups = 1;
  // Number of groups that participate in channel-series ownership
  // (protocol-0 warmup). Capped by number of channels; remaining groups can
  // join in channel-point ownership.
  size_t channel_owner_groups = 1;
  // Number of subworlds created for execution.
  size_t execution_groups = 1;
  // First protocol index that may switch from channel-series ownership to
  // channel-point ownership.
  size_t point_parallel_start_protocol_index = 1;
  // Human-readable description of the protocol ownership policy.
  std::string frequency_partition_policy = "channel_series_then_point";
  // Auto-mode threshold: minimum state count before enabling parallel mapping.
  size_t min_states = 1;
  // Number of generated perturbation channels considered by the planner.
  size_t num_channels = 0;
  // Number of (state,frequency) points across all linear states.
  size_t num_points = 0;
  // Size of the universe communicator.
  size_t world_size = 1;
  // Point-ownership lane count after capping by available points.
  size_t effective_point_groups = 1;
  // True when the solve stage is allowed to run via planner output.
  bool execution_enabled = false;
  // True when execution uses macrotask subworlds instead of single-world serial.
  bool subgroup_parallel_enabled = false;
  // Diagnostic reason describing why a mode was chosen.
  std::string reason;
  // Per-channel ownership map used by solve_all_states.
  std::vector<PerturbationChannelAssignment> channel_assignments;

  [[nodiscard]] bool
  use_channel_series_ownership_for_protocol(size_t protocol_index) const {
    if (mapping_groups <= 1) {
      return true;
    }
    return protocol_index < point_parallel_start_protocol_index;
  }

  // Runtime override used by the solver: on restart, protocol-0 can promote to
  // channel-point ownership when protocol-0 points are already saved.
  [[nodiscard]] size_t effective_point_parallel_start_protocol_index(
      bool owner_group_schedule, bool has_protocol_thresholds,
      bool restart_protocol0_saved_complete) const {
    if (mapping_groups <= 1) {
      return 0;
    }
    if (owner_group_schedule && has_protocol_thresholds &&
        point_parallel_start_protocol_index == 1 &&
        restart_protocol0_saved_complete) {
      return 0;
    }
    return point_parallel_start_protocol_index;
  }

  [[nodiscard]] json to_json() const {
    json rows = json::array();
    for (const auto &a : channel_assignments)
      rows.push_back(a.to_json());
    return {{"requested_mode", requested_mode},
            {"effective_mode", effective_mode},
            {"requested_groups", requested_groups},
            {"mapping_groups", mapping_groups},
            {"channel_owner_groups", channel_owner_groups},
            // Legacy key kept for downstream tooling compatibility.
            {"state_owner_groups", channel_owner_groups},
            {"execution_groups", execution_groups},
            {"point_parallel_start_protocol_index",
             point_parallel_start_protocol_index},
            {"frequency_partition_policy", frequency_partition_policy},
            {"min_states", min_states},
            {"num_channels", num_channels},
            // Legacy key kept for downstream tooling compatibility.
            {"num_states", num_channels},
            {"num_points", num_points},
            {"world_size", world_size},
            {"effective_point_groups", effective_point_groups},
            {"execution_enabled", execution_enabled},
            {"subgroup_parallel_enabled", subgroup_parallel_enabled},
            {"reason", reason},
            {"channel_assignments", rows},
            // Legacy key kept for downstream tooling compatibility.
            {"assignments", rows}};
  }
};

class PointOwnershipScheduler {
public:
  PointOwnershipScheduler(const std::vector<LinearResponseDescriptor> &states,
                          size_t owner_groups)
      : owner_groups_(std::max<size_t>(1, owner_groups)),
        state_point_offsets_(states.size() + 1, 0) {
    for (size_t i = 0; i < states.size(); ++i) {
      state_point_offsets_[i + 1] =
          state_point_offsets_[i] + states[i].num_frequencies();
    }
  }

  [[nodiscard]] size_t owner_groups() const { return owner_groups_; }

  [[nodiscard]] size_t num_points() const {
    return state_point_offsets_.empty() ? 0 : state_point_offsets_.back();
  }

  [[nodiscard]] size_t owner_group(size_t channel_index,
                                   size_t freq_index) const {
    if (owner_groups_ <= 1) {
      return 0;
    }
    const size_t linear_point_index =
        state_point_offsets_[channel_index] + freq_index;
    return linear_point_index % owner_groups_;
  }

private:
  size_t owner_groups_ = 1;
  std::vector<size_t> state_point_offsets_;
};

class StateParallelPlanner {
public:
  static StateParallelPlan
  build(const ResponseParameters &params, size_t world_size,
        const std::vector<LinearResponseDescriptor> &states) {
    StateParallelPlan plan;
    plan.requested_mode = params.state_parallel();
    plan.requested_groups = std::max<size_t>(1, params.state_parallel_groups());
    plan.min_states = params.state_parallel_min_states();
    plan.num_channels = states.size();
    for (const auto &state : states) {
      plan.num_points += state.num_frequencies();
    }
    plan.world_size = world_size;

    const bool requested_on = (plan.requested_mode == "on");
    const bool requested_auto = (plan.requested_mode == "auto");
    const bool enough_channels = (plan.num_channels >= plan.min_states);
    const bool groups_valid = (plan.requested_groups > 1) &&
                              (plan.requested_groups <= plan.world_size);

    bool should_plan_parallel_mapping = false;

    if (plan.requested_mode == "off") {
      plan.reason = "state_parallel=off";
    } else if (!groups_valid) {
      plan.reason =
          "invalid group request (must be >1 and <= world size); using serial";
    } else if (requested_auto && !enough_channels) {
      plan.reason = "auto mode disabled: channel count below min_states";
    } else if (requested_on || (requested_auto && enough_channels)) {
      should_plan_parallel_mapping = true;
      plan.reason = "ownership mapping planned with subgroup execution";
    }

    if (should_plan_parallel_mapping) {
      plan.mapping_groups = plan.requested_groups;
      plan.channel_owner_groups =
          std::min(plan.mapping_groups, std::max<size_t>(1, plan.num_channels));
      plan.execution_groups = plan.requested_groups;
      plan.execution_enabled = true;
      plan.subgroup_parallel_enabled = true;
      plan.effective_mode = "owner_group_subworld";
      plan.point_parallel_start_protocol_index =
          params.state_parallel_point_start_protocol();
      plan.frequency_partition_policy =
          plan.point_parallel_start_protocol_index == 0
              ? "channel_point_from_t0"
              : "channel_series_then_point";
      plan.effective_point_groups =
          plan.num_points == 0
              ? 1
              : std::min(plan.mapping_groups,
                         std::max<size_t>(1, plan.num_points));
      if (plan.effective_point_groups < plan.mapping_groups) {
        plan.reason +=
            "; channel-point ownership lanes capped by available points";
      }
      if (plan.channel_owner_groups < plan.mapping_groups) {
        plan.reason +=
            "; protocol-0 warmup uses a subset of groups until channel-point ownership";
      }
    } else {
      plan.mapping_groups = 1;
      plan.channel_owner_groups = 1;
      plan.execution_groups = 1;
      plan.execution_enabled = true;
      plan.subgroup_parallel_enabled = false;
      plan.effective_mode = "serial";
      plan.point_parallel_start_protocol_index = 0;
      plan.frequency_partition_policy = "channel_series_only";
      plan.effective_point_groups = 1;
    }

    plan.channel_assignments.reserve(states.size());
    for (size_t i = 0; i < states.size(); ++i) {
      const auto owner = i % plan.channel_owner_groups;
      plan.channel_assignments.push_back(
          PerturbationChannelAssignment{i,
                                        states[i].perturbationDescription(),
                                        owner});
    }
    return plan;
  }
};
