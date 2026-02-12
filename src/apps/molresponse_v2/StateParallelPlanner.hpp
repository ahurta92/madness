#pragma once

#include "ResponseState.hpp"
#include "../../madness/chem/ResponseParameters.hpp"

#include <algorithm>
#include <string>
#include <vector>

#include <madness/external/nlohmann_json/json.hpp>

using json = nlohmann::json;

struct StateAssignment {
  size_t state_index = 0;
  std::string perturbation;
  size_t owner_group = 0;

  [[nodiscard]] json to_json() const {
    return {{"state_index", state_index},
            {"perturbation", perturbation},
            {"owner_group", owner_group}};
  }
};

struct StateParallelPlan {
  std::string requested_mode = "off";
  std::string effective_mode = "serial";
  size_t requested_groups = 1;
  size_t mapping_groups = 1;
  size_t execution_groups = 1;
  size_t min_states = 1;
  size_t num_states = 0;
  size_t world_size = 1;
  bool execution_enabled = false;
  std::string reason;
  std::vector<StateAssignment> assignments;

  [[nodiscard]] json to_json() const {
    json rows = json::array();
    for (const auto &a : assignments)
      rows.push_back(a.to_json());
    return {{"requested_mode", requested_mode},
            {"effective_mode", effective_mode},
            {"requested_groups", requested_groups},
            {"mapping_groups", mapping_groups},
            {"execution_groups", execution_groups},
            {"min_states", min_states},
            {"num_states", num_states},
            {"world_size", world_size},
            {"execution_enabled", execution_enabled},
            {"reason", reason},
            {"assignments", rows}};
  }
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
    plan.num_states = states.size();
    plan.world_size = world_size;

    const bool requested_on = (plan.requested_mode == "on");
    const bool requested_auto = (plan.requested_mode == "auto");
    const bool enough_states = (plan.num_states >= plan.min_states);
    const bool groups_valid = (plan.requested_groups > 1) &&
                              (plan.requested_groups <= plan.world_size);

    bool should_plan_parallel_mapping = false;

    if (plan.requested_mode == "off") {
      plan.reason = "state_parallel=off";
    } else if (!groups_valid) {
      plan.reason =
          "invalid group request (must be >1 and <= world size); using serial";
    } else if (requested_auto && !enough_states) {
      plan.reason = "auto mode disabled: state count below min_states";
    } else if (requested_on || (requested_auto && enough_states)) {
      should_plan_parallel_mapping = true;
      plan.reason = "ownership mapping planned; subgroup execution not enabled yet";
    }

    if (should_plan_parallel_mapping) {
      plan.mapping_groups = plan.requested_groups;
      plan.execution_groups = 1;
      plan.execution_enabled = false;
      plan.effective_mode = "planned_serial";
    } else {
      plan.mapping_groups = 1;
      plan.execution_groups = 1;
      plan.execution_enabled = false;
      plan.effective_mode = "serial";
    }

    plan.assignments.reserve(states.size());
    for (size_t i = 0; i < states.size(); ++i) {
      const auto owner = i % plan.mapping_groups;
      plan.assignments.push_back(
          StateAssignment{i, states[i].perturbationDescription(), owner});
    }
    return plan;
  }
};
