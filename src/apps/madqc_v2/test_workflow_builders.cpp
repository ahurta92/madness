#include <WorkflowBuilders.hpp>
#include <apps/molresponse_v2/StateParallelPlanner.hpp>
#include <madness/chem/ResponseParameters.hpp>

#include <array>
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace {

LinearResponseDescriptor make_dipole_state(char axis) {
  return LinearResponseDescriptor(DipolePerturbation{axis}, {0.0, 0.1},
                                  {1.0e-3, 1.0e-4}, true);
}

} // namespace

int main() {
  using madness::workflow_builders::WorkflowKind;
  using madness::workflow_builders::workflow_kind_from_name;

  struct Case {
    const char *name;
    WorkflowKind kind;
  };

  constexpr std::array<Case, 8> expected = {{
      {"scf", WorkflowKind::Scf},
      {"nemo", WorkflowKind::Nemo},
      {"response", WorkflowKind::Response},
      {"mp2", WorkflowKind::Mp2Cc2},
      {"cc2", WorkflowKind::Mp2Cc2},
      {"cis", WorkflowKind::Cis},
      {"oep", WorkflowKind::Oep},
      {"optimize", WorkflowKind::Optimize},
  }};

  bool ok = true;
  for (const auto &entry : expected) {
    const auto actual = workflow_kind_from_name(entry.name);
    if (actual != entry.kind) {
      ok = false;
      std::cerr << "workflow_kind_from_name mismatch for '" << entry.name
                << "'\n";
    }
  }

  if (workflow_kind_from_name("unknown_workflow") != WorkflowKind::Unknown) {
    ok = false;
    std::cerr << "workflow_kind_from_name should return Unknown for invalid "
                 "workflow\n";
  }

  const std::string workflow_list =
      madness::workflow_builders::runnable_workflow_list();
  for (const char *name : madness::workflow_builders::runnable_workflows) {
    if (workflow_list.find(name) == std::string::npos) {
      ok = false;
      std::cerr << "runnable_workflow_list missing '" << name << "'\n";
    }
  }

  {
    ResponseParameters params;
    params.set_user_defined_value<std::string>("state_parallel", "on");
    params.set_user_defined_value<size_t>("state_parallel_groups", 8);
    params.set_user_defined_value<size_t>("state_parallel_min_states", 1);
    params.set_user_defined_value<size_t>("state_parallel_point_start_protocol",
                                          1);

    const std::vector<LinearResponseDescriptor> states = {
        make_dipole_state('x'),
        make_dipole_state('y'),
        make_dipole_state('z'),
    };
    const auto plan = StateParallelPlanner::build(params, 32, states);

    if (plan.point_parallel_start_protocol_index != 1) {
      ok = false;
      std::cerr << "state_parallel_point_start_protocol=1 not reflected in "
                   "planner output\n";
    }
    if (plan.frequency_partition_policy != "state_then_point") {
      ok = false;
      std::cerr << "expected state_then_point frequency policy for "
                   "point_start=1\n";
    }
    if (plan.num_points != 6) {
      ok = false;
      std::cerr << "expected 6 linear response points for 3 states x 2 "
                   "frequencies\n";
    }
    if (plan.effective_point_groups != 6) {
      ok = false;
      std::cerr << "expected effective point groups to cap at number of points "
                   "(6)\n";
    }
    if (!plan.use_state_ownership_for_protocol(0) ||
        plan.use_state_ownership_for_protocol(1)) {
      ok = false;
      std::cerr << "planner protocol ownership boundary incorrect for "
                   "point_start=1\n";
    }

    PointOwnershipScheduler point_scheduler(states, plan.effective_point_groups);
    std::set<size_t> owners;
    for (size_t state_index = 0; state_index < states.size(); ++state_index) {
      for (size_t freq_index = 0;
           freq_index < states[state_index].num_frequencies(); ++freq_index) {
        owners.insert(point_scheduler.owner_group(state_index, freq_index));
      }
    }
    if (owners.size() != 6) {
      ok = false;
      std::cerr << "expected all 6 point-owner lanes to be used\n";
    }

    const size_t restart_effective =
        plan.effective_point_parallel_start_protocol_index(
            /*owner_group_schedule=*/true,
            /*has_protocol_thresholds=*/true,
            /*restart_protocol0_saved_complete=*/true);
    if (restart_effective != 0) {
      ok = false;
      std::cerr << "restart-ready runtime policy should promote point "
                   "ownership to protocol 0\n";
    }

    const size_t fresh_effective =
        plan.effective_point_parallel_start_protocol_index(
            /*owner_group_schedule=*/true,
            /*has_protocol_thresholds=*/true,
            /*restart_protocol0_saved_complete=*/false);
    if (fresh_effective != 1) {
      ok = false;
      std::cerr << "fresh-run runtime policy should preserve configured "
                   "point-start protocol\n";
    }
  }

  {
    ResponseParameters params;
    params.set_user_defined_value<std::string>("state_parallel", "on");
    params.set_user_defined_value<size_t>("state_parallel_groups", 8);
    params.set_user_defined_value<size_t>("state_parallel_min_states", 1);
    params.set_user_defined_value<size_t>("state_parallel_point_start_protocol",
                                          0);

    const std::vector<LinearResponseDescriptor> states = {
        make_dipole_state('x'),
        make_dipole_state('y'),
    };
    const auto plan = StateParallelPlanner::build(params, 32, states);

    if (plan.point_parallel_start_protocol_index != 0) {
      ok = false;
      std::cerr << "state_parallel_point_start_protocol=0 not reflected in "
                   "planner output\n";
    }
    if (plan.frequency_partition_policy != "point_from_t0") {
      ok = false;
      std::cerr << "expected point_from_t0 frequency policy for point_start=0\n";
    }
    if (plan.num_points != 4) {
      ok = false;
      std::cerr << "expected 4 linear response points for 2 states x 2 "
                   "frequencies\n";
    }
    if (plan.effective_point_groups != 4) {
      ok = false;
      std::cerr << "expected point group cap to 4 points for point_start=0 "
                   "scenario\n";
    }
    if (plan.use_state_ownership_for_protocol(0)) {
      ok = false;
      std::cerr << "point_start=0 should disable state ownership at protocol 0\n";
    }
  }

  {
    ResponseParameters params;
    params.set_user_defined_value<std::string>("state_parallel", "on");
    params.set_user_defined_value<size_t>("state_parallel_groups", 8);
    params.set_user_defined_value<size_t>("state_parallel_min_states", 1);
    params.set_user_defined_value<size_t>("state_parallel_point_start_protocol",
                                          1);

    const std::vector<LinearResponseDescriptor> states = {
        make_dipole_state('x'),
    };
    const auto plan = StateParallelPlanner::build(params, 32, states);
    if (plan.effective_point_groups != 2) {
      ok = false;
      std::cerr << "single-state point-group cap should be 2 for two "
                   "frequencies\n";
    }
  }

  return ok ? 0 : 1;
}
