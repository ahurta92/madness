#pragma once

#include <madness/world/world.h>

#include <string>
#include <vector>

struct ExcitedBundleProtocolInput {
  double threshold = 0.0;
  size_t protocol_index = 0;
  bool restart_saved = false;
  bool restart_converged = false;
  size_t owner_group = 0;
  bool tda = false;
  size_t num_states = 1;
  size_t guess_max_iter = 5;
  size_t maxiter = 20;
  size_t maxsub = 8;
};

struct ExcitedBundleProtocolResult {
  bool attempted = false;
  bool saved = false;
  bool converged = false;
  bool failed = false;
  bool skipped = false;
  bool restart_reused = false;
  std::string stage_status = "placeholder_pending_solver";
  std::vector<double> energies;
};

class ExcitedStateBundleSolver {
public:
  virtual ~ExcitedStateBundleSolver() = default;
  [[nodiscard]] virtual std::string name() const = 0;
  [[nodiscard]] virtual ExcitedBundleProtocolResult
  solve_protocol(madness::World &world,
                 const ExcitedBundleProtocolInput &input) const = 0;
};

// Phase-4 adapter seam default: behavior-preserving no-op implementation.
class ExcitedStateBundleNoopSolver final : public ExcitedStateBundleSolver {
public:
  [[nodiscard]] std::string name() const override {
    return "noop_scaffold_adapter";
  }

  [[nodiscard]] ExcitedBundleProtocolResult
  solve_protocol(madness::World &,
                 const ExcitedBundleProtocolInput &input) const override {
    ExcitedBundleProtocolResult result;
    result.attempted = true;
    result.saved = input.restart_saved;
    result.converged = input.restart_converged;
    result.failed = false;
    result.skipped = false;
    result.restart_reused = false;
    result.stage_status = "placeholder_pending_solver";
    return result;
  }
};
