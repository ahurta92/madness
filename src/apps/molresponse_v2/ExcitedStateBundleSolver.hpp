#pragma once
// Backward-compat shim. Types have moved to ResponseRecord.hpp.
// This file will be removed once ExcitedResponse replaces ExcitedStateBundleSolver.
#include "ResponseRecord.hpp"

#include <madness/world/world.h>
#include <memory>
#include <string>
#include <vector>

// Legacy aliases — keep the old names compiling during the transition.
using ExcitedBundleProtocolInput  = ExcitedProtocolInput;
using ExcitedBundleProtocolResult = ExcitedProtocolResult;

struct ExcitedBundleSolverConfig {
    std::string          archive_file;
    std::string          output_prefix = "excited_bundle";
    std::vector<double>  protocols;
    int                  print_level   = 0;
};

class ExcitedStateBundleSolver {
public:
    virtual ~ExcitedStateBundleSolver() = default;
    [[nodiscard]] virtual std::string name() const = 0;
    [[nodiscard]] virtual ExcitedProtocolResult
    solve_protocol(madness::World &world,
                   const ExcitedProtocolInput &input) const = 0;
};

class ExcitedStateBundleNoopSolver final : public ExcitedStateBundleSolver {
public:
    [[nodiscard]] std::string name() const override {
        return "noop_scaffold_adapter";
    }
    [[nodiscard]] ExcitedProtocolResult
    solve_protocol(madness::World &,
                   const ExcitedProtocolInput &input) const override {
        ExcitedProtocolResult result;
        result.attempted      = true;
        result.saved          = input.restart_saved;
        result.converged      = input.restart_converged;
        result.stage_status   = "placeholder_pending_solver";
        return result;
    }
};

[[nodiscard]] std::unique_ptr<ExcitedStateBundleSolver>
make_excited_state_bundle_solver_adapter(const ExcitedBundleSolverConfig &config);
