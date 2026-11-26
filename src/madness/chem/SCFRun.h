#pragma once
#include <madness/chem/SCFInput.h>
#include <madness/chem/IterationMetrics.h>
#include <nlohmann/json.hpp>

namespace madness {
class SCF;
std::filesystem::path scf_checkpoint_directory(const SCF& scf);
/// Shared electronic execution for workflow and standalone SCF callers.
double run_scf(World& world, SCF& scf);
/// Not available in this version: it returns as a loop over a run's node descriptions.
void inspect_scf(World& world, const std::filesystem::path& root,
                 const CalculationParameters& plotting);
/// The orbital archive of a run's saved ground-state (or virtuals) node.
std::filesystem::path scf_checkpoint_archive(World& world, const std::filesystem::path& root,
                                             bool frozen);
std::filesystem::path scf_checkpoint_archive(World& world, const SCF& scf);
/// A completed ground state as later calculations see it: where its orbital
/// archive is, the level it was saved at, and an identity (its settings, that
/// level and its save time) that changes whenever it is saved again.
struct GroundReference {
    std::filesystem::path archive;
    iteration::Protocol protocol;
    nlohmann::json identity;
};
GroundReference scf_ground_reference(World& world, const SCF& scf);
/// Collective. Add the dipole, the gradient and the populations, each when it
/// was computed (not empty), to the ground-state entry of results.json that
/// run_scf wrote.
void report_scf_properties(World& world, const SCF& scf, const Tensor<double>& dipole,
                           const Tensor<double>& gradient, const nlohmann::json& populations);
/// Collective: whether the run holds a saved ground-state node.
bool has_saved_ground_state(World& world, const std::filesystem::path& root);
void require_scf_restart(World& world, const std::filesystem::path& root);
}
