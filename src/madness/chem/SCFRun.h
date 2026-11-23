#pragma once
#include <madness/chem/SCFInput.h>
#include <madness/chem/IterationMetrics.h>
#include <nlohmann/json.hpp>

namespace madness {
class SCF;
std::filesystem::path scf_checkpoint_directory(const SCF& scf);
/// Shared electronic execution for workflow and standalone SCF callers.
double run_scf(World& world, SCF& scf);
/// Inspect saved evidence and optionally plot its native functions. Never solve or publish.
void inspect_scf(World& world, const std::filesystem::path& root,
                 const CalculationParameters& plotting);
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
void require_scf_restart(World& world, const std::filesystem::path& root, bool frozen);
}
