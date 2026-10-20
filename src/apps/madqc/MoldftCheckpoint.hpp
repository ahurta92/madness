#pragma once

#include "NamedStages.hpp"
#include <filesystem>
#include <optional>
#include <madness/chem/StageCheckpoint.h>

namespace named_stages {
namespace fs = std::filesystem;
using json = nlohmann::json;

// Application policy lives here; StageCheckpoint only checks physical storage.
class MoldftCheckpoint {
  madness::World &world_;
  const Stage &stage_;
  fs::path store_;
  std::string run_id_;
  std::size_t sequence_ = 0;
  json physical_;
  madness::checkpoint::Manifest latest_;
  bool has_latest_ = false;
  std::optional<madness::SCF::IterationCheckpoint> comparison_;

  static void require(bool condition, const std::string &message) {
    if (!condition)
      throw std::runtime_error("moldft checkpoint: " + message);
  }

  static void validate_measurements(const json& metrics) {
    if (!metrics.contains("convergence_report")) return;
    const auto& report = metrics.at("convergence_report");
    auto validate = [&](const json& measurements) {
      for (const auto& metric : measurements.items())
        require(metric.value().value("availability", "") != "invalid",
                "invalid ground-state measurement: " + metric.key());
    };
    validate(report.at("measurements"));
    if (report.contains("saved_state_validation"))
      validate(report.at("saved_state_validation").at("measurements"));
  }

  static json physical(const madness::SCF &scf) {
    auto molecule = scf.molecule.to_json();
    // File provenance changes each attempt and is not molecular identity.
    molecule["parameters"].erase("source_name");
    molecule["parameters"].erase("source_type");
    return {{"schema", "madqc.moldft.v1"},
            {"molecule", molecule},
            {"xc", scf.param.xc()},
            {"charge", scf.param.charge()},
            {"spin_restricted", scf.param.spin_restricted()},
            {"nalpha", scf.param.nalpha()},
            {"nbeta", scf.param.nbeta()},
            {"localize", scf.param.localize_method()},
            {"L", scf.param.L()},
            {"representation", "mo"}};
  }

public:
  std::size_t iterations = 0;
  std::string work = "solved";

  MoldftCheckpoint(madness::World &world, const Stage &stage,
                   const fs::path &directory, std::string run_id)
      : world_(world), stage_(stage), store_(directory / "checkpoints"),
        run_id_(std::move(run_id)) {}

  // Validate the complete selected generation before the engine opens an
  // archive. A private copy lets the existing restart planner use its ordinary
  // prefix without ever changing an immutable generation or using an AO
  // sidecar.
  void prepare(madness::SCF &scf, bool continuing) {
    physical_ = physical(scf);
    if (!continuing)
      return;
    latest_ = madness::checkpoint::read_published(world_, store_);
    require(latest_.stage == stage_.name && latest_.method == "HF",
            "stage/method mismatch");
    require(latest_.physical == physical_,
            "geometry, model or cell changed; use a new run");
    require(latest_.archives.size() == 1 &&
                latest_.archives[0].role == "orbitals" &&
                latest_.archives[0].basename == "state.restartdata",
            "required orbital inventory is missing");
    require(latest_.numerical.at("representation") == "mo",
            "unsupported representation");
    const auto &cp = scf.param;
    const bool satisfied =
        latest_.stop_reason == "converged" &&
        latest_.numerical.at("threshold").get<double>() <=
            cp.protocol().back() &&
        latest_.numerical.at("k").get<int>() ==
            (cp.k() == -1 ? protocol_k(cp.protocol().back()) : cp.k()) &&
        latest_.targets.at("effective_dconv").get<double>() <=
            std::max(cp.dconv(), cp.protocol().back()) &&
        latest_.targets.at("econv").get<double>() <= cp.econv();
    work = satisfied ? "loaded" : "resumed";
    madness::checkpoint::detail::on_root(world_, [&] {
      validate_measurements(latest_.metrics);
      const auto directory = store_ / latest_.generation;
      for (const auto &shard : latest_.archives[0].shards) {
        const auto suffix = shard.path.substr(std::string("state").size());
        fs::copy_file(directory / shard.path, cp.prefix() + suffix,
                      fs::copy_options::none);
      }
      return std::string{};
    });
    // Explicit iteration cannot silently fall back to an atomic guess.
    scf.param.set_user_defined_value(
        "restart", satisfied ? std::string("auto") : std::string("iterate"));
    has_latest_ = true;
  }

  const madness::checkpoint::Manifest &latest() const { return latest_; }
  fs::path prefix() const { return store_ / latest_.generation / "state"; }

  void observe(const madness::SCF &scf,
               const madness::SCF::IterationCheckpoint &step) {
    ++iterations;
    comparison_ = step;
    json metrics = {
        {"energy",
         {{"value", step.input_energy},
          {"units", "hartree"},
          {"phase", "iteration_input_before_rotation_and_update"}}}};
    if (step.has_previous) {
      metrics["density_change"] = {
          {"alpha", step.density_alpha},
          {"beta", step.density_beta},
          {"phase", "successive_iteration_input_densities"}};
      metrics["bsh_residual"] = {{"value", step.previous_bsh_residual},
                                 {"phase", "previous_iteration_update"}};
    }
    // Measurements describe iteration evaluations, not the archived output
    // after update/rotation. Keep that distinction even in a completed report.
    metrics["kinetic_energy"] = {{"value", step.input_kinetic}, {"units", "hartree"},
        {"phase", "iteration_kinetic_evaluation_before_update"}};
    const double effective = std::max(scf.param.dconv(), step.threshold);
    auto measured = [&](double value, const std::string& units,
                        const std::string& phase, bool active, json target,
                        json limit, int sample_iteration) {
      const bool available = step.has_previous;
      const bool finite = std::isfinite(value);
      return json{{"value", available && finite ? json(value) : json(nullptr)},
                  {"availability", !available ? "unavailable" : finite ? "measured" : "invalid"},
                  {"reason", available ? "" : "no_comparable_previous_sample"},
                  {"units", units}, {"phase", phase}, {"active_gate", active},
                  {"target", target}, {"effective_limit", limit},
                  {"satisfied", available && finite && !limit.is_null()
                      ? json(value < limit.get<double>()) : json(nullptr)},
                  {"valid_for_saved_state", false},
                  {"sample", {{"run_id", run_id_}, {"stage", stage_.name},
                              {"threshold", step.threshold}, {"k", scf.amo.at(0).k()},
                              {"iteration", sample_iteration}}}};
    };
    json evidence;
    evidence["density_change"] = measured(std::max(step.density_alpha, step.density_beta),
        "a0^(-3/2)", "successive_iteration_input_densities", true, scf.param.dconv(),
        effective * std::max(std::size_t(5), scf.molecule.natom()), step.iteration);
    evidence["density_change"]["reduction"] = "maximum_over_spin_density_L2_norms";
    evidence["density_change"]["components"] = step.has_previous
        ? json{{"alpha", step.density_alpha}, {"beta", step.density_beta}} : json(nullptr);
    evidence["bsh_residual"] = measured(step.previous_bsh_residual, "dimensionless",
        "raw_fixed_point_difference_before_previous_update", !scf.param.get<bool>("conv_only_dens"),
        nullptr, 5.0 * effective, step.iteration - 1);
    evidence["bsh_residual"]["reduction"] = "maximum_over_orbital_and_spin_L2_norms";
    evidence["energy_change"] = measured(step.energy_change, "hartree",
        "successive_iteration_energy_evaluations", false, scf.param.econv(), scf.param.econv(), step.iteration);
    evidence["delta_kinetic_energy"] = measured(step.kinetic_change, "hartree",
        "successive_iteration_kinetic_evaluations", false, nullptr, nullptr, step.iteration);
    evidence["delta_kinetic_energy"]["reduction"] = "absolute_change_of_occupation_weighted_spin_sum";
    metrics["convergence_report"] = {{"schema", "madqc.convergence.v1"},
        {"gate_policy", "existing_engine"}, {"measurements", evidence}};
    publish(scf, false, json::array(), metrics,
            {{"iteration", step.iteration},
             {"iterations_this_run", iterations},
             {"phase", step.final_rotation ? "after_final_rotation"
                                           : "after_orbital_update"}},
            false);
  }

  // Only these final properties have actually been evaluated on the saved
  // state.
  void finish(madness::SCF &scf, bool converged) {
    if (work == "loaded" && iterations == 0)
      return;
    const auto evaluation = scf.evaluate_restricted_hf(world_);
    const double energy = evaluation.energy;
    require(std::isfinite(energy), "nonfinite final energy");
    const auto& dipole = evaluation.dipole;
    std::vector<double> dip(3);
    for (int i = 0; i < 3; ++i) {
      dip[i] = dipole(i);
      require(std::isfinite(dip[i]), "nonfinite final dipole");
    }
    scf.current_energy = energy;
    json properties = {
        {{"name", "energy"}, {"value", energy}, {"units", "hartree"}},
        {{"name", "dipole"}, {"value", dip}, {"units", "e*a0"}}};
    auto metrics = has_latest_ ? latest_.metrics : json::object();
    // History remains intact. These additional checks evaluate the saved
    // orbitals and use a same-rung iteration input as the comparison sample.
    json evidence = metrics.at("convergence_report").at("measurements");
    const bool comparable = comparison_ && comparison_->comparison_iteration >= 0;
    double density = 0.0;
    if (comparable)
      density = (evaluation.alpha_density - comparison_->comparison_alpha).norm2();
    const std::map<std::string, double> values = {
        {"bsh_residual", evaluation.bsh_residual}, {"density_change", density},
        {"energy_change", comparable ? std::abs(energy - comparison_->comparison_energy) : 0.0},
        {"delta_kinetic_energy", comparable ? std::abs(evaluation.kinetic - comparison_->comparison_kinetic) : 0.0}};
    for (auto& item : evidence.items()) {
      auto& m = item.value();
      const double value = values.at(item.key());
      const bool available = item.key() == "bsh_residual" || comparable;
      const bool finite = std::isfinite(value);
      m["value"] = available && finite ? json(value) : json(nullptr);
      m["availability"] = !available ? "unavailable" : finite ? "measured" : "invalid";
      m["reason"] = available ? "" : "no_comparable_previous_sample";
      m["valid_for_saved_state"] = available && finite;
      m["phase"] = "saved_orbitals_before_archive";
      m["active_gate"] = false; // diagnostic; the engine verdict is unchanged
      m["satisfied"] = available && finite && !m.at("effective_limit").is_null()
          ? json(value < m.at("effective_limit").get<double>()) : json(nullptr);
      m.erase("sample");
      m.erase("components");
    }
    metrics["convergence_report"]["saved_state_validation"] = {
        {"affects_stopping", false}, {"measurements", evidence},
        {"energy", energy}, {"kinetic_energy", evaluation.kinetic},
        {"comparison", comparable ? json{{"run_id", run_id_}, {"stage", stage_.name},
            {"iteration", comparison_->comparison_iteration}, {"threshold", comparison_->threshold},
            {"k", scf.amo.at(0).k()}, {"energy", comparison_->comparison_energy},
            {"kinetic_energy", comparison_->comparison_kinetic}} : json(nullptr)}};
    publish(scf, converged, properties, metrics,
            {{"phase", "final_properties_on_saved_orbitals"},
             {"iterations_this_run", iterations}},
            true);
  }

private:
  void publish(const madness::SCF &scf, bool converged, json properties,
               const json &metrics, const json &position, bool evaluated) {
    namespace checkpoint = madness::checkpoint;
    checkpoint::Manifest m;
    m.stage = stage_.name;
    m.method = "HF";
    m.generation = "g-" + run_id_ + "-" + std::to_string(++sequence_);
    m.physical = physical_;
    m.numerical = {{"threshold", madness::FunctionDefaults<3>::get_thresh()},
                   {"k", madness::FunctionDefaults<3>::get_k()},
                   {"representation", "mo"},
                   {"protocol", scf.param.protocol()},
                   {"parameters", scf.param.to_json()}};
    m.working_orbitals = {{"state_id", m.generation},
                          {"ordering", "archive_order"}};
    m.position = position;
    if (has_latest_)
      m.position["previous_generation"] = latest_.generation;
    m.metrics = metrics;
    m.targets = {{"dconv", scf.param.dconv()},
                 {"econv", scf.param.econv()},
                 {"effective_dconv",
                  std::max(scf.param.dconv(),
                           madness::FunctionDefaults<3>::get_thresh())}};
    m.stop_reason = converged   ? "converged"
                    : evaluated ? "iteration_limit"
                                : "in_progress";
    for (auto &p : properties)
      p["generation"] = m.generation;
    m.properties = properties;
    const auto directory = store_ / m.generation;
    checkpoint::detail::on_root(world_, [&] {
      validate_measurements(m.metrics);
      fs::create_directories(store_);
      require(fs::create_directory(directory), "generation already exists");
      return std::string{};
    });
    auto meta = scf.orbital_metadata();
    if (!evaluated)
      meta.current_energy =
          1.e10; // input energy is not an output-state property
    if (!converged) {
      meta.converged_for_thresh = 1.e10;
      meta.converged_for_dconv = 1.e10;
    }
    std::string error =
        scf.write_mos(world_, (directory / "state.restartdata").string(), meta);
    // Inventory is built only after every writer has closed. The header records
    // the actual writer count, including the parallel archive's clamp.
    checkpoint::detail::collective_error(world_, error);
    world_.gop.fence();
    if (world_.rank() == 0) {
      try {
        const int writers = checkpoint::detail::header(
            directory / "state.restartdata.00000", true);
        m.archives.push_back(checkpoint::inspect_archive(
            directory, "orbitals", "state.restartdata", {}, writers));
      } catch (const std::exception &e) {
        error = e.what();
      }
    }
    checkpoint::publish(world_, store_, m, error);
    // Reading back also gives every rank the exact, validated published
    // inventory.
    latest_ = checkpoint::read_published(world_, store_);
    has_latest_ = true;
  }
};
} // namespace named_stages
