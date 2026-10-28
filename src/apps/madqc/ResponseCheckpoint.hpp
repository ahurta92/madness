#pragma once

#include "../molresponse/solvers/function_kinetics.hpp"

#include "NamedStages.hpp"
#include <apps/molresponse/solvers/es_save_load.hpp>
#include <apps/molresponse/solvers/fd_save_load.hpp>
#include <apps/molresponse/solvers/vbc_save_load.hpp>
#include <madness/chem/StageCheckpoint.h>

#include <apps/molresponse/calc/calc_manager.hpp>
#include <madness/chem/LadderPlanner.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <functional>
#include <map>
#include <set>

#include <iomanip>
#include <sstream>

namespace named_stages {

// B9: one response stage's convergence settings, built once. The engine rows
// come from the engine's own targets (ConvergencePolicy::effective_for_thresh);
// its gates are per leg and relative, so they are shown as formulas, with the
// per-channel verdicts in the channel records. The written targets come from
// fd_requested_accuracy_targets, the same source the solver stops on.
[[nodiscard]] inline madness::convergence::ConvergenceSettings
response_convergence_settings(const Stage& stage, double final_thresh) {
  namespace cv = madness::convergence;
  using cv::Group;
  using cv::Relation;
  using cv::SetBy;
  const auto& rp = stage.params.get<ResponseParameters>();
  const auto text = [](double value) { return cv::detail::scientific(value); };
  cv::ConvergenceSettings settings;
  settings.stage = stage.name;
  settings.application = "molresponse";
  settings.decider = rp.get<bool>(cv::use_default_criteria_key)
                         ? cv::Decider::engine_check
                         : cv::Decider::requested_criteria;
  settings.decider_set_by_user
      = stage.written.contains(cv::use_default_criteria_key);
  molresponse_v3::ConvergencePolicy policy;
  policy.dconv_user = rp.dconv();
  const auto limits = policy.effective_for_thresh(final_thresh);
  const auto dconv_by
      = stage.written.contains("dconv") ? SetBy::user : SetBy::engine_default;
  const std::string effective = "max(dconv " + text(rp.dconv()) + ", thresh "
                                + text(final_thresh) + ")";
  settings.metrics.push_back(
      {Group::engine_check, "BSH residual per leg", "dconv",
       "engine_bsh_residual", "", "response_amplitude_L2_atomic_units",
       rp.dconv(), dconv_by, Relation::less_than_or_equal, std::nullopt,
       text(limits.bsh_residual) + " + " + text(limits.rtol)
           + " x ||x|| (relative), from 5 x " + effective,
       "per channel; verdicts in the channel records", ""});
  settings.metrics.push_back(
      {Group::engine_check, "density change per leg", "dconv",
       "engine_density_change", "", "kernel_induced_density_L2_atomic_units",
       rp.dconv(), dconv_by, Relation::less_than_or_equal, std::nullopt,
       text(limits.density_residual) + " + " + text(limits.rtol)
           + " x ||rho|| (relative), from 5 x " + effective,
       "per channel; verdicts in the channel records", ""});
  // D-025: an excited-state bundle's own check, per root, absolute and
  // strict (es_root_within).
  if (std::find(stage.properties.begin(), stage.properties.end(),
                "excited_states")
      != stage.properties.end()) {
    settings.metrics.push_back({Group::engine_check, "density change per root",
                                "dconv", "engine_es_density_change", "",
                                "kernel_induced_density_L2_atomic_units",
                                rp.dconv(), dconv_by, Relation::less_than,
                                limits.density_residual, "5 x " + effective,
                                "per root; verdicts in the bundle record", ""});
    settings.metrics.push_back(
        {Group::engine_check, "excitation energy change per root", "dconv",
         "engine_es_omega_change", "", "hartree", rp.dconv(), dconv_by,
         Relation::less_than, limits.omega_residual, "5 x " + effective,
         "per root; verdicts in the bundle record", ""});
  }

  struct Requested {
    const char *name, *keyword, *measurement, *quantity, *units;
    std::string engine_instead;
  };
  const std::vector<Requested> requested = {
      {"induced density change", "dconv", "charge_density_change",
       "kernel_induced_charge_density_change_L2",
       "kernel_induced_density_L2_atomic_units",
       "compares each leg's density change against "
           + text(limits.density_residual) + " + " + text(limits.rtol)
           + " x ||rho||"},
      {"BSH residual (max over functions)", "accuracy.bsh_residual",
       "bsh_residual", "maximum_raw_fixed_point_function_L2",
       "response_amplitude_L2_atomic_units",
       "compares each leg's BSH residual against " + text(limits.bsh_residual)
           + " + " + text(limits.rtol) + " x ||x||"},
      {"kinetic trace change", "accuracy.kinetic_trace_change",
       "delta_kinetic_energy", "absolute_kinetic_trace_change",
       "hartree_times_response_amplitude_scale_squared",
       "has no kinetic term"}};
  const auto targets = molresponse_v3::fd_requested_accuracy_targets(rp);
  for (const auto& row : requested) {
    const bool written = targets.contains(row.measurement);
    settings.metrics.push_back(
        {Group::requested_criteria, row.name, row.keyword, row.measurement,
         row.quantity, row.units,
         written ? std::optional<double>(
                       targets.at(row.measurement).at("value").get<double>())
                 : std::nullopt,
         written ? SetBy::user : SetBy::not_set, Relation::less_than_or_equal,
         std::nullopt, "", "", row.engine_instead});
  }
  settings.metrics.push_back({Group::requested_criteria, "energy change", "",
                              "energy_change", "", "hartree", std::nullopt,
                              SetBy::not_set, Relation::less_than_or_equal,
                              std::nullopt, "",
                              "not applicable: a driven response has no local "
                              "energy; see the reference stage",
                              ""});
  settings.iteration_caps = {
      {static_cast<int>(rp.maxiter())},
      stage.written.contains("maxiter") ? SetBy::user : SetBy::engine_default};
  // A3: no save boundary inside a solve yet (DESIGN §1.4).
  settings.checkpoints.fixed
      = "at every solve end (molresponse has no save boundary inside a "
        "solve); every generation kept";
  settings.warnings = cv::unenforced_target_warnings(settings);
  for (auto& warning : cv::incomplete_request_warnings(
           settings, {"charge_density_change", "bsh_residual"}))
    settings.warnings.push_back(std::move(warning));
  return settings;
}

// The end-of-stage outcome of a response stage (B9): for each written-target
// measurement, the largest saved-state value over the channels and
// frequencies of the final rung, so the limiting channel decides the row.
[[nodiscard]] inline std::vector<madness::convergence::MetricOutcome>
response_saved_outcomes(const json& metadata,
                        const std::string& final_protocol) {
  std::map<std::string, madness::convergence::MetricOutcome> largest;
  // Held in a local: iterating items() of the temporary value() returns
  // would dangle.
  const auto states = metadata.value("fd_states", json::object());
  for (const auto& channel : states.items()) {
    if (!channel.value().contains(final_protocol))
      continue;
    for (const auto& point : channel.value().at(final_protocol).items()) {
      const auto measurements = molresponse_v3::fd_accuracy_measurements(
          point.value()
              .value("saved_state_validation", json::object())
              .value("measurements", json::object()));
      for (const auto& item : measurements.items()) {
        const auto& record = item.value();
        auto& outcome = largest[item.key()];
        outcome.measurement = item.key();
        if (record.value("availability", "") != "measured"
            || !record.at("value").is_number())
          continue;
        const double value = record.at("value").get<double>();
        if (!outcome.value || value > *outcome.value)
          outcome.value = value;
        outcome.availability = "measured";
      }
    }
  }
  std::vector<madness::convergence::MetricOutcome> outcomes;
  outcomes.reserve(largest.size());
  for (auto& item : largest)
    outcomes.push_back(std::move(item.second));
  return outcomes;
}

// B5, B6: how an FD state's rung ended. Entries written since 4f74ffe0c carry
// rung_exit's reason as `stop`; older ones are read through a documented
// mapping: best-effort acceptance and an exhausted rung are an iteration
// limit (or a stall), never convergence.
[[nodiscard]] inline std::string response_stop(const json& entry) {
  if (entry.contains("stop"))
    return entry.at("stop").get<std::string>();
  if (entry.value("converged", false) && !entry.value("accepted", false))
    return "converged";
  if (entry.value("stalled", false))
    return "stalled";
  if (entry.value("accepted", false) || entry.value("rung_exhausted", false))
    return "iteration_limit";
  return "in_progress";
}

// B6: today's engine check as planner rows, judged at the final rung: per
// leg, BSH residual and density change each at most absolute + rtol x norm,
// the solver's own gate (fd_leg_within), from the same policy the solve uses.
[[nodiscard]] inline std::vector<madness::planning::EngineLimit>
response_engine_limits(const Stage& stage, double final_thresh) {
  molresponse_v3::ConvergencePolicy policy;
  policy.dconv_user = stage.params.get<ResponseParameters>().dconv();
  const auto limits = policy.effective_for_thresh(final_thresh);
  return {{"bsh_residual", limits.bsh_residual, limits.rtol, true},
          {"density_change", limits.density_residual, limits.rtol, true}};
}

// B6: a state record (state_records) as the planner reads it. The written
// targets are judged on the saved state's measurements, under the names the
// targets use (fd_accuracy_measurements).
[[nodiscard]] inline madness::planning::SavedState
response_saved_state(const json& record) {
  madness::planning::SavedState state;
  state.identity = record.at("identity").get<std::string>();
  state.reference_generation
      = record.at("reference_generation").get<std::string>();
  state.perturbation = record.at("perturbation").get<std::string>();
  // An excited-state bundle has no frequency (null).
  state.frequency = record.at("frequency").is_number()
                        ? record.at("frequency").get<double>()
                        : 0.0;
  state.rung = {record.at("rung").at("thresh").get<double>(),
                record.at("rung").at("k").get<int>()};
  state.iterations_at_rung = record.at("iterations_at_rung").get<int>();
  state.stop = madness::convergence::stop_reason_from_string(
      record.at("stop").get<std::string>());
  state.diverged = record.at("diverged").get<bool>();
  state.engine_measurements = record.at("engine_measurements");
  state.saved_measurements = molresponse_v3::fd_accuracy_measurements(
      record.at("saved_measurements"));
  return state;
}

// The response states a request's selected values are built from (D-017,
// D-019): FD legs (perturbation, frequency) and quadratic sources (VBC node
// ids), in request order, each once. The engine may solve more (it takes one
// grid and one axis set for every property); those extras feed no selected
// value.
//   polarizability alpha_AB(w): the A and B legs at w.
//   hyperpolarizability beta_ABC(-2w; w, w): the A leg at 2w, the B and C
//     legs at w, the source of B and C at (w, w).
//   raman d(alpha_AB)/dQ(w): the A and B legs at w, the nuclear leg Q at 0,
//     the source of B and Q at (w, 0); Q spans atoms x directions, every atom
//     of the reference and x, y, z when one is empty (raman_coordinates).
// A request's properties. One saved before the request named its properties
// (polarizability only) carried its frequencies and components at the top.
[[nodiscard]] inline json requested_properties(const json& request) {
  if (request.contains("properties"))
    return request.at("properties");
  return {{"polarizability",
           {{"frequencies", request.at("frequencies")},
            {"components", request.at("components")}}}};
}

struct ResponseNeeds {
  std::vector<std::pair<molresponse_v3::Perturbation, double>> legs;
  std::vector<std::string> sources;
};
[[nodiscard]] inline ResponseNeeds response_needs(const json& request,
                                                  int natom) {
  using molresponse_v3::Perturbation;
  ResponseNeeds needs;
  std::set<std::string> seen;
  const auto leg = [&](const Perturbation& p, double frequency) {
    if (seen.insert(molresponse_v3::fd_node_id(p, frequency)).second)
      needs.legs.emplace_back(p, frequency);
  };
  const auto source = [&](const Perturbation& b, const Perturbation& c,
                          double fb, double fc) {
    auto id = molresponse_v3::vbc_node_id(b, c, fb, fc);
    if (seen.insert(id).second)
      needs.sources.push_back(std::move(id));
  };
  const auto dipole
      = [](char axis) { return Perturbation::dipole(axis - 'x'); };
  const auto properties = requested_properties(request);
  if (properties.contains("polarizability")) {
    const auto& r = properties.at("polarizability");
    std::set<char> axes;
    for (const auto& component : r.at("components"))
      for (char axis : component.get<std::string>())
        axes.insert(axis);
    for (double frequency : r.at("frequencies"))
      for (char axis : axes)
        leg(dipole(axis), frequency);
  }
  if (properties.contains("hyperpolarizability")) {
    const auto& r = properties.at("hyperpolarizability");
    for (double frequency : r.at("frequencies"))
      for (const auto& component : r.at("components")) {
        const auto abc = component.get<std::string>();
        leg(dipole(abc[0]), 2.0 * frequency);
        leg(dipole(abc[1]), frequency);
        leg(dipole(abc[2]), frequency);
        source(dipole(abc[1]), dipole(abc[2]), frequency, frequency);
      }
  }
  if (properties.contains("raman")) {
    const auto& r = properties.at("raman");
    auto atoms = r.at("atoms").get<std::vector<int>>();
    if (atoms.empty())
      for (int atom = 0; atom < natom; ++atom)
        atoms.push_back(atom);
    auto directions = r.at("directions").get<std::string>();
    if (directions.empty())
      directions = "xyz";
    for (double frequency : r.at("frequencies"))
      for (int atom : atoms)
        for (char axis : directions) {
          const auto nuclear = Perturbation::nuclear(atom, axis - 'x');
          for (const auto& component : r.at("components")) {
            const auto ab = component.get<std::string>();
            leg(dipole(ab[0]), frequency);
            leg(dipole(ab[1]), frequency);
            leg(nuclear, 0.0);
            source(dipole(ab[1]), nuclear, frequency, 0.0);
          }
        }
  }
  return needs;
}

// The excited-state bundle a request names (D-025): one state, keyed by the
// engine's node id, which carries the method and the root count.
struct BundleNeed {
  std::string id;
  bool tda = true;
  int roots = 0;
};
[[nodiscard]] inline std::optional<BundleNeed>
requested_bundle(const json& request) {
  const auto properties = requested_properties(request);
  if (!properties.contains("excited_states"))
    return std::nullopt;
  const auto& r = properties.at("excited_states");
  const bool tda = r.at("method").get<std::string>() == "tda";
  const int roots = r.at("roots").get<int>();
  return BundleNeed{molresponse_v3::es_node_id(tda, roots), tda, roots};
}

// The engine's name for each named response property's rows.
[[nodiscard]] inline std::string engine_property(const std::string& name) {
  static const std::map<std::string, std::string> names
      = {{"polarizability", "alpha"},
         {"hyperpolarizability", "beta"},
         {"raman", "raman"}};
  return names.at(name);
}

// A requested property's value as the engine recorded it: its rows, or for
// excited states the bundle at the final rung `final_key` (no rows of its own).
[[nodiscard]] inline json property_value(const json& metadata,
                                         const std::string& name,
                                         const std::string& final_key) {
  if (name == "excited_states")
    return metadata.value("excited_states", json::object())
        .value(final_key, json());
  return metadata.at("properties").value(engine_property(name), json());
}

// D-025: the bundle's engine check as planner rows: every root's density
// change and eigenvalue change below the engine's targets (es_root_within,
// strict, absolute), judged on the largest value over the roots
// (bundle_measurements).
[[nodiscard]] inline std::vector<madness::planning::EngineLimit>
bundle_engine_limits(const Stage& stage, double final_thresh) {
  molresponse_v3::ConvergencePolicy policy;
  policy.dconv_user = stage.params.get<ResponseParameters>().dconv();
  const auto limits = policy.effective_for_thresh(final_thresh);
  return {{"density_residual", limits.density_residual, 0.0, false},
          {"omega_residual", limits.omega_residual, 0.0, false}};
}

// A saved bundle's engine-check sample: per measurement the largest value
// over the roots, and the root it came from; null when a root has none yet.
[[nodiscard]] inline json bundle_measurements(const json& bundle) {
  json measurements = json::object();
  for (const char* name : {"density_residual", "omega_residual"}) {
    json worst = {{"value", nullptr}, {"availability", "unavailable"}};
    bool every_root = true;
    for (const auto& root : bundle.value("roots", json::array())) {
      const auto value = root.value(name, json());
      if (!value.is_number()) {
        every_root = false;
        break;
      }
      if (worst.at("value").is_null()
          || value.get<double>() > worst.at("value").get<double>())
        worst = {{"value", value},
                 {"availability", "measured"},
                 {"root_id", root.value("root_id", "")}};
    }
    if (!every_root)
      worst = {{"value", nullptr}, {"availability", "unavailable"}};
    measurements[name] = std::move(worst);
  }
  return measurements;
}

// B6: iterations recorded per state and rung before this attempt (keyed by
// the engine's channel|protocol|frequency keys), the rungs the planner
// resumes, and the rungs this attempt solved, so the total across runs is
// kept.
struct IterationTotals {
  std::map<std::string, int> previous;
  std::set<std::string> resumed, solved;
};

// Iterations at one state's rung across runs, from `iter`, the count of
// the attempt that last solved it: a rung solved in this attempt adds to
// what it had when the planner resumed it, else counts this attempt alone;
// an untouched one keeps its recorded total.
[[nodiscard]] inline int iterations_across_runs(const IterationTotals& totals,
                                                const std::string& key,
                                                int iter) {
  const auto previous = totals.previous.find(key);
  const int before = previous != totals.previous.end() ? previous->second : 0;
  if (totals.solved.contains(key))
    return (totals.resumed.contains(key) ? before : 0) + iter;
  if (previous != totals.previous.end())
    return before;
  return iter;
}

// B6, D-016: the planner's decision per requested state, for the engine (by
// its node id), its record, for the manifest and results.json, and the state
// keys it resumes at their own saved rung (IterationTotals::resumed).
struct ResponsePlans {
  molresponse_v3::StatePlans engine;
  json records = json::array();
  std::set<std::string> resumed;
};

// S3b (D-048): the read side of --continue, as functions of their arguments,
// defined after ResponseCheckpoint (whose state keys they use).
[[nodiscard]] inline json response_request(const Stage& stage);
[[nodiscard]] inline json
response_state_records(const json& metadata,
                       const std::string& reference_generation,
                       const IterationTotals& totals);
[[nodiscard]] inline ResponsePlans
plan_response_states(const Stage& stage, const json& stage_request, int natom,
                     const std::string& reference_generation,
                     const json& records);

// D-045, the response AppSpec's observe (S3c1): an engine save event as a
// work item, the event with the work it did. A bundle restarted from its own
// saved rungs is resumed, as an FD state restarted from its own is. Nothing
// is carried across events.
[[nodiscard]] inline json response_observe(const Stage&, int, const json& event,
                                           json&) {
  const auto& seed = event.at("seed");
  auto item = event;
  item["work"] = seed == "fd_restart" || seed == "es_restart" ? "resumed"
                 : seed == "nearest_frequency"                ? "seeded"
                                                              : "solved";
  return item;
}

// Scientific policy for the named response frontend. The engine still owns
// state serialization, restart guesses, protocol scheduling and analysis.
class ResponseCheckpoint {
  using Manifest = madness::checkpoint::Manifest;
  madness::World& world_;
  const Stage& stage_;
  std::filesystem::path store_, attempt_;
  std::string run_id_;
  json reference_, request_;
  int natom_;
  Manifest latest_;
  bool restored_ = false, loaded_ = false;
  std::size_t sequence_ = 0;
  // B6: this attempt's iteration totals and plans (S3b: the functions that
  // read them take them as arguments).
  IterationTotals totals_;
  ResponsePlans plans_;
  // The skeleton's publish (publish_generation), from the run binding.
  GenerationPublisher publish_;

  friend json response_state_records(const json& metadata,
                                     const std::string& reference_generation,
                                     const IterationTotals& totals);
  friend ResponsePlans
  plan_response_states(const Stage& stage, const json& stage_request, int natom,
                       const std::string& reference_generation,
                       const json& records);

  static std::string state_key(const std::string& channel,
                               const std::string& protocol,
                               const std::string& frequency) {
    return channel + "|" + protocol + "|" + frequency;
  }
  // An excited-state bundle's state key at one rung (no frequency).
  static std::string bundle_key(const std::string& node_id,
                                const std::string& protocol) {
    return node_id + "|" + protocol;
  }
  // The engine key of a request's final rung.
  static std::string request_final_key(const json& request) {
    const double threshold = request.at("protocol").back();
    return molresponse_v3::protocol_key(threshold, protocol_k(threshold));
  }
  // Where the engine saves the bundle of rung `protocol` (relative to the
  // calc directory): the per-root archives and their index, roots.json.
  static std::string bundle_directory(const std::string& protocol) {
    return "es__" + protocol;
  }
  // The engine's transition analysis of a converged bundle at rung `protocol`.
  static std::string analysis_file(const std::string& protocol) {
    return "es_analysis__" + protocol + ".json";
  }

public:
  // The generation checks below are functions of their arguments: the
  // response AppSpec binds members and validate_generation (S3a, D-036).
  static void require(bool ok, const std::string& why) {
    if (!ok) throw std::runtime_error("response checkpoint: " + why);
  }

  // Derive the required archive inventory from state definitions, rather than
  // trusting whatever subset happens to be listed in the storage manifest.
  // Dipole legs and nuclear-displacement legs (nuc_<atom>_<axis>) are FD
  // states; quadratic sources (vbc_states) are archives of their own; an
  // excited-state bundle (D-025) is one archive per root under es__<key>/.
  static std::vector<std::string> members(const json& metadata) {
    std::vector<std::string> names;
    const auto nuclear_channel = [](const std::string& channel) {
      const auto axis = channel.rfind('_');
      return channel.starts_with("nuc_") && axis > 4
             && axis + 2 == channel.size()
             && std::string("xyz").find(channel.back()) != std::string::npos
             && channel.find_first_not_of("0123456789", 4) == axis;
    };
    for (const auto& channel : metadata.at("fd_states").items()) {
      require(channel.key() == "dipole_x" || channel.key() == "dipole_y"
                  || channel.key() == "dipole_z"
                  || nuclear_channel(channel.key()),
              "unsupported channel");
      for (const auto& rung : channel.value().items()) {
        const auto& protocol = metadata.at("protocols").at(rung.key());
        require(molresponse_v3::protocol_key(protocol.at("thresh"), protocol.at("k")) == rung.key(),
                "inconsistent protocol definition");
        for (const auto& point : rung.value().items()) {
          const auto& e = point.value();
          // Named continuation must not silently turn a failed state into a
          // fresh solve. This check runs both before loading and publication,
          // so a diverged solve cannot replace the previous usable generation.
          require(!e.value("diverged", false), "invalid response state: diverged");
          // Older executor verdicts could certify one update without a density
          // comparison. Do not trust that flag on reuse, or certify a stall.
          // Explicit best-effort acceptance remains a separate engine outcome.
          if (e.value("converged", false) && !e.value("accepted", false))
            require(e.value("iter", 0) >= 2 && !e.value("stalled", false),
                    "invalid response state: ineligible strict convergence");
          auto validate_bsh = [&](const json& metric) {
            if (metric.contains("per_function"))
              require(madness::valid_function_residual_record(metric.at("per_function")),
                      "invalid response state: per-function BSH measurement");
          };
          auto validate_kinetic = [&](const json& measurements) {
            if (measurements.contains("bsh_residual")) validate_bsh(measurements.at("bsh_residual"));
            if (measurements.contains("kinetic_energy"))
              require(molresponse_v3::valid_function_kinetic_record(measurements.at("kinetic_energy")),
                      "invalid response state: per-function kinetic measurement");
          };
          if (e.contains("convergence_measurements")) validate_kinetic(e.at("convergence_measurements"));
          if (e.contains("saved_state_validation")) validate_kinetic(e.at("saved_state_validation").at("measurements"));
          if (e.contains("function_bsh_history")) {
            const auto& history = e.at("function_bsh_history");
            require(history.is_array(), "invalid response state: per-function BSH history");
            for (const auto& sample : history) {
              require(sample.contains("per_function") && sample.contains("value") &&
                      sample.at("value").is_number() && sample.value("availability", "") == "measured",
                      "invalid response state: per-function BSH history");
              const double value = sample.at("value").get<double>();
              require(std::isfinite(value) && value >= 0.0,
                      "invalid response state: per-function BSH history aggregate");
              validate_bsh(sample);
            }
          }
          if (e.contains("function_kinetic_history"))
            for (const auto& sample : e.at("function_kinetic_history"))
              require(molresponse_v3::valid_function_kinetic_record(sample),
                      "invalid response state: per-function kinetic history");
          if (e.contains("saved_state_validation"))
            for (const auto& metric : e.at("saved_state_validation").at("measurements").items())
              require(metric.value().value("availability", "") != "invalid",
                      "invalid response state: nonfinite saved-state measurement");
          if (e.contains("convergence_measurements"))
            for (const auto& metric : e.at("convergence_measurements").items())
              require(metric.value().value("availability", "") != "invalid",
                      "invalid response state: nonfinite measurement");
          // Bind progress evidence to the request that produced this attempt,
          // including historical/pending entries retained by later requests.
          const auto origins = metadata.value("checkpoint_requests", json::object());
          const auto history = e.value("function_kinetic_history", json::array());
          const std::string origin = history.empty() ? "" :
              history.at(0).at("sample").value("attempt_directory", "");
          require(!e.contains("checkpoint_progress") || origins.contains(origin),
                  "invalid response state: missing originating progress request");
          if (origins.contains(origin)) {
            const auto& request = origins.at(origin);
            json progress = request.value("progress", json::object());
            if (!progress.empty()) progress["protocol"] = request.at("protocol");
            require(progress.empty() ? !e.contains("checkpoint_progress") :
                        e.value("checkpoint_progress", json::object()) == progress,
                    "invalid response state: originating progress request mismatch");
          }
          require(molresponse_v3::valid_fd_explicit_accuracy(e, protocol.at("thresh"), protocol.at("k")),
                  "invalid response state: explicit accuracy evidence");
          require(molresponse_v3::valid_metric_diagnostics(e),
                  "invalid response state: metric comparison");
          const double frequency = e.at("freq");
          require(std::isfinite(frequency) &&
                  molresponse_v3::ResponseMetadata::freq_key(frequency) == point.key(),
                  "inconsistent frequency definition");
          const auto name = molresponse_v3::response_filename(channel.key(), rung.key(), frequency);
          require(e.at("archive") == name && e.at("backend") == "native" &&
                  e.at("shell") == "closed_shell" &&
                  e.at("type") == (std::abs(frequency) < 1e-12 ? "static" : "full"),
                  "inconsistent archive/channel/model definition");
          require(e.at("checkpoint_dconv").get<double>() > 0.0,
                  "missing response convergence target");
          names.push_back(name);
        }
      }
    }
    // A source is built once from two converged legs at its rung and saved
    // converged; the engine drops it when a leg is solved again there.
    for (const auto& source : metadata.at("vbc_states").items())
      for (const auto& rung : source.value().items()) {
        const auto& protocol = metadata.at("protocols").at(rung.key());
        require(molresponse_v3::protocol_key(protocol.at("thresh"),
                                             protocol.at("k"))
                    == rung.key(),
                "inconsistent protocol definition");
        const auto& e = rung.value();
        const auto name
            = molresponse_v3::vbc_archive_basename(source.key(), rung.key());
        require(e.value("converged", false) && !e.value("diverged", false),
                "invalid quadratic source: not converged");
        require(e.at("archive") == name && e.at("backend") == "native",
                "inconsistent quadratic source archive");
        names.push_back(name);
      }
    // A bundle at a rung: the FD rules (never diverged; converged only after
    // two iterations and not stalled; finite measurements), the engine's main
    // bundle only (not a seed twin or a warmup cache), one archive per root.
    for (const auto& rung : metadata.at("excited_states").items()) {
      const auto& protocol = metadata.at("protocols").at(rung.key());
      require(
          molresponse_v3::protocol_key(protocol.at("thresh"), protocol.at("k"))
              == rung.key(),
          "inconsistent protocol definition");
      const auto& e = rung.value();
      const auto dir = bundle_directory(rung.key());
      require(!e.value("diverged", false),
              "invalid excited-state bundle: diverged");
      if (e.value("converged", false))
        require(e.value("iter", 0) >= 2 && !e.value("stalled", false),
                "invalid excited-state bundle: ineligible convergence");
      const auto type = e.value("type", std::string());
      require((type == "tda" || type == "full")
                  && e.value("shell", std::string()) == "closed_shell"
                  && e.value("backend", std::string()) == "native"
                  && e.value("bundle_dir", std::string()) == dir,
              "inconsistent excited-state bundle definition");
      const int count = e.at("n_roots").get<int>();
      const auto& roots = e.at("roots");
      require(count >= 1 && roots.is_array()
                  && roots.size() == static_cast<std::size_t>(count),
              "inconsistent excited-state root count");
      for (int slot = 0; slot < count; ++slot) {
        const auto& root = roots.at(slot);
        const auto omega = root.value("omega", json());
        require(root.value("slot", -1) == slot && omega.is_number()
                    && std::isfinite(omega.get<double>()),
                "invalid excited-state root: eigenvalue");
        // A residual is a finite number, or null before the first step;
        // a converged bundle has measured every one.
        for (const char* name : {"density_residual", "omega_residual"}) {
          const auto value = root.value(name, json());
          require(value.is_number()
                      ? std::isfinite(value.get<double>())
                            && value.get<double>() >= 0.0
                      : value.is_null() && !e.value("converged", false),
                  "invalid excited-state root: nonfinite measurement");
        }
        names.push_back(dir + "/"
                        + molresponse_v3::detail_save_load::root_file(slot));
      }
    }
    std::sort(names.begin(), names.end());
    require(std::adjacent_find(names.begin(), names.end()) == names.end(), "duplicate state member");
    return names;
  }

private:
  // Every leg and source the requested values are built from
  // (response_needs) is at the final rung and converged. `loaded`: the state
  // keys the planner loaded (D-016); their verdict is the planner's, so only
  // their identity is checked here.
  static bool request_complete(const json& metadata, const json& request,
                               int natom,
                               const std::set<std::string>& loaded = {}) {
    const double threshold = request.at("protocol").back();
    const auto key = molresponse_v3::protocol_key(threshold, protocol_k(threshold));
    if (!metadata.at("protocols").contains(key) ||
        metadata.at("protocols").at(key).at("thresh").get<double>() != threshold)
      return false;
    const auto needs = response_needs(request, natom);
    for (const auto& [perturbation, frequency] : needs.legs) {
      const auto fkey = molresponse_v3::ResponseMetadata::freq_key(frequency);
      const auto channel = perturbation.description();
      const auto& states = metadata.at("fd_states");
      if (!states.contains(channel) || !states.at(channel).contains(key)
          || !states.at(channel).at(key).contains(fkey))
        return false;
      const auto& e = states.at(channel).at(key).at(fkey);
      if (loaded.contains(state_key(channel, key, fkey))) {
        if (e.at("freq").get<double>() != frequency
            || e.value("diverged", false))
          return false;
        continue;
      }
      // Written targets block completion only when they decide (B8).
      if (request.value("decided_by", "") == "requested_criteria"
          && request.contains("accuracy_targets")) {
        const auto measurements
            = e.value("saved_state_validation", json::object())
                  .value("measurements", json::object());
        if (madness::evaluate_explicit_criteria(
                molresponse_v3::fd_accuracy_measurements(measurements),
                request.at("accuracy_targets"))
                .at("status")
            != "satisfied")
          return false;
      }
      if (e.at("freq").get<double>() != frequency
          || !e.at("converged").get<bool>() || e.value("accepted", false)
          || e.value("diverged", false)
          || e.at("checkpoint_dconv").get<double>()
                 > request.at("dconv").get<double>())
        return false;
    }
    const auto& sources = metadata.at("vbc_states");
    for (const auto& id : needs.sources)
      if (!sources.contains(id) || !sources.at(id).contains(key)
          || !sources.at(id).at(key).value("converged", false))
        return false;
    // The bundle at the final rung, the requested method and root count.
    if (const auto bundle = requested_bundle(request)) {
      const auto& bundles = metadata.at("excited_states");
      if (!bundles.contains(key))
        return false;
      const auto& e = bundles.at(key);
      if (e.value("type", std::string()) != (bundle->tda ? "tda" : "full")
          || e.value("n_roots", 0) != bundle->roots
          || e.value("diverged", false))
        return false;
      if (!loaded.contains(bundle_key(bundle->id, key))
          && (!e.value("converged", false) || e.value("stalled", false)))
        return false;
    }
    return true;
  }

public:
  // A2: the science checks on a complete generation of stage `stage`, whose
  // reference has `natom` atoms.
  static void validate_generation(const std::string& stage, int natom,
                                  const Manifest& m) {
    require(m.stage == stage && m.method == "TDHF"
                && m.physical
                       == json{{"schema", "madqc.response.v1"},
                               {"shell", "closed_shell"}},
            "stage or response model mismatch");
    if (m.numerical.at("request").contains("progress")) {
      const auto origins = m.numerical.at("engine_metadata").value("checkpoint_requests", json::object());
      const auto attempt = m.numerical.value("attempt_directory", "");
      require(origins.contains(attempt) && origins.at(attempt) == m.numerical.at("request"),
              "invalid response state: missing or inconsistent originating request registry");
    }
    auto expected = members(m.numerical.at("engine_metadata"));
    std::vector<std::string> actual;
    for (const auto& a : m.archives) {
      require(a.role == "response" && a.indices.empty(), "invalid response archive role");
      actual.push_back(a.basename);
    }
    std::sort(actual.begin(), actual.end());
    require(!expected.empty() && expected == actual, "required response inventory is missing or inconsistent");
    // Each bundle's index (its roots.json, kept in the manifest because it
    // is not an archive) agrees with the bundle and names its archives.
    const auto& engine_metadata = m.numerical.at("engine_metadata");
    const auto indices = m.numerical.value("bundle_indices", json::object());
    require(indices.size() == engine_metadata.at("excited_states").size(),
            "excited-state bundle index is missing or extra");
    for (const auto& rung : engine_metadata.at("excited_states").items()) {
      require(indices.contains(rung.key()),
              "excited-state bundle index is missing");
      const auto& index = indices.at(rung.key());
      const auto& e = rung.value();
      const auto writers = index.value("writer_nproc", json());
      bool files = index.value("roots", json()).is_array()
                   && index.at("roots").size() == e.at("roots").size();
      for (std::size_t slot = 0; files && slot < index.at("roots").size();
           ++slot)
        files = index.at("roots").at(slot).value("file", std::string())
                == molresponse_v3::detail_save_load::root_file(
                    static_cast<int>(slot));
      require(files && index.value("type", json()) == e.at("type")
                  && index.value("shell", json()) == e.at("shell")
                  && index.value("n_roots", json()) == e.at("n_roots")
                  && index.value("backend", json()) == e.at("backend")
                  && writers.is_number_integer() && writers.get<int>() >= 1,
              "excited-state bundle index is inconsistent with the bundle");
    }
    // One record per requested property, its value the engine's rows (for
    // excited states the bundle at the final rung).
    const auto requested = requested_properties(m.numerical.at("request"));
    const auto final_key = request_final_key(m.numerical.at("request"));
    std::set<std::string> recorded;
    for (const auto& p : m.properties) {
      const auto name = p.at("name").get<std::string>();
      require(p.at("generation") == m.generation && requested.contains(name)
                  && recorded.insert(name).second
                  && p.at("value")
                         == property_value(engine_metadata, name, final_key),
              "property/state record is inconsistent");
    }
    if (m.stop_reason == "converged")
      require(recorded.size() == requested.size()
                  && request_complete(
                      m.numerical.at("engine_metadata"),
                      m.numerical.at("request"), natom,
                      loaded_states(m.numerical.value("plans", json::array()))),
              "completed property lacks converged prerequisites");
  }

private:
  // The state keys a generation's plan records say were loaded and met.
  static std::set<std::string> loaded_states(const json& plan_records) {
    std::set<std::string> loaded;
    for (const auto& record : plan_records)
      if (record.contains("loaded_state"))
        loaded.insert(record.at("loaded_state").get<std::string>());
    return loaded;
  }

  // R4/R7: the same perturbation's saved state nearest in frequency (on a
  // tie the lower), at the finest rung the start rung can load (coarser or
  // equal, not aliased): the record of the state the engine's picker loads.
  [[nodiscard]] static std::optional<madness::planning::SavedState>
  frequency_seed(const json& records, const std::string& perturbation,
                 double frequency, const std::string& reference_generation,
                 const std::map<std::string, double>& requested_rungs,
                 const madness::planning::Rung& start) {
    std::vector<madness::planning::SavedState> usable;
    for (const auto& record : records) {
      if (record.at("perturbation") != perturbation)
        continue;
      auto saved = response_saved_state(record);
      const auto rung = requested_rungs.find(
          molresponse_v3::protocol_key(saved.rung.thresh, saved.rung.k));
      const bool aliased
          = rung != requested_rungs.end() && rung->second != saved.rung.thresh;
      if (!aliased && saved.rung.thresh >= start.thresh
          && saved.rung.k <= start.k)
        usable.push_back(std::move(saved));
    }
    auto seed = madness::planning::nearest_seed(
        usable, perturbation, reference_generation, frequency);
    if (!seed)
      return seed;
    for (const auto& state : usable)
      if (!state.diverged && state.frequency == seed->frequency
          && (state.rung.k > seed->rung.k
              || (state.rung.k == seed->rung.k
                  && state.rung.thresh < seed->rung.thresh)))
        seed = state;
    return seed;
  }

  // The request every state of this stage is planned against, less its
  // identity; `requested_rungs` maps each requested rung's engine key to its
  // threshold.
  [[nodiscard]] static madness::planning::LadderRequest
  ladder_request(const Stage& stage, const std::string& reference_generation,
                 std::map<std::string, double>& requested_rungs) {
    madness::planning::LadderRequest request;
    for (double threshold :
         stage.params.get<CalculationParameters>().protocol()) {
      request.ladder.push_back({threshold, protocol_k(threshold)});
      requested_rungs.emplace(
          molresponse_v3::protocol_key(threshold, protocol_k(threshold)),
          threshold);
    }
    const double final_thresh = request.ladder.back().thresh;
    request.reference_generation = reference_generation;
    request.settings = response_convergence_settings(stage, final_thresh);
    request.engine_limits = response_engine_limits(stage, final_thresh);
    request.forced_iterations = stage.params.get<ResponseParameters>().get<int>(
                                    "diagnostic.iterations_per_rung")
                                > 0;
    return request;
  }

  // The record a state is planned from: its finest saved rung the engine can
  // use (see plan_response_states).
  [[nodiscard]] static std::optional<madness::planning::SavedState>
  own_state(const json& records, const std::string& identity,
            const std::map<std::string, double>& requested_rungs,
            const madness::planning::Rung& final_rung) {
    std::optional<madness::planning::SavedState> own;
    for (const auto& record : records) {
      if (record.at("identity") != identity)
        continue;
      auto saved = response_saved_state(record);
      const auto rung = requested_rungs.find(
          molresponse_v3::protocol_key(saved.rung.thresh, saved.rung.k));
      const bool aliased
          = rung != requested_rungs.end() && rung->second != saved.rung.thresh;
      const bool tighter_than_final
          = saved.rung.thresh < final_rung.thresh
            && !madness::scf_continuation::same_rung(saved.rung, final_rung);
      if (aliased || tighter_than_final)
        continue;
      if (!own || saved.rung.thresh < own->rung.thresh)
        own = std::move(saved);
    }
    return own;
  }

  // One state's plan, added to `plans`: the engine's input (by the engine's
  // node id), the iteration bookkeeping (resumed) and the record for
  // results.json and the manifest. `fields` name the state in the record;
  // `key_at` gives its state key at a rung (by engine protocol key).
  static void
  adopt(ResponsePlans& plans, const madness::planning::LadderPlan& plan,
        const std::optional<madness::planning::SavedState>& own,
        const std::optional<madness::planning::SavedState>& seed,
        const madness::planning::LadderRequest& request,
        const std::string& node_id, const json& fields,
        const std::function<std::string(const std::string&)>& key_at) {
    namespace planning = madness::planning;
    using Kind = molresponse_v3::StatePlan::Kind;
    const std::map<planning::Action, const char*> words
        = {{planning::Action::load, "loaded"},
           {planning::Action::resume, "resumed"},
           {planning::Action::fresh, "fresh"},
           {planning::Action::seed, "seeded"}};
    const auto& start = request.ladder.at(plan.start_rung);
    molresponse_v3::StatePlan engine_plan;
    engine_plan.start_thresh = start.thresh;
    json record = {{"state", request.identity},
                   {"action", words.at(plan.action)},
                   {"rung",
                    {{"index", plan.start_rung + 1},
                     {"thresh", start.thresh},
                     {"k", start.k}}},
                   {"iterations_done_at_rung", plan.iterations_done},
                   {"why", plan.why},
                   {"unmet", plan.unmet},
                   {"unmet_not_enforceable", plan.unmet_not_enforceable},
                   {"summary", planning::describe(plan, request.ladder)}};
    record.update(fields);
    if (plan.action == planning::Action::load) {
      engine_plan.kind = Kind::load;
      engine_plan.loaded_key
          = molresponse_v3::protocol_key(own->rung.thresh, own->rung.k);
      // R2 loads a stalled final rung with the request unmet: not met.
      engine_plan.met = plan.unmet.empty();
      if (engine_plan.met)
        record["loaded_state"] = key_at(engine_plan.loaded_key);
    } else if (plan.action == planning::Action::resume) {
      engine_plan.kind = Kind::resume;
      if (madness::scf_continuation::same_rung(own->rung, start))
        plans.resumed.insert(
            key_at(molresponse_v3::protocol_key(start.thresh, start.k)));
    } else if (plan.action == planning::Action::seed) {
      engine_plan.kind = Kind::seed;
      engine_plan.seed_frequency = seed->frequency;
      record["seed_from"] = plan.seed_from;
      record["seed_frequency"] = seed->frequency;
    } else {
      engine_plan.kind = Kind::fresh;
    }
    plans.engine[node_id] = engine_plan;
    plans.records.push_back(std::move(record));
  }

public:
  std::size_t iterations = 0;
  json work_items = json::array();
  std::string work = "solved";
  bool reference_changed = false;

  // `natom`: the reference molecule's atom count, which "every atom" in a
  // Raman request means. `publish`: the skeleton's publish
  // (publish_generation, S3c2).
  ResponseCheckpoint(madness::World& world, const Stage& stage,
                     const std::filesystem::path& stage_directory,
                     const std::filesystem::path& attempt, std::string run_id,
                     json reference, int natom, GenerationPublisher publish)
      : world_(world), stage_(stage), store_(stage_directory / "checkpoints"),
        attempt_(attempt), run_id_(std::move(run_id)),
        reference_(std::move(reference)), natom_(natom),
        publish_(std::move(publish)) {
    request_ = response_request(stage);
  }

  [[nodiscard]] bool targets_decide() const {
    return request_.value("decided_by", "") == "requested_criteria";
  }

  // D-049: --continue reads a generation only when the stage has one; a
  // stage with none (newly added, or no completed first save) starts fresh.
  // Collective.
  [[nodiscard]] bool has_generation() const {
    bool exists = false;
    madness::checkpoint::detail::on_root(world_, [&] {
      exists = std::filesystem::exists(store_ / "latest");
      return std::string{};
    });
    world_.gop.broadcast(exists, 0);
    return exists;
  }

  // The generation --continue read, validated (validate_generation): return
  // true only when the requested analysis itself is already complete.
  // Otherwise restore prerequisites privately and let the engine schedule work.
  bool load(const Manifest& latest) {
    namespace checkpoint = madness::checkpoint;
    latest_ = latest;
    reference_changed = latest_.reference != reference_;
    if (reference_changed) return false; // no old descendants enter this attempt
    restored_ = true;
    loaded_ = stage_.params.get<ResponseParameters>().get<int>("diagnostic.iterations_per_rung") == 0 &&
              latest_.stop_reason == "converged" &&
              latest_.numerical.at("request") == request_ && !latest_.properties.empty();
    if (loaded_) { work = "loaded"; return true; }

    // A bundle's per-root archives load only with the rank count that wrote
    // them (DESIGN B6, until the N-to-M question of §5 is decided).
    const auto indices
        = latest_.numerical.value("bundle_indices", json::object());
    for (const auto& index : indices.items()) {
      const int writers = index.value().at("writer_nproc").get<int>();
      require(writers == world_.size(),
              "ES bundle in " + bundle_directory(index.key())
                  + " was written with " + std::to_string(writers)
                  + " process(es) but loaded with "
                  + std::to_string(world_.size()) + "; re-run with "
                  + std::to_string(writers) + " rank(s)");
    }
    auto metadata = latest_.numerical.at("engine_metadata");
    // B6: the generation's state records (or, for one written before them,
    // records derived from its entries): totals so far, and the plans.
    const auto reference_generation = reference_.value("generation", "");
    auto records = latest_.numerical.value("states", json::array());
    if (records.empty())
      records = response_state_records(metadata, reference_generation, totals_);
    for (const auto& record : records)
      totals_.previous[record.at("key").get<std::string>()]
          = record.at("iterations_at_rung").get<int>();
    plans_ = plan_response_states(stage_, request_, natom_,
                                  reference_generation, records);
    totals_.resumed = plans_.resumed;
    // The engine's state entries are copied as saved (D7, D-016); the plans
    // above steer it. Analysis is derived: old property rows and the Raman
    // normal-mode analysis are left out so they cannot certify this request.
    clear_analysis(metadata);
    // A bundle the planner loads keeps its transition analysis: the engine
    // writes one only when it solves the bundle.
    const auto analyses
        = latest_.numerical.value("transition_analysis", json::object());
    std::vector<std::string> keep;
    if (const auto bundle = requested_bundle(request_)) {
      const auto plan = plans_.engine.find(bundle->id);
      if (plan != plans_.engine.end()
          && plan->second.kind == molresponse_v3::StatePlan::Kind::load
          && analyses.contains(plan->second.loaded_key))
        keep.push_back(plan->second.loaded_key);
    }
    checkpoint::detail::on_root(world_, [&] {
      for (const auto& a : latest_.archives)
        for (const auto& shard : a.shards) {
          std::filesystem::create_directories(
              (attempt_ / shard.path).parent_path());
          std::filesystem::copy_file(store_ / latest_.generation / shard.path,
                                     attempt_ / shard.path);
        }
      for (const auto& index : indices.items())
        checkpoint::detail::write_exclusive(
            attempt_ / bundle_directory(index.key()) / "roots.json",
            index.value().dump(2));
      for (const auto& key : keep)
        checkpoint::detail::write_exclusive(attempt_ / analysis_file(key),
                                            analyses.at(key).dump(2));
      checkpoint::detail::write_exclusive(attempt_ / "response_metadata.json", metadata.dump(2));
      return std::string{};
    });
    return false;
  }

  const Manifest& latest() const { return latest_; }
  [[nodiscard]] const molresponse_v3::StatePlans& state_plans() const {
    return plans_.engine;
  }
  /// What --continue did per state and why (B6), for results.json.
  [[nodiscard]] json continuation() const {
    return {{"planner", "plan_ladder"}, {"states", plans_.records}};
  }
  json metadata() const { return latest_.numerical.at("engine_metadata"); }
  json properties() const { return metadata().at("properties"); }
  // The full Raman set's normal-mode analysis, where the engine keeps it.
  [[nodiscard]] json analysis(const std::string& key) const {
    return metadata().value(key, json());
  }
  // The engine's transition analysis of the bundle at the final rung, as the
  // last generation kept it (null when there is none).
  [[nodiscard]] json transition() const {
    return latest_.numerical.value("transition_analysis", json::object())
        .value(request_final_key(request_), json());
  }
  // Who decided each bundle's verdict (D-016): the planner for a bundle it
  // loaded and found met, else the engine check (no field). Rank 0's
  // metadata, before the selection reads it.
  void annotate_bundles(json& metadata) const {
    // A loaded generation is used as it was published.
    if (loaded_ || !metadata.contains("excited_states"))
      return;
    for (auto& entry : metadata.at("excited_states"))
      entry.erase("decided_by");
    const auto bundle = requested_bundle(request_);
    if (!bundle)
      return;
    const auto plan = plans_.engine.find(bundle->id);
    if (plan != plans_.engine.end()
        && plan->second.kind == molresponse_v3::StatePlan::Kind::load
        && plan->second.met
        && metadata.at("excited_states").contains(plan->second.loaded_key))
      metadata["excited_states"][plan->second.loaded_key]["decided_by"]
          = "planner";
  }

  // Called collectively after a native FD or bundle save and the metadata
  // update have finished. A bundle's event (D-025) has no frequency.
  // `work_item` is the AppSpec's observe of the event (response_observe).
  void observe(madness::World&, const json& event, const json& work_item) {
    iterations += event.at("iterations").get<std::size_t>();
    const auto& perturbation
        = event.at("perturbation").get_ref<const std::string&>();
    const auto& protocol = event.at("protocol").get_ref<const std::string&>();
    if (event.value("kind", "") == "excited_states")
      totals_.solved.insert(bundle_key(perturbation, protocol));
    else
      totals_.solved.insert(
          state_key(perturbation, protocol,
                    molresponse_v3::ResponseMetadata::freq_key(
                        event.at("frequency").get<double>())));
    work_items.push_back(work_item);
    const auto metadata = read_metadata();
    publish(metadata, false, false, event);
  }

  void finish(const json& root_metadata, bool converged) {
    if (loaded_) return;
    std::string text = world_.rank() == 0 ? root_metadata.dump() : std::string{};
    world_.gop.broadcast_serializable(text, 0);
    work = iterations == 0 ? "analyzed" : restored_ ? "resumed" : "solved";
    publish(json::parse(text), true, converged, {{"phase", "property_analysis"}});
  }

private:
  // The engine's derived analysis: property rows, and the Raman normal-mode
  // analysis and spectra of a full set.
  static void clear_analysis(json& metadata) {
    metadata["properties"] = json::object();
    metadata.erase("vibrational_analysis");
    metadata.erase("raman_spectra");
  }

  json read_metadata() const {
    const auto text = madness::checkpoint::detail::on_root(world_, [&] {
      return madness::checkpoint::detail::read_json(attempt_ / "response_metadata.json").dump();
    });
    return json::parse(text);
  }

  void publish(json metadata, bool analyzed, bool converged, const json& position) {
    namespace checkpoint = madness::checkpoint;
    if (!analyzed)
      clear_analysis(metadata);
    annotate_bundles(metadata);
    // Native saves replace entries but preserve their sample attempt identity.
    // Keep the corresponding request registry across intermediate publications.
    auto origins = metadata.value("checkpoint_requests", json::object());
    if (!latest_.generation.empty() && latest_.reference == reference_) {
      const auto previous = latest_.numerical.at("engine_metadata").value("checkpoint_requests", json::object());
      for (const auto& origin : previous.items()) origins[origin.key()] = origin.value();
    }
    origins[attempt_.string()] = request_;
    metadata["checkpoint_requests"] = std::move(origins);
    Manifest m;
    m.stage = stage_.name;
    m.method = "TDHF";
    m.generation = "g-" + stage_.name + "-" + run_id_ + "-" + std::to_string(++sequence_);
    m.physical = {{"schema", "madqc.response.v1"}, {"shell", "closed_shell"}};
    m.reference = reference_;
    m.numerical = {
        {"request", request_},
        {"attempt_directory", attempt_.string()},
        {"engine_metadata", metadata},
        {"states", response_state_records(
                       metadata, reference_.value("generation", ""), totals_)},
        {"plans", plans_.records}};
    m.targets = {{"dconv", request_.at("dconv")}};
    m.position = position;
    m.position["iterations_this_run"] = iterations;
    if (!latest_.generation.empty()) m.position["previous_generation"] = latest_.generation;
    m.stop_reason = converged ? "converged" : analyzed ? "unfinished" : "in_progress";
    // Keep analyzed values as diagnostics in engine_metadata, but an unmet
    // explicit request cannot publish them as eligible stage properties.
    if (analyzed && (converged || !targets_decide()))
      for (const auto& name : stage_.properties)
        m.properties.push_back(
            {{"name", name},
             {"generation", m.generation},
             {"converged", converged},
             {"value",
              property_value(metadata, name, request_final_key(request_))}});
    GenerationWrite write;
    write.manifest = std::move(m);
    write.engine_metadata = metadata;
    write.check_before_write = [&](const Manifest&) {
      (void)members(metadata);
      if (converged)
        require(request_complete(metadata, request_, natom_,
                                 loaded_states(plans_.records)),
                "analysis lacks converged prerequisites");
    };
    // Writers have already closed collectively. Copy the complete prerequisite
    // set; immutable generations never share writable files with an attempt.
    write.write_archives = [&](const std::filesystem::path& directory,
                               Manifest& m) {
      return checkpoint::detail::on_root(world_, [&] {
        for (const auto& name : members(metadata)) {
          const int writers
              = checkpoint::detail::header(attempt_ / (name + ".00000"), true);
          const auto inventory = checkpoint::inspect_archive(
              attempt_, "response", name, {}, writers);
          for (const auto& shard : inventory.shards) {
            std::filesystem::create_directories(
                (directory / shard.path).parent_path());
            std::filesystem::copy_file(attempt_ / shard.path,
                                       directory / shard.path);
          }
          m.archives.push_back(checkpoint::inspect_archive(
              directory, "response", name, {}, writers));
        }
        // A bundle's index (roots.json) is not an archive: the manifest keeps
        // it, and the transition analysis of a converged bundle with it.
        json indices = json::object(), analyses = json::object();
        for (const auto& rung : metadata.at("excited_states").items()) {
          indices[rung.key()] = checkpoint::detail::read_json(
              attempt_ / bundle_directory(rung.key()) / "roots.json");
          const auto analysis = attempt_ / analysis_file(rung.key());
          if (rung.value().value("converged", false)
              && std::filesystem::exists(analysis))
            analyses[rung.key()] = checkpoint::detail::read_json(analysis);
        }
        m.numerical["bundle_indices"] = std::move(indices);
        if (!analyses.empty())
          m.numerical["transition_analysis"] = std::move(analyses);
        return std::string{};
      });
    };
    latest_ = publish_(std::move(write));
  }
};

// The stage's request (D-019): what its states are planned against and what
// a generation's request must equal to be loaded whole.
inline json response_request(const Stage& stage) {
  const auto& cp = stage.params.get<CalculationParameters>();
  // The request is keyed on the properties the stage names (D-019): each
  // with its frequencies and components, Raman with its atoms and
  // directions as written (empty: every atom, x y z).
  json properties = json::object();
  for (const auto& name : stage.properties) {
    if (name == "polarizability")
      properties[name] = {{"frequencies", stage.frequencies},
                          {"components", stage.components}};
    else if (name == "hyperpolarizability")
      properties[name]
          = {{"frequencies", stage.hyperpolarizability.frequencies},
             {"components", stage.hyperpolarizability.components}};
    else if (name == "raman")
      properties[name] = {{"frequencies", stage.raman.frequencies},
                          {"components", stage.raman.components},
                          {"atoms", stage.raman_atoms},
                          {"directions", stage.raman_directions},
                          {"full_set", stage.raman_full_set}};
    else if (name == "excited_states")
      properties[name] = {{"roots", stage.excited_states.roots},
                          {"method", stage.excited_states.method},
                          {"guess", stage.excited_states.guess},
                          {"guess_basis", stage.excited_states.guess_basis}};
    else
      ResponseCheckpoint::require(false, "unsupported property " + name);
  }
  json request = {{"protocol", cp.protocol()},
                  {"properties", properties},
                  {"dconv", stage.params.get<ResponseParameters>().dconv()}};
  const auto targets = molresponse_v3::fd_requested_accuracy_targets(
      stage.params.get<ResponseParameters>());
  if (!targets.empty())
    request["accuracy_targets"] = targets;
  // Who decides convergence (B8); reuse and property eligibility follow it.
  request["decided_by"] = stage.params.get<ResponseParameters>().get<bool>(
                              madness::convergence::use_default_criteria_key)
                              ? "engine_check"
                              : "requested_criteria";
  const auto progress = molresponse_v3::fd_requested_progress(
      stage.params.get<ResponseParameters>());
  if (!progress.empty())
    request["progress"] = progress;
  const int diagnostic = stage.params.get<ResponseParameters>().get<int>(
      "diagnostic.iterations_per_rung");
  if (diagnostic > 0)
    request["diagnostic_iterations_per_rung"] = diagnostic;
  return request;
}

// B6, the planner's record of every response state in `metadata`, one per
// (perturbation, frequency, rung): identity, reference generation, rung,
// iterations at the rung across runs (recorded, never compared, D-007), how
// it ended, whether it diverged, the deciding mode, and the raw
// measurements (the deciding sample and the saved state, each with the norm
// its relative limit uses). An excited-state bundle is one state per rung
// (D-025), judged on its worst root, with the per-root rows inside.
inline json response_state_records(const json& metadata,
                                   const std::string& reference_generation,
                                   const IterationTotals& totals) {
  json records = json::array();
  // Held in a local: iterating items() of the temporary value() returns
  // would dangle.
  const auto states = metadata.value("fd_states", json::object());
  for (const auto& channel : states.items())
    for (const auto& rung : channel.value().items()) {
      const auto& protocol = metadata.at("protocols").at(rung.key());
      for (const auto& point : rung.value().items()) {
        const auto& e = point.value();
        const auto key = ResponseCheckpoint::state_key(channel.key(),
                                                       rung.key(), point.key());
        const int total
            = iterations_across_runs(totals, key, e.value("iter", 0));
        std::ostringstream frequency;
        frequency << std::setprecision(17) << e.at("freq").get<double>();
        records.push_back(
            {{"schema", "madqc.response.state.v1"},
             {"key", key},
             {"identity", channel.key() + "|" + frequency.str()},
             {"perturbation", channel.key()},
             {"frequency", e.at("freq")},
             {"reference_generation", reference_generation},
             {"rung",
              {{"thresh", protocol.at("thresh")}, {"k", protocol.at("k")}}},
             {"iterations_at_rung", total},
             {"stop", response_stop(e)},
             {"diverged", e.value("diverged", false)},
             {"decided_by", e.value("checkpoint_decided_by", "engine_check")},
             {"engine_measurements",
              e.value("convergence_measurements", json::object())},
             {"saved_measurements",
              e.value("saved_state_validation", json::object())
                  .value("measurements", json::object())}});
      }
    }
  const auto bundles = metadata.value("excited_states", json::object());
  for (const auto& rung : bundles.items()) {
    const auto& e = rung.value();
    const auto& protocol = metadata.at("protocols").at(rung.key());
    const auto id = molresponse_v3::es_node_id(
        e.at("type").get<std::string>() == "tda", e.at("n_roots").get<int>());
    const auto key = ResponseCheckpoint::bundle_key(id, rung.key());
    json roots = json::array();
    for (const auto& root : e.at("roots"))
      roots.push_back(
          {{"root_id", root.value("root_id", "")},
           {"omega", root.value("omega", json())},
           {"converged", root.value("converged", false)},
           {"locked", root.value("locked", false)},
           {"stalled", root.value("stalled", false)},
           {"density_residual", root.value("density_residual", json())},
           {"omega_residual", root.value("omega_residual", json())}});
    records.push_back(
        {{"schema", "madqc.response.state.v1"},
         {"key", key},
         {"kind", "excited_states"},
         {"identity", id},
         {"perturbation", id},
         {"frequency", nullptr},
         {"reference_generation", reference_generation},
         {"rung", {{"thresh", protocol.at("thresh")}, {"k", protocol.at("k")}}},
         {"iterations_at_rung",
          iterations_across_runs(totals, key, e.value("iter", 0))},
         {"stop", e.contains("stop_reason")
                      ? e.at("stop_reason").get<std::string>()
                      : response_stop(e)},
         {"diverged", e.value("diverged", false)},
         {"decided_by", e.value("decided_by", "engine_check")},
         {"engine_measurements", bundle_measurements(e)},
         {"saved_measurements", json::object()},
         {"roots", std::move(roots)}});
  }
  return records;
}

// B6: one plan per state the requested values need (response_needs), from
// the generation's state records. `own` is the state's finest saved rung
// the engine can use. Left out: a rung whose engine key equals a requested
// rung's but whose threshold differs (keys round it; the engine holds one
// definition per key), and a rung tighter than the final one (R3, deferred
// for response, D-018): the analysis reads only states at or coarser than
// the final rung, so such a state is solved again. A state with no usable
// record starts from the nearest saved frequency (R4/R7, frequency_seed).
// Quadratic sources are not planned: the engine builds them from their
// legs and rebuilds one whose leg is solved again. R6 needs nothing here:
// members() never publishes a diverged state, so the generation holds the
// state's converged coarser rung and the next rung restarts from it.
inline ResponsePlans
plan_response_states(const Stage& stage, const json& stage_request, int natom,
                     const std::string& reference_generation,
                     const json& records) {
  ResponsePlans plans;
  std::map<std::string, double> requested_rungs;
  auto request = ResponseCheckpoint::ladder_request(stage, reference_generation,
                                                    requested_rungs);
  for (const auto& [perturbation, frequency] :
       response_needs(stage_request, natom).legs) {
    std::ostringstream identity;
    identity << perturbation.description() << '|' << std::setprecision(17)
             << frequency;
    request.identity = identity.str();
    const auto own = ResponseCheckpoint::own_state(
        records, request.identity, requested_rungs, request.ladder.back());
    const auto seed = own ? std::nullopt
                          : ResponseCheckpoint::frequency_seed(
                                records, perturbation.description(), frequency,
                                request.reference_generation, requested_rungs,
                                request.ladder.front());
    const auto channel = perturbation.description();
    const auto fkey = molresponse_v3::ResponseMetadata::freq_key(frequency);
    ResponseCheckpoint::adopt(
        plans, madness::planning::plan_ladder(own, seed, request), own, seed,
        request, molresponse_v3::fd_node_id(perturbation, frequency),
        {{"perturbation", channel}, {"frequency", frequency}},
        [&](const std::string& protocol) {
          return ResponseCheckpoint::state_key(channel, protocol, fkey);
        });
  }
  // The excited-state bundle (D-025): one state, judged on its worst root
  // against the engine's own check; no seed (the engine cannot start a
  // bundle from another root count or method, DESIGN §5).
  if (const auto bundle = requested_bundle(stage_request)) {
    request.identity = bundle->id;
    request.engine_limits
        = bundle_engine_limits(stage, request.ladder.back().thresh);
    // The written targets are the dipole legs' (B3: excited states have
    // none yet; the engine check decides a bundle).
    for (auto& metric : request.settings.metrics)
      if (metric.group == madness::convergence::Group::requested_criteria)
        metric.target.reset();
    const auto own = ResponseCheckpoint::own_state(
        records, bundle->id, requested_rungs, request.ladder.back());
    auto plan = madness::planning::plan_ladder(own, std::nullopt, request);
    if (!own)
      for (const auto& record : records)
        if (record.value("kind", "") == "excited_states") {
          plan.why += "; the saved bundle is "
                      + record.at("identity").get<std::string>()
                      + ", another method or root count, so this one "
                        "starts fresh";
          break;
        }
    ResponseCheckpoint::adopt(
        plans, plan, own, std::nullopt, request, bundle->id,
        {{"kind", "excited_states"},
         {"method", bundle->tda ? "tda" : "rpa"},
         {"roots", bundle->roots}},
        [&](const std::string& protocol) {
          return ResponseCheckpoint::bundle_key(bundle->id, protocol);
        });
  }
  return plans;
}
} // namespace named_stages
