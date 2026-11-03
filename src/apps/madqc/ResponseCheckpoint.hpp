#pragma once

#include "../molresponse/solvers/function_kinetics.hpp"

#include "ConvergenceReport.hpp"
#include "NamedStages.hpp"
#include "ResponseProperties.hpp"
#include "DaltonSeed.hpp"
#include "Sha256.hpp"
#include <apps/molresponse/solvers/dalton_import.hpp>
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
    // D-099: the bundle's stall detector, B's rows.
    for (auto& row : bundle_stall_rows(stage))
      settings.metrics.push_back(std::move(row));
  }

  const std::vector<cv::RequestedMetric> requested = {
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
       "absolute_kinetic_trace_change", "absolute_kinetic_trace_change",
       "hartree_times_response_amplitude_scale_squared",
       "has no kinetic term"}};
  const auto targets = molresponse_v3::fd_requested_accuracy_targets(rp);
  for (const auto& row : requested)
    settings.metrics.push_back(cv::requested_metric(row, targets));
  // D-061: on an excited-state stage the written econv is the bundle's
  // omega target, in place of the energy change a driven response lacks.
  if (std::find(stage.properties.begin(), stage.properties.end(),
                "excited_states")
      != stage.properties.end())
    settings.metrics.push_back(
        excited_state_omega_target_row(stage, limits.omega_residual));
  else
    settings.metrics.push_back(
        {Group::requested_criteria, "energy change", "", "energy_change", "",
         "hartree", std::nullopt, SetBy::not_set, Relation::less_than_or_equal,
         std::nullopt, "",
         "not applicable: a driven response has no local energy; see the "
         "reference stage",
         ""});
  settings.scope
      = cv::protocol_scope(stage.params.get<CalculationParameters>().protocol(),
                           protocol_k(final_thresh));
  settings.iteration_caps = {
      {static_cast<int>(rp.maxiter())},
      stage.written.contains("maxiter") ? SetBy::user : SetBy::engine_default};
  // A3: no save boundary inside a solve yet (DESIGN §1.4); what it keeps is
  // the one retention's (D-070), checkpoint.keep_all the request for more.
  settings.checkpoints.fixed
      = "at every solve end (molresponse has no save boundary inside a "
        "solve)";
  settings.checkpoints.keep_all = stage.checkpoint_keep_all;
  settings.checkpoints.keep_all_set_by
      = stage.written.contains("checkpoint.keep_all") ? SetBy::user
                                                      : SetBy::engine_default;
  // S5: the Dalton seed, as written (a starting guess for states with no
  // saved state).
  if (!stage.seed_dalton_dir.empty()) {
    std::ostringstream seed;
    seed << "Dalton directory " << stage.seed_dalton_dir << " (user)";
    if (stage.seed_freq_tol > 0.0)
      seed << "; driven legs within " << stage.seed_freq_tol
           << " au of an RSPVEC frequency (seed.freq_tol, user)";
    settings.seed = seed.str();
  }
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
        outcome.units = record.value("units", "");
        if (record.value("phase", json()).is_string())
          outcome.phase = record.at("phase").get<std::string>();
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
  // D-095: the engine's per-leg checks, each leg judged on its own effective
  // limit as the engine judged it, from its last solver step.
  const std::map<std::string, std::string> engine_rows
      = {{"engine_bsh_residual", "bsh_residual"},
         {"engine_density_change", "density_change"}};
  for (const auto& channel : states.items()) {
    if (!channel.value().contains(final_protocol))
      continue;
    for (const auto& point : channel.value().at(final_protocol).items()) {
      const auto checks
          = point.value().value("convergence_measurements", json::object());
      for (const auto& [row, name] : engine_rows) {
        if (!checks.contains(name))
          continue;
        const auto& check = checks.at(name);
        madness::convergence::MetricOutcome leg;
        leg.measurement = row;
        leg.subject = channel.key() + " " + point.key();
        leg.units = check.value("units", "");
        if (check.value("phase", json()).is_string())
          leg.phase = check.at("phase").get<std::string>();
        if (check.value("effective_limit", json()).is_number())
          leg.limit = check.at("effective_limit").get<double>();
        if (check.value("availability", "") == "measured"
            && check.value("value", json()).is_number()) {
          leg.value = check.at("value").get<double>();
          leg.availability = "measured";
        }
        outcomes.push_back(std::move(leg));
      }
    }
  }
  // D-095: the excited-state bundle's checks, per root (B's rows).
  if (metadata.contains("excited_states")
      && metadata.at("excited_states").contains(final_protocol))
    for (auto& leg :
         bundle_root_outcomes(metadata.at("excited_states").at(final_protocol)))
      outcomes.push_back(std::move(leg));
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
// D-113: two_photon's legs are the dipole legs at half each root's
// frequency, so they are known once the bundle is: `bundle_omegas` are the
// roots' omegas at the final rung (empty before the bundle is saved).
[[nodiscard]] inline ResponseNeeds
response_needs(const json& request, int natom,
               const std::vector<double>& bundle_omegas = {}) {
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
    const auto coordinates
        = raman_coordinates(r.at("atoms").get<std::vector<int>>(),
                            r.at("directions").get<std::string>(), natom);
    for (double frequency : r.at("frequencies"))
      for (const auto& [atom, axis] : coordinates) {
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
  // The engine's derived legs (expand_converged_es): es_freq_factor 0.5 times
  // the root's omega, per axis.
  if (properties.contains("two_photon"))
    for (double omega : bundle_omegas)
      for (char axis : {'x', 'y', 'z'})
        leg(dipole(axis), 0.5 * omega);
  return needs;
}

// D-116: what each frequency-dependent state is for, by the requests in
// their written order (as response_needs walks them), so a reader keys a leg
// on its request, not on its frequency. Keyed "<perturbation>@<freq key>" (a
// state record's perturbation and the engine's frequency key); each leg has
// every role it serves (dipole_z at w can serve alpha and beta):
//   driven           alpha's or Raman's leg at w (omega_index), beta's legs at
//                    omega_B and omega_C (n_b, n_c);
//   shg_2w           beta's leg at 2w (n_b, n_c of its driver);
//   nuclear          a Raman nuclear leg (atom, direction);
//   two_photon_half  a two-photon leg at Omega_f/2 (parent_root: the bundle
//                    root's stable root_id, not its slot).
// `bundle_roots`: the bundle's roots at the final rung (root_id, omega).
[[nodiscard]] inline std::map<std::string, json>
state_roles(const json& request, int natom, const json& bundle_roots = json::array()) {
  using molresponse_v3::Perturbation;
  std::map<std::string, json> roles;
  const auto add = [&](const Perturbation& p, double frequency, json role) {
    auto& list = roles[p.description() + "@"
                       + molresponse_v3::ResponseMetadata::freq_key(frequency)];
    if (list.is_null())
      list = json::array();
    if (std::find(list.begin(), list.end(), role) == list.end())
      list.push_back(std::move(role));
  };
  const auto dipole = [](char axis) { return Perturbation::dipole(axis - 'x'); };
  const auto properties = requested_properties(request);
  if (properties.contains("polarizability")) {
    const auto& r = properties.at("polarizability");
    std::set<char> axes;
    for (const auto& component : r.at("components"))
      for (char axis : component.get<std::string>())
        axes.insert(axis);
    const auto& frequencies = r.at("frequencies");
    for (std::size_t i = 0; i < frequencies.size(); ++i)
      for (char axis : axes)
        add(dipole(axis), frequencies[i].get<double>(),
            {{"property", "polarizability"}, {"role", "driven"}, {"omega_index", i}});
  }
  if (properties.contains("hyperpolarizability")) {
    const auto& r = properties.at("hyperpolarizability");
    const auto& frequencies = r.at("frequencies");
    for (std::size_t n = 0; n < frequencies.size(); ++n) {
      const double f = frequencies[n].get<double>();
      for (const auto& component : r.at("components")) {
        const auto abc = component.get<std::string>();
        add(dipole(abc[0]), 2.0 * f,
            {{"property", "hyperpolarizability"},
             {"role", f == 0.0 ? "driven" : "shg_2w"}, {"n_b", n}, {"n_c", n}});
        for (char axis : {abc[1], abc[2]})
          add(dipole(axis), f,
              {{"property", "hyperpolarizability"}, {"role", "driven"},
               {"n_b", n}, {"n_c", n}});
      }
    }
  }
  if (properties.contains("raman")) {
    const auto& r = properties.at("raman");
    const auto coordinates
        = raman_coordinates(r.at("atoms").get<std::vector<int>>(),
                            r.at("directions").get<std::string>(), natom);
    const auto& frequencies = r.at("frequencies");
    for (std::size_t i = 0; i < frequencies.size(); ++i)
      for (const auto& [atom, axis] : coordinates) {
        add(Perturbation::nuclear(atom, axis - 'x'), 0.0,
            {{"property", "raman"}, {"role", "nuclear"}, {"atom", atom},
             {"direction", std::string(1, axis)}});
        for (const auto& component : r.at("components"))
          for (char a : component.get<std::string>())
            add(dipole(a), frequencies[i].get<double>(),
                {{"property", "raman"}, {"role", "driven"}, {"omega_index", i}});
      }
  }
  if (properties.contains("two_photon"))
    for (const auto& root : bundle_roots)
      if (root.contains("omega") && root.at("omega").is_number())
        for (char axis : {'x', 'y', 'z'})
          add(dipole(axis), 0.5 * root.at("omega").get<double>(),
              {{"property", "two_photon"}, {"role", "two_photon_half"},
               {"parent_root", root.value("root_id", json())}});
  return roles;
}

// D-118: what the Dalton seed of a response stage cost, from its engine
// metadata: the import (reading the Dalton files, building the MRA guesses
// and their gauge rotation, once per stage), the bundle's seed guard per rung,
// and each leg's nearest-record seed; the total and how many states each
// part seeded. Null when nothing was seeded. Never part of wall_s.
[[nodiscard]] inline json response_seed_cost(const json& metadata) {
  double imported = -1.0, guards = 0.0, legs = 0.0;
  int guard_count = 0, leg_count = 0;
  if (metadata.contains("seeded_from")
      && metadata.at("seeded_from").contains("import_wall_s"))
    imported = metadata.at("seeded_from").at("import_wall_s").get<double>();
  const auto bundles = metadata.value("excited_states", json::object());
  const auto legs_metadata = metadata.value("fd_states", json::object());
  for (const auto& rung : bundles.items())
    if (rung.value().contains("seed_guard")
        && rung.value().at("seed_guard").contains("wall_s")) {
      guards += rung.value().at("seed_guard").at("wall_s").get<double>();
      ++guard_count;
    }
  for (const auto& channel : legs_metadata.items())
    for (const auto& rung : channel.value().items())
      for (const auto& point : rung.value().items())
        if (point.value().contains("seed_wall_s")) {
          legs += point.value().at("seed_wall_s").get<double>();
          ++leg_count;
        }
  if (imported < 0.0 && guard_count == 0 && leg_count == 0)
    return nullptr;
  return {{"wall_s", std::max(imported, 0.0) + guards + legs},
          {"import_wall_s", imported < 0.0 ? json() : json(imported)},
          {"seed_guard_wall_s", guards},
          {"legs_wall_s", legs},
          {"seeded_states", {{"bundle_rungs", guard_count}, {"legs", leg_count}}}};
}

// D-113: the omegas of the bundle's roots in a state record of the bundle
// (response_state_records), or in the engine metadata's entry for it.
[[nodiscard]] inline std::vector<double> root_omegas(const json& bundle) {
  std::vector<double> omegas;
  for (const auto& root : bundle.value("roots", json::array()))
    if (root.contains("omega") && root.at("omega").is_number())
      omegas.push_back(root.at("omega").get<double>());
  return omegas;
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

// D-082: the written targets a response state's kind has rows for, by the
// state record's kind. The verdict (requested_states) and the planner
// (plan_response_states) judge each state on these only; a row of its kind
// the record lacks stays unavailable.
[[nodiscard]] inline json kind_measurements(const std::string& kind) {
  if (kind == "excited_states")
    return json::array({"excitation_energy_change"});
  return json::array({"charge_density_change", "bsh_residual",
                      "absolute_kinetic_trace_change"});
}

// D-082: a state record's saved-state measurements under the names of the
// written targets, by its kind: an FD state's through fd_accuracy_measurements,
// a bundle's computed from its engine row (bundle_saved_measurements; a
// record without engine_measurements gives none).
[[nodiscard]] inline json record_measurements(const json& record) {
  if (record.value("kind", "") == "excited_states")
    return bundle_saved_measurements(
        record.value("engine_measurements", json::object()));
  return molresponse_v3::fd_accuracy_measurements(
      record.value("saved_measurements", json::object()));
}

// D-082: the request's settings with only the written targets in
// `measurements` kept (the others not set), for the planner.
[[nodiscard]] inline madness::convergence::ConvergenceSettings
targets_of_kind(madness::convergence::ConvergenceSettings settings,
                const json& measurements) {
  for (auto& metric : settings.metrics)
    if (metric.group == madness::convergence::Group::requested_criteria
        && std::find(measurements.begin(), measurements.end(),
                     metric.measurement)
               == measurements.end())
      metric.target.reset();
  return settings;
}

// The engine's name for each named response property's rows.
[[nodiscard]] inline std::string engine_property(const std::string& name) {
  static const std::map<std::string, std::string> names
      = {{"polarizability", "alpha"},
         {"hyperpolarizability", "beta"},
         {"raman", "raman"},
         {"two_photon", "tpa"}};
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
          {"excitation_energy_change", limits.omega_residual, 0.0, false}};
}

// A saved bundle's engine-check sample: per measurement the largest value
// over the roots, and the root it came from; null when a root has none yet.
// The phase is the engine's (D-059): the density change between the step's
// input densities, the eigenvalue change across the step; null with no
// value (D-063). The rows take the record's names (D-106): the engine's
// omega_residual is the excitation_energy_change.
[[nodiscard]] inline json bundle_measurements(const json& bundle) {
  json measurements = json::object();
  struct Row {
    const char *name, *engine, *phase;
  };
  static const std::array<Row, 2> rows = {
      {{"density_residual", "density_residual",
        "solver_step_input_before_acceleration"},
       {"excitation_energy_change", "omega_residual", "solver_step_output"}}};
  for (const auto& [name, engine, phase] : rows) {
    json worst = {{"value", nullptr}, {"availability", "unavailable"}};
    bool every_root = true;
    for (const auto& root : bundle.value("roots", json::array())) {
      const auto value = root.value(engine, json());
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
    worst["phase"] = worst.at("value").is_null() ? json(nullptr) : json(phase);
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

// B6, D-016: the plans as this attempt uses them (response_engine_plans):
// the engine's input per requested state (by its node id), the plan records,
// for the manifest and results.json, and the state keys it resumes at their
// own saved rung (IterationTotals::resumed).
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
[[nodiscard]] inline StagePlans
plan_response_states(const Stage& stage, const json& stage_request, int natom,
                     const std::string& reference_generation,
                     const json& records);
[[nodiscard]] inline ResponsePlans
response_engine_plans(const StagePlans& plans);

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
  json reference_, request_;
  int natom_;
  Manifest latest_;
  std::vector<std::string> attempt_copies_;
  json whole_load_;
  json reference_evidence_;
  // S5: the seed block of a seeded stage (null when unseeded).
  json seed_ = nullptr;
  bool restored_ = false, loaded_ = false;
  // B6: this attempt's iteration totals and plans (S3b: the functions that
  // read them take them as arguments).
  IterationTotals totals_;
  ResponsePlans plans_;
  // The skeleton's publish (publish_generation), from the run binding.
  GenerationPublisher publish_;

  friend json response_state_records(const json& metadata,
                                     const std::string& reference_generation,
                                     const IterationTotals& totals);
  friend StagePlans
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
    // D-113: two_photon's legs from the bundle's roots at the final rung,
    // as the engine derived them; the bundle itself is checked below.
    std::vector<double> omegas;
    if (metadata.contains("excited_states")
        && metadata.at("excited_states").contains(key))
      omegas = root_omegas(metadata.at("excited_states").at(key));
    const auto needs = response_needs(request, natom, omegas);
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
      // D-113: two_photon's values, one row per root at the final rung.
      if (requested_properties(request).contains("two_photon")) {
        const auto rows = metadata.at("properties")
                              .value("tpa", json::object())
                              .value(key, json::array());
        if (!rows.is_array() || rows.size() != omegas.size())
          return false;
      }
      // Its written omega target blocks completion only when targets decide
      // (B8, D-061).
      if (!loaded.contains(bundle_key(bundle->id, key))
          && request.value("decided_by", "") == "requested_criteria"
          && request.contains("excited_state_targets")
          && madness::evaluate_explicit_criteria(
                 bundle_saved_measurements(bundle_measurements(e)),
                 request.at("excited_state_targets"))
                     .at("status")
                 != "satisfied")
        return false;
    }
    return true;
  }

public:
  // D-016, D-025: the bundle verdicts in `metadata` as a generation records
  // them: a bundle the planner loaded as meeting today's request
  // (`planner_key`, null when it loaded none) says decided_by planner, and
  // no other bundle carries a decided_by. A stage that loads as published
  // (`loaded`) keeps its metadata as published.
  static void mark_bundle_verdicts(json& metadata, bool loaded,
                                   const json& planner_key) {
    if (loaded || !metadata.contains("excited_states"))
      return;
    for (auto& entry : metadata.at("excited_states"))
      entry.erase("decided_by");
    if (!planner_key.is_null()
        && metadata.at("excited_states")
               .contains(planner_key.get<std::string>()))
      metadata["excited_states"][planner_key.get<std::string>()]["decided_by"]
          = "planner";
  }

  // A3, the response AppSpec's manifest (S4c3, D-054): every field of a
  // response generation but the shared ones (set_shared_fields) and the
  // archives (with the bundle indices and transition analyses the binding
  // adds while copying them), from the publish sample: the engine's
  // metadata as read, whether it holds the stage's analysis (`analyzed`) and
  // its verdict, the position (the save event, or the analysis), the
  // previous generation (D-066: any generation read or published), the
  // previous generation's request registry when its reference is today's,
  // the attempt directory, the reference, this attempt's iteration totals,
  // the plan records, and the bundle verdicts' inputs. The request is the
  // stage's (response_request).
  static Manifest manifest(const Stage& stage, int natom, const json& sample) {
    const auto request = response_request(stage);
    const bool analyzed = sample.at("analyzed").get<bool>();
    const bool converged = sample.at("converged").get<bool>();
    auto metadata = sample.at("metadata");
    if (!analyzed)
      clear_analysis(metadata);
    const auto& verdicts = sample.at("bundle_verdicts");
    mark_bundle_verdicts(metadata, verdicts.at("loaded").get<bool>(),
                         verdicts.at("planner_key"));
    // Native saves replace entries but preserve their sample attempt
    // identity. Keep the corresponding request registry across
    // intermediate publications.
    auto origins = metadata.value("checkpoint_requests", json::object());
    const auto& previous = sample.at("previous_checkpoint_requests");
    if (!previous.is_null())
      for (const auto& origin : previous.items())
        origins[origin.key()] = origin.value();
    origins[sample.at("attempt_directory").get<std::string>()] = request;
    metadata["checkpoint_requests"] = std::move(origins);
    const auto& reference = sample.at("reference");
    const auto& totals_record = sample.at("totals");
    IterationTotals totals;
    totals.previous
        = totals_record.at("previous").get<std::map<std::string, int>>();
    totals.resumed = totals_record.at("resumed").get<std::set<std::string>>();
    totals.solved = totals_record.at("solved").get<std::set<std::string>>();
    Manifest m;
    m.method = "TDHF";
    m.physical = {{"schema", "madqc.response.v1"}, {"shell", "closed_shell"}};
    m.reference = reference;
    m.numerical
        = {{"request", request},
           {"attempt_directory", sample.at("attempt_directory")},
           {"engine_metadata", metadata},
           {"states", response_state_records(
                          metadata, reference.value("generation", ""), totals)},
           {"plans", sample.at("plans")}};
    // D-116: each frequency-dependent state's roles in the request.
    {
      const auto pkey = request_final_key(request);
      const auto bundle_roots
          = metadata.contains("excited_states")
                    && metadata.at("excited_states").contains(pkey)
                ? metadata.at("excited_states").at(pkey).value("roots", json::array())
                : json::array();
      const auto roles = state_roles(request, natom, bundle_roots);
      for (auto& state : m.numerical.at("states"))
        if (state.value("kind", "") == "frequency_dependent_response") {
          const auto found = roles.find(
              state.at("perturbation").get<std::string>() + "@"
              + molresponse_v3::ResponseMetadata::freq_key(
                  state.at("frequency").get<double>()));
          if (found != roles.end())
            state["roles"] = found->second;
        }
    }
    m.targets = {{"dconv", request.at("dconv")}};
    // S5: a seeded stage's seed block: the files its generations hold, and
    // what the engine recorded seeding from.
    if (sample.contains("seed")) {
      m.seed = manifest_seed_block(sample.at("seed"));
      if (metadata.contains("seeded_from"))
        m.seed["seeded_from"] = metadata.at("seeded_from");
    }
    m.position = sample.at("position");
    if (!sample.at("previous_generation").is_null())
      m.position["previous_generation"] = sample.at("previous_generation");
    m.stop_reason = converged  ? "converged"
                    : analyzed ? "unfinished"
                               : "in_progress";
    // D-098 step 1: the scientific report, in the form moldft's manifest
    // keeps it (metrics.convergence_report): the settings, and for the
    // analysed generation the outcome rows and each needed state's written
    // targets judged (the shared verdict's, when the publish carried it).
    namespace cv = madness::convergence;
    const auto final_thresh = request.at("protocol").back().get<double>();
    const auto settings = response_convergence_settings(stage, final_thresh);
    json report = {{"schema", "madqc.convergence.v1"},
                   {"decided_by", cv::to_string(settings.decider)},
                   {"convergence_settings", cv::to_json(settings)}};
    if (analyzed)
      report["convergence_outcome"] = cv::outcome_record(
          settings,
          response_saved_outcomes(metadata, request_final_key(request)));
    // D-105: the assembled report body, once, here: what the stage knew of
    // its reference, and each leg's measurements and histories (the
    // generation names itself, so no state_id is copied).
    if (sample.contains("reference_evidence"))
      report["reference"] = sample.at("reference_evidence");
    auto body = response_convergence_report(metadata, "", json());
    for (auto& channel : body["channels"]) {
      channel.erase("state_id");
      if (channel.contains("saved_state_validation"))
        channel["saved_state_validation"].erase("state_id");
    }
    report["channels"] = std::move(body["channels"]);
    if (sample.contains("verdict")) {
      json accuracy = json::object();
      for (const auto& state : sample.at("verdict").at("states"))
        if (state.contains("requested_accuracy"))
          accuracy[state.at("key").get<std::string>()]
              = state.at("requested_accuracy");
      report["requested_accuracy"] = std::move(accuracy);
    }
    m.metrics["convergence_report"] = std::move(report);
    // Keep analyzed values as diagnostics in engine_metadata, but an unmet
    // explicit request cannot publish them as eligible stage properties.
    const bool targets_decide
        = request.value("decided_by", "") == "requested_criteria";
    if (analyzed && (converged || !targets_decide))
      for (const auto& name : stage.properties)
        m.properties.push_back(
            {{"name", name},
             {"converged", converged},
             {"value",
              property_value(metadata, name, request_final_key(request))}});
    return m;
  }

  // B8 (S4d1), the response AppSpec's requested_states: the FD states the
  // request needs (response_needs) and its excited-state bundle, at the
  // request's final protocol, as request_complete lists them; each with the
  // written targets its kind has rows for (kind_measurements, D-082).
  // D-113: two_photon's legs are the bundle's: their frequencies are half
  // the roots' omegas in the bundle's record at the final rung, so every leg
  // the engine solved for a root is judged.
  static json requested_states(const Stage& stage, int natom,
                               const json& records) {
    const auto request = response_request(stage);
    const double threshold = request.at("protocol").back();
    const int k = protocol_k(threshold);
    const auto key = molresponse_v3::protocol_key(threshold, k);
    const json rung = {{"thresh", threshold}, {"k", k}};
    std::vector<double> omegas;
    if (const auto bundle = requested_bundle(request))
      for (const auto& record : records)
        if (record.at("key") == bundle_key(bundle->id, key))
          omegas = root_omegas(record);
    json states = json::array();
    for (const auto& [perturbation, frequency] :
         response_needs(request, natom, omegas).legs)
      states.push_back(
          {{"key",
            state_key(perturbation.description(), key,
                      molresponse_v3::ResponseMetadata::freq_key(frequency))},
           {"protocol_step", rung},
           {"measurements",
            kind_measurements("frequency_dependent_response")}});
    if (const auto bundle = requested_bundle(request))
      states.push_back({{"key", bundle_key(bundle->id, key)},
                        {"protocol_step", rung},
                        {"measurements", kind_measurements("excited_states")}});
    return states;
  }

  // B1, B8 (S4d1), the response AppSpec's saved_measurements: a state's
  // saved-state measurements under the names of the written targets, by its
  // kind (record_measurements, as the planner reads them; D-082).
  static json saved_measurements(const Stage&, int, const json& record) {
    return record_measurements(record);
  }

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
    // S5: the seed block's form (seed_block_form_problem).
    if (m.seed.is_object()) {
      const auto problem = seed_block_form_problem(m.seed);
      require(problem.empty(), problem);
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
    // S4a: state records, when the generation has them (D-052: the same
    // check on both sides; publish_generation requires them).
    if (m.numerical.contains("states"))
      validate_state_records(m.numerical.at("states"));
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
      auto saved = saved_state_from_record(record);
      saved.saved_measurements = record_measurements(record);
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
      auto saved = saved_state_from_record(record);
      saved.saved_measurements = record_measurements(record);
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
  // B6: one state's plan, as the field returns it (PlannedState) and as
  // results.json records it; the engine's input is built from it by
  // response_engine_plans. `key_at` names the state key at a rung.
  static void
  adopt(StagePlans& plans, const madness::planning::LadderPlan& plan,
        const std::optional<madness::planning::SavedState>& own,
        const std::optional<madness::planning::SavedState>& seed,
        const madness::planning::LadderRequest& request,
        const std::string& node_id, const json& fields,
        const std::function<std::string(const std::string&)>& key_at) {
    namespace planning = madness::planning;
    PlannedState state;
    state.identity = request.identity;
    state.plan = plan;
    state.start = request.ladder.at(plan.start_rung);
    if (own)
      state.own_rung = own->rung;
    if (seed)
      state.seed_frequency = seed->frequency;
    state.fields = {{"node", node_id},
                    {"key_at_start", key_at(molresponse_v3::protocol_key(
                                         state.start.thresh, state.start.k))}};
    json record = continuation_record(state, request.ladder);
    record["state"] = request.identity;
    record.update(fields);
    // R2 loads a stalled final rung with the request unmet: not met.
    if (plan.action == planning::Action::load && plan.unmet.empty())
      record["loaded_state"]
          = key_at(molresponse_v3::protocol_key(own->rung.thresh, own->rung.k));
    else if (plan.action == planning::Action::seed) {
      record["seed_from"] = plan.seed_from;
      record["seed_frequency"] = seed->frequency;
    }
    plans.continuation["states"].push_back(std::move(record));
    plans.states.push_back(std::move(state));
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
                     const std::filesystem::path& attempt, json reference,
                     int natom, GenerationPublisher publish)
      : world_(world), stage_(stage), store_(stage_directory / "checkpoints"),
        attempt_(attempt), reference_(std::move(reference)), natom_(natom),
        publish_(std::move(publish)) {
    request_ = response_request(stage);
    if (!stage.seed_dalton_dir.empty())
      seed_ = json::parse(madness::checkpoint::detail::on_root(world_, [&] {
        return dalton_seed_block(stage.seed_dalton_dir).dump();
      }));
  }

  [[nodiscard]] bool targets_decide() const {
    return request_.value("decided_by", "") == "requested_criteria";
  }

  // The generation --continue read, validated (validate_generation): return
  // Whether its states are planned (the plan's context): not when the
  // reference changed, and not when the requested analysis itself is already
  // complete (the stage loads as published).
  bool accept(const Manifest& latest) {
    latest_ = latest;
    // S5: a generation's seed copies are checked before it is used.
    if (latest_.seed.is_object()) {
      const auto problem = madness::checkpoint::detail::on_root(world_, [&] {
        return check_seed_files(store_ / latest_.generation, latest_.seed);
      });
      require(problem.empty(), problem);
    }
    reference_changed = latest_.reference != reference_;
    if (reference_changed) return false; // no old descendants enter this attempt
    restored_ = true;
    loaded_ = stage_.params.get<ResponseParameters>().get<int>("diagnostic.iterations_per_rung") == 0 &&
              latest_.stop_reason == "converged" &&
              latest_.numerical.at("request") == request_ && !latest_.properties.empty();
    if (loaded_) {
      work = "loaded";
      // D-094: what --continue does with a met request, and why.
      whole_load_
          = {{"planner", "request_complete"},
             {"action", "loaded"},
             {"generation", latest_.generation},
             {"why", "it holds today's request: the request is equal, it "
                     "stopped converged and its properties are present"},
             {"summary", "loaded whole from " + latest_.generation}};
    }
    return !loaded_;
  }

  // The accepted generation, with its state records and their plans when
  // they were planned: returns true only when the requested analysis itself
  // is already complete. Otherwise restore prerequisites privately and let
  // the engine schedule work.
  bool load(const json& records, const std::optional<StagePlans>& plans) {
    namespace checkpoint = madness::checkpoint;
    if (!plans)
      return loaded_;

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
    // B6: the generation's state records (the AppSpec's state_records): the
    // totals so far; and the plans, as the engine takes them.
    for (const auto& record : records)
      totals_.previous[record.at("key").get<std::string>()]
          = record.at("iterations_at_protocol_step").get<int>();
    plans_ = response_engine_plans(*plans);
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
    // D-090: every shard copied in is deleted at the end of the stage, a
    // root or state the request no longer names included.
    for (const auto& a : latest_.archives)
      for (const auto& shard : a.shards)
        attempt_copies_.push_back(shard.path);
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
  // D-094: the record of a stage loaded whole (null otherwise).
  [[nodiscard]] const json& whole_load() const { return whole_load_; }
  // D-105: what this stage records about its reference, kept in each
  // manifest's report.
  void set_reference_evidence(json evidence) {
    reference_evidence_ = std::move(evidence);
  }
  // D-090: the archive shards the load copied into the attempt.
  [[nodiscard]] const std::vector<std::string>& attempt_copies() const {
    return attempt_copies_;
  }
  [[nodiscard]] const molresponse_v3::StatePlans& state_plans() const {
    return plans_.engine;
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
    mark_bundle_verdicts(metadata, loaded_, planner_bundle_key());
  }

  // Called collectively after a native FD or bundle save and the metadata
  // update have finished. A bundle's event (D-025) has no frequency.
  // `work_item` is the AppSpec's observe of the event (response_observe);
  // `decision`, its publish_now (S4c2: every event).
  void observe(madness::World&, const json& event, const json& work_item,
               const PublishDecision& decision) {
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
    if (!decision.publish)
      return;
    const auto metadata = read_metadata();
    publish(metadata, false, false, event);
  }

  // The analysed generation; the skeleton forms its verdict from its state
  // records and the publish takes it (D-077).
  void finish(const json& root_metadata) {
    if (loaded_) return;
    std::string text = world_.rank() == 0 ? root_metadata.dump() : std::string{};
    world_.gop.broadcast_serializable(text, 0);
    work = iterations == 0 ? "analyzed" : restored_ ? "resumed" : "solved";
    publish(json::parse(text), true, false, {{"phase", "property_analysis"}},
            true);
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

  // The bundle the planner loaded as meeting today's request, as a key of
  // the engine metadata's excited_states; null when it loaded none.
  [[nodiscard]] json planner_bundle_key() const {
    const auto bundle = requested_bundle(request_);
    if (!bundle)
      return nullptr;
    const auto plan = plans_.engine.find(bundle->id);
    if (plan != plans_.engine.end()
        && plan->second.kind == molresponse_v3::StatePlan::Kind::load
        && plan->second.met)
      return plan->second.loaded_key;
    return nullptr;
  }

  void publish(json metadata, bool analyzed, bool converged,
               const json& position, bool final = false) {
    namespace checkpoint = madness::checkpoint;
    GenerationWrite write;
    // The publish sample (the AppSpec's manifest,
    // ResponseCheckpoint::manifest).
    const bool same_reference
        = !latest_.generation.empty() && latest_.reference == reference_;
    write.publish_sample = {
        {"metadata", std::move(metadata)},
        {"analyzed", analyzed},
        {"converged", converged},
        {"position", position},
        {"previous_generation",
         latest_.generation.empty() ? json(nullptr) : json(latest_.generation)},
        {"previous_checkpoint_requests",
         same_reference ? latest_.numerical.at("engine_metadata")
                              .value("checkpoint_requests", json::object())
                        : json(nullptr)},
        {"attempt_directory", attempt_.string()},
        {"reference", reference_},
        {"totals",
         {{"previous", totals_.previous},
          {"resumed", totals_.resumed},
          {"solved", totals_.solved}}},
        {"plans", plans_.records},
        {"bundle_verdicts",
         {{"loaded", loaded_}, {"planner_key", planner_bundle_key()}}}};
    if (!reference_evidence_.is_null())
      write.publish_sample["reference_evidence"] = reference_evidence_;
    if (!seed_.is_null())
      write.publish_sample["seed"] = seed_;
    write.iterations_this_run = iterations;
    write.final = final;
    write.check_before_write = [&](const Manifest& m) {
      const auto& metadata = m.numerical.at("engine_metadata");
      (void)members(metadata);
      if (m.stop_reason == "converged")
        require(request_complete(metadata, request_, natom_,
                                 loaded_states(plans_.records)),
                "analysis lacks converged prerequisites");
    };
    // Writers have already closed collectively. Copy the complete prerequisite
    // set; immutable generations never share writable files with an attempt.
    write.write_archives = [&](const std::filesystem::path& directory,
                               Manifest& m) {
      return checkpoint::detail::on_root(world_, [&] {
        const auto metadata = m.numerical.at("engine_metadata");
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
        // S5: every generation of a seeded stage holds its seed files; the
        // copies are checked here, so a generation without them stays
        // unpublished (D-052).
        if (!seed_.is_null())
          return copy_seed_files(directory, seed_);
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
    else if (name == "two_photon")
      properties[name] = json::object(); // D-113: every root of the bundle
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
  // D-061: the bundle's written omega target, beside the legs' (each leg is
  // judged on every accuracy_targets entry, request_complete). Only when
  // written, so a request without econv is the one stored before.
  if (const auto& omega = stage.excited_states.omega_target)
    request["excited_state_targets"]
        = {{"excitation_energy_change",
            {{"value", *omega}, {"units", "hartree"}}}};
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

// D-059, D-063: a state record's measurement rows, each with its phase; the
// engine writes one on every row it sampled, and a row with no sample (the
// not-applicable energy change, the placeholder kinetic change) gets null.
[[nodiscard]] inline json with_phases(json rows) {
  for (auto& row : rows)
    if (row.is_object() && !row.contains("phase") && !row.contains("sample"))
      row["phase"] = nullptr;
  return rows;
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
        records.push_back(
            {{"schema", state_schema},
             {"key", key},
             {"kind", "frequency_dependent_response"},
             {"identity", frequency_dependent_identity(
                              channel.key(), e.at("freq").get<double>())},
             {"perturbation", channel.key()},
             {"frequency", e.at("freq")},
             {"reference_generation", reference_generation},
             {"protocol_step",
              {{"thresh", protocol.at("thresh")}, {"k", protocol.at("k")}}},
             {"iterations_at_protocol_step", total},
             {"stop_reason", response_stop(e)},
             {"diverged", e.value("diverged", false)},
             {"decided_by", e.value("checkpoint_decided_by", "engine_check")},
             {"engine_measurements",
              with_phases(e.value("convergence_measurements", json::object()))},
             {"saved_measurements",
              with_phases(e.value("saved_state_validation", json::object())
                              .value("measurements", json::object()))}});
        // D-118: a leg seeded from Dalton carries its seed's own time.
        if (e.contains("seed_wall_s"))
          records.back()["seed_wall_s"] = e.at("seed_wall_s");
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
           {"excitation_energy_change", root.value("omega_residual", json())}});
    records.push_back(
        {{"schema", state_schema},
         {"key", key},
         {"kind", "excited_states"},
         {"identity", id},
         {"perturbation", id},
         {"frequency", nullptr},
         {"reference_generation", reference_generation},
         {"protocol_step",
          {{"thresh", protocol.at("thresh")}, {"k", protocol.at("k")}}},
         {"iterations_at_protocol_step",
          iterations_across_runs(totals, key, e.value("iter", 0))},
         {"stop_reason", e.contains("stop_reason")
                             ? e.at("stop_reason").get<std::string>()
                             : response_stop(e)},
         {"diverged", e.value("diverged", false)},
         {"decided_by", e.value("decided_by", "engine_check")},
         {"engine_measurements", bundle_measurements(e)},
         {"saved_measurements", json::object()},
         {"roots", std::move(roots)}});
    // D-118: a seeded bundle carries its seed guard's time at this rung.
    if (e.contains("seed_guard") && e.at("seed_guard").contains("wall_s"))
      records.back()["seed_wall_s"] = e.at("seed_guard").at("wall_s");
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
inline StagePlans plan_response_states(const Stage& stage,
                                       const json& stage_request, int natom,
                                       const std::string& reference_generation,
                                       const json& records) {
  StagePlans plans;
  plans.continuation = {{"planner", "plan_ladder"}, {"states", json::array()}};
  std::map<std::string, double> requested_rungs;
  auto request = ResponseCheckpoint::ladder_request(stage, reference_generation,
                                                    requested_rungs);
  // D-082: each state is planned on the written targets its kind has rows
  // for, the set the shared verdict judges (requested_states).
  const auto stage_settings = request.settings;
  request.settings = targets_of_kind(
      stage_settings, kind_measurements("frequency_dependent_response"));
  // S5 (R4's wording): a state that starts fresh in a seeded stage takes its
  // starting guess from the Dalton seed; the plan says so.
  const auto seeded = [&](madness::planning::LadderPlan plan) {
    if (!stage.seed_dalton_dir.empty()
        && plan.action == madness::planning::Action::fresh)
      plan.why += "; its starting guess from the Dalton seed (seed_source: "
                  "dalton)";
    return plan;
  };
  // One frequency-dependent leg's plan; `leg_request` is the legs' request
  // (the bundle's section below changes `request`).
  const auto leg_request = request;
  std::set<std::string> planned_legs;
  const auto plan_leg = [&](const molresponse_v3::Perturbation& perturbation,
                            double frequency) {
    if (!planned_legs.insert(molresponse_v3::fd_node_id(perturbation, frequency))
             .second)
      return;
    auto request = leg_request;
    request.identity
        = frequency_dependent_identity(perturbation.description(), frequency);
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
        plans, seeded(madness::planning::plan_ladder(own, seed, request)), own,
        seed, request, molresponse_v3::fd_node_id(perturbation, frequency),
        {{"perturbation", channel}, {"frequency", frequency}},
        [&](const std::string& protocol) {
          return ResponseCheckpoint::state_key(channel, protocol, fkey);
        });
  };
  for (const auto& [perturbation, frequency] :
       response_needs(stage_request, natom).legs)
    plan_leg(perturbation, frequency);
  // The excited-state bundle (D-025): one state, judged on its worst root
  // against the engine's own check; no seed (the engine cannot start a
  // bundle from another root count or method, DESIGN §5).
  if (const auto bundle = requested_bundle(stage_request)) {
    request.identity = bundle->id;
    request.engine_limits
        = bundle_engine_limits(stage, request.ladder.back().thresh);
    // The bundle's written target is its omega row (econv, D-061).
    request.settings
        = targets_of_kind(stage_settings, kind_measurements("excited_states"));
    const auto own = ResponseCheckpoint::own_state(
        records, bundle->id, requested_rungs, request.ladder.back());
    auto plan
        = seeded(madness::planning::plan_ladder(own, std::nullopt, request));
    // D-099 (R-063): a bundle loaded stalled at its final rung (R2) is done at
    // that rung for the engine; the protocol and dconv do not move it.
    if (plan.action == madness::planning::Action::load && !plan.unmet.empty()) {
      std::string stalled;
      for (const auto& record : records)
        if (record.value("identity", "") == bundle->id
            && record.value("kind", "") == "excited_states"
            && madness::scf_continuation::same_rung(
                {record.at("protocol_step").at("thresh").get<double>(),
                 record.at("protocol_step").at("k").get<int>()},
                request.ladder.back()))
          for (const auto& root : record.value("roots", json::array()))
            if (root.value("stalled", false))
              stalled += (stalled.empty() ? "" : ", ")
                         + root.value("root_id", std::string("?"));
      plan.why = "the bundle stalled on the engine check ("
                 + (stalled.empty() ? std::string("a root") : "root " + stalled)
                 + ": its density change and its eigenvalue change both flat "
                   "over the stall detector's window); to go on, request more "
                   "roots when the last root borders a cluster, or write "
                   "stall.window 0 or a larger window, or raise maxiter; the "
                   "protocol and dconv do not move a stalled bundle";
    }
    if (!own)
      for (const auto& record : records)
        if (record.value("kind", "") == "excited_states") {
          plan.why += "; the saved bundle is "
                      + excited_state_words(
                          record.at("identity").get<std::string>())
                      + ", another method or root count, so this one "
                        "starts fresh";
          break;
        }
    // D-113: two_photon's legs are at half the roots' omegas. A loaded
    // bundle keeps its saved omegas, so its legs are planned here (loaded,
    // resumed or seeded from the nearest saved frequency). A bundle solved in
    // this run moves its omegas; the engine adds its legs once it converges
    // (expand_converged_es) and plans them by its own restart rule, and the
    // shared verdict judges them from the saved bundle's roots.
    const bool two_photon
        = requested_properties(stage_request).contains("two_photon");
    if (two_photon && plan.action != madness::planning::Action::load)
      plan.why += "; its two-photon legs are planned by the engine once it "
                  "converges";
    ResponseCheckpoint::adopt(
        plans, plan, own, std::nullopt, request, bundle->id,
        {{"kind", "excited_states"},
         {"method", bundle->tda ? "tda" : "rpa"},
         {"roots", bundle->roots}},
        [&](const std::string& protocol) {
          return ResponseCheckpoint::bundle_key(bundle->id, protocol);
        });
    if (two_photon && plan.action == madness::planning::Action::load)
      for (const auto& record : records)
        if (record.value("identity", "") == bundle->id
            && record.value("kind", "") == "excited_states"
            && madness::scf_continuation::same_rung(
                {record.at("protocol_step").at("thresh").get<double>(),
                 record.at("protocol_step").at("k").get<int>()},
                request.ladder.back()))
          for (double omega : root_omegas(record))
            for (int axis = 0; axis < 3; ++axis)
              plan_leg(molresponse_v3::Perturbation::dipole(axis), 0.5 * omega);
  }
  return plans;
}

// The response binding's half of the plans (S4b): each state's plan as the
// engine's input (molresponse_v3::StatePlan), the state keys it resumes, and
// the plan records.
[[nodiscard]] inline ResponsePlans
response_engine_plans(const StagePlans& plans) {
  namespace planning = madness::planning;
  using Kind = molresponse_v3::StatePlan::Kind;
  ResponsePlans out;
  out.records = plans.continuation.at("states");
  for (const auto& state : plans.states) {
    molresponse_v3::StatePlan engine_plan;
    engine_plan.start_thresh = state.start.thresh;
    switch (state.plan.action) {
    case planning::Action::load:
      engine_plan.kind = Kind::load;
      engine_plan.loaded_key = molresponse_v3::protocol_key(
          state.own_rung->thresh, state.own_rung->k);
      // R2 loads a stalled final rung with the request unmet: not met.
      engine_plan.met = state.plan.unmet.empty();
      break;
    case planning::Action::resume:
      engine_plan.kind = Kind::resume;
      if (madness::scf_continuation::same_rung(*state.own_rung, state.start))
        out.resumed.insert(state.fields.at("key_at_start").get<std::string>());
      break;
    case planning::Action::seed:
      engine_plan.kind = Kind::seed;
      engine_plan.seed_frequency = *state.seed_frequency;
      break;
    case planning::Action::fresh:
      engine_plan.kind = Kind::fresh;
      break;
    }
    out.engine[state.fields.at("node").get<std::string>()] = engine_plan;
  }
  return out;
}

// B6, the response AppSpec's state_records (S4b): the generation's records,
// or, for one written before them, the records its engine metadata gives,
// against its own reference (the one this attempt uses: a generation of
// another reference is not planned, ResponseCheckpoint::accept) and with no
// iterations of this attempt yet.
[[nodiscard]] inline json
response_generation_records(const Stage&, int,
                            const madness::checkpoint::Manifest& manifest) {
  auto records = manifest.numerical.value("states", json::array());
  if (records.empty())
    records = response_state_records(manifest.numerical.at("engine_metadata"),
                                     manifest.reference.value("generation", ""),
                                     {});
  // D-106: older generations' records in the current keys.
  return upgrade_state_records(records);
}
} // namespace named_stages
