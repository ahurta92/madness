#pragma once
#include "NamedStages.hpp"
#include <apps/molresponse/calc/calc_manager.hpp>
#include <madness/chem/ConvergenceSettings.h>
#include <nlohmann/json.hpp>

#include <array>

// Select what a molresponse stage requested out of the engine's property rows.
// The response calculation has already run; nothing here computes a number.
// One function per property (D-017): the engine writes one row per tensor
// element it assembled under properties.<alpha|beta|raman>.<protocol key>,
// and the stage asks for frequencies and Cartesian components. A requested
// element the engine did not write is never a silent gap: alpha treats it as
// an error; a quadratic element (beta, Raman) is published without a value
// and with the engine's reason, since its source needs two converged legs
// (D-024). The row an element came from travels along as engine_evidence.
namespace named_stages {
using json = nlohmann::json;

// Select requested frequencies/components from the engine's final-rung tensor.
// The response calculation has already run; this function only formats results.
inline json selected_alpha(const Stage& stage, const json& properties) {
  const auto& cp = stage.params.get<CalculationParameters>();
  const auto pkey = molresponse_v3::protocol_key(
      cp.protocol().back(), protocol_k(cp.protocol().back()));
  const auto& rows = properties.at("alpha").at(pkey);
  json selected = json::array();
  for (double omega : stage.frequencies) {
    const json* found = nullptr;
    for (const auto& row : rows) {
      if (std::abs(row.at("omega").get<double>() - omega) < 1e-12) {
        if (found)
          throw std::runtime_error(
              "duplicate alpha frequency in engine results");
        found = &row;
      }
    }
    if (!found)
      throw std::runtime_error("missing requested polarizability frequency");
    const std::string dirs = found->at("directions");
    json values = json::object();
    for (const auto& component : stage.components) {
      const auto i = dirs.find(component[0]), j = dirs.find(component[1]);
      if (i == std::string::npos || j == std::string::npos)
        throw std::runtime_error("missing requested alpha component");
      const double value = found->at("alpha").at(i).at(j).get<double>();
      if (!std::isfinite(value))
        throw std::runtime_error("nonfinite polarizability");
      values[component] = value;
    }
    selected.push_back({{"frequency", omega},
                        {"frequency_units", "au"},
                        {"units", "a0^3"},
                        {"components", values},
                        {"converged", found->value("converged", false)},
                        {"engine_evidence", *found}});
  }
  return selected;
}

// The protocol key of the stage's final rung: the engine files its rows by it.
inline std::string final_protocol_key(const Stage& stage) {
  const auto& cp = stage.params.get<CalculationParameters>();
  return molresponse_v3::protocol_key(cp.protocol().back(),
                                      protocol_k(cp.protocol().back()));
}

// The engine's rows of one quadratic property at the stage's final rung, or
// no rows when it assembled none there (every source undelivered: D-024).
inline const json& quadratic_rows(const json& properties,
                                  const std::string& property,
                                  const std::string& pkey) {
  static const json none = json::array();
  if (!properties.is_object() || !properties.contains(property))
    return none;
  const auto& by_rung = properties.at(property);
  if (!by_rung.is_object() || !by_rung.contains(pkey))
    return none;
  return by_rung.at(pkey);
}

// The one engine row of a quadratic property (beta or d(alpha)/dQ) with output
// axis A and inputs B and C at (freq_b, freq_c), or nullptr when the engine
// wrote none (D-024: the caller publishes the element as unavailable). The
// engine assembles each element once; two rows would mean two assemblies,
// which is an error here, not a choice.
inline const json* quadratic_row(const json& rows, const std::string& what,
                                 char a, const std::string& b,
                                 const std::string& c, double freq_b,
                                 double freq_c) {
  const json* found = nullptr;
  for (const auto& row : rows) {
    if (row.at("A").get<std::string>() != std::string(1, a)
        || row.at("B").get<std::string>() != b
        || row.at("C").get<std::string>() != c
        || std::abs(row.at("freq_b").get<double>() - freq_b) >= 1e-12
        || std::abs(row.at("freq_c").get<double>() - freq_c) >= 1e-12)
      continue;
    if (found)
      throw std::runtime_error("duplicate " + what + " row in engine results");
    found = &row;
  }
  if (!found)
    return nullptr;
  if (!std::isfinite(found->at("beta").get<double>()))
    throw std::runtime_error("nonfinite " + what);
  return found;
}

// Why the engine wrote no row for an element built from the quadratic source
// (B, C) at (freq_b, freq_c). The run audits the planned sources it never
// delivered into run_summary.dropped_work (calc_executor.hpp, one entry per
// source with the engine's reason), keyed by the source's node id
// (vbc_node_id, calc_manager.hpp): its reason if it names the source, else
// the fact that it does not.
inline std::string quadratic_unavailable(const json& metadata,
                                         const std::string& b,
                                         const std::string& c, double freq_b,
                                         double freq_c) {
  using molresponse_v3::ResponseMetadata;
  const std::string source = "vbc:" + b + "__" + c + "@"
                             + ResponseMetadata::freq_key(freq_b) + "_"
                             + ResponseMetadata::freq_key(freq_c);
  if (metadata.is_object() && metadata.contains("run_summary")) {
    const auto& summary = metadata.at("run_summary");
    if (summary.is_object() && summary.contains("dropped_work"))
      for (const auto& dropped : summary.at("dropped_work"))
        if (dropped.value("id", "") == source)
          return "the engine built no quadratic source " + source + ": "
                 + dropped.value("reason", "no reason recorded");
  }
  return "the engine wrote no row and run_summary.dropped_work does not name "
         "its quadratic source "
         + source;
}

// The requested elements of one quadratic record under construction: each
// element's value and engine row, the reason for each element without a row,
// and whether every element's legs met their targets.
struct QuadraticElements {
  json values = json::object();
  json evidence = json::object();
  json unavailable = json::object();
  bool converged = true;
};
// One requested element: its value and row, or (D-024) a null value, the
// engine's reason under components_unavailable, and an unconverged record.
inline void select_quadratic_element(const json* row, const std::string& name,
                                     const std::string& reason_if_missing,
                                     QuadraticElements& elements) {
  if (row != nullptr) {
    elements.values[name] = row->at("beta").get<double>();
    elements.converged = elements.converged && row->value("converged", false);
    elements.evidence[name] = *row;
    return;
  }
  elements.values[name] = nullptr;
  elements.unavailable[name] = reason_if_missing;
  elements.converged = false;
}

// The hyperpolarizability: for each driver frequency w and triple ABC, the
// engine's SHG row beta_ABC(-2w; w, w) (static beta at w = 0). The value is
// in atomic units (e^3 a0^3 Eh^-2), the benchmark's convention. `metadata`
// is the engine's response_metadata.json, read for the reason behind a
// missing row.
inline json selected_beta(const Stage& stage, const json& properties,
                          const json& metadata) {
  const auto& rows
      = quadratic_rows(properties, "beta", final_protocol_key(stage));
  json selected = json::array();
  for (double omega : stage.hyperpolarizability.frequencies) {
    QuadraticElements elements;
    for (const auto& component : stage.hyperpolarizability.components) {
      const std::string b = std::string("dipole_") + component[1];
      const std::string c = std::string("dipole_") + component[2];
      select_quadratic_element(
          quadratic_row(rows, "hyperpolarizability", component[0], b, c, omega,
                        omega),
          component, quadratic_unavailable(metadata, b, c, omega, omega),
          elements);
    }
    json record = {{"frequency", omega},
                   {"frequency_units", "au"},
                   {"process", omega == 0.0 ? "static" : "shg"},
                   {"output_frequency", 2.0 * omega},
                   {"units", "e^3 a0^3 Eh^-2"},
                   {"components", elements.values},
                   {"converged", elements.converged},
                   {"engine_evidence", elements.evidence}};
    if (!elements.unavailable.empty())
      record["components_unavailable"] = elements.unavailable;
    selected.push_back(std::move(record));
  }
  return selected;
}

// The nuclear coordinates a Raman request covers, as the engine names them
// (nuc_<atom>_<axis>): the written atoms or every atom of the reference, times
// the written directions or x, y, z.
inline std::vector<std::pair<int, char>> raman_coordinates(const Stage& stage,
                                                           int natom) {
  std::vector<int> atoms = stage.raman_atoms;
  if (atoms.empty())
    for (int atom = 0; atom < natom; ++atom)
      atoms.push_back(atom);
  for (int atom : atoms)
    if (atom >= natom)
      throw std::runtime_error("raman.atoms names atom " + std::to_string(atom)
                               + " but the reference has "
                               + std::to_string(natom) + " atom(s)");
  const std::string directions
      = stage.raman_directions.empty() ? "xyz" : stage.raman_directions;
  std::vector<std::pair<int, char>> coordinates;
  for (int atom : atoms)
    for (char axis : directions)
      coordinates.emplace_back(atom, axis);
  return coordinates;
}

// The polarizability gradient: for each optical frequency w and nuclear
// coordinate Q, the engine's rows d(alpha_AB)/dQ(w) for the requested pairs
// AB, in the engine's units, a.u. (e^2 bohr^2 / Eh per bohr of displacement).
// With the full set the
// engine's normal-mode analysis is passed through unchanged, and the stage is
// converged only if its completeness block says every row was.
inline json selected_raman(const Stage& stage, const json& properties,
                           const json& metadata,
                           const json& vibrational_analysis,
                           const json& raman_spectra, int natom) {
  const auto pkey = final_protocol_key(stage);
  const auto& rows = quadratic_rows(properties, "raman", pkey);
  json values_out = json::array();
  for (double omega : stage.raman.frequencies) {
    for (const auto& [atom, axis] : raman_coordinates(stage, natom)) {
      const std::string coordinate
          = "nuc_" + std::to_string(atom) + "_" + std::string(1, axis);
      QuadraticElements elements;
      for (const auto& component : stage.raman.components) {
        const std::string b = std::string("dipole_") + component[1];
        select_quadratic_element(
            quadratic_row(rows, "raman", component[0], b, coordinate, omega,
                          0.0),
            component,
            quadratic_unavailable(metadata, b, coordinate, omega, 0.0),
            elements);
      }
      json record = {{"frequency", omega},
                     {"frequency_units", "au"},
                     {"coordinate", coordinate},
                     {"atom", atom},
                     {"axis", std::string(1, axis)},
                     {"units", "a.u. (e^2 bohr^2 / Eh per bohr)"},
                     {"components", elements.values},
                     {"converged", elements.converged},
                     {"engine_evidence", elements.evidence}};
      if (!elements.unavailable.empty())
        record["components_unavailable"] = elements.unavailable;
      values_out.push_back(std::move(record));
    }
  }
  json record = {{"values", values_out},
                 {"coordinates", stage.raman_full_set ? "all" : "written"}};
  if (!stage.raman_full_set) {
    record["normal_mode_analysis"] = false;
    record["normal_mode_analysis_unavailable"]
        = "a written subset of coordinates gives derivative rows only; omit "
          "raman.atoms and raman.directions for the full set";
    return record;
  }
  // The engine assembled the spectrum from every 3N x 9 row at this rung; its
  // completeness block, per frequency, is the record of whether it could.
  if (!raman_spectra.is_object() || !raman_spectra.contains("completeness"))
    throw std::runtime_error("the full Raman set has no completeness record");
  bool complete = true;
  for (double omega : stage.raman.frequencies) {
    const auto key = molresponse_v3::ResponseMetadata::freq_key(omega);
    const auto& completeness = raman_spectra.at("completeness");
    if (!completeness.contains(key))
      throw std::runtime_error("the Raman completeness record has no entry at "
                               "frequency "
                               + key);
    const auto& entry = completeness.at(key);
    complete = complete && entry.at("rows_found") == entry.at("rows_expected")
               && entry.value("all_rows_converged", false);
  }
  record["normal_mode_analysis"] = true;
  record["spectrum_complete"] = complete;
  record["vibrational_analysis"] = vibrational_analysis;
  record["raman_spectra"] = raman_spectra;
  return record;
}

// The transition analysis of one root out of the engine's es_analysis file
// (es_analysis.hpp: roots[] with slot and, when the bundle carries stable
// identities, root_id), matched by root_id, else by slot; nullptr when the
// file has none for it.
inline const json* transition_of_root(const json& transition,
                                      const json& root) {
  if (!transition.is_object() || !transition.contains("roots"))
    return nullptr;
  for (const auto& candidate : transition.at("roots")) {
    if (root.contains("root_id") && candidate.contains("root_id")) {
      if (candidate.at("root_id") == root.at("root_id"))
        return &candidate;
      continue;
    }
    if (candidate.value("slot", -1) == root.value("slot", -2))
      return &candidate;
  }
  return nullptr;
}

// D-061: the written econv beside the engine's gate, in the record: the
// request and the evidence, not the verdict. `target` is the written value
// or null (D8: the engine's default is not a request); `worst` the largest
// last eigenvalue change over the roots with its root, the same reduction as
// the state record's engine_measurements.omega_residual (null when a root
// has none yet); with a target each root says whether its own change is
// within it (null without a value). The status (satisfied, unmet,
// unavailable, not_requested) and, in engine mode, the "unmet, not
// enforceable" outcome are the shared verdict's, one judge (B1):
// convergence_report.requested_accuracy under the bundle's state key, where
// `verdict_in` points. No status, enforced or outcome field here.
inline json requested_omega_target(const Stage& stage, json& roots,
                                   const std::string& state_key) {
  namespace cv = madness::convergence;
  const auto& target = stage.excited_states.omega_target;
  json worst;
  for (const auto& root : roots) {
    const auto value = root.value("omega_residual", json());
    if (!value.is_number()) {
      worst = nullptr;
      break;
    }
    if (worst.is_null()
        || value.get<double>() > worst.at("value").get<double>())
      worst = {{"value", value}, {"root_id", root.value("root_id", "")}};
  }
  if (target)
    for (auto& root : roots) {
      const auto value = root.value("omega_residual", json());
      root["omega_residual_within_requested"]
          = value.is_number() ? json(value.get<double>() <= *target)
                              : json(nullptr);
    }
  return {
      {"target", target ? json(*target) : json(nullptr)},
      {"set_by", cv::to_string(target ? cv::SetBy::user : cv::SetBy::not_set)},
      {"units", "hartree"},
      {"relation", cv::to_string(cv::Relation::less_than_or_equal)},
      {"worst", worst},
      {"verdict_in",
       {{"path", "convergence_report.requested_accuracy"},
        {"state_key", state_key},
        {"criterion", "omega_residual"}}}};
}

// The excited-state bundle the stage requested (D-025). The engine saves one
// bundle per rung under excited_states.<protocol key> of its
// response_metadata.json (es_save_load.hpp: the bundle's verdict, stop
// reason and iterations, one entry per root with its verdict, residuals and
// targets since B1), and writes the transition analysis of a converged
// bundle to es_analysis__<key>.json (es_analysis.hpp), passed here as
// `transition`. The record is the bundle at the stage's final rung: the
// method, root count and guess the stage asked for (the guess named even
// when the engine's default), the bundle's exit, and its roots keyed by the
// engine's stable root_id (the slot and display name can change between
// rungs, both are kept beside it), each with the engine's fields and the
// transition properties or the reason they are missing. The bundle is
// converged when it ended through `converged` and every root is converged
// and not stalled. A bundle of another method or root count at that rung is
// an error (the engine keys bundles by rung alone); no bundle at the final
// rung is a record without roots, unconverged, with the reason.
inline json selected_excited_states(const Stage& stage, const json& metadata,
                                    const json& transition) {
  const auto& request = stage.excited_states;
  const auto& rp = stage.params.get<ResponseParameters>();
  const auto pkey = final_protocol_key(stage);
  // The bundle's state key at the final rung, as the state record carries it
  // (ResponseCheckpoint::bundle_key joins the node id and the protocol key).
  const auto state_key
      = molresponse_v3::es_node_id(request.method == "tda", request.roots) + "|"
        + pkey;
  json record = {{"method", request.method},
                 {"roots_requested", request.roots},
                 {"guess", rp.excited_guess()},
                 {"guess_basis", rp.excited_guess_basis()},
                 {"guess_written", !request.guess.empty()},
                 {"protocol_key", pkey},
                 {"omega_units", "au"},
                 {"converged", false},
                 {"roots", json::array()}};
  const json* bundle = nullptr;
  if (metadata.is_object() && metadata.contains("excited_states")) {
    const auto& bundles = metadata.at("excited_states");
    if (bundles.is_object() && bundles.contains(pkey))
      bundle = &bundles.at(pkey);
  }
  if (bundle == nullptr) {
    record["unavailable"]
        = "the engine saved no excited-state bundle at the final rung " + pkey;
    record["requested_omega_target"]
        = requested_omega_target(stage, record["roots"], state_key);
    return record;
  }
  const std::string type = request.method == "tda" ? "tda" : "full";
  if (bundle->value("type", type) != type
      || bundle->value("n_roots", request.roots) != request.roots)
    throw std::runtime_error("the excited-state bundle at " + pkey + " is "
                             + bundle->value("type", std::string("?"))
                             + " with "
                             + std::to_string(bundle->value("n_roots", 0))
                             + " root(s), not the requested " + request.method
                             + " with " + std::to_string(request.roots));
  const std::string stop = bundle->value("stop_reason", std::string());
  const bool bundle_converged = bundle->value("converged", false)
                                && stop == "converged"
                                && !bundle->value("diverged", false)
                                && !bundle->value("stalled", false);
  record["stop_reason"] = stop;
  record["iterations_at_rung"] = bundle->value("iter", 0);
  record["diverged"] = bundle->value("diverged", false);
  record["stalled"] = bundle->value("stalled", false);
  // Who decided the bundle's verdict: the engine's check, or the planner
  // when a met load certified it (D-016); the checkpoint writes the field.
  const std::string decided_by
      = bundle->value("decided_by", std::string("engine_check"));
  record["decided_by"] = decided_by;
  bool every_root = true;
  json roots = json::array();
  for (const auto& root : bundle->value("roots", json::array())) {
    json entry = root;
    // Archive bookkeeping, not the record.
    for (const auto* key : {"file", "coeffs", "bytes"})
      entry.erase(key);
    entry["omega_units"] = "au";
    entry["decided_by"] = root.value("decided_by", decided_by);
    every_root = every_root && root.value("converged", false)
                 && !root.value("stalled", false);
    if (const json* analysis = transition_of_root(transition, root))
      entry["transition"] = *analysis;
    else
      entry["transition_unavailable"]
          = bundle_converged ? "no transition analysis for this root"
                             : "bundle not converged (the engine writes the "
                               "transition analysis of a converged bundle "
                               "only)";
    roots.push_back(std::move(entry));
  }
  record["roots"] = std::move(roots);
  // D-104: the engine's units of the per-root transition blocks, once.
  if (transition.is_object() && transition.contains("units"))
    record["transition_units"] = transition.at("units");
  record["requested_omega_target"]
      = requested_omega_target(stage, record["roots"], state_key);
  record["converged"] = bundle_converged && every_root;
  json evidence = *bundle;
  for (const auto* key : {"roots", "kinetic_history"})
    evidence.erase(key);
  record["engine_evidence"] = std::move(evidence);
  return record;
}

// D-061: a written econv on a stage that requests excited_states is the
// target of the excitation energy change per root (D8: the engine's default
// is never a request), judged on the worst root with <=, in hartree. It is
// never the loop's target: the bundle's solve stops on the engine's own gate
// on the same quantity, 5 x max(dconv, thresh), strict, per root (the dconv
// engine row, D-025), in both modes (D-092); in criteria mode the written
// econv is enforced as the stage verdict, by the shared judge (B1, B8).
// This is the row's place in the settings table (B9); `engine_omega_limit`
// is that gate at the stage's final rung. The response settings table
// (ResponseCheckpoint.hpp) places it in place of its "energy change" row.
[[nodiscard]] inline madness::convergence::MetricSetting
excited_state_omega_target_row(const Stage& stage, double engine_omega_limit) {
  namespace cv = madness::convergence;
  const auto& target = stage.excited_states.omega_target;
  cv::MetricSetting row;
  row.group = cv::Group::requested_criteria;
  row.name = "excitation energy change per root";
  row.keyword = "econv";
  row.measurement = "omega_residual";
  row.units = "hartree";
  row.target = target;
  row.set_by = target ? cv::SetBy::user : cv::SetBy::not_set;
  row.relation = cv::Relation::less_than_or_equal;
  if (!target)
    row.note = "not set; the bundle stops on the engine's omega gate in both "
               "modes (the dconv row above); the engine's default is never a "
               "target";
  row.engine_check_instead
      = "stops the bundle when each root's eigenvalue change is below "
        + cv::detail::scientific(engine_omega_limit)
        + " (5 x max(dconv, thresh)), strictly, in both modes; a written "
          "econv is judged by the stage verdict";
  return row;
}

// D-099 (B3): the bundle's plateau detector as two engine rows of the
// settings table: a root is stalled when neither its density change nor its
// eigenvalue change has improved by stall.ratio over stall.window
// iterations while the other active roots are converged or likewise flat;
// the solve then stops as stalled. The keys are the engine's own; unwritten,
// its defaults apply (6 iterations, 0.1). stall.window 0 disables the
// detector, so the bundle runs to maxiter and --continue resumes it. The
// legs' detector takes progress.* instead. Not a target: no relation to a
// measurement, the limit is the setting itself. The response settings table
// (ResponseCheckpoint.hpp) places them beside the bundle's engine rows.
[[nodiscard]] inline std::vector<madness::convergence::MetricSetting>
bundle_stall_rows(const Stage& stage) {
  namespace cv = madness::convergence;
  const auto& rp = stage.params.get<ResponseParameters>();
  auto row = [&](const char* name, const char* keyword, const char* units,
                 double value, bool written, std::string note) {
    cv::MetricSetting r;
    r.group = cv::Group::engine_check;
    r.name = name;
    r.keyword = keyword;
    r.measurement = keyword;
    r.units = units;
    r.target = value;
    r.set_by = written ? cv::SetBy::user : cv::SetBy::engine_default;
    r.effective_limit = value;
    r.note = std::move(note);
    return r;
  };
  const bool window_written = rp.is_user_defined("stall.window");
  const bool ratio_written = rp.is_user_defined("stall.ratio");
  const int window = rp.stall_window();
  return {row("bundle stall window", "stall.window", "iterations",
              static_cast<double>(window), window_written,
              window == 0 ? "detector off: the bundle runs to maxiter"
                          : "a root is stalled when neither its density "
                            "change nor its eigenvalue change improves by "
                            "stall.ratio over this many iterations"),
          row("bundle stall ratio", "stall.ratio", "fraction", rp.stall_ratio(),
              ratio_written,
              "the relative improvement over stall.window that counts as "
              "progress; the legs take progress.* instead")};
}

// D-061, B1: the saved bundle's measurement under the name of its written
// target, for the shared judge (evaluate_explicit_criteria over
// convergence::requested_targets, which needs value, availability and
// units): the engine's worst-root eigenvalue change, the omega_residual row
// of the bundle's engine measurements (value, root_id, availability, phase;
// ResponseCheckpoint::bundle_measurements) with its units. The change across
// the last step is the saved bundle's evidence: the bundle saved is that
// step's output. A row without a value (a root with no eigenvalue step yet)
// stays unavailable. No engine rows, no judged row.
[[nodiscard]] inline json
bundle_saved_measurements(const json& engine_measurements) {
  json saved = json::object();
  if (!engine_measurements.contains("omega_residual"))
    return saved;
  json row = engine_measurements.at("omega_residual");
  row["units"] = "hartree";
  saved["omega_residual"] = std::move(row);
  return saved;
}

// D-095: the bundle's engine checks as per-root outcomes under the two
// engine rows (ResponseCheckpoint.hpp: "density change per root",
// "excitation energy change per root"), each root judged on its own
// recorded target, so the table and the record show the evidence the engine
// judged on (ConvergenceSettings.h: leg_outcomes, row_verdict). The phases
// are the engine's, as bundle_measurements names them: the density change
// between a step's input densities, the eigenvalue change across the step.
// A root without a value (no eigenvalue step yet) is unavailable; the
// subject is the root id, the record's key.
[[nodiscard]] inline std::vector<madness::convergence::MetricOutcome>
bundle_root_outcomes(const json& bundle) {
  struct Row {
    const char* measurement;
    const char* value;
    const char* limit;
    const char* units;
    const char* phase;
  };
  static const std::array<Row, 2> rows
      = {{{"engine_es_density_change", "density_residual", "density_target",
           "kernel_induced_density_L2_atomic_units",
           "solver_step_input_before_acceleration"},
          {"engine_es_omega_change", "omega_residual", "omega_target",
           "hartree", "solver_step_output"}}};
  std::vector<madness::convergence::MetricOutcome> outcomes;
  for (const auto& row : rows)
    for (const auto& root : bundle.value("roots", json::array())) {
      madness::convergence::MetricOutcome leg;
      leg.measurement = row.measurement;
      leg.subject = root.value("root_id", std::string());
      leg.units = row.units;
      leg.phase = row.phase;
      if (root.value(row.limit, json()).is_number())
        leg.limit = root.at(row.limit).get<double>();
      if (root.value(row.value, json()).is_number()) {
        leg.value = root.at(row.value).get<double>();
        leg.availability = "measured";
      }
      outcomes.push_back(std::move(leg));
    }
  return outcomes;
}

// Every property the stage requested, keyed by its name, each with its values
// (the runner adds the state identities). `converged` is the stage's verdict:
// every selected row converged (an element without a value is not), for the
// full Raman set the spectrum complete, and the excited-state bundle
// converged. `metadata` is the engine's response_metadata.json of the run
// (or the saved copy on a load); `transition` its transition analysis of the
// bundle at the final rung, when one was written.
inline json select_response_properties(const Stage& stage,
                                       const json& properties,
                                       const json& metadata,
                                       const json& vibrational_analysis,
                                       const json& raman_spectra, int natom,
                                       const json& transition = json()) {
  json selected = json::object();
  for (const auto& property : stage.properties) {
    if (property == "polarizability")
      selected[property] = {{"values", selected_alpha(stage, properties)}};
    else if (property == "hyperpolarizability")
      selected[property]
          = {{"values", selected_beta(stage, properties, metadata)}};
    else if (property == "raman")
      selected[property]
          = selected_raman(stage, properties, metadata, vibrational_analysis,
                           raman_spectra, natom);
    else if (property == "excited_states")
      selected[property] = selected_excited_states(stage, metadata, transition);
    else
      throw std::runtime_error("no selector for property " + property);
  }
  return selected;
}
} // namespace named_stages
