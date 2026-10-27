#pragma once
#include "NamedStages.hpp"
#include <nlohmann/json.hpp>

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
// AB, in atomic units (a0^2 per bohr of displacement). With the full set the
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
                     {"units", "a0^2 per bohr"},
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

// Every property the stage requested, keyed by its name, each with its values
// (the runner adds the state identities). `converged` is the stage's verdict:
// every selected row converged (an element without a value is not), and for
// the full Raman set the spectrum complete. `metadata` is the engine's
// response_metadata.json of the run (or the saved copy on a load).
inline json select_response_properties(const Stage& stage,
                                       const json& properties,
                                       const json& metadata,
                                       const json& vibrational_analysis,
                                       const json& raman_spectra, int natom) {
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
    else
      throw std::runtime_error("no selector for property " + property);
  }
  return selected;
}
inline bool all_requested_converged(const json& selected) {
  for (const auto& [name, record] : selected.items()) {
    for (const auto& row : record.at("values"))
      if (!row.at("converged").get<bool>())
        return false;
    if (record.contains("spectrum_complete")
        && !record.at("spectrum_complete").get<bool>())
      return false;
  }
  return true;
}
} // namespace named_stages
