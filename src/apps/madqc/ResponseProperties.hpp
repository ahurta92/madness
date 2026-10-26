#pragma once
#include "NamedStages.hpp"
#include <nlohmann/json.hpp>

// Select what a molresponse stage requested out of the engine's property rows.
// The response calculation has already run; nothing here computes a number.
// One function per property (D-017): the engine writes one row per tensor
// element it assembled under properties.<alpha|beta|raman>.<protocol key>,
// and the stage asks for frequencies and Cartesian components. A requested
// element the engine did not write is an error, never a silent gap, and the
// row it came from travels along as engine_evidence.
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

// The one engine row of a quadratic property (beta or d(alpha)/dQ) with output
// axis A and inputs B and C at (freq_b, freq_c). The engine assembles each
// element once; two rows would mean two assemblies, which is an error here,
// not a choice.
inline const json& quadratic_row(const json& rows, const std::string& what,
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
    throw std::runtime_error("missing requested " + what + " element " + a + " "
                             + b + " " + c);
  if (!std::isfinite(found->at("beta").get<double>()))
    throw std::runtime_error("nonfinite " + what);
  return *found;
}

// The hyperpolarizability: for each driver frequency w and triple ABC, the
// engine's SHG row beta_ABC(-2w; w, w) (static beta at w = 0). The value is
// in atomic units (e^3 a0^3 Eh^-2), the benchmark's convention.
inline json selected_beta(const Stage& stage, const json& properties) {
  const auto& rows = properties.at("beta").at(final_protocol_key(stage));
  json selected = json::array();
  for (double omega : stage.hyperpolarizability.frequencies) {
    json values = json::object(), evidence = json::object();
    bool converged = true;
    for (const auto& component : stage.hyperpolarizability.components) {
      const auto& row
          = quadratic_row(rows, "hyperpolarizability", component[0],
                          std::string("dipole_") + component[1],
                          std::string("dipole_") + component[2], omega, omega);
      values[component] = row.at("beta").get<double>();
      converged = converged && row.value("converged", false);
      evidence[component] = row;
    }
    selected.push_back({{"frequency", omega},
                        {"frequency_units", "au"},
                        {"process", omega == 0.0 ? "static" : "shg"},
                        {"output_frequency", 2.0 * omega},
                        {"units", "e^3 a0^3 Eh^-2"},
                        {"components", values},
                        {"converged", converged},
                        {"engine_evidence", evidence}});
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
                           const json& vibrational_analysis,
                           const json& raman_spectra, int natom) {
  const auto pkey = final_protocol_key(stage);
  const auto& rows = properties.at("raman").at(pkey);
  json values_out = json::array();
  for (double omega : stage.raman.frequencies) {
    for (const auto& [atom, axis] : raman_coordinates(stage, natom)) {
      const std::string coordinate
          = "nuc_" + std::to_string(atom) + "_" + std::string(1, axis);
      json values = json::object(), evidence = json::object();
      bool converged = true;
      for (const auto& component : stage.raman.components) {
        const auto& row = quadratic_row(rows, "raman", component[0],
                                        std::string("dipole_") + component[1],
                                        coordinate, omega, 0.0);
        values[component] = row.at("beta").get<double>();
        converged = converged && row.value("converged", false);
        evidence[component] = row;
      }
      values_out.push_back({{"frequency", omega},
                            {"frequency_units", "au"},
                            {"coordinate", coordinate},
                            {"atom", atom},
                            {"axis", std::string(1, axis)},
                            {"units", "a0^2 per bohr"},
                            {"components", values},
                            {"converged", converged},
                            {"engine_evidence", evidence}});
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
// every selected row converged, and for the full Raman set the spectrum
// complete.
inline json select_response_properties(const Stage& stage,
                                       const json& properties,
                                       const json& vibrational_analysis,
                                       const json& raman_spectra, int natom) {
  json selected = json::object();
  for (const auto& property : stage.properties) {
    if (property == "polarizability")
      selected[property] = {{"values", selected_alpha(stage, properties)}};
    else if (property == "hyperpolarizability")
      selected[property] = {{"values", selected_beta(stage, properties)}};
    else if (property == "raman")
      selected[property] = selected_raman(
          stage, properties, vibrational_analysis, raman_spectra, natom);
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
