#pragma once

// The derived quantities of results.json: for each item (or single-answer
// property) whose inputs are complete, a `derived` object, and on its property
// `definitions` (the formula of each, in plain text) and `derived_units`. The
// numbers come from DerivedQuantities.h; this file only finds the inputs in
// results.json and writes the outputs beside them. It reads nothing else, so a
// unit test can feed it a hand-written results.json.
//
// Two-photon absorption is derived by its producer while it assembles (the
// same functions), and Raman's per-mode quantities are rows of the spectrum;
// neither is touched here.
#include <madness/chem/DerivedQuantities.h>

#include <nlohmann/json.hpp>

#include <optional>
#include <string>
#include <vector>

namespace madness::iteration {

namespace derived_results {

using json = nlohmann::json;
namespace dq = madness::derived;

/// A dipole smaller than this in norm (atomic units) has no direction.
inline constexpr double dipole_floor = 1e-8;

inline const char* const axes = "xyz";

inline bool complete(const json& value) { return value.value("status", std::string()) == "complete"; }

inline std::optional<double> number(const json& object, const std::string& key) {
    if (!object.is_object() || !object.contains(key) || !object.at(key).is_number()) return std::nullopt;
    return object.at(key).get<double>();
}

/// The quantities this file derives, by property: the derived names it owns
/// and their definitions. Writing them first removes any earlier ones, so a
/// value that is no longer complete loses its derived quantities.
struct Definition {
    const char* name;
    const char* text;
};

inline const std::vector<Definition>& polarizability_definitions() {
    static const std::vector<Definition> d{
        {"isotropic", "(xx + yy + zz) / 3"},
        {"anisotropy", "sqrt( 1/2 [ (xx-yy)^2 + (yy-zz)^2 + (zz-xx)^2 ] + 3 [ xy^2 + xz^2 + yz^2 ] ); an "
                       "off-diagonal pair ij, ji enters as the mean of the two components"}};
    return d;
}

inline const std::vector<Definition>& hyperpolarizability_definitions() {
    static const std::vector<Definition> d{
        {"vector", "beta_i = (1/3) sum_j ( beta_ijj + beta_jij + beta_jji ), with beta_ijk the component named ijk"},
        {"vector_magnitude", "|vector|"},
        {"parallel", "(3/5) sum_i beta_i mu_i / |mu|, with mu the ground-state dipole; absent when the dipole was "
                     "not computed or |mu| < 1e-8"},
        {"hrs", "sqrt( <beta_ZZZ^2> + <beta_ZXX^2> ), the hyper-Rayleigh orientation averages of the laboratory "
                "components (Bersohn, Pao and Frisch 1966), no Kleinman symmetry assumed"},
        {"hrs_depolarization", "<beta_ZZZ^2> / <beta_ZXX^2>; absent when hrs < 1e-8, where the tensor is "
                               "numerical noise"}};
    return d;
}

inline const std::vector<Definition>& excited_state_definitions() {
    static const std::vector<Definition> d{
        {"energy_ev", "excitation_energy x 27.211386245981 (CODATA hartree-electron volt relationship)"},
        {"wavelength_nm", "1e7 / ( excitation_energy x 219474.63136314 ), the vacuum wavelength in nm of a photon of "
                          "the excitation energy"}};
    return d;
}

inline const std::vector<Definition>& orbital_energy_definitions() {
    static const std::vector<Definition> d{
        {"homo_energy", "the highest of the occupied orbital energies, per spin when the property has beta "
                        "energies; present only when the energies are canonical (not localized) or each spin has "
                        "at most one occupied orbital"}};
    return d;
}

inline const std::vector<Definition>& dipole_definitions() {
    static const std::vector<Definition> d{
        {"dipole_magnitude", "|dipole|, in atomic units (e a0)"},
        {"dipole_magnitude_debye", "|dipole| in Debye (e a0 = 2.541746 D, from constants.h)"}};
    return d;
}

/// Remove this file's earlier output from `property` and its items.
inline void clear(json& property, const std::vector<Definition>& definitions) {
    const auto strip = [&](const char* key) {
        if (!property.contains(key) || !property.at(key).is_object()) return;
        for (const auto& d : definitions) property.at(key).erase(d.name);
        if (property.at(key).empty()) property.erase(key);
    };
    strip("definitions");
    strip("derived_units");
    property.erase("derived");
    if (property.contains("values") && property.at("values").is_array())
        for (auto& item : property.at("values")) item.erase("derived");
}

/// Write `derived` into `target` and, once any value of the property has
/// one, the definitions and units of the names it holds.
inline void put(json& property, json& target, const json& derived, const std::vector<Definition>& definitions,
                const json& units) {
    if (derived.empty()) return;
    target["derived"] = derived;
    for (const auto& d : definitions) {
        if (!derived.contains(d.name)) continue;
        property["definitions"][d.name] = d.text;
        if (units.contains(d.name)) property["derived_units"][d.name] = units.at(d.name);
    }
}

/// The unit of a property whose `units` is one string, or `fallback`.
inline json property_unit(const json& property, const char* fallback) {
    const auto units = property.value("units", json());
    return units.is_string() ? units : json(fallback);
}

inline void derive_polarizability(json& property) {
    clear(property, polarizability_definitions());
    const auto unit = property_unit(property, "atomic_units");
    const json units{{"isotropic", unit}, {"anisotropy", unit}};
    if (!property.contains("values") || !property.at("values").is_array()) return;
    for (auto& item : property.at("values")) {
        if (!complete(item)) continue;
        const auto components = item.value("components", json::object());
        const auto xx = number(components, "xx"), yy = number(components, "yy"), zz = number(components, "zz");
        if (!xx || !yy || !zz) continue;
        json derived{{"isotropic", dq::polarizability_isotropic(*xx, *yy, *zz)}};
        // An off-diagonal pair enters as the mean of the components present.
        const auto pair = [&](const char* ij, const char* ji) -> std::optional<double> {
            const auto a = number(components, ij), b = number(components, ji);
            if (a && b) return 0.5 * (*a + *b);
            return a ? a : b;
        };
        const auto xy = pair("xy", "yx"), xz = pair("xz", "zx"), yz = pair("yz", "zy");
        if (xy && xz && yz) derived["anisotropy"] = dq::polarizability_anisotropy(*xx, *yy, *zz, *xy, *xz, *yz);
        put(property, item, derived, polarizability_definitions(), units);
    }
}

/// The full 27-component tensor of an item, if every component is present.
inline std::optional<dq::Tensor3> full_hyperpolarizability(const json& components) {
    dq::Tensor3 b{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) {
                const auto value = number(components, std::string{axes[i], axes[j], axes[k]});
                if (!value) return std::nullopt;
                b[i][j][k] = *value;
            }
    return b;
}

inline void derive_hyperpolarizability(json& property, const std::optional<dq::Vector3>& dipole) {
    clear(property, hyperpolarizability_definitions());
    const auto unit = property_unit(property, "atomic_units");
    const json units{{"vector", unit}, {"vector_magnitude", unit}, {"parallel", unit}, {"hrs", unit},
                     {"hrs_depolarization", "dimensionless"}};
    if (!property.contains("values") || !property.at("values").is_array()) return;
    for (auto& item : property.at("values")) {
        if (!complete(item)) continue;
        const auto b = full_hyperpolarizability(item.value("components", json::object()));
        if (!b) continue;
        const auto vector = dq::hyperpolarizability_vector(*b);
        const auto averages = dq::hyper_rayleigh_averages(*b);
        json derived{{"vector", {{"x", vector[0]}, {"y", vector[1]}, {"z", vector[2]}}},
                     {"vector_magnitude", dq::norm(vector)},
                     {"hrs", dq::hyper_rayleigh(averages)}};
        if (dq::hyper_rayleigh(averages) >= dq::hyper_rayleigh_floor && averages.zxx != 0.0)
            derived["hrs_depolarization"] = dq::hyper_rayleigh_depolarization(averages);
        if (dipole && dq::norm(*dipole) >= dipole_floor)
            derived["parallel"] = dq::hyperpolarizability_parallel(vector, *dipole);
        put(property, item, derived, hyperpolarizability_definitions(), units);
    }
}

inline void derive_excited_states(json& property) {
    clear(property, excited_state_definitions());
    const json units{{"energy_ev", "eV"}, {"wavelength_nm", "nm"}};
    if (!property.contains("values") || !property.at("values").is_array()) return;
    for (auto& item : property.at("values")) {
        if (!complete(item)) continue;
        const auto energy = number(item, "excitation_energy");
        if (!energy || *energy <= 0.0) continue;
        put(property, item,
            {{"energy_ev", dq::excitation_energy_ev(*energy)}, {"wavelength_nm", dq::wavelength_nm(*energy)}},
            excited_state_definitions(), units);
    }
}

/// The highest of the first `occupied` energies of a spin, if there are that many.
inline std::optional<double> highest_occupied(const json& energies, const json& occupied) {
    if (!energies.is_array() || !occupied.is_number_integer()) return std::nullopt;
    const auto n = occupied.get<long>();
    if (n <= 0 || n > static_cast<long>(energies.size())) return std::nullopt;
    std::optional<double> highest;
    for (long i = 0; i < n; ++i) {
        if (!energies.at(i).is_number()) return std::nullopt;
        const double e = energies.at(i).get<double>();
        if (!highest || e > *highest) highest = e;
    }
    return highest;
}

inline void derive_orbital_energies(json& property) {
    clear(property, orbital_energy_definitions());
    if (!complete(property) || !property.contains("value") || !property.at("value").is_object()) return;
    const auto& value = property.at("value");
    const auto occupied = property.value("occupied", json::object());
    // Localized orbitals carry diagonal Fock elements, which are not the
    // orbital energies unless a spin has a single occupied orbital.
    const bool canonical = property.value("orbitals", std::string()) == "canonical";
    json homo = json::object();
    for (const auto& [spin, energies] : value.items()) {
        const auto n = occupied.value(spin, json());
        if (!canonical && !(n.is_number_integer() && n.get<long>() <= 1)) return;
        if (n.is_number_integer() && n.get<long>() == 0) continue;
        const auto highest = highest_occupied(energies, n);
        if (!highest) return;
        homo[spin] = *highest;
    }
    if (homo.empty()) return;
    const json derived{{"homo_energy", value.contains("beta") ? homo : homo.begin().value()}};
    put(property, property, derived, orbital_energy_definitions(), {{"homo_energy", property_unit(property, "hartree")}});
}

/// The dipole of a complete `dipole` property: its three components.
inline std::optional<dq::Vector3> dipole_vector(const json& property) {
    if (!complete(property) || !property.contains("value")) return std::nullopt;
    const auto& value = property.at("value");
    if (!value.is_array() || value.size() != 3) return std::nullopt;
    dq::Vector3 v{};
    for (std::size_t i = 0; i < 3; ++i) {
        if (!value.at(i).is_number()) return std::nullopt;
        v[i] = value.at(i).get<double>();
    }
    return v;
}

inline void derive_dipole(json& property) {
    clear(property, dipole_definitions());
    const auto dipole = dipole_vector(property);
    if (!dipole) return;
    const double magnitude = dq::norm(*dipole);
    put(property, property,
        {{"dipole_magnitude", magnitude}, {"dipole_magnitude_debye", dq::debye_from_atomic_units(magnitude)}},
        dipole_definitions(), {{"dipole_magnitude", property_unit(property, "atomic_units")},
                               {"dipole_magnitude_debye", "Debye"}});
}

} // namespace derived_results

/// Write the derived quantities of every entry of `results` (a parsed or
/// assembled results.json), replacing any written before. The
/// hyperpolarizability's `parallel` takes the dipole from the `ground_state`
/// entry, when that entry has a complete one.
inline void add_derived_quantities(nlohmann::json& results) {
    using namespace derived_results;
    if (!results.contains("calculations") || !results.at("calculations").is_array()) return;
    std::optional<dq::Vector3> dipole;
    for (const auto& entry : results.at("calculations"))
        if (entry.value("calculation_id", std::string()) == "ground_state" && entry.contains("properties") &&
            entry.at("properties").contains("dipole"))
            dipole = dipole_vector(entry.at("properties").at("dipole"));
    for (auto& entry : results.at("calculations")) {
        if (!entry.contains("properties") || !entry.at("properties").is_object()) continue;
        for (auto& [name, property] : entry.at("properties").items()) {
            if (name == "polarizability") derive_polarizability(property);
            else if (name == "hyperpolarizability") derive_hyperpolarizability(property, dipole);
            else if (name == "excited_states") derive_excited_states(property);
            else if (name == "orbital_energies") derive_orbital_energies(property);
            else if (name == "dipole") derive_dipole(property);
        }
    }
}

} // namespace madness::iteration
