#pragma once

#include <madness/chem/ParameterManager.hpp>
#include <apps/molresponse/ResponseProtocol.hpp>
#include <apps/molresponse/solvers/response_metadata.hpp>
#include <cmath>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

/// Deliberately small, flat input surface for fresh moldft -> molresponse runs.
/// Scientific values are parsed by the existing typed parameter groups.
namespace named_stages {

struct Line { std::size_t number; std::string text; };
struct Block { std::string kind, name; std::vector<Line> lines; };

inline std::string clean(std::string s) {
    s = s.substr(0, s.find('#'));
    const auto first = s.find_first_not_of(" \t\r");
    return first == std::string::npos ? "" : s.substr(first, s.find_last_not_of(" \t\r") - first + 1);
}
inline std::string key(const std::string& s) {
    auto word = s.substr(0, s.find_first_of(" \t="));
    std::transform(word.begin(), word.end(), word.begin(), [](unsigned char c) { return std::tolower(c); });
    return word;
}
inline void fail(const Line& line, const std::string& why) {
    throw std::runtime_error("named stages line " + std::to_string(line.number) + ": " + why);
}
inline bool named_input(const std::string& text) {
    std::istringstream input(text);
    std::string line;
    while (std::getline(input, line)) {
        line = clean(line);
        if (line.empty()) continue;
        std::istringstream words(line);
        std::string kind, name;
        words >> kind >> name;
        return !name.empty() && line.front() != '{';
    }
    return false;
}
inline std::vector<Block> blocks(const std::string& text) {
    std::vector<Block> result;
    std::set<std::string> names;
    std::istringstream input(text);
    std::string textline;
    bool inside = false;
    std::size_t n = 0;
    while (std::getline(input, textline)) {
        Line line{++n, clean(textline)};
        if (line.text.empty()) continue;
        if (!inside) {
            std::istringstream words(line.text);
            Block b;
            std::string extra;
            words >> b.kind >> b.name;
            if (b.kind != "molecule" && b.kind != "moldft" && b.kind != "molresponse")
                fail(line, "unsupported stage type " + b.kind);
            if (b.name.empty() || words >> extra ||
                !std::isalpha(static_cast<unsigned char>(b.name.front())) ||
                b.name.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-") != std::string::npos)
                fail(line, "expected a safe stage name after " + b.kind);
            if (!names.insert(b.name).second) fail(line, "duplicate name " + b.name);
            result.push_back(std::move(b));
            inside = true;
        } else if (line.text == "end") {
            inside = false;
        } else {
            if (key(line.text) == "end") fail(line, "unexpected text after end");
            result.back().lines.push_back(std::move(line));
        }
    }
    if (inside) throw std::runtime_error("unterminated named block " + result.back().name);
    return result;
}

struct Controls : madness::QCCalculationParametersBase {
    Controls() {
        initialize<std::string>("molecule", "", "named molecule");
        initialize<std::string>("reference", "", "earlier moldft stage");
        initialize<std::vector<std::string>>("requested_properties", {}, "requested properties");
        initialize<std::vector<double>>("polarizability.frequencies", {}, "explicit frequencies");
        initialize<std::string>("polarizability.frequency_units", "", "frequency units", {"", "au"});
        initialize<std::vector<std::string>>("polarizability.components", {}, "Cartesian pairs");
    }
    [[nodiscard]] std::string get_tag() const override { return "stage"; }
    [[nodiscard]] std::set<std::string> case_sensitive_keys() const override { return {"molecule", "reference"}; }
};
struct Stage {
    std::string kind, name, reference;
    Params params;
    std::vector<std::string> properties, components;
    std::vector<double> frequencies;
    Stage(madness::World& world, const Block& b) : kind(b.kind), name(b.name), params(world) {}
};
inline void read_value(madness::World& world, madness::QCCalculationParametersBase& group, const Line& line) {
    std::string text = "stage\n" + line.text + "\nend\n";
    try { group.read_internal(world, text, "stage"); }
    catch (const std::exception& e) { fail(line, e.what()); }
}
inline void valid_protocol(const std::vector<double>& protocol) {
    if (protocol.empty()) throw std::runtime_error("protocol must not be empty");
    double previous = std::numeric_limits<double>::infinity();
    for (double t : protocol) {
        if (!std::isfinite(t) || t <= 0 || t >= previous)
            throw std::runtime_error("protocol must contain positive, strictly decreasing thresholds");
        previous = t;
    }
}
inline int protocol_k(double t) {
    return molresponse_v3::default_k_for_thresh(t);
}

inline std::vector<Stage> parse(madness::World& world, const std::string& text) {
    std::map<std::string, madness::Molecule> molecules;
    std::map<std::string, std::size_t> indices;
    std::vector<Stage> stages;
    // This bounded allowlist prevents accepting settings whose outputs or
    // scientific semantics have not yet been wired through the named adapter.
    const std::set<std::string> scf_keys = {"xc", "protocol", "dconv", "econv", "maxiter", "l", "k",
        "charge", "localize", "print_level", "maxsub", "maxrotn"};
    const std::set<std::string> response_keys = {"dconv", "maxiter", "kain", "print_level"};
    for (const auto& b : blocks(text)) {
        std::set<std::string> seen;
        if (b.kind == "molecule") {
            madness::Molecule mol;
            std::string body = "molecule\n";
            const std::set<std::string> geometry_keys = {"units", "eprec", "no_orient"};
            for (const auto& line : b.lines) {
                const auto k = key(line.text);
                if (geometry_keys.count(k)) {
                    if (!seen.insert(k).second) fail(line, "duplicate parameter " + k);
                    read_value(world, mol.parameters, line);
                } else {
                    std::istringstream atom(line.text);
                    std::string symbol, extra;
                    double x = 0, y = 0, z = 0;
                    if (!(atom >> symbol >> x >> y >> z) || atom >> extra ||
                        !std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
                        fail(line, "expected atom and three finite coordinates");
                    if (madness::symbol_to_atomic_number(symbol) <= 0)
                        fail(line, "expected an ordinary atom in coordinates");
                    body += line.text + "\n";
                }
            }
            body += "end\n";
            std::istringstream stream(body);
            mol.read(stream);
            if (!mol.natom()) throw std::runtime_error("empty molecule " + b.name);
            molecules.emplace(b.name, std::move(mol));
            continue;
        }
        Stage stage(world, b);
        Controls controls;
        auto& cp = stage.params.get<CalculationParameters>();
        auto& rp = stage.params.get<ResponseParameters>();
        for (const auto& line : b.lines) {
            const auto k = key(line.text);
            if (!seen.insert(k).second) fail(line, "duplicate parameter " + k);
            if (controls.parameter_exists(k)) read_value(world, controls, line);
            else if (k == "protocol" || (b.kind == "moldft" && scf_keys.count(k))) read_value(world, cp, line);
            else if (b.kind == "molresponse" && response_keys.count(k) && rp.parameter_exists(k)) read_value(world, rp, line);
            else fail(line, "unsupported parameter " + k + " in " + b.kind);
        }
        stage.properties = controls.get<std::vector<std::string>>("requested_properties");
        if (b.kind == "moldft") {
            const auto name = controls.get<std::string>("molecule");
            if (!molecules.count(name)) throw std::runtime_error("unresolved molecule " + name + " for " + b.name);
            if (seen.count("reference")) throw std::runtime_error("moldft takes molecule, not reference");
            for (const auto& k : seen) if (k.find("polarizability.") == 0) throw std::runtime_error("polarizability belongs to molresponse");
            stage.params.get<madness::Molecule>() = molecules.at(name);
            if (cp.get<std::string>("xc") != "hf") throw std::runtime_error("named moldft currently supports xc hf");
            if (stage.properties.empty()) throw std::runtime_error("requested_properties must be explicit for " + b.name);
            for (const auto& p : stage.properties)
                if (p != "energy" && p != "dipole") throw std::runtime_error("unsupported moldft property " + p);
            cp.set_user_defined_value("dipole", std::find(stage.properties.begin(), stage.properties.end(), "dipole") != stage.properties.end());
            cp.set_user_defined_value("restart", std::string("none"));
            cp.set_user_defined_value("save", true);
        } else {
            stage.reference = controls.get<std::string>("reference");
            if (!indices.count(stage.reference) || stages.at(indices.at(stage.reference)).kind != "moldft")
                throw std::runtime_error("unresolved earlier moldft reference " + stage.reference + " for " + b.name);
            if (seen.count("molecule")) throw std::runtime_error("molresponse takes reference, not molecule");
            const auto& source = stages.at(indices.at(stage.reference));
            const auto own_protocol = cp.protocol();
            cp = source.params.get<CalculationParameters>();
            cp.set_user_defined_value("protocol", own_protocol);
            stage.params.get<madness::Molecule>() = source.params.get<madness::Molecule>();
            if (stage.properties != std::vector<std::string>{"polarizability"})
                throw std::runtime_error("molresponse requested_properties must be [polarizability]");
            stage.frequencies = controls.get<std::vector<double>>("polarizability.frequencies");
            stage.components = controls.get<std::vector<std::string>>("polarizability.components");
            if (stage.frequencies.empty() || stage.components.empty() || controls.get<std::string>("polarizability.frequency_units") != "au")
                throw std::runtime_error("polarizability needs explicit frequencies, frequency_units au and components");
            std::set<std::string> frequencies;
            for (double f : stage.frequencies)
                if (!std::isfinite(f) || f < 0 || !frequencies.insert(molresponse_v3::ResponseMetadata::freq_key(f)).second)
                    throw std::runtime_error("polarizability frequencies must be finite, nonnegative and distinct at engine frequency-key precision");
            std::set<char> axes;
            std::set<std::string> components;
            for (const auto& c : stage.components) {
                if (c.size() != 2 || c.find_first_not_of("xyz") != std::string::npos || !components.insert(c).second)
                    throw std::runtime_error("polarizability components must be distinct Cartesian pairs");
                axes.insert(c.begin(), c.end());
            }
            rp.set_user_defined_value("requested_properties", stage.properties);
            rp.set_user_defined_value("dipole.frequencies", stage.frequencies);
            rp.set_user_defined_value("dipole.directions", std::string(axes.begin(), axes.end()));
            // Defaults must come from this stage, not the upstream SCF maxiter.
            rp.set_user_defined_value("maxiter", rp.get<std::size_t>("maxiter"));
            valid_protocol(own_protocol);
            std::set<std::string> protocol_keys;
            for (double t : own_protocol)
                if (!protocol_keys.insert(molresponse_v3::protocol_key(t, protocol_k(t))).second)
                    throw std::runtime_error("response protocol rungs collide at engine protocol-key precision");
            const auto& reference_cp = source.params.get<CalculationParameters>();
            const int source_k = reference_cp.k() == -1 ? protocol_k(reference_cp.protocol().back()) : reference_cp.k();
            if (own_protocol.back() < reference_cp.protocol().back() || protocol_k(own_protocol.back()) > source_k)
                throw std::runtime_error("response protocol is tighter than reference " + stage.reference);
        }
        valid_protocol(cp.protocol());
        if (!std::isfinite(cp.dconv()) || cp.dconv() <= 0 || cp.maxiter() <= 0)
            throw std::runtime_error("positive finite dconv and positive maxiter required");
        stage.params.set_derived_values();
        if (!cp.spin_restricted() || cp.nalpha() <= 0 || cp.nalpha() != cp.nbeta())
            throw std::runtime_error("named stages currently require a closed-shell reference");
        if (!std::isfinite(cp.econv()) || cp.econv() <= 0 || !std::isfinite(cp.L()) || cp.L() <= 0)
            throw std::runtime_error("positive finite econv and L required");
        if (!std::isfinite(rp.dconv()) || rp.dconv() <= 0 || rp.maxiter() == 0)
            throw std::runtime_error("positive finite response dconv and positive maxiter required");
        indices.emplace(b.name, stages.size());
        stages.push_back(std::move(stage));
    }
    if (stages.empty()) throw std::runtime_error("named input has no application stages");
    return stages;
}
} // namespace named_stages
