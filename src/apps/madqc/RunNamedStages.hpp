#pragma once
#include "NamedStages.hpp"
#include <madness/chem/Applications.hpp>
#include <apps/molresponse/madqc_adapter.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <unistd.h>

namespace named_stages {
namespace fs = std::filesystem;
using json = nlohmann::json;

template<class F> void on_root(madness::World& world, F operation) {
    std::string error;
    if (world.rank() == 0) {
        try { operation(); }
        catch (const std::exception& e) { error = e.what(); }
    }
    world.gop.broadcast_serializable(error, 0);
    if (!error.empty()) throw std::runtime_error(error);
}
inline void write_result(madness::World& world, const fs::path& file, const json& value) {
    on_root(world, [&] {
        std::ofstream output(file);
        output.exceptions(std::ios::failbit | std::ios::badbit);
        output << value.dump(2) << '\n';
        output.close();
    });
}
inline json selected_alpha(const Stage& stage, const json& properties) {
    const auto& cp = stage.params.get<CalculationParameters>();
    const auto pkey = molresponse_v3::protocol_key(cp.protocol().back(), protocol_k(cp.protocol().back()));
    const auto& rows = properties.at("alpha").at(pkey);
    json selected = json::array();
    for (double omega : stage.frequencies) {
        const json* found = nullptr;
        for (const auto& row : rows) {
            if (std::abs(row.at("omega").get<double>() - omega) < 1e-12) {
                if (found) throw std::runtime_error("duplicate alpha frequency in engine results");
                found = &row;
            }
        }
        if (!found) throw std::runtime_error("missing requested polarizability frequency");
        const std::string dirs = found->at("directions");
        json values = json::object();
        for (const auto& component : stage.components) {
            const auto i = dirs.find(component[0]), j = dirs.find(component[1]);
            if (i == std::string::npos || j == std::string::npos) throw std::runtime_error("missing requested alpha component");
            const double value = found->at("alpha").at(i).at(j).get<double>();
            if (!std::isfinite(value)) throw std::runtime_error("nonfinite polarizability");
            values[component] = value;
        }
        selected.push_back({{"frequency", omega}, {"frequency_units", "au"},
            {"units", "a0^3"}, {"components", values}, {"converged", found->value("converged", false)},
            {"engine_evidence", *found}});
    }
    return selected;
}

/// Returns false for the existing unnamed input surface. Missing default input
/// is left to that path, which also accepts geometry supplied on the command line.
inline bool run_if_named(madness::World& world, const madness::commandlineparser& parser) {
    std::string text;
    on_root(world, [&] {
        std::ifstream input(parser.value("input"));
        if (!input) return;
        std::ostringstream contents;
        contents << input.rdbuf();
        if (input.bad()) throw std::runtime_error("reading named input failed");
        text = contents.str();
    });
    world.gop.broadcast_serializable(text, 0);
    if (!named_input(text)) return false;
    for (const auto& option : parser.keyval) {
        if (option.first != "input" && option.first != "prefix" && option.first != "user_defined_input_file" && option.first != "user_defined_prefix")
            throw std::runtime_error("named stages: put settings inside their stage; unsupported option --" + option.first);
    }
    auto stages = parse(world, text); // validate the complete deck before creating outputs
    fs::path root;
    std::string root_string, run_id;
    on_root(world, [&] {
        root = fs::absolute(parser.key_exists("user_defined_prefix") ? fs::path(parser.value_raw("prefix"))
                                  : fs::path(parser.value("input")).stem()).lexically_normal();
        fs::create_directories(root.parent_path());
        if (!fs::create_directory(root)) throw std::runtime_error("fresh named-stage run requires a new output directory: " + root.string());
        root_string = root.string();
        run_id = std::to_string(std::chrono::system_clock::now().time_since_epoch().count()) + "-" + std::to_string(getpid());
        for (const auto& stage : stages)
            if (!fs::create_directory(root / stage.name))
                throw std::runtime_error("stage output directory collision: " + stage.name);
    });
    world.gop.broadcast_serializable(root_string, 0);
    world.gop.broadcast_serializable(run_id, 0);
    root = root_string;
    std::map<std::string, std::shared_ptr<madness::SCF>> references;
    std::map<std::string, std::string> state_ids;
    json index = {{"schema", "madqc.named-stages.v1"}, {"run_id", run_id},
                  {"mode", "fresh"}, {"mpi_ranks", world.size()}, {"stages", json::array()}};
    for (auto& stage : stages) {
        const fs::path directory = root / stage.name;
        auto& cp = stage.params.get<CalculationParameters>();
        cp.set_user_defined_value("prefix", (directory / "state").string());
        const std::string state_id = run_id + ":" + stage.name;
        json result = {{"stage", stage.name}, {"application", stage.kind}, {"state_id", state_id},
                       {"identity_kind", "fresh_run_state"}, {"checkpoint_generation", nullptr},
                       {"requested_properties", stage.properties}, {"protocol", cp.protocol()},
                       {"output_directory", directory.string()}, {"properties", json::object()}};
        bool converged = false;
        if (stage.kind == "moldft") {
            moldft_lib library((directory / "resolved.in").string());
            auto values = library.run(world, stage.params, NextAction::Redo);
            const auto& properties = std::get<1>(values);
            const auto& evidence = std::get<2>(values);
            auto scf = library.calc(world, stage.params);
            scf->work_dir = directory;
            references.emplace(stage.name, scf);
            state_ids.emplace(stage.name, state_id);
            converged = evidence.converged_for_thresh <= cp.protocol().back() && evidence.converged_for_dconv <= cp.dconv();
            result["method"] = "HF";
            result["convergence"] = evidence.to_json();
            result["targets"] = {{"dconv", cp.dconv()}, {"econv", cp.econv()}};
            result["parameters"] = cp.to_json_if_precedence("defined");
            result["geometry"] = scf->molecule.to_json();
            result["representation"] = {{"L", scf->param.L()}, {"k", scf->amo.at(0).k()},
                                         {"threshold", cp.protocol().back()}};
            for (const auto& property : stage.properties) {
                if (property == "energy") result["properties"][property] = {{"value", properties.energy}, {"units", "hartree"}, {"source_state_id", state_id}};
                else {
                    std::vector<double> dipole;
                    dipole.reserve(3);
                    if (!properties.dipole || properties.dipole->size() != 3) throw std::runtime_error("missing ground-state dipole");
                    for (int i = 0; i < 3; ++i) dipole.push_back((*properties.dipole)(i));
                    result["properties"][property] = {{"value", dipole}, {"units", "e*a0"}, {"source_state_id", state_id}};
                }
            }
        } else {
            // SCF may orient/center the input geometry and rederive the cell.
            // Bind to the completed reference, not the pre-SCF input estimate.
            const auto& source = references.at(stage.reference);
            cp.set_user_defined_value("l", source->param.L());
            result["reference_representation"] = {{"L", source->param.L()}, {"k", source->amo.at(0).k()}};
            auto response = molresponse_v3_lib::run_response(world, stage.params, references.at(stage.reference),
                                                            directory, root / stage.reference / "state");
            result["method"] = "TDHF";
            result["reference"] = {{"stage", stage.reference}, {"state_id", state_ids.at(stage.reference)}};
            result["targets"] = {{"dconv", stage.params.get<ResponseParameters>().dconv()}};
            result["parameters"] = stage.params.get<ResponseParameters>().to_json_if_precedence("defined");
            // The engine's property/evidence JSON exists on root only.
            on_root(world, [&] {
                auto alpha = selected_alpha(stage, response.properties);
                converged = std::all_of(alpha.begin(), alpha.end(), [](const json& row) { return row.at("converged").get<bool>(); });
                result["properties"]["polarizability"] = {{"values", alpha}, {"source_state_id", state_id}, {"reference_state_id", state_ids.at(stage.reference)}};
                result["engine_metadata"] = response.metadata;
            });
            world.gop.broadcast(converged, 0);
        }
        result["status"] = converged ? "completed" : "unconverged";
        write_result(world, directory / "results.json", result);
        index["stages"].push_back(result);
        write_result(world, root / "results.json", index);
        if (!converged) throw std::runtime_error("stage " + stage.name + " did not satisfy the existing engine convergence test; see results.json");
    }
    if (world.rank() == 0) madness::print("Named-stage results:", (root / "results.json").string());
    return true;
}
} // namespace named_stages
