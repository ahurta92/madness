#include <madness/chem/SCF.h>
#include <madness/chem/CheckpointRecords.h>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>

using namespace madness;
namespace it = madness::iteration;
namespace fs = std::filesystem;

// Numerical calibration and regression against curated, unchanged QC references.
// Optional explicit tolerances are experiment controls for this test executable,
// not new madqc input aliases. Normal CTest execution uses the model defaults.
int main(int argc, char** argv) {
    initialize(argc, argv);
    int failures = 0;
    {
        World world(SafeMPI::COMM_WORLD);
        startup(world, argc, argv, true);
        commandlineparser options(argc, argv);
        const std::string name = options.key_exists("case") ? options.value("case") : "scf_he_hf";
        const fs::path source = fs::path(MAD_ROOT_DIR) / "src/examples/qc" / name;
        const fs::path original = fs::current_path();
        const fs::path output = options.key_exists("evidence") ? fs::absolute(options.value_raw("evidence")) : fs::path();
        std::string directory;
        if (world.rank() == 0) {
            char pattern[] = "/tmp/madness-scf-defaults-XXXXXX";
            if (const char* made = mkdtemp(pattern)) directory = made;
        }
        world.gop.broadcast_serializable(directory, 0);
        MADNESS_CHECK_THROW(!directory.empty(), "cannot create calibration directory");
        fs::current_path(directory);
        nlohmann::json evidence{{"case", name}, {"history", nlohmann::json::array()},
                                {"protocol_decisions", nlohmann::json::array()}};
        try {
            commandlineparser input;
            input.set_keyval_keep_case("input", (source / (name + ".in")).string());
            CalculationParameters parameters;
            parameters.enable_scf_inputs();
            parameters.read_input_and_commandline_options(world, input, "dft");
            parameters.set_user_defined_value("prefix", name);
            parameters.set_user_defined_value("print_level", 0);
            parameters.set_user_defined_value("maxiter", 25);
            Molecule molecule(world, input);
            SCF scf(world, parameters, molecule);
            auto protocol = scf.param.protocol();
            if (options.key_exists("protocol")) {
                protocol = {std::stod(options.value("protocol"))};
                scf.param.set_user_defined_value("protocol", protocol);
            }
            it::Request request;
            request.criteria = {it::Metric::bsh_residual, it::Metric::density_change, it::Metric::energy_change};
            if (options.key_exists("kinetic-required")) request.criteria.push_back(it::Metric::kinetic_change);
            request.max_updates = 25;
            for (auto metric : {it::Metric::bsh_residual, it::Metric::density_change,
                                it::Metric::energy_change, it::Metric::kinetic_change}) {
                const std::string key = it::metric_name(metric);
                if (options.key_exists(key)) request.tolerances[metric] = {std::stod(options.value(key)), 0.0};
            }
            scf.set_protocol<3>(world, protocol.front());
            // Fail before costly guess construction if a default is missing.
            it::resolve(request, scf.iteration_model(), {protocol.front(), FunctionDefaults<3>::get_k()}, protocol.size() == 1);
            scf.ao = scf.project_ao_basis(world, scf.aobasis);
            scf.make_nuclear_potential(world);
            scf.initial_guess(world);
            it::Decision last;
            std::size_t total_updates = 0;
            const auto started = std::chrono::steady_clock::now();
            for (std::size_t stage = 0; stage < protocol.size(); ++stage) {
                scf.set_protocol<3>(world, protocol[stage]);
                if (stage) scf.project(world);
                scf.make_nuclear_potential(world);
                scf.reset_aobasis("sto-3g");
                scf.ao.clear(); world.gop.fence();
                scf.ao = scf.project_ao_basis(world, scf.aobasis);
                auto settings = it::resolve(request, scf.iteration_model(),
                    {protocol[stage], FunctionDefaults<3>::get_k()}, stage + 1 == protocol.size());
                SCF::IterationControl control{settings, {"scf", name, {"ground"}}, "calibration"};
                it::IterationMetrics latest;
                control.observer = [&](World&, const SCF&, const it::IterationMetrics& metrics,
                                       const it::Decision& decision, const it::Progress&) {
                    latest = metrics;
                    evidence["history"].push_back({{"threshold", metrics.info.protocol.thresh},
                        {"k", metrics.info.protocol.k}, {"iteration", metrics.info.iteration},
                        {"reason", decision.reason},
                        {"metrics", checkpoint::records::measurements(metrics, SCF::iteration_metric_set())}});
                };
                last = scf.solve(world, control);
                total_updates += last.info.iteration;
                // Harness fault injection: a genuine controller stop must not
                // be hidden by a successful solve at the following protocol.
                if (stage == 0 && options.key_exists("invalid-stage")) {
                    evidence["injected_invalid_state"] = true;
                    last = it::decide(settings, control.task, {last.info.iteration, last.info.iteration, false}, latest);
                }
                evidence["protocol_decisions"].push_back({{"threshold", protocol[stage]},
                    {"action", static_cast<int>(last.action)}, {"reason", last.reason}});
                if (stage + 1 < protocol.size() && last.action != it::Action::finish_protocol) break;
            }
            evidence["solve_seconds"] = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
            evidence["updates"] = total_updates;
            evidence["complete"] = last.action == it::Action::complete;
            evidence["energy"] = scf.current_energy;
            evidence["metric_definitions"] = SCF::iteration_metric_definitions();
            evidence["tolerances"] = nlohmann::json::object();
            for (const auto& check : last.checks)
                evidence["tolerances"][it::metric_name(check.metric)] = {{"absolute", check.absolute}, {"relative", check.relative}};
            nlohmann::json reference;
            std::ifstream(source / "reference" / (name + ".calc_info.json")) >> reference;
            const double reference_energy = reference.at("tasks").at(0).at("properties").at("energy");
            evidence["reference_energy"] = reference_energy;
            evidence["reference_difference"] = scf.current_energy - reference_energy;
            // Overrides of the representation are calibration experiments, not
            // claims that a different discretization matches a coarse oracle.
            const bool same_protocol = !options.key_exists("protocol") && !options.key_exists("invalid-stage");
            const bool expected_limit = options.key_exists("expect-limit");
            bool expected_outcome = last.action == it::Action::complete;
            if (options.key_exists("invalid-stage"))
                expected_outcome = last.action == it::Action::stop && last.reason == "invalid_state" &&
                                   evidence["protocol_decisions"].size() == 1;
            if (expected_limit) {
                expected_outcome = last.action == it::Action::stop && last.reason == "iteration_limit";
                bool residual_failed = false;
                for (const auto& check : last.checks)
                    if (check.metric == it::Metric::bsh_residual)
                        residual_failed = check.result == it::CheckResult::unmet;
                expected_outcome = expected_outcome && residual_failed;
            }
            if (!expected_outcome ||
                (same_protocol && std::abs(scf.current_energy - reference_energy) > 1e-5)) {
                ++failures;
                std::cerr << "SCF default/reference check failed: " << evidence["complete"]
                          << ", energy difference " << evidence["reference_difference"] << '\n';
            }
        } catch (const std::exception& error) {
            ++failures; evidence["error"] = error.what(); std::cerr << error.what() << '\n';
        }
        if (world.rank() == 0 && !output.empty()) {
            try {
                std::ofstream file(output);
                file.exceptions(std::ios::failbit | std::ios::badbit);
                file << evidence.dump(2) << '\n';
            } catch (const std::exception& error) {
                ++failures;
                std::cerr << "cannot write calibration evidence: " << error.what() << '\n';
            }
        }
        world.gop.fence();
        fs::current_path(original);
        if (world.rank() == 0) {
            std::error_code error;
            fs::remove_all(directory, error);
            if (error) { ++failures; std::cerr << "cannot remove calibration directory: " << error.message() << '\n'; }
        }
        world.gop.max(failures);
    }
    finalize();
    return failures ? 1 : 0;
}
