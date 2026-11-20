#include <madness/chem/SCFInput.h>
#include <madness/chem/ResponseParameters.hpp>
#include <madness/world/worldgop.h>
#include <iostream>
#include <sstream>

using namespace madness;
namespace {
int failures = 0;
void check(bool value, const char* message) { if (!value) { ++failures; std::cerr << message << '\n'; } }
template<class F> void rejects(F action, const char* message) {
    bool rejected = false;
    try { action(); } catch (const std::exception&) { rejected = true; }
    check(rejected, message);
}
}
int main(int argc, char** argv) {
    initialize(argc, argv);
    {
        World world(SafeMPI::COMM_WORLD);
        CalculationParameters p;
        p.enable_scf_inputs();
        const auto printed = p.print_to_string();
        for (const auto* key : {"dconv", "econv", "accuracy.bsh_residual", "progress.window",
                               "convergence.use_default_criteria", "restartao", "no_compute"})
            check(printed.find(std::string(key) + " ") == std::string::npos,
                  "retired SCF input is still advertised");
        ResponseParameters response;
        const auto response_listing = response.print_to_string();
        for (const auto* key : {"accuracy.bsh_residual", "accuracy.kinetic_trace_change",
                 "convergence.use_default_criteria", "progress.window", "progress.minimum_relative_improvement",
                 "diagnostic.iterations_per_rung", "stall.window", "stall.ratio",
                 "seed.start_rung", "seed.es_warmup", "excited.maxiter"}) {
            check(response_listing.find(std::string(key) + " ") == std::string::npos,
                  "retired response input is still advertised");
            rejects([&] { response.validate_input_value(key, "0"); }, "retired response input accepted");
            std::vector<std::string> args{"test_scf_inputs", std::string("--response=") + key + "=0",
                                           "--response=maxiter=2"};
            std::vector<char*> pointers;
            pointers.reserve(args.size());
            for (auto& arg : args) pointers.push_back(arg.data());
            commandlineparser input(static_cast<int>(pointers.size()), pointers.data());
            bool rejected = false;
            try { response.read_commandline_options(world, input, "response"); }
            catch (const std::exception& error) {
                rejected = std::string(error.what()).find(std::string("unknown input: ") + key) != std::string::npos;
            }
            check(rejected, "discarded retired response occurrence must fail collectively as unknown input");
        }
        // A complete parameter listing must preserve default/derived versus
        // explicit controls when copied into a deck, including unset tolerances.
        struct ScfDeck : CalculationParameters {
            ScfDeck() { enable_scf_inputs(); }
            void read(World& w, std::string text) { read_internal(w, text, "dft"); }
        };
        struct ResponseDeck : ResponseParameters {
            void read(World& w, std::string text) { read_internal(w, text, "response"); }
        };
        ScfDeck scf_deck, scf_restored;
        scf_deck.set_user_defined_value("maxiter", 7);
        scf_deck.set_derived_value("prefix", std::string("derived-prefix"));
        try {
            scf_restored.read(world, "dft\n" + scf_deck.print_to_string() + "\nend\n");
            check(scf_restored.maxiter() == 7 && scf_restored.is_user_defined("maxiter"), "explicit SCF setting lost");
            check(!scf_restored.is_user_defined("convergence.criteria") && !scf_restored.is_user_defined("prefix"),
                  "default or derived SCF setting became an explicit request");
        } catch (const std::exception&) { check(false, "actual SCF parameter listing cannot be read"); }
        ResponseDeck response_deck, response_restored;
        response_deck.set_user_defined_value("maxiter", std::size_t(7));
        try {
            response_restored.read(world, "response\n" + response_deck.print_to_string() + "\nend\n");
            check(response_restored.maxiter() == 7 && response_restored.is_user_defined("maxiter"), "explicit response setting lost");
            check(!response_restored.is_user_defined("convergence.criteria") &&
                  !response_restored.is_user_defined("checkpoint.every_iterations"),
                  "default listing enabled the shared response path");
        } catch (const std::exception&) { check(false, "actual response parameter listing cannot be read"); }
        auto request = scf_iteration_request(p);
        check(request.criteria.size() == 3 && request.max_updates == 25 && request.checkpoint_interval == 0,
              "ordinary defaults select three metrics and endpoint checkpoints");
        p.set_user_defined_value("freeze_occupied", true);
        check(scf_iteration_request(p).criteria == std::vector<iteration::Metric>{iteration::Metric::bsh_residual},
              "frozen default only selects the applicable residual");
        p.set_user_defined_value("convergence.criteria", std::vector<std::string>{"bsh_residual", "kinetic_change"});
        p.set_user_defined_value("convergence.bsh_residual.relative", 1e-4);
        check(!scf_iteration_request(p).tolerances.at(iteration::Metric::bsh_residual).absolute,
              "relative-only override retains the model absolute default");
        rejects([&] { p.validate_input_value("dconv", "1e-5"); }, "dconv is retired for SCF");
        rejects([&] { p.validate_input_value("save", "true"); }, "save is retired for SCF");
        rejects([&] { p.validate_input_value("convergence.energy_change", "nan"); }, "nonfinite tolerance rejected");
        rejects([&] { p.validate_input_value("checkpoint.every_iterations", "-1"); }, "negative cadence rejected");
        std::vector<std::string> arguments{"test_scf_inputs", "--dft=dconv=1e-4", "--dft=maxiter=2"};
        std::vector<char*> argument_pointers;
        argument_pointers.reserve(arguments.size());
        for (auto& argument : arguments) argument_pointers.push_back(argument.data());
        commandlineparser parser(static_cast<int>(argument_pointers.size()), argument_pointers.data());
        rejects([&] { p.read_commandline_options(world, parser, "dft"); }, "discarded retired occurrence must be rejected");
        CalculationParameters legacy;
        legacy.validate_input_value("dconv", "1e-4");
        rejects([&] { legacy.validate_input_value("convergence.energy_change", "1e-5"); },
                "unadopted engines cannot silently accept shared SCF controls");
        world.gop.max(failures);
    }
    finalize();
    return failures != 0;
}
