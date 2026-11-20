#include <madness/mra/QCCalculationParametersBase.h>
#include <madness/world/worldinit.h>
#include <madness/world/worldgop.h>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <cstdlib>

using namespace madness;
namespace {
struct Parameters : QCCalculationParametersBase {
    Parameters() {
        initialize<int>("maxiter", 7, "iteration allowance");
        initialize<std::string>("kind", "first", "algorithm", {"first", "second"});
        initialize<std::string>("path", "/Original/Directory", "case-sensitive path");
        initialize<double>("legacy_tol", 1e-4, "retired key fixture");
        initialize<double>("accuracy", 1e-4, "type-check fixture");
        ignore_unknown_keys = false;
    }
    [[nodiscard]] std::string get_tag() const override { return "calc"; }
    using QCCalculationParametersBase::read_input_and_commandline_options;
    [[nodiscard]] std::set<std::string> case_sensitive_keys() const override { return {"path"}; }
    void validate_input_value(const std::string& key, const std::string& value) const override {
        if (key == "legacy_tol") throw std::invalid_argument("legacy_tol is retired; use the explicit metric name");
        if (key == "maxiter" && std::stoi(value) < 0) throw std::invalid_argument("maxiter must be nonnegative");
    }
};
commandlineparser arguments(std::vector<std::string> values) {
    values.insert(values.begin(), "test_cli_occurrences");
    std::vector<char*> argv;
    argv.reserve(values.size());
    for (auto& value : values) argv.push_back(value.data());
    return {static_cast<int>(argv.size()), argv.data()};
}
}
int main(int argc, char** argv) {
    initialize(argc, argv);
    int failures = 0;
    {
        World world(SafeMPI::COMM_WORLD);
        auto check = [&](bool ok, const char* why) {
            if (!ok) { ++failures; std::cerr << why << '\n'; }
        };
        auto rejects = [&](const commandlineparser& parser, const char* why) {
            Parameters p;
            int type = 0;
            std::string error;
            try { p.read_commandline_options(world, parser, "calc"); }
            catch (const std::invalid_argument& e) { type = 1; error = e.what(); }
            catch (const std::runtime_error& e) { type = 2; error = e.what(); }
            int root_type = type;
            std::string root_error = error;
            world.gop.broadcast(root_type, 0);
            world.gop.broadcast_serializable(root_error, 0);
            check(type != 0 && type == root_type && !error.empty() && error == root_error, why);
            check(p.get<int>("maxiter") == 7, "rejected CLI input leaves the current settings intact");
        };
        rejects(arguments({"--calc=maxiter bad", "--calc=maxiter 5"}), "invalid discarded integer must fail on every rank");
        rejects(arguments({"--calc=accuracy 1e-4oops", "--calc=accuracy 1e-5"}), "discarded malformed float must fail through registered type checking");
        rejects(arguments({"--calc=typo 2", "--calc=maxiter 5"}), "unknown discarded key must fail on every rank");
        rejects(arguments({"--calc=kind bad", "--calc=kind second"}), "invalid discarded enum must fail on every rank");
        rejects(arguments({"--calc=legacy_tol 1e-3", "--calc=maxiter 5"}), "discarded retired key must fail through its owning hook");
        rejects(arguments({"--calc=maxiter -2", "--calc=maxiter 5"}), "discarded invalid bound must fail through its owning hook");
        rejects(arguments({"--calc=maxiter bad"}), "winning invalid input must fail collectively");
        rejects(arguments({"--calc=maxiter 2; kind bad"}), "partially read invalid group must not change settings");
        {
            auto parser = arguments({"--calc=maxiter 2", "--calc=maxiter 5", "--flag"});
            parser.set_keyval("calc", "maxiter 9");
            const auto& original = parser.original_arguments();
            check(original.size() == 3 && original[0].value == "maxiter 2" &&
                  original[1].value == "maxiter 5" && original[2].key == "flag" &&
                  original[2].value.empty(), "original arguments survive programmatic overrides without invented occurrences");
            Parameters p;
            p.read_commandline_options(world, parser, "calc");
            check(p.get<int>("maxiter") == 9, "programmatic final value retains precedence");
        }
        {
            auto parser = arguments({"--calc=maxiter bad"});
            parser.set_keyval("calc", "maxiter 5");
            rejects(parser, "programmatic overrides cannot hide invalid human input");
        }
        {
            Parameters p;
            p.set_user_defined_value("maxiter", 9); // existing file value
            p.set_user_defined_value("path", std::string("/Original/Directory"));
            auto parser = arguments({"--calc=maxiter 2; path /Discarded/Path", "--calc=kind second"});
            p.read_commandline_options(world, parser, "calc");
            check(p.get<int>("maxiter") == 9 && p.get<std::string>("path") == "/Original/Directory" &&
                  p.get<std::string>("kind") == "second", "last group replaces earlier groups without leaking discarded values");
        }
        {
            Parameters p;
            auto parser = arguments({"--calc=maxiter 2", "--calc=maxiter 5; path /Case/Sensitive"});
            p.read_commandline_options(world, parser, "calc");
            check(p.get<int>("maxiter") == 5 && p.get<std::string>("path") == "/Case/Sensitive",
                  "winning CLI overrides preserve path case");
        }
        // Rank zero owns parsing. Other ranks must receive the same failure even
        // when they have no local copy of the group, not exit before a collective.
        if (world.size() > 1) rejects(world.rank() == 0 ? arguments({"--calc=maxiter bad"}) : commandlineparser(),
            "rank-zero-only argument errors reach every rank");
        {
            std::string filename;
            if (world.rank() == 0) {
                char pattern[] = "/tmp/madness-cli-occurrences-XXXXXX";
                if (const char* directory = mkdtemp(pattern)) {
                    filename = std::string(directory) + "/input";
                    std::ofstream(filename) << "calc\nmaxiter 2\nmaxiter -1\nend\n";
                }
            }
            world.gop.broadcast_serializable(filename, 0);
            MADNESS_CHECK_THROW(!filename.empty(), "cannot create CLI test input");
            Parameters p;
            commandlineparser parser;
            parser.set_keyval_keep_case("input", filename);
            std::string error;
            try { p.read_input_and_commandline_options(world, parser, "calc"); }
            catch (const std::invalid_argument& e) { error = e.what(); }
            check(error == "maxiter must be nonnegative" && p.get<int>("maxiter") == 7,
                  "file input uses the same collective validation and leaves settings intact on failure");
            if (world.rank() == 0) std::ofstream(filename) << "calc\nmaxiter 2\nend\n";
            world.gop.fence();
            error.clear();
            try {
                auto combined = arguments({"--input=" + filename, "--calc=maxiter bad"});
                p.read_input_and_commandline_options(world, combined, "calc");
            } catch (const std::invalid_argument& e) { error = e.what(); }
            check(!error.empty() && p.get<int>("maxiter") == 7,
                  "combined file and CLI failure rolls back the entire read on every rank");
            error.clear();
            const std::string missing = filename + ".missing";
            try {
                auto root_only = world.rank() == 0 ? arguments({"--input=" + missing}) : commandlineparser();
                p.read_input_and_commandline_options(world, root_only, "calc");
            } catch (const std::invalid_argument& e) { error = e.what(); }
            check(error == "could not find user-defined input file: " + missing + "\n",
                  "combined read uses rank-zero input selection and identical missing-file errors");
            auto root_only_valid = world.rank() == 0
                ? arguments({"--input=" + filename, "--calc=maxiter 5; path /Root/Chosen"})
                : commandlineparser();
            p.read_input_and_commandline_options(world, root_only_valid, "calc");
            check(p.get<int>("maxiter") == 5 && p.get<std::string>("path") == "/Root/Chosen",
                  "successful combined input uses rank-zero sources and broadcasts resolved settings");
            world.gop.fence();
            if (world.rank() == 0) std::filesystem::remove_all(std::filesystem::path(filename).parent_path());
        }
        world.gop.max(failures);
    }
    finalize();
    return failures ? 1 : 0;
}
