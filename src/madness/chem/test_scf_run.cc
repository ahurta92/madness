#include <madness/chem/SCFRun.h>
#include <madness/chem/SCFCheckpoint.h>
#include <madness/chem/RootCoordination.h>
#include <iostream>
using namespace madness;
namespace {
// This process's launch entry. Every SCF run in one process is part of one
// launch; each run adds its own top-level steps, tagged with its part number.
checkpoint::json timing_records(World& world, const std::filesystem::path& root) {
    const auto launches = checkpoint::json::parse(checkpoint::detail::on_root(world, [&] {
        MADNESS_CHECK_THROW(std::filesystem::is_regular_file(root / "timings.json"),
                            "SCF must persist shared timings even with print_level zero");
        return checkpoint::detail::read_json(root / "timings.json").dump();
    })).at("launches");
    MADNESS_CHECK_THROW(launches.size() == 1 && launches[0].contains("started_at") &&
                        !launches[0].contains("attempt_id"), "one process writes one launch entry");
    return launches[0];
}
checkpoint::json steps(const checkpoint::json& launch) { return launch.at("timing").at("children"); }
// Recursively inspect the finite timing tree saved by the application.
// NOLINTNEXTLINE(misc-no-recursion)
void collect_iterations(const checkpoint::json& node, checkpoint::json& iterations) {
    if (node.at("name") == "iteration") iterations.push_back(node);
    for (const auto& child : node.value("children", checkpoint::json::array())) collect_iterations(child, iterations);
}
}
int main(int argc, char** argv) {
    initialize(argc, argv);
    int failed = 0;
    {
        World world(SafeMPI::COMM_WORLD);
        startup(world, argc, argv, true);
        try {
            const auto scratch = checkpoint::detail::on_root(world, [] {
                char pattern[] = "/tmp/test-scf-run-XXXXXX";
                const auto directory = mkdtemp(pattern);
                if (!directory) throw std::runtime_error("cannot create SCF test directory");
                return std::string(directory);
            });
            Molecule molecule; molecule.add_atom(0.1, 0, 0, 2, 2);
            CalculationParameters p; p.enable_scf_inputs();
            p.set_user_defined_value("protocol", std::vector<double>{1e-4});
            p.set_user_defined_value("localize", std::string("canon"));
            p.set_user_defined_value("print_level", 0);
            // One update with loose change targets, so a later loose BSH target
            // is judged on saved evidence for every default criterion.
            p.set_user_defined_value("maxiter", 1);
            p.set_user_defined_value("convergence.bsh_residual", 0.0);
            p.set_user_defined_value("convergence.density_change", 1000.0);
            p.set_user_defined_value("convergence.energy_change", 1000.0);
            SCF scf(world, p, molecule); scf.work_dir = scratch;
            MolecularEnergy energy(world, scf);
            energy.value(molecule.get_all_coords().flat());
            MADNESS_CHECK_THROW(!scf.last_calculation_outcome->complete, "zero cap must remain incomplete");
            const auto ground = ground_state_node(scf.work_dir);
            const auto before = read_node(world, ground);
            MADNESS_CHECK_THROW(before.has_value(), "endpoint must save its metrics and native functions");
            const auto initial_times = steps(timing_records(world, scratch));
            MADNESS_CHECK_THROW(initial_times.size() == 1 && initial_times[0].at("part") == 0,
                                "fresh SCF timing record is missing");
            auto iterations = checkpoint::json::array();
            collect_iterations(initial_times.back(), iterations);
            MADNESS_CHECK_THROW(iterations.size() == 2 && iterations[0].at("iteration") == 0 &&
                                iterations[1].at("iteration") == 1,
                                "one-update endpoint has two evaluation timings");
            bool advanced = false;
            for (const auto& step : iterations[0].at("children")) advanced |= step.at("name") == "advance";
            MADNESS_CHECK_THROW(advanced, "the update after the first evaluation must be timed");
            for (const auto& step : iterations[1].at("children"))
                MADNESS_CHECK_THROW(step.at("name") != "advance", "terminal check must not report an update");
            scf.param.set_user_defined_value("convergence.bsh_residual", 1000.0);
            energy.value(molecule.get_all_coords().flat());
            MADNESS_CHECK_THROW(scf.last_calculation_outcome->complete, "satisfied saved metrics must be reusable");
            MADNESS_CHECK_THROW(read_node(world, ground) == before, "reuse must not save again");
            const auto reused_times = steps(timing_records(world, scratch));
            MADNESS_CHECK_THROW(reused_times.size() == 2 && reused_times[0] == initial_times[0] &&
                                reused_times[1].at("part") == 1, "reuse must retain the original timing record");
            iterations = checkpoint::json::array();
            collect_iterations(reused_times.back(), iterations);
            MADNESS_CHECK_THROW(iterations.empty(), "reuse must not copy historical iterations");
            auto moved = molecule.get_all_coords(); moved(0, 0) = 0; moved(0, 1) = 0.1;
            bool rejected = false;
            std::string model_failure;
            try { energy.value(moved.flat()); }
            catch (const std::exception& error) { rejected = true; model_failure = error.what(); }
            MADNESS_CHECK_THROW(rejected, "ordinary restart must reject a changed geometry");
            MADNESS_CHECK_THROW(read_node(world, ground) == before, "model rejection must preserve the saved state");
            const auto failed_launch = timing_records(world, scratch);
            const auto failed_times = steps(failed_launch);
            MADNESS_CHECK_THROW(failed_times.size() == 3 && failed_times[0] == initial_times[0] &&
                                failed_times[1] == reused_times[1] && failed_launch.contains("failure"),
                                "failed SCF must retain previous records and record its failure");
            // A secondary timing-file error must not hide the collective model error.
            checkpoint::detail::on_root(world, [&] {
                std::ofstream(std::filesystem::path(scratch) / "timings.json") << "{}";
                return std::string{};
            });
            std::string preserved_failure;
            try { energy.value(moved.flat()); }
            catch (const std::exception& error) { preserved_failure = error.what(); }
            checkpoint::detail::on_root(world, [&] {
                const auto file = std::filesystem::path(scratch) / "timings.json";
                MADNESS_CHECK_THROW(checkpoint::detail::read_json(file).empty(),
                                    "failed timing publication replaced the existing record");
                std::ofstream(file) << checkpoint::json{{"version", 2}, {"launches", {failed_launch}}}.dump();
                return std::string{};
            });
            MADNESS_CHECK_THROW(preserved_failure == model_failure,
                                "timing publication masked the original model rejection");
            scf.param.set_user_defined_value("protocol", std::vector<double>{1e-6});
            energy.value(molecule.get_all_coords().flat());
            MADNESS_CHECK_THROW(scf.amo.front().k() == 8 && scf.amo.front().thresh() == 1e-6,
                                "native restart converts to the requested final representation");
            // Raw archives are unqualified guesses. Their occupation holes
            // must not survive a new request for default occupations.
            scf.aocc[0] = 0.0;
            const auto seed = (std::filesystem::path(scratch) / "seed").string();
            MADNESS_CHECK_THROW(scf.write_mos(world, seed + ".restartdata", scf.orbital_metadata()).empty(),
                                "seed fixture write failed");
            auto seeded_parameters = p;
            seeded_parameters.set_user_defined_value("prefix", seed);
            seeded_parameters.set_user_defined_value("convergence.bsh_residual", 1000.0);
            SCF seeded(world, seeded_parameters, molecule);
            seeded.work_dir = std::filesystem::path(scratch) / "seeded";
            run_scf(world, seeded);
            MADNESS_CHECK_THROW(seeded.aocc[0] == 1.0, "raw seed retained a removed occupation hole");
            if (world.rank() == 0) std::filesystem::remove_all(scratch);
        } catch (const std::exception& error) { std::cerr << error.what() << '\n'; failed = 1; }
        world.gop.max(failed);
    }
    finalize(); return failed;
}
