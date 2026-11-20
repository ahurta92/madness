#include <madness/chem/SCFRun.h>
#include <madness/chem/SCFCheckpoint.h>
#include <apps/molresponse/GroundState.hpp>
#include <iostream>

using namespace madness;
int main(int argc, char** argv) {
    initialize(argc, argv);
    int failed = 0;
    {
        World world(SafeMPI::COMM_WORLD);
        startup(world, argc, argv, true);
        try {
            const auto scratch = checkpoint::detail::on_root(world, [] {
                char pattern[] = "/tmp/shared-reference-XXXXXX";
                const auto directory = mkdtemp(pattern);
                if (!directory) throw std::runtime_error("cannot create reference test directory");
                return std::string(directory);
            });
            for (bool charged : {false, true}) {
                CalculationParameters p; p.enable_scf_inputs();
                p.set_user_defined_value("protocol", std::vector<double>{1e-4});
                p.set_user_defined_value("localize", std::string("canon"));
                p.set_user_defined_value("aobasis", std::string("6-31gss"));
                p.set_user_defined_value("nvalpha", charged ? 0 : 1);
                p.set_user_defined_value("charge", charged ? 1.0 : 0.0);
                p.set_user_defined_value("nopen", charged ? 1 : 0);
                p.set_user_defined_value("spin_restricted", !charged);
                p.set_user_defined_value("maxiter", 0);
                p.set_user_defined_value("print_level", 0);
                p.set_user_defined_value("convergence.criteria", std::vector<std::string>{"bsh_residual"});
                p.set_user_defined_value("convergence.bsh_residual", 1000.0);
                Molecule molecule; molecule.add_atom(0, 0, 0, 2, 2);
                SCF scf(world, p, molecule);
                scf.work_dir = std::filesystem::path(scratch) / (charged ? "charged" : "virtuals");
                run_scf(world, scf);
                const auto archive = scf_checkpoint_archive(world, scf.work_dir, false);
                const auto reference = molresponse_v3::GroundState::from_archive(
                    world, archive.string(), scf.molecule, &scf.param, &scf.pcm_param);
                MADNESS_CHECK_THROW(reference.num_alpha() == 1, "response must use the occupied prefix");
                MADNESS_CHECK_THROW(reference.params().charge() == (charged ? 1.0 : 0.0),
                                    "response must retain the verified charge, not infer it from bundle size");
                MADNESS_CHECK_THROW(reference.is_spin_restricted() == !charged, "response spin layout changed");
                if (charged) MADNESS_CHECK_THROW(reference.num_beta() == 0, "empty beta channel changed");
                if (charged) {
                    Molecule core; core.add_atom(0, 0, 0, 4, 4);
                    const char* old = std::getenv("MRA_DATA_DIR");
                    const auto old_path = old ? std::optional<std::string>(old) : std::nullopt;
                    setenv("MRA_DATA_DIR", MAD_ROOT_DIR "/src/madness/chem", 1);
                    core.read_core_file("mcp");
                    if (old_path) setenv("MRA_DATA_DIR", old_path->c_str(), 1);
                    else unsetenv("MRA_DATA_DIR");
                    auto beta_virtuals = scf.param;
                    beta_virtuals.set_user_defined_value("nbeta", 1);
                    beta_virtuals.set_user_defined_value("nmo_beta", 2);
                    bool guarded = false;
                    try {
                        (void)molresponse_v3::GroundState::from_archive(world, archive.string(), core,
                            &beta_virtuals, &scf.pcm_param);
                    } catch (const MadnessException& error) {
                        guarded = std::string(error.msg).find("native core-orbital trimming") != std::string::npos;
                    }
                    MADNESS_CHECK_THROW(guarded, "beta-only core trimming must fail before native reading");
                }
            }
            if (world.rank() == 0) std::filesystem::remove_all(scratch);
        } catch (const std::exception& error) { std::cerr << error.what() << '\n'; failed = 1; }
        world.gop.max(failed);
    }
    finalize(); return failed;
}
