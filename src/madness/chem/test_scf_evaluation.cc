#include <madness/chem/SCF.h>
#include <madness/mra/mra.h>
#include <filesystem>
#include <cstdlib>
#include <fstream>
#include <iostream>

using namespace madness;
namespace {
int failures = 0;
void check(bool condition, const char* message) {
    if (!condition) { ++failures; std::cerr << "FAILED: " << message << '\n'; }
}
void current_state_case(World& world, const std::string& name, const std::string& xc,
                        int atom, const std::string& localization, bool immediate = false) {
    Molecule molecule;
    molecule.add_atom(0.0, 0.0, 0.0, double(atom), atom);
    if (atom == 8) {
        molecule.add_atom(0.0, 1.421178283, -1.064386617, 1.0, 1);
        molecule.add_atom(0.0, -1.421178283, -1.064386617, 1.0, 1);
    }
    CalculationParameters p;
    p.set_user_defined_value("xc", xc);
    if (atom == 3) p.set_user_defined_value("nopen", 1);
    p.set_user_defined_value("protocol", std::vector<double>{1e-4});
    p.set_user_defined_value("k", 6);
    p.set_user_defined_value("l", 20.0);
    p.set_user_defined_value("localize", localization);
    p.set_user_defined_value("maxiter", 1);
    p.set_user_defined_value("print_level", name == "water_localized" ? 2 : 0);
    p.set_user_defined_value("save", true);
    p.set_user_defined_value("prefix", name);
    if (immediate) {
        p.set_user_defined_value(convergence::use_default_criteria_key, false);
        p.set_user_defined_value("accuracy.bsh_residual", 10.0);
    }
    SCF scf(world, p, molecule);
    std::vector<SCF::IterationCheckpoint> samples;
    scf.checkpoint_observer = [&](World& w, const SCF&, const SCF::IterationCheckpoint& step) {
        const auto before = copy(w, scf.amo);
        auto measured = scf.evaluate_orbitals(w, true);
        const auto raw = scf.evaluate_residual(w, measured);
        check(norm2(w, sub(w, scf.amo, before)) < 1e-10, "evaluation preserves orbitals");
        check(std::abs(step.input_energy - measured.energy) < 1e-9,
              "reported energy belongs to observed orbitals");
        check(std::abs(step.input_kinetic - measured.kinetic) < 1e-9,
              "reported kinetic energy belongs to observed orbitals");
        check(std::abs(step.bsh_residual - raw.maximum) < 1e-9,
              "reported residual belongs to observed orbitals");
        if (xc == "hf" && atom != 3) {
            // Existing independent closed-shell evaluator guards the extracted
            // energy, kinetic and fixed-point definitions as well as timing.
            const auto reference = scf.evaluate_restricted_hf(w);
            check(std::abs(reference.energy - measured.energy) < 1e-9, "HF evaluator agrees with existing energy definition");
            check(std::abs(reference.bsh_residual - raw.maximum) < 1e-9, "HF evaluator agrees with existing residual definition");
        }
        if (samples.empty()) check(!step.has_previous, "initial changes have no fabricated comparison");
        else {
            check(step.has_previous && step.comparison_iteration == samples.back().iteration,
                  "endpoint comparison is the preceding iteration");
            check(std::abs(step.energy_change - std::abs(measured.energy - samples.back().input_energy)) < 1e-9,
                  "energy change retains the preceding iteration through final rotations");
            check(std::abs(step.kinetic_change - std::abs(measured.kinetic - samples.back().input_kinetic)) < 1e-9,
                  "kinetic change retains the preceding iteration through final rotations");
        }
        samples.push_back(step);
    };
    MolecularEnergy energy(world, scf);
    energy.value(molecule.get_all_coords().flat());
    check(samples.size() == (immediate ? 1u : 2u), "evaluations count only allowed updates");
    check(!samples.empty() && samples.back().final_rotation, "capped state is prepared and evaluated");
    check(!samples.empty() && samples.back().iteration == (immediate ? 0 : 1), "terminal evaluation consumes no update");
    if (immediate) check(scf.last_calculation_outcome->complete, "satisfied initial residual needs no KAIN update");
    const auto final_energy = scf.current_energy;
    scf.load_mos(world);
    auto saved = scf.evaluate_orbitals(world);
    const auto saved_residual = scf.evaluate_residual(world, saved);
    check(std::abs(final_energy - saved.energy) < 1e-9, "returned energy belongs to native saved orbitals");
    check(!samples.empty() && std::abs(samples.back().input_energy - saved.energy) < 1e-9,
          "terminal metrics describe reloaded archive");
    check(std::abs(samples.back().bsh_residual - saved_residual.maximum) < 1e-9,
          "terminal residual describes reloaded archive");
    if (atom == 3) {
        for (long i = 0; i < scf.aeps.size(); ++i)
            check(std::abs(scf.aeps[i] - saved.fock_alpha(i, i)) < 1e-9,
                  "saved alpha energies describe alpha orbitals independently of beta");
        for (long i = 0; i < scf.beps.size(); ++i)
            check(std::abs(scf.beps[i] - saved.fock_beta(i, i)) < 1e-9,
                  "saved beta energies describe beta orbitals");
    }
    if (name == "water_localized") {
        tensorT rotation, spectrum;
        sygvp(world, saved.fock_alpha, matrix_inner(world, scf.amo, scf.amo, true),
              1, rotation, spectrum);
        std::ifstream input(name + ".energies.dat");
        std::string header;
        std::getline(input, header);
        for (long i = 0; i < spectrum.size(); ++i) {
            double value = 0.0;
            check(bool(input >> value) && std::abs(value - spectrum[i]) < 6e-8,
                  "localized DOS file contains the final Fock eigenvalue spectrum");
        }
    }
    if (name == "he_hf") {
        // An intentionally modified cached residual proves advancement uses
        // that supplied work rather than evaluating the BSH equation again.
        auto cached = saved_residual;
        cached.functions = copy(world, cached.functions);
        for (auto& function : cached.functions) function.scale(0.5);
        auto expected = sub(world, scf.amo, cached.functions);
        scf.do_step_restriction(world, scf.amo, expected, "test");
        scf.orthonormalize(world, expected, scf.param.nalpha());
        subspaceT history;
        tensorT coefficients;
        scf.update_subspace(world, cached, history, coefficients);
        check(norm2(world, sub(world, scf.amo, expected)) < 1e-9,
              "KAIN consumes cached residual work without another BSH application");
    }
}
}

int main(int argc, char** argv) {
    initialize(argc, argv);
    {
        World world(SafeMPI::COMM_WORLD);
        startup(world, argc, argv, true);
        const auto original = std::filesystem::current_path();
        std::string directory;
        if (world.rank() == 0) {
            char pattern[] = "/tmp/madness-scf-evaluation-XXXXXX";
            const char* made = mkdtemp(pattern);
            if (made) directory = made;
        }
        world.gop.broadcast_serializable(directory, 0);
        MADNESS_CHECK_THROW(!directory.empty(), "cannot create SCF test directory");
        std::filesystem::current_path(directory);
        try {
            if (argc > 1 && std::string(argv[1]) == "--localized-spectrum") {
                current_state_case(world, "water_localized", "hf", 8, "new");
            } else {
            current_state_case(world, "he_hf", "hf", 2, "canon");
            current_state_case(world, "he_lda", "lda", 2, "new");
            current_state_case(world, "li_hf", "hf", 3, "new");
            current_state_case(world, "water_hf", "hf", 8, "canon");
            current_state_case(world, "he_immediate", "hf", 2, "canon", true);
            current_state_case(world, "water_localized", "hf", 8, "new");
            }
        }
        catch (const std::exception& error) { ++failures; std::cerr << error.what() << '\n'; }
        world.gop.fence();
        std::filesystem::current_path(original);
        if (world.rank() == 0) std::filesystem::remove_all(directory);
        world.gop.max(failures);
    }
    finalize();
    return failures ? 1 : 0;
}
