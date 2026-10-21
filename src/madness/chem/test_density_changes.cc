#include <madness/chem/DensityChanges.h>
#include <madness/chem/SCF.h>
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace madness;

static double gaussian(const coord_3d& r) {
    return std::exp(-(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]));
}
static double odd_gaussian(const coord_3d& r) { return r[0]*gaussian(r); }

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    {
        startup(world, argc, argv, true);
        FunctionDefaults<3>::set_cubic_cell(-12.0, 12.0);
        FunctionDefaults<3>::set_k(8);
        FunctionDefaults<3>::set_thresh(1e-8);
        // Explicit Function handles avoid retaining a factory's lazy state.
        real_function_3d even = real_factory_3d(world).f(gaussian);
        real_function_3d odd = real_factory_3d(world).f(odd_gaussian);
        even.scale(1.0/even.norm2()); odd.scale(1.0/odd.norm2());
        auto scaled = [&](double factor) { auto x = copy(even); x.scale(factor); return x; };
        auto before = scaled(2.0), higher = scaled(3.0), lower = scaled(1.0), zero = scaled(0.0);
        auto check = [&](const real_function_3d& a, const real_function_3d& b,
                         const real_function_3d& pa, const real_function_3d& pb,
                         double charge, double spin) {
            vector_real_function_3d inputs{a, b, pa, pb};
            auto snapshots = copy(world, inputs);
            std::vector<TreeState> states;
            for (const auto& input : inputs) states.push_back(get_tree_state(input));
            const auto result = density_change_norms(world, a, b, pa, pb);
            for (std::size_t i = 0; i < inputs.size(); ++i)
                MADNESS_CHECK(get_tree_state(inputs[i]) == states[i]);
            // Compare independent copies so the check itself preserves inputs.
            auto after = copy(world, inputs);
            for (double norm : norm2s(world, sub(world, after, snapshots)))
                MADNESS_CHECK(norm < 1e-12);
            MADNESS_CHECK(std::abs(result.charge-charge) < 1e-10);
            MADNESS_CHECK(std::abs(result.spin-spin) < 1e-10);
        };
        // Each spin changes by one normalized field: charge=2, spin=0.
        check(higher, higher, before, before, 2.0, 0.0);
        // Opposing changes cancel in charge, while spin changes by two.
        check(higher, lower, before, before, 0.0, 2.0);
        check(higher, before, before, before, 1.0, 1.0);
        check(higher, zero, before, zero, 1.0, 1.0);
        // Orthogonal even/odd changes distinguish field norms from sums or
        // differences of scalar spin norms (both would give wrong answers).
        check(even, odd, zero, zero, std::sqrt(2.0), std::sqrt(2.0));
        MADNESS_CHECK(std::abs(before.norm2()-2.0) < 1e-12);
        MADNESS_CHECK(std::abs(higher.norm2()-3.0) < 1e-12);
        MADNESS_CHECK(std::abs(lower.norm2()-1.0) < 1e-12);
        MADNESS_CHECK(zero.norm2() < 1e-12);
        // Exercise the real unrestricted observer path, including zero beta
        // population. This short solve checks reporting, not energy accuracy.
        CalculationParameters parameters;
        parameters.set_user_defined_value<std::string>("xc", "hf");
        parameters.set_user_defined_value<std::string>("prefix", "density-observer-h");
        parameters.set_user_defined_value("l", 12.0);
        parameters.set_user_defined_value("nopen", 1);
        parameters.set_user_defined_value("protocol", std::vector<double>{1e-4});
        parameters.set_user_defined_value("dconv", 1e-4);
        parameters.set_user_defined_value("econv", 1e-6);
        parameters.set_user_defined_value("maxiter", 3);
        parameters.set_user_defined_value("save", false);
        parameters.set_user_defined_value("print_level", 0);
        Molecule molecule;
        molecule.add_atom(0.0, 0.0, 0.0, 1.0, 1);
        SCF scf(world, parameters, molecule);
        MADNESS_CHECK(!scf.param.spin_restricted() && scf.param.nbeta() == 0);
        int measured_samples = 0;
        double largest_alpha_change = 0.0;
        scf.checkpoint_observer = [&](World&, const SCF&, const SCF::IterationCheckpoint& step) {
            if (!step.has_previous) return;
            ++measured_samples;
            largest_alpha_change = std::max(largest_alpha_change, step.density_alpha);
            if (world.rank() == 0)
                std::cout << std::setprecision(17) << "SCF_DENSITY_SAMPLE " << step.iteration
                          << " " << step.density_alpha << " " << step.density_beta
                          << " " << step.density_charge << " " << step.density_spin << "\n";
            MADNESS_CHECK(std::abs(step.density_charge-step.density_alpha) < 1e-10);
            MADNESS_CHECK(std::abs(step.density_spin-step.density_alpha) < 1e-10);
        };
        MolecularEnergy energy(world, scf);
        MADNESS_CHECK(std::isfinite(energy.value(molecule.get_all_coords().flat())));
        MADNESS_CHECK(measured_samples > 0 && largest_alpha_change > 1e-12);
        if (world.rank() == 0) std::cout << "charge/spin density changes passed\n";
    }
    finalize();
}
