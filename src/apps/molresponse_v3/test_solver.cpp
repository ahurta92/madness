#include "FDSolver.hpp"
#include "GroundState.hpp"
#include "Perturbations.hpp"
#include "ResponseFunctions.hpp"
#include "ResponseKernel.hpp"

#include <madness/mra/mra.h>
#include <madness/misc/info.h>
#include <madness/world/MADworld.h>
#include <madness/world/worldmem.h>

#include <cmath>
#include <string>

using namespace madness;
using namespace molresponse_v3;

/// Reference data from legacy molresponse (k=8, HF, protocol [1e-4,1e-6,1e-6])
/// Source: /gpfs/projects/rjh/adrian/post_watoc/september/hf/
struct ReferenceAlpha {
    double xx, yy, zz;
    double isotropic() const { return (xx + yy + zz) / 3.0; }
};

struct ReferenceExcited {
    int num_roots;
    std::vector<double> omega_au;
};

struct MoleculeReference {
    std::string name;
    double scf_energy;
    ReferenceAlpha alpha;
    ReferenceExcited excited;
};

// H2 legacy reference
static const MoleculeReference H2_REF = {
    "H2", -1.1336222987,
    {4.614102, 4.614102, 6.456580},
    {4, {0.46535391, 0.47696774, 0.48065960, 0.48069615}}
};

// H2O legacy reference
static const MoleculeReference H2O_REF = {
    "H2O", -76.0674405361,
    {7.903662, 9.190789, 8.532181},
    {4, {0.31778600, 0.39966806, 0.49630384, 0.65994464}}
};

// BeH2 legacy reference
static const MoleculeReference BeH2_REF = {
    "BeH2", -15.7732181961,
    {18.809856, 18.809856, 19.412931},
    {4, {0.25948557, 0.25948661, 0.33376130, 0.33377215}}
};

// Li atom (UHF, nopen=1) — no legacy MRA reference.
// Literature HF alpha ~ 164 au (spherical: xx=yy=zz).
// We use 0.0 as expected to signal "no reference, symmetry check only".
static const MoleculeReference Li_REF = {
    "Li (UHF)", 0.0,
    {0.0, 0.0, 0.0},  // no reference — will run symmetry-only checks
    {0, {}}
};

struct TestResult {
    std::string test_name;
    bool passed;
    double expected;
    double actual;
    double tolerance;
    double error;
};

TestResult check_value(const std::string& name, double expected, double actual,
                        double tolerance) {
    double error = std::abs(actual - expected);
    bool passed = error < tolerance;
    return {name, passed, expected, actual, tolerance, error};
}

void print_result(const TestResult& r) {
    const char* status = r.passed ? "PASS" : "FAIL";
    print("  [", status, "] ", r.test_name,
          ": expected=", r.expected, " actual=", r.actual,
          " error=", r.error, " tol=", r.tolerance);
}

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);

    try {
        startup(world, argc, argv, true);

        if (world.rank() == 0) {
            print("\n============================================");
            print("  molresponse_v3 — Solver Validation Tests");
            print("============================================\n");
        }

        commandlineparser parser(argc, argv);

        std::string archive_path;
        if (parser.key_exists("archive")) {
            archive_path = parser.value("archive");
        } else {
            if (world.rank() == 0) {
                print("Usage: test_solver --archive=<path> [--molecule=<json>]");
                print("  Tests static polarizability against legacy reference.");
            }
            finalize();
            return 1;
        }

        // Tolerance for comparison (medium tier: dconv=1e-4, expect ~1e-3 agreement)
        double tol = 0.01;  // 0.01 au tolerance for alpha components
        if (parser.key_exists("tol")) {
            tol = std::stod(parser.value("tol"));
        }

        // Load molecule from calc_info.json
        std::string molecule_json_path;
        if (parser.key_exists("molecule")) {
            molecule_json_path = parser.value("molecule");
        } else {
            auto archive_dir = std::filesystem::path(archive_path).parent_path();
            for (const auto& name : {"moldft.calc_info.json", "mad.calc_info.json"}) {
                auto candidate = archive_dir / name;
                if (std::filesystem::exists(candidate)) {
                    molecule_json_path = candidate.string();
                    break;
                }
            }
        }

        Molecule molecule;
        if (!molecule_json_path.empty() && std::filesystem::exists(molecule_json_path)) {
            std::ifstream ifs(molecule_json_path);
            nlohmann::json j;
            ifs >> j;
            nlohmann::json mol_json;
            if (j.contains("tasks") && j["tasks"].is_array() && !j["tasks"].empty()) {
                mol_json = j["tasks"][0]["molecule"];
            } else if (j.contains("molecule")) {
                mol_json = j["molecule"];
            }
            // Legacy flat format
            if (mol_json.is_null() && j.contains("parameters")) {
                // Molecule loaded from archive directly
            }
            if (!mol_json.is_null()) {
                molecule.from_json(mol_json);
            }
        }

        // Load ground state
        auto gs = GroundState::from_archive(world, archive_path, molecule);

        if (world.rank() == 0) {
            gs.print_info();
        }

        // Prepare ground state
        double thresh = FunctionDefaults<3>::get_thresh();
        auto coulop = poperatorT(CoulombOperatorPtr(world, gs.params().lo(), 0.001 * thresh));
        gs.prepare(world, 0.001 * thresh, coulop);

        // Determine which reference to use based on orbital count + spin
        const MoleculeReference* ref = nullptr;
        bool symmetry_only = false;  // no reference values, just check symmetry
        if (gs.num_alpha() == 1 && gs.is_spin_restricted()) {
            ref = &H2_REF;
        } else if (gs.num_alpha() == 5 && gs.is_spin_restricted()) {
            ref = &H2O_REF;
        } else if (gs.num_alpha() == 3 && gs.is_spin_restricted()) {
            ref = &BeH2_REF;
        } else if (!gs.is_spin_restricted()) {
            ref = &Li_REF;  // generic unrestricted — symmetry checks only
            symmetry_only = true;
        }

        if (!ref) {
            if (world.rank() == 0) {
                print("No reference data for molecule with ", gs.num_alpha(),
                      " alpha orbitals, restricted=", gs.is_spin_restricted());
                print("Supported: H2 (1), BeH2 (3), H2O (5), any unrestricted (symmetry only)");
            }
            finalize();
            return 1;
        }

        if (world.rank() == 0) {
            if (symmetry_only) {
                print("\nReference: ", ref->name, " (NO reference values — symmetry checks only)");
                print("  spin_restricted=", gs.is_spin_restricted());
                print("  num_alpha=", gs.num_alpha(), " num_beta=", gs.num_beta());
            } else {
                print("\nReference: ", ref->name, " (legacy k=8)");
                print("  alpha_xx=", ref->alpha.xx, " alpha_yy=", ref->alpha.yy,
                      " alpha_zz=", ref->alpha.zz);
                print("  isotropic=", ref->alpha.isotropic());
            }
            print("  tolerance=", tol, " au\n");
        }

        // Solve static polarizability in all 3 directions
        std::vector<TestResult> results;
        double alpha_tensor[3] = {};
        const char* dir_names[] = {"x", "y", "z"};
        double ref_vals[] = {ref->alpha.xx, ref->alpha.yy, ref->alpha.zz};

        for (int d = 0; d < 3; d++) {
            if (world.rank() == 0) {
                print("--- Solving dipole_", dir_names[d], " ---");
            }

            RealResponseState pert;
            pert.x_alpha = dipole_perturbation(world, gs, d);
            if (!gs.is_spin_restricted()) {
                pert.x_beta = dipole_perturbation_beta(world, gs, d);
            }

            auto solve_result = fd_solve(
                world, ResponseType::Static, pert, gs,
                /*omega=*/0.0, /*maxiter=*/15,
                /*dconv=*/thresh);

            alpha_tensor[d] = solve_result.alpha;

            if (!symmetry_only) {
                auto tr = check_value(
                    std::string("alpha_") + dir_names[d] + dir_names[d],
                    ref_vals[d], solve_result.alpha, tol);
                results.push_back(tr);
            }

            if (world.rank() == 0) {
                print("");
            }
        }

        // Isotropic
        double iso = (alpha_tensor[0] + alpha_tensor[1] + alpha_tensor[2]) / 3.0;
        if (!symmetry_only) {
            results.push_back(check_value("isotropic", ref->alpha.isotropic(), iso, tol));
        }

        // Symmetry check: xx == yy for linear/spherical molecules
        results.push_back(check_value("symmetry_xx_yy",
            0.0, std::abs(alpha_tensor[0] - alpha_tensor[1]), 1e-6));

        // For atoms (1 atom, spherical): xx == yy == zz
        if (gs.molecule().natom() == 1) {
            results.push_back(check_value("spherical_xx_zz",
                0.0, std::abs(alpha_tensor[0] - alpha_tensor[2]), tol * 0.1));
        }

        // Sanity: alpha should be positive and reasonable
        results.push_back(check_value("alpha_positive",
            1.0, (iso > 0.0) ? 1.0 : 0.0, 0.5));

        if (symmetry_only && world.rank() == 0) {
            print("  Computed alpha_xx=", alpha_tensor[0],
                  " alpha_yy=", alpha_tensor[1],
                  " alpha_zz=", alpha_tensor[2]);
            print("  isotropic=", iso);
        }

        // Print summary
        if (world.rank() == 0) {
            print("\n============================================");
            print("  Test Results vs Legacy Reference");
            print("============================================");
            int passed = 0, failed = 0;
            for (const auto& r : results) {
                print_result(r);
                if (r.passed) passed++;
                else failed++;
            }
            print("\n  ", passed, " passed, ", failed, " failed");
            print("============================================\n");
        }

        world.gop.fence();

        // Return nonzero if any test failed
        bool all_passed = true;
        for (const auto& r : results) {
            if (!r.passed) all_passed = false;
        }

        world.gop.fence();
        print_stats(world);
        finalize();
        return all_passed ? 0 : 1;

    } catch (const std::exception& e) {
        if (world.rank() == 0) {
            print("ERROR:", e.what());
        }
        finalize();
        return 2;
    }
}
