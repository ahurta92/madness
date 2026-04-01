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

// Li atom (UHF, nopen=1)
// Literature HF alpha = 170.02019 au (spherical: xx=yy=zz)
// Source: Phys. Rev. A 91, 022501 (2015)
static const MoleculeReference Li_REF = {
    "Li (UHF)", -7.432672458,
    {170.02019, 170.02019, 170.02019},
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
        // Try to use fock.json from the archive directory for correct Fock matrix
        double thresh = FunctionDefaults<3>::get_thresh();
        auto coulop = poperatorT(CoulombOperatorPtr(world, gs.params().lo(), 0.001 * thresh));
        auto archive_dir = std::filesystem::path(archive_path).parent_path();
        std::string fock_json;
        for (const auto& name : {"moldft.fock.json", "mad.fock.json"}) {
            auto candidate = archive_dir / name;
            if (std::filesystem::exists(candidate)) {
                fock_json = candidate.string();
                if (world.rank() == 0) print("Using Fock from:", fock_json);
                break;
            }
        }
        gs.prepare(world, 0.001 * thresh, coulop, fock_json);

        // Determine which reference to use based on orbital count + spin
        const MoleculeReference* ref = nullptr;
        bool symmetry_only = false;  // no reference values, just check symmetry
        if (gs.num_alpha() == 1 && gs.is_spin_restricted()) {
            ref = &H2_REF;
        } else if (gs.num_alpha() == 5 && gs.is_spin_restricted()) {
            ref = &H2O_REF;
        } else if (gs.num_alpha() == 3 && gs.is_spin_restricted()) {
            ref = &BeH2_REF;
        } else if (!gs.is_spin_restricted() && gs.molecule().natom() == 1) {
            ref = &Li_REF;  // Li atom (or any UHF atom)
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
                /*omega=*/0.0, /*maxiter=*/25,
                /*dconv=*/thresh,
                /*maxrotn=*/10.0, /*maxsub=*/10,
                PrintLevel::Verbose);

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

        // ============================================================
        // Dynamic (frequency-dependent) test
        // For H2: test at omega=0.029 and compare alpha_zz against
        // legacy dynamic data (should be slightly larger than static)
        // ============================================================
        bool run_dynamic = parser.key_exists("dynamic") ||
                           (!symmetry_only && gs.num_orbitals() == 1);

        if (run_dynamic && !symmetry_only) {
            double omega = 0.029;  // close to legacy 0.029092
            if (parser.key_exists("omega")) {
                omega = std::stod(parser.value("omega"));
            }

            if (world.rank() == 0) {
                print("\n--- Dynamic Polarizability (omega=", omega, ") ---");
            }

            // Solve dipole_z at omega (Full response: x + y channels)
            RealResponseState pert_dyn;
            pert_dyn.x_alpha = dipole_perturbation(world, gs, 2); // z
            pert_dyn.y_alpha = pert_dyn.x_alpha; // y perturbation = x perturbation for dipole

            auto dyn_result = fd_solve(
                world, ResponseType::Full, pert_dyn, gs,
                omega, /*maxiter=*/15, /*dconv=*/thresh,
                /*maxrotn=*/0.5, /*maxsub=*/10,
                PrintLevel::Verbose);

            if (world.rank() == 0) {
                print("  alpha_zz(omega=", omega, ") =", dyn_result.alpha,
                      " converged=", dyn_result.converged);
            }

            // Dynamic alpha should be larger than static (dispersion)
            results.push_back(check_value("dynamic_gt_static",
                1.0, (dyn_result.alpha > alpha_tensor[2]) ? 1.0 : 0.0, 0.5));

            // Should be positive
            results.push_back(check_value("dynamic_positive",
                1.0, (dyn_result.alpha > 0.0) ? 1.0 : 0.0, 0.5));

            // For H2 at omega~0.029: legacy gives alpha_zz ~ 6.58 (vs static 6.46)
            if (ref == &H2_REF) {
                results.push_back(check_value("dynamic_zz_h2",
                    6.58, dyn_result.alpha, 0.2));  // loose tolerance
            }
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
    finalize();
}
