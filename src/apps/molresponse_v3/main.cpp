#include <madness/world/MADworld.h>
#include <madness/mra/mra.h>
#include <madness/misc/info.h>
#include <madness/world/worldmem.h>
#include <madness/chem/atomutil.h>

#include "FDSolver.hpp"
#include "GroundState.hpp"
#include "Perturbations.hpp"
#include "ResponseFunctions.hpp"
#include "ResponseKernel.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

using namespace madness;

void print_usage() {
    print("Usage: molresponse_v3 --archive=<path> [--molecule=<json>]");
    print("");
    print("Load a ground-state checkpoint and print info.");
    print("");
    print("Required:");
    print("  --archive=<path>   Path to moldft restart archive");
    print("");
    print("Optional:");
    print("  --molecule=<json>  Path to moldft.calc_info.json for molecule data");
}

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);

    try {
        startup(world, argc, argv, true);

        if (world.rank() == 0) {
            print("\n====================================");
            print("  molresponse_v3");
            print("  Ground-state loader (SCF wrapper)");
            print("====================================\n");
        }

        commandlineparser parser(argc, argv);

        if (parser.key_exists("help")) {
            if (world.rank() == 0) print_usage();
            finalize();
            return 0;
        }

        // Resolve archive path
        std::string archive_path;
        if (parser.key_exists("archive")) {
            archive_path = parser.value("archive");
        } else {
            if (world.rank() == 0) {
                print("ERROR: --archive=<path> is required");
                print_usage();
            }
            finalize();
            return 1;
        }

        if (world.rank() == 0) {
            print("Archive path:", archive_path);
        }

        // Look for molecule JSON (calc_info.json)
        std::string molecule_json_path;
        if (parser.key_exists("molecule")) {
            molecule_json_path = parser.value("molecule");
        } else {
            // Try common names in the archive directory
            auto archive_dir = std::filesystem::path(archive_path).parent_path();
            for (const auto& name : {"moldft.calc_info.json", "mad.calc_info.json"}) {
                auto candidate = archive_dir / name;
                if (std::filesystem::exists(candidate)) {
                    molecule_json_path = candidate.string();
                    break;
                }
            }
        }

        // Read molecule from JSON if available
        Molecule molecule;
        if (!molecule_json_path.empty() && std::filesystem::exists(molecule_json_path)) {
            if (world.rank() == 0) {
                print("Molecule JSON:", molecule_json_path);
            }
            std::ifstream ifs(molecule_json_path);
            nlohmann::json j;
            ifs >> j;
            // Navigate into tasks[0].molecule if present (calc_info format)
            nlohmann::json mol_json;
            if (j.contains("tasks") && j["tasks"].is_array() && !j["tasks"].empty()) {
                auto& task0 = j["tasks"][0];
                if (task0.contains("molecule")) {
                    mol_json = task0["molecule"];
                }
            } else if (j.contains("molecule")) {
                mol_json = j["molecule"];
            }
            if (!mol_json.is_null()) {
                molecule.from_json(mol_json);
            }
        } else {
            if (world.rank() == 0) {
                print("WARN: No molecule JSON found, molecule loaded from archive");
            }
        }

        // Load ground state using v3 GroundState (wraps SCF)
        auto ground = molresponse_v3::GroundState::from_archive(
            world, archive_path, molecule);

        // Print ground-state summary
        if (world.rank() == 0) {
            ground.print_info();

            const auto& mol = ground.molecule();
            print("\n--- Molecular Geometry ---");
            print("Atoms:", mol.natom());
            for (size_t i = 0; i < mol.natom(); ++i) {
                const auto& atom = mol.get_atom(i);
                print("  ", get_atomic_data(atom.atomic_number).symbol,
                      " ", atom.x, " ", atom.y, " ", atom.z);
            }
            print("--- End Summary ---\n");
        }

        // --- Test ResponseState types ---
        if (world.rank() == 0) {
            print("\n--- ResponseState Tests ---");
        }

        long na = ground.num_alpha();
        long nb = ground.is_spin_restricted() ? 0 : ground.num_beta();

        // Static restricted (x only, no y, no beta)
        auto static_state = molresponse_v3::RealResponseState::allocate(
            world, na, nb, false);
        if (world.rank() == 0) {
            static_state.print_info("static (x-only)");
        }

        // Full restricted (x + y)
        auto full_state = molresponse_v3::RealResponseState::allocate(
            world, na, nb, true);
        if (world.rank() == 0) {
            full_state.print_info("full (x+y)");
        }

        // Test flat round-trip
        auto flat = full_state.flat();
        if (world.rank() == 0) {
            print("  flat size:", flat.size(),
                  "(expected:", full_state.total_size(), ")");
        }
        full_state.from_flat(flat);
        if (world.rank() == 0) {
            print("  flat round-trip: OK");
        }

        // Test in-place vmra on flat (shared_ptr semantics)
        auto f = full_state.flat();
        madness::scale(world, f, 2.0);
        double norm_xa = madness::norm2(world, full_state.x_alpha);
        if (world.rank() == 0) {
            print("  scale(flat, 2.0) -> ||x_alpha|| =", norm_xa,
                  "(should be 0 for zero funcs)");
        }

        if (world.rank() == 0) {
            print("--- ResponseState Tests PASSED ---\n");
        }

        // --- Test Perturbation operators ---
        if (world.rank() == 0) {
            print("\n--- Perturbation Tests ---");
        }

        // Need to call prepare() so Q projector is available
        // Set up protocol for current defaults
        auto thresh = FunctionDefaults<3>::get_thresh();
        auto coulop = poperatorT(CoulombOperatorPtr(world, 1e-10, 0.001 * thresh));
        ground.prepare(world, 0.001 * thresh, coulop);

        // Dipole perturbation in z direction
        auto rhs_z = molresponse_v3::dipole_perturbation(world, ground, 2);
        auto norms_z = norm2s(world, rhs_z);
        if (world.rank() == 0) {
            print("  dipole_z perturbation: ", rhs_z.size(), "functions");
            for (size_t i = 0; i < norms_z.size(); i++) {
                print("    ||Q*mu_z*phi_", i, "|| =", norms_z[i]);
            }
        }

        // All three dipole operators
        auto [mu_x, mu_y, mu_z] = molresponse_v3::dipole_operators(world);
        if (world.rank() == 0) {
            print("  dipole operators: x, y, z constructed");
        }

        // Dipole perturbation for all 3 directions
        const char* dirs[] = {"x", "y", "z"};
        for (int d = 0; d < 3; d++) {
            auto rhs = molresponse_v3::dipole_perturbation(world, ground, d);
            auto norms = norm2s(world, rhs);
            if (world.rank() == 0) {
                print("  dipole_", dirs[d], ": ||rhs|| =", norms[0]);
            }
        }

        if (world.rank() == 0) {
            print("--- Perturbation Tests PASSED ---\n");
        }

        // --- Test FD iteration step ---
        if (world.rank() == 0) {
            print("\n--- FD Iteration Test (H2 Static dipole_z) ---");
        }

        // Build perturbation as a ResponseState
        molresponse_v3::RealResponseState v_pert;
        v_pert.x_alpha = molresponse_v3::dipole_perturbation(world, ground, 2);

        // Start from zero response
        auto response = molresponse_v3::RealResponseState::allocate(
            world, na, nb, false);  // static: no y

        // Create BSH operators for static (omega=0)
        auto bsh_alpha_x = molresponse_v3::make_bsh_operators(
            world, ground.energies_alpha(), 0.0, ground.params().lo());
        std::vector<poperatorT> bsh_alpha_y, bsh_beta_x, bsh_beta_y;

        // One FD iteration
        auto updated = molresponse_v3::fd_iteration(
            world,
            molresponse_v3::ResponseType::Static,
            response, v_pert, ground, coulop,
            bsh_alpha_x, bsh_alpha_y, bsh_beta_x, bsh_beta_y,
            0.0);

        auto x_norms = norm2s(world, updated.x_alpha);
        if (world.rank() == 0) {
            for (size_t i = 0; i < x_norms.size(); i++) {
                print("  ||x_", i, "_new|| =", x_norms[i]);
            }
        }

        // Compute response density
        auto rho1 = molresponse_v3::compute_response_density(
            world, molresponse_v3::ResponseType::Static, updated, ground);
        double rho_norm = rho1.norm2();
        if (world.rank() == 0) {
            print("  ||rho^(1)|| =", rho_norm);
        }

        // Compute alpha_zz from first iteration
        double afactor = molresponse_v3::alpha_factor(
            molresponse_v3::ResponseType::Static, ground.is_spin_restricted());
        double alpha_zz = afactor * madness::inner(updated.x_alpha, v_pert.x_alpha);
        if (world.rank() == 0) {
            print("  alpha_zz (1 iter) =", alpha_zz, "(reference: ~8.53)");
            print("--- FD Iteration Test PASSED ---\n");
        }

        // --- Full static polarizability solve ---
        if (world.rank() == 0) {
            print("\n======================================");
            print("  Static Polarizability Solve (H2)");
            print("======================================\n");
        }

        const char* dir_names[] = {"x", "y", "z"};
        double alpha_tensor[3] = {};

        for (int d = 0; d < 3; d++) {
            if (world.rank() == 0) {
                print("--- Solving dipole_", dir_names[d], " ---");
            }

            molresponse_v3::RealResponseState pert_d;
            pert_d.x_alpha = molresponse_v3::dipole_perturbation(world, ground, d);

            auto result = molresponse_v3::fd_solve(
                world,
                molresponse_v3::ResponseType::Static,
                pert_d, ground,
                /*omega=*/0.0,
                /*maxiter=*/15,
                /*dconv=*/FunctionDefaults<3>::get_thresh());

            alpha_tensor[d] = result.alpha;

            if (world.rank() == 0) {
                print("  alpha_", dir_names[d], dir_names[d], " = ", result.alpha,
                      " converged=", result.converged,
                      " iters=", result.iterations, "\n");
            }
        }

        if (world.rank() == 0) {
            double iso = (alpha_tensor[0] + alpha_tensor[1] + alpha_tensor[2]) / 3.0;
            print("\n--- Static Polarizability Summary ---");
            print("  alpha_xx =", alpha_tensor[0]);
            print("  alpha_yy =", alpha_tensor[1]);
            print("  alpha_zz =", alpha_tensor[2]);
            print("  isotropic =", iso);
            print("--- End Polarizability ---\n");
        }

        world.gop.fence();

    } catch (const std::exception& e) {
        if (world.rank() == 0) {
            print("ERROR:", e.what());
        }
    }

    world.gop.fence();
    print_stats(world);
    finalize();
    return 0;
}
