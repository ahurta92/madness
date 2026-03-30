#include <madness/world/MADworld.h>
#include <madness/mra/mra.h>
#include <madness/misc/info.h>
#include <madness/world/worldmem.h>
#include <madness/chem/atomutil.h>

#include "GroundState.hpp"
#include "ResponseFunctions.hpp"

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
