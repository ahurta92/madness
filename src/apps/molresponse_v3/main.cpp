#include <madness/world/MADworld.h>
#include <madness/mra/mra.h>
#include <madness/misc/info.h>
#include <madness/world/worldmem.h>
#include <madness/chem/CalculationParameters.h>
#include <madness/chem/ResponseParameters.hpp>

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
using namespace molresponse_v3;

void print_usage() {
    print("Usage: molresponse_v3 [input_file] [options]");
    print("");
    print("Reads a standard MADNESS input file with dft/response/molecule sections.");
    print("Loads the ground-state checkpoint and solves frequency-dependent response.");
    print("");
    print("The input file defaults to 'input' if not specified.");
    print("");
    print("Input file format:");
    print("  dft");
    print("    xc \"hf\"");
    print("    protocol [1e-4,1e-6]");
    print("    l 200.0");
    print("  end");
    print("  response");
    print("    archive ../moldft.restartdata");
    print("    dipole true");
    print("    dipole.directions xyz");
    print("    dipole.frequencies [0.0,0.02,0.04]");
    print("    dconv 1e-4");
    print("    maxiter 25");
    print("    kain true");
    print("  end");
    print("  molecule");
    print("    eprec 1e-4");
    print("    units bohr");
    print("    Li 0.0 0.0 0.0");
    print("  end");
    print("");
    print("Options:");
    print("  --help           Print this help message");
    print("  --print_level=N  Override print level (0-3)");
}

/// Parse dipole.directions string ("xyz", "z", "xz", etc.) into axis indices.
std::vector<int> parse_directions(const std::string& dirs) {
    std::vector<int> axes;
    for (char c : dirs) {
        if (c == 'x' || c == 'X') axes.push_back(0);
        else if (c == 'y' || c == 'Y') axes.push_back(1);
        else if (c == 'z' || c == 'Z') axes.push_back(2);
    }
    if (axes.empty()) axes = {0, 1, 2};  // default: all three
    return axes;
}

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);

    try {
        startup(world, argc, argv, true);
        commandlineparser parser(argc, argv);

        if (parser.key_exists("help")) {
            if (world.rank() == 0) print_usage();
            finalize();
            return 0;
        }

        if (world.rank() == 0) {
            print("\n====================================");
            print("  molresponse_v3");
            print("  Frequency-dependent response solver");
            print("====================================\n");
            print(info::print_revision_information());
            print("input file:", parser.value("input"));
        }

        // ----------------------------------------------------------------
        // 1. Read parameters from input file
        // ----------------------------------------------------------------
        CalculationParameters dft_params(world, parser);
        ResponseParameters response_params(world, parser);
        Molecule molecule(world, parser);

        // Extract what we need
        auto protocol = dft_params.protocol();
        double L = dft_params.L();
        std::string archive_path = response_params.archive();
        auto frequencies = response_params.dipole_frequencies();
        auto directions = parse_directions(response_params.dipole_directions());
        size_t maxiter = response_params.maxiter();
        double dconv = response_params.dconv();
        double maxrotn = response_params.maxrotn();
        size_t maxsub = response_params.maxsub();
        bool use_kain = response_params.kain();

        // Print level: from input or command-line override
        int print_level_int = response_params.print_level();
        if (parser.key_exists("print_level")) {
            print_level_int = std::stoi(parser.value("print_level"));
        }
        auto plevel = static_cast<PrintLevel>(
            std::min(print_level_int, static_cast<int>(PrintLevel::Debug)));

        // Resolve archive path relative to input file directory
        auto input_dir = std::filesystem::path(parser.value("input")).parent_path();
        if (!input_dir.empty() && !std::filesystem::path(archive_path).is_absolute()) {
            auto candidate = input_dir / archive_path;
            if (std::filesystem::exists(candidate.string() + ".00000")) {
                archive_path = candidate.string();
            }
        }

        // Look for Fock JSON near the archive
        std::string fock_json_file;
        auto archive_dir = std::filesystem::path(archive_path).parent_path();
        for (const auto& name : {"moldft.fock.json", "mad.fock.json"}) {
            auto candidate = archive_dir / name;
            if (std::filesystem::exists(candidate)) {
                fock_json_file = candidate.string();
                break;
            }
        }

        if (world.rank() == 0) {
            print("\nPARAMETERS:");
            print("  archive    =", archive_path);
            print("  fock_json  =", fock_json_file.empty() ? "(none)" : fock_json_file);
            print("  L          =", L);
            print("  protocol   =", protocol);
            print("  dipole     =", response_params.dipole());
            print("  directions =", response_params.dipole_directions());
            print("  frequencies=", frequencies);
            print("  maxiter    =", maxiter);
            print("  dconv      =", dconv);
            print("  maxrotn    =", maxrotn);
            print("  kain       =", use_kain);
            print("  print_level=", print_level_int);
            print("");
        }

        // ----------------------------------------------------------------
        // 2. Load ground state
        // ----------------------------------------------------------------
        // Set L before loading (from_archive checks L)
        FunctionDefaults<3>::set_cubic_cell(-L, L);

        auto gs = GroundState::from_archive(world, archive_path, molecule);

        if (world.rank() == 0) {
            print("\n--- Molecular Geometry ---");
            const auto& mol = gs.molecule();
            for (size_t i = 0; i < mol.natom(); ++i) {
                const auto& atom = mol.get_atom(i);
                print("  ", atomic_number_to_symbol(atom.atomic_number),
                      " ", atom.x, " ", atom.y, " ", atom.z);
            }
            print("");
        }

        // ----------------------------------------------------------------
        // 3. Protocol loop
        // ----------------------------------------------------------------
        // Storage for response functions carried between protocols.
        // Key: (direction, frequency_index) -> converged ResponseState
        struct SolveKey {
            int direction;
            size_t freq_idx;
            bool operator<(const SolveKey& o) const {
                return std::tie(direction, freq_idx) < std::tie(o.direction, o.freq_idx);
            }
        };
        std::map<SolveKey, RealResponseState> prev_solutions;

        // Final results storage: alpha_tensor[direction][freq_idx]
        const char* dir_names[] = {"x", "y", "z"};

        for (size_t proto_idx = 0; proto_idx < protocol.size(); proto_idx++) {
            double thresh = protocol[proto_idx];

            if (world.rank() == 0) {
                print("\n========================================================");
                print("  PROTOCOL ", proto_idx, ": thresh=", thresh);
                print("========================================================\n");
            }

            // Set numerical parameters for this protocol step
            int k;
            if      (thresh >= 0.9e-2) k = 4;
            else if (thresh >= 0.9e-4) k = 6;
            else if (thresh >= 0.9e-6) k = 8;
            else if (thresh >= 0.9e-8) k = 10;
            else                        k = 12;

            FunctionDefaults<3>::set_k(k);
            FunctionDefaults<3>::set_thresh(thresh);
            FunctionDefaults<3>::set_refine(true);
            FunctionDefaults<3>::set_initial_level(2);
            FunctionDefaults<3>::set_autorefine(false);
            FunctionDefaults<3>::set_apply_randomize(false);
            FunctionDefaults<3>::set_project_randomize(false);
            FunctionDefaults<3>::set_cubic_cell(-L, L);
            GaussianConvolution1DCache<double>::map.clear();

            double vtol = thresh * 0.1;
            auto coulop = poperatorT(CoulombOperatorPtr(world, gs.params().lo(), 0.001 * thresh));

            // Prepare ground state for this protocol
            gs.prepare(world, vtol, coulop, fock_json_file);

            // Effective dconv for this protocol (don't demand more than thresh allows)
            double eff_dconv = std::max(dconv, thresh);

            // ----------------------------------------------------------------
            // 4. Solve for each direction and frequency
            // ----------------------------------------------------------------
            for (int d : directions) {
                for (size_t fi = 0; fi < frequencies.size(); fi++) {
                    double omega = frequencies[fi];
                    bool is_static = (omega == 0.0);
                    auto rtype = is_static ? ResponseType::Static : ResponseType::Full;

                    if (world.rank() == 0) {
                        print("--- dipole_", dir_names[d],
                              " omega=", omega,
                              " type=", is_static ? "static" : "full",
                              " ---");
                    }

                    // Build perturbation
                    RealResponseState pert;
                    pert.x_alpha = dipole_perturbation(world, gs, d);
                    if (!is_static) {
                        pert.y_alpha = pert.x_alpha;
                    }
                    if (!gs.is_spin_restricted()) {
                        pert.x_beta = dipole_perturbation_beta(world, gs, d);
                        if (!is_static) {
                            pert.y_beta = pert.x_beta;
                        }
                    }

                    // Initial guess from previous protocol (if available)
                    SolveKey key{d, fi};
                    const RealResponseState* guess = nullptr;
                    if (prev_solutions.count(key)) {
                        guess = &prev_solutions[key];
                    }

                    // Solve
                    auto result = fd_solve(
                        world, rtype, pert, gs, omega,
                        static_cast<int>(maxiter),
                        eff_dconv,
                        maxrotn,
                        static_cast<int>(maxsub),
                        plevel,
                        guess);

                    if (world.rank() == 0) {
                        print("  alpha_", dir_names[d], dir_names[d],
                              "(omega=", omega, ") =", result.alpha,
                              " converged=", result.converged,
                              " iters=", result.iterations, "\n");
                    }

                    // Store for next protocol
                    prev_solutions[key] = std::move(result.response);
                }
            }
        }

        // ----------------------------------------------------------------
        // 5. Final summary
        // ----------------------------------------------------------------
        if (world.rank() == 0) {
            print("\n========================================================");
            print("  FINAL RESULTS (last protocol)");
            print("========================================================\n");

            for (size_t fi = 0; fi < frequencies.size(); fi++) {
                double omega = frequencies[fi];
                print("  omega =", omega, " au");

                // Re-solve is not needed — we already printed per-solve.
                // But let's print a summary table by re-computing alpha
                // from the stored response functions.
                for (int d : directions) {
                    SolveKey key{d, fi};
                    if (prev_solutions.count(key)) {
                        auto& x = prev_solutions[key];
                        bool is_static = (omega == 0.0);
                        auto rtype = is_static ? ResponseType::Static : ResponseType::Full;

                        // Recompute alpha from stored response
                        RealResponseState pert;
                        pert.x_alpha = dipole_perturbation(world, gs, d);
                        if (!is_static) pert.y_alpha = pert.x_alpha;
                        if (!gs.is_spin_restricted()) {
                            pert.x_beta = dipole_perturbation_beta(world, gs, d);
                            if (!is_static) pert.y_beta = pert.x_beta;
                        }

                        double ip_xa = inner(x.x_alpha, pert.x_alpha);
                        double ip_ya = (!is_static && !x.y_alpha.empty())
                            ? inner(x.y_alpha, pert.y_alpha) : 0.0;
                        double ip_xb = (!gs.is_spin_restricted() && !x.x_beta.empty())
                            ? inner(x.x_beta, pert.x_beta) : 0.0;
                        double ip_yb = (!gs.is_spin_restricted() && !is_static && !x.y_beta.empty())
                            ? inner(x.y_beta, pert.y_beta) : 0.0;

                        double af = alpha_factor(rtype, gs.is_spin_restricted());
                        double alpha = af * (ip_xa + ip_ya + ip_xb + ip_yb);

                        print("    alpha_", dir_names[d], dir_names[d], " =", alpha);
                    }
                }
                print("");
            }
        }

        world.gop.fence();
        print_stats(world);

    } catch (const std::exception& e) {
        if (world.rank() == 0) {
            print("ERROR:", e.what());
        }
        finalize();
        return 1;
    }

    finalize();
    return 0;
}
