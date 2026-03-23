// test_tda_h2.cpp — End-to-end TDA excited-state test for H2.
//
// Usage:
//   test_tda_h2 [options] <ground_state_dir>
//
// Positional:
//   ground_state_dir   Directory containing moldft.restartdata (and
//                      optionally moldft.fock.json) from a completed
//                      restricted HF/DFT moldft run.
//
// Options:
//   --num-states N     Number of excitation roots (default: 4)
//   --max-iter N       Max KAIN iterations (default: 30)
//   --dconv X          Convergence threshold on residual norm (default: 1e-4)
//   --print-level N    Verbosity: 0=quiet, 1=per-iter summary, 2=per-state (default: 2)
//   --thresh X         MRA threshold (default: 1e-4)
//   --k N              Wavelet order (default: 7)
//
// What this test does:
//   1. Load GroundStateData from <ground_state_dir>/moldft.restartdata.
//   2. Build localized-Gaussian guess x-states.
//   3. Pack into ResponseBundle<TDARestrictedResponse> (TDA).
//   4. Estimate initial omega from orbital energy weighted average.
//   5. Call iterate_excited(world, bundle, omega, gs, params) — exercises the
//      bundle-level KAIN path added in ResponseBundle.hpp / ExcitedResponse.cpp.
//   6. Print excitation energies in Ha and eV for human inspection.
//
// Expected results for H2 at R=1.4 bohr, restricted HF, thresh~1e-4:
//   - First few excitation energies should be positive (>0) and physically
//     plausible (typically in the range 0.2–0.6 Ha for low-lying states).
//   - The test does NOT hard-fail on exact numerical values; it prints them
//     for comparison against legacy molresponse output.  A separate Python
//     smoke test (test_tda_h2_smoke.py) performs the numerical regression.

#include "ExcitedResponse.hpp"
#include "GroundStateData.hpp"
#include "ResponseBundle.hpp"
#include "ResponseDebugLogger.hpp"
#include "ResponseVector.hpp"

#include <madness/chem/molecule.h>
#include <madness/mra/funcdefaults.h>
#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/world/world.h>

#include <algorithm>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace madness;

// ── Argument parsing ──────────────────────────────────────────────────────────

struct TestArgs {
    std::string gs_dir;
    size_t      num_states   = 4;
    size_t      max_iter     = 30;
    double      dconv        = 1.0e-4;
    int         print_level  = 2;
    double      thresh       = 1.0e-4;
    int         k            = 7;
};

static void print_usage(const char *prog) {
    std::printf("Usage: %s [options] <ground_state_dir>\n", prog);
    std::printf("Options:\n");
    std::printf("  --num-states N   (default 4)\n");
    std::printf("  --max-iter N     (default 30)\n");
    std::printf("  --dconv X        (default 1e-4)\n");
    std::printf("  --print-level N  (default 2)\n");
    std::printf("  --thresh X       (default 1e-4)\n");
    std::printf("  --k N            (default 7)\n");
}

static TestArgs parse_args(int argc, char **argv) {
    TestArgs args;
    bool gs_dir_set = false;
    for (int i = 1; i < argc; ++i) {
        std::string a(argv[i]);
        if (a == "--help" || a == "-h") {
            print_usage(argv[0]);
            std::exit(0);
        } else if (a.rfind("--num-states=", 0) == 0) {
            args.num_states = static_cast<size_t>(std::stoul(a.substr(13)));
        } else if (a == "--num-states" && i + 1 < argc) {
            args.num_states = static_cast<size_t>(std::stoul(argv[++i]));
        } else if (a.rfind("--max-iter=", 0) == 0) {
            args.max_iter = static_cast<size_t>(std::stoul(a.substr(11)));
        } else if (a == "--max-iter" && i + 1 < argc) {
            args.max_iter = static_cast<size_t>(std::stoul(argv[++i]));
        } else if (a.rfind("--dconv=", 0) == 0) {
            args.dconv = std::stod(a.substr(8));
        } else if (a == "--dconv" && i + 1 < argc) {
            args.dconv = std::stod(argv[++i]);
        } else if (a.rfind("--print-level=", 0) == 0) {
            args.print_level = std::stoi(a.substr(14));
        } else if (a == "--print-level" && i + 1 < argc) {
            args.print_level = std::stoi(argv[++i]);
        } else if (a.rfind("--thresh=", 0) == 0) {
            args.thresh = std::stod(a.substr(9));
        } else if (a == "--thresh" && i + 1 < argc) {
            args.thresh = std::stod(argv[++i]);
        } else if (a.rfind("--k=", 0) == 0) {
            args.k = std::stoi(a.substr(4));
        } else if (a == "--k" && i + 1 < argc) {
            args.k = std::stoi(argv[++i]);
        } else if (a[0] != '-') {
            args.gs_dir   = a;
            gs_dir_set    = true;
        } else {
            std::fprintf(stderr, "Unknown option: %s\n", a.c_str());
            print_usage(argv[0]);
            std::exit(1);
        }
    }
    if (!gs_dir_set) {
        std::fprintf(stderr, "Error: ground_state_dir argument is required.\n");
        print_usage(argv[0]);
        std::exit(1);
    }
    // Normalize trailing slash.
    if (!args.gs_dir.empty() && args.gs_dir.back() == '/')
        args.gs_dir.pop_back();
    return args;
}

// ── Main ──────────────────────────────────────────────────────────────────────

int main(int argc, char **argv) {
    World &world = initialize(argc, argv);

    // Parse arguments before startup() so --help works without MPI.
    // argv is still intact after initialize().
    TestArgs args = parse_args(argc, argv);
    const bool enable_debug_log = args.print_level >= 3;
    const std::string debug_log_file =
        args.gs_dir.empty() ? "tda_excited_debug_log.json"
                            : args.gs_dir + "/tda_excited_debug_log.json";
    std::unique_ptr<ResponseDebugLogger> debug_logger;
    if (enable_debug_log)
        debug_logger = std::make_unique<ResponseDebugLogger>(debug_log_file, true);

    try {
        startup(world, argc, argv, /*print_banner=*/true);

        // ── MRA defaults ───────────────────────────────────────────────────
        FunctionDefaults<3>::set_k(args.k);
        FunctionDefaults<3>::set_thresh(args.thresh);
        FunctionDefaults<3>::set_refine(true);
        FunctionDefaults<3>::set_autorefine(false);
        FunctionDefaults<3>::set_apply_randomize(true);
        FunctionDefaults<3>::set_project_randomize(true);
        FunctionDefaults<3>::set_truncate_mode(1);

        if (world.rank() == 0) {
            std::printf("\ntest_tda_h2: TDA excited-state end-to-end test\n");
            std::printf("  ground_state_dir = %s\n", args.gs_dir.c_str());
            std::printf("  num_states       = %zu\n", args.num_states);
            std::printf("  max_iter         = %zu\n", args.max_iter);
            std::printf("  dconv            = %.2e\n", args.dconv);
            std::printf("  thresh           = %.2e\n", args.thresh);
            std::printf("  k                = %d\n", args.k);
            if (enable_debug_log)
                std::printf("  debug_log        = %s\n", debug_log_file.c_str());
            std::printf("\n");
        }

        // ── Load ground state ──────────────────────────────────────────────
        const std::string archive  = args.gs_dir + "/moldft.restartdata";
        const std::string fock_json = args.gs_dir + "/moldft.fock.json";

        if (world.rank() == 0)
            std::printf("Loading ground state from: %s\n", archive.c_str());

        GroundStateData gs(world, archive, Molecule());
        gs.prepareOrbitals(world, args.k, args.thresh);

        const double vtol   = std::max(1.0e-12, 0.1 * args.thresh);
        const auto   coulop = CoulombOperator(world, 1.0e-8, args.thresh);
        gs.computePreliminaries(world, coulop, vtol, fock_json);
        world.gop.fence();

        if (world.rank() == 0) {
            std::printf("  num_orbitals = %ld\n", gs.getNumOrbitals());
            std::printf("  spin_restricted = %s\n",
                        gs.isSpinRestricted() ? "true" : "false");
        }

        if (!gs.isSpinRestricted()) {
            if (world.rank() == 0)
                std::printf("ERROR: test_tda_h2 requires a spin-restricted ground state.\n");
            finalize();
            return 1;
        }

        const size_t N = static_cast<size_t>(gs.getNumOrbitals());
        const size_t M = args.num_states;

        // ── Build guess states ─────────────────────────────────────────────
        const size_t n_gen = std::max<size_t>(2 * M, 8);
        if (world.rank() == 0)
            std::printf("\nBuilding %zu localized Gaussian guess states...\n", n_gen);

        auto x_guesses = build_fresh_guess_x_states(
            world, n_gen, N, gs.getMolecule(), gs.Qhat, 42u);

        if (x_guesses.size() < M) {
            if (world.rank() == 0)
                std::printf("ERROR: only %zu linearly independent guesses (need %zu).\n",
                            x_guesses.size(), M);
            finalize();
            return 1;
        }
        x_guesses.resize(M);

        // ── Pack into ResponseBundle<TDARestrictedResponse> ───────────────
        auto packed = pack_guess_states<TDARestrictedResponse>(world, x_guesses, N);
        ResponseBundle<TDARestrictedResponse> bundle(std::move(packed));

        // ── Estimate initial omega ─────────────────────────────────────────
        std::vector<double> omega =
            estimate_initial_omega(world, x_guesses, gs.getEnergies());

        if (world.rank() == 0) {
            std::printf("Initial omega estimates:\n");
            for (size_t i = 0; i < M; ++i)
                std::printf("  state %2zu: %.6f Ha  (%.4f eV)\n",
                            i, omega[i], omega[i] * 27.2114);
            std::printf("\n");
        }

        // ── Run TDA iteration ──────────────────────────────────────────────
        ExcitedSolverParams params;
        params.maxiter      = args.max_iter;
        params.maxsub       = 8;
        params.dconv        = args.dconv;
        params.print_level  = args.print_level;
        params.tda          = true;
        params.max_rotation = 0.25;
        params.debug_logger = debug_logger.get();

        if (debug_logger) {
            TimedValueLogger::set_console_enabled(true);
            debug_logger->start_named_state("excited_tda_bundle", args.thresh, 0.0);
        }

        if (world.rank() == 0)
            std::printf("Running TDA iterate_excited (bundle-level KAIN)...\n\n");

        auto diag = iterate_excited(world, bundle, omega, gs, params);
        if (debug_logger && world.rank() == 0) {
            debug_logger->finalize_state();
            debug_logger->write_to_disk();
        }

        // ── Report results ─────────────────────────────────────────────────
        if (world.rank() == 0) {
            std::printf("\n=== TDA H2 excitation energies ===\n");
            for (size_t i = 0; i < M; ++i)
                std::printf("  State %2zu:  omega = %12.8f Ha  (%8.4f eV)\n",
                            i, omega[i], omega[i] * 27.2114);
            std::printf("\nConverged: %s  (iterations used: %zu)\n",
                        diag.converged ? "YES" : "NO",
                        diag.iterations_used);
            if (!diag.iteration_max_residuals.empty())
                std::printf("Final max residual: %.4e\n",
                            diag.iteration_max_residuals.back());
            std::printf("\n");

            // Basic sanity checks (non-fatal).
            bool sane = true;
            for (size_t i = 0; i < M; ++i) {
                if (omega[i] <= 0.0) {
                    std::printf("WARN: state %zu has non-positive omega = %g\n",
                                i, omega[i]);
                    sane = false;
                }
            }
            // Check that energies are monotonically non-decreasing.
            for (size_t i = 1; i < M; ++i) {
                if (omega[i] < omega[i-1] - 1.0e-6) {
                    std::printf("WARN: omega not sorted: omega[%zu]=%.6f > omega[%zu]=%.6f\n",
                                i-1, omega[i-1], i, omega[i]);
                    sane = false;
                }
            }
            if (sane && diag.converged)
                std::printf("PASS: all energies positive and monotone.\n");
        }
    } catch (const MadnessException &e) {
        if (debug_logger && world.rank() == 0) {
            debug_logger->finalize_state();
            debug_logger->write_to_disk();
        }
        if (world.rank() == 0)
            std::printf("MADNESS exception: %s  (%s:%d)\n",
                        e.msg, e.filename, e.line);
        finalize();
        return 1;
    } catch (const std::exception &e) {
        if (debug_logger && world.rank() == 0) {
            debug_logger->finalize_state();
            debug_logger->write_to_disk();
        }
        if (world.rank() == 0)
            std::printf("Exception: %s\n", e.what());
        finalize();
        return 1;
    }

    world.gop.fence();
    world.gop.fence();
    finalize();
    return 0;
}
