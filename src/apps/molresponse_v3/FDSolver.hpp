#ifndef MOLRESPONSE_V3_FDSOLVER_HPP
#define MOLRESPONSE_V3_FDSOLVER_HPP

#include "GroundState.hpp"
#include "ResponseFunctions.hpp"
#include "ResponseKernel.hpp"

#include <madness/mra/mra.h>
#include <madness/mra/nonlinsol.h>

namespace molresponse_v3 {

using namespace madness;

/// Print level for solver output.
/// 0: silent, 1: convergence only, 2: per-iteration, 3: debug (inner products, Fock)
enum class PrintLevel { Silent = 0, Normal = 1, Verbose = 2, Debug = 3 };

/// Result of one frequency-dependent response solve.
struct FDSolveResult {
    RealResponseState response;     // converged response functions
    double alpha;                   // polarizability component
    double density_residual;        // final density change
    double function_residual;       // final function residual norm
    int iterations;                 // iterations used
    bool converged;
};

/// Frequency-dependent response solver.
///
/// Solves (A - omega*B) x = v_p for one perturbation at one frequency
/// using BSH iteration with KAIN acceleration. Handles Static (omega=0),
/// Full (omega!=0), and TDA response types.
/// Project a response state to the current FunctionDefaults k/thresh.
/// Used when carrying over response functions between protocol steps.
inline void project_response_state(World& world, RealResponseState& state) {
    int k = FunctionDefaults<3>::get_k();
    double thresh = FunctionDefaults<3>::get_thresh();

    auto reproject = [&](vector_real_function_3d& v) {
        if (v.empty()) return;
        reconstruct(world, v);
        for (auto& f : v) {
            f = project(f, k, thresh, true);
        }
        truncate(world, v, thresh);
    };

    reproject(state.x_alpha);
    reproject(state.y_alpha);
    reproject(state.x_beta);
    reproject(state.y_beta);
}

inline FDSolveResult fd_solve(
    World& world,
    ResponseType type,
    const RealResponseState& perturbation,
    GroundState& gs,
    double omega,
    int maxiter = 25,
    double dconv = 1e-4,
    double maxrotn = 0.5,
    int maxsub = 10,
    PrintLevel print_level = PrintLevel::Normal,
    const RealResponseState* initial_guess = nullptr) {

    double lo = gs.params().lo();
    double thresh = FunctionDefaults<3>::get_thresh();
    long na = gs.num_alpha();
    long nb = gs.is_spin_restricted() ? 0 : gs.num_beta();
    bool include_y = (type == ResponseType::Full);

    // Convergence targets (matching v2 pattern from FrequencyLoop)
    double density_target = std::max(thresh * 10.0, dconv)
                          * std::max(5.0, static_cast<double>(gs.molecule().natom()));
    double residual_target = density_target * 10.0;

    // Coulomb operator
    auto coulop = poperatorT(CoulombOperatorPtr(world, lo, 0.001 * thresh));

    // BSH operators
    auto bsh_alpha_x = make_bsh_operators(world, gs.energies_alpha(), omega, lo);
    std::vector<poperatorT> bsh_alpha_y;
    if (include_y) {
        bsh_alpha_y = make_bsh_operators_y(world, gs.energies_alpha(), omega, lo);
    }

    std::vector<poperatorT> bsh_beta_x, bsh_beta_y;
    if (!gs.is_spin_restricted()) {
        bsh_beta_x = make_bsh_operators(world, gs.energies_beta(), omega, lo);
        if (include_y) {
            bsh_beta_y = make_bsh_operators_y(world, gs.energies_beta(), omega, lo);
        }
    }

    // Initialize response state: from initial guess or zero
    RealResponseState x;
    if (initial_guess && initial_guess->total_size() > 0) {
        x = *initial_guess;
        // Project to current k/thresh in case protocol changed
        project_response_state(world, x);
        if (print_level >= PrintLevel::Verbose && world.rank() == 0) {
            print("  Starting from initial guess (projected to current protocol)");
        }
    } else {
        x = RealResponseState::allocate(world, na, nb, include_y);
    }

    // KAIN solver operates on the flat vector
    long flat_size = x.total_size();
    auto allocator = vector_function_allocator<double, 3>(world, flat_size);
    XNonlinearSolver<vector_real_function_3d, double,
                      vector_function_allocator<double, 3>>
        kain(allocator, false);
    kain.set_maxsub(maxsub);

    // Previous density for convergence tracking
    real_function_3d rho_old = FunctionFactory<double, 3>(world);

    FDSolveResult result;
    result.converged = false;
    result.iterations = 0;

    // Print alpha equation being solved
    if (print_level >= PrintLevel::Verbose && world.rank() == 0) {
        double af = alpha_factor(type, gs.is_spin_restricted());
        if (gs.is_spin_restricted() && type == ResponseType::Static) {
            print("  Equation: alpha = ", af, " * <x_a|r|phi_a>");
        } else if (gs.is_spin_restricted() && type == ResponseType::Full) {
            print("  Equation: alpha = ", af, " * [<x_a|r|phi_a> + <y_a|r|phi_a>]");
        } else if (!gs.is_spin_restricted() && type == ResponseType::Static) {
            print("  Equation: alpha = ", af, " * [<x_a|r|phi_a> + <x_b|r|phi_b>]");
        } else {
            print("  Equation: alpha = ", af, " * [<x_a|r|phi_a> + <y_a|r|phi_a> + <x_b|r|phi_b> + <y_b|r|phi_b>]");
        }
        print("  n_alpha=", na, " n_beta=", nb, " omega=", omega,
              " include_y=", include_y);
    }

    // Debug: print Fock matrices once
    // Fock is a replicated Tensor — safe to print from rank 0 only
    if (print_level >= PrintLevel::Debug && world.rank() == 0) {
        const auto& Fa = gs.focka();
        print("  DIAG focka:");
        print(Fa);
        if (!gs.is_spin_restricted()) {
            const auto& Fb = gs.fockb();
            print("  DIAG fockb:");
            print(Fb);
        }
    }

    for (int iter = 0; iter < maxiter; iter++) {

        // Debug: <x|x> — collective operation, ALL ranks must participate
        if (print_level >= PrintLevel::Debug) {
            double xx = madness::inner(x.flat(), x.flat());
            if (world.rank() == 0) {
                print("  <x|x> =", xx, " (iter", iter, "start)");
            }
        }

        // 1. One FD iteration step
        int debug_level = static_cast<int>(print_level);
        auto x_new = fd_iteration(world, type, x, perturbation, gs, coulop,
                                   bsh_alpha_x, bsh_alpha_y,
                                   bsh_beta_x, bsh_beta_y, omega,
                                   debug_level);

        // 2. Flatten for KAIN
        auto x_new_flat = x_new.flat();
        auto x_old_flat = x.flat();

        // 3. Residual (total and per-channel)
        auto residual_flat = sub(world, x_new_flat, x_old_flat);
        std::vector<double> rnorms = norm2s(world, residual_flat);
        double res_norm = 0.0;
        for (auto n : rnorms) res_norm += n * n;
        res_norm = std::sqrt(res_norm);

        // Per-channel residuals
        if (print_level >= PrintLevel::Debug) {
            auto res_xa = sub(world, x_new.x_alpha, x.x_alpha);
            double rn_xa = norm2(world, res_xa);

            double rn_ya = 0.0;
            if (include_y) {
                auto res_ya = sub(world, x_new.y_alpha, x.y_alpha);
                rn_ya = norm2(world, res_ya);
            }

            double rn_xb = 0.0, rn_yb = 0.0;
            if (!gs.is_spin_restricted()) {
                auto res_xb = sub(world, x_new.x_beta, x.x_beta);
                rn_xb = norm2(world, res_xb);
                if (include_y) {
                    auto res_yb = sub(world, x_new.y_beta, x.y_beta);
                    rn_yb = norm2(world, res_yb);
                }
            }

            if (world.rank() == 0) {
                print("    res: x_a=", rn_xa, " y_a=", rn_ya,
                      " x_b=", rn_xb, " y_b=", rn_yb, " total=", res_norm);
            }
        }

        // 4. KAIN update
        auto kain_flat = kain.update(x_old_flat, residual_flat);

        // 5. Step restriction
        if (res_norm > maxrotn) {
            double s = maxrotn / res_norm;
            gaxpy(world, s, kain_flat, 1.0 - s, x_old_flat);
            if (print_level >= PrintLevel::Verbose && world.rank() == 0) {
                print("  STEP_RESTRICT iter=", iter, " norm=", res_norm,
                      " scale=", s);
            }
        }

        // 6. Update x from KAIN result
        x.from_flat(kain_flat);
        truncate(world, x.x_alpha, thresh);
        if (include_y) truncate(world, x.y_alpha, thresh);
        if (!gs.is_spin_restricted()) {
            truncate(world, x.x_beta, thresh);
            if (include_y) truncate(world, x.y_beta, thresh);
        }

        // 7. Density change — collective operation (norm2 involves MPI)
        auto rho = compute_response_density(world, type, x, gs);
        double drho = (rho - rho_old).norm2();
        rho_old = rho;

        // 8. Compute alpha with per-channel decomposition
        // All inner products are collective — all ranks participate
        double ip_xa = madness::inner(x.x_alpha, perturbation.x_alpha);
        double ip_ya = include_y ? madness::inner(x.y_alpha, perturbation.y_alpha) : 0.0;
        double ip_xb = (!gs.is_spin_restricted()) ?
            madness::inner(x.x_beta, perturbation.x_beta) : 0.0;
        double ip_yb = (!gs.is_spin_restricted() && include_y) ?
            madness::inner(x.y_beta, perturbation.y_beta) : 0.0;

        double afactor = alpha_factor(type, gs.is_spin_restricted());
        double alpha = afactor * (ip_xa + ip_ya + ip_xb + ip_yb);

        // 9. Report
        if (print_level >= PrintLevel::Verbose && world.rank() == 0) {
            print("  iter=", iter, " res=", res_norm,
                  " drho=", drho, " alpha=", alpha);
        }

        // Debug: per-channel norms and alpha decomposition
        if (print_level >= PrintLevel::Debug) {
            // Norms — collective operations
            auto norm_xa = norm2(world, x.x_alpha);
            auto norm_ya = include_y ? norm2(world, x.y_alpha) : 0.0;
            auto norm_xb = (!gs.is_spin_restricted()) ? norm2(world, x.x_beta) : 0.0;
            auto norm_yb = (!gs.is_spin_restricted() && include_y) ? norm2(world, x.y_beta) : 0.0;

            if (world.rank() == 0) {
                print("    ||x_a||=", norm_xa, " ||y_a||=", norm_ya,
                      " ||x_b||=", norm_xb, " ||y_b||=", norm_yb);
                // Alpha decomposition
                // Static restricted:  alpha = -4 * <x_a|r|phi_a>
                // Dynamic restricted: alpha = -2 * [<x_a|r|phi_a> + <y_a|r|phi_a>]
                // Static UHF:         alpha = -2 * [<x_a|r|phi_a> + <x_b|r|phi_b>]
                // Dynamic UHF:        alpha = -1 * [<x_a|r|phi_a> + <y_a|r|phi_a> + <x_b|r|phi_b> + <y_b|r|phi_b>]
                print("    <x_a|vp_a>=", ip_xa, " <y_a|vp_a>=", ip_ya,
                      " <x_b|vp_b>=", ip_xb, " <y_b|vp_b>=", ip_yb);
                print("    factor=", afactor,
                      " sum=", ip_xa + ip_ya + ip_xb + ip_yb,
                      " alpha=", alpha);
            }
        }

        // 10. Convergence
        result.iterations = iter + 1;
        result.density_residual = drho;
        result.function_residual = res_norm;
        result.alpha = alpha;

        if (drho < density_target && res_norm < residual_target && iter > 0) {
            result.converged = true;
            result.response = std::move(x);
            if (print_level >= PrintLevel::Normal && world.rank() == 0) {
                print("  CONVERGED at iter=", iter, " alpha=", alpha);
            }
            return result;
        }
    }

    // Did not converge
    result.response = std::move(x);
    if (print_level >= PrintLevel::Normal && world.rank() == 0) {
        print("  NOT CONVERGED after ", maxiter, " iterations",
              " alpha=", result.alpha);
    }
    return result;
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_FDSOLVER_HPP
