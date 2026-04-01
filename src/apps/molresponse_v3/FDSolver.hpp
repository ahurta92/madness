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
    PrintLevel print_level = PrintLevel::Normal) {

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

    // Initialize response state (zero)
    auto x = RealResponseState::allocate(world, na, nb, include_y);

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

    // Debug: print Fock matrix once
    // Fock is a replicated Tensor — safe to print from rank 0 only
    if (print_level >= PrintLevel::Debug && world.rank() == 0) {
        const auto& F = gs.focka();
        print("  DIAG focka:");
        for (long i = 0; i < F.dim(0); i++) {
            std::string row = "    [" + std::to_string(i) + "]";
            for (long j = 0; j < F.dim(1); j++) {
                char buf[16];
                std::snprintf(buf, sizeof(buf), " %10.6f", F(i, j));
                row += buf;
            }
            print(row);
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
        auto x_new = fd_iteration(world, type, x, perturbation, gs, coulop,
                                   bsh_alpha_x, bsh_alpha_y,
                                   bsh_beta_x, bsh_beta_y, omega);

        // 2. Flatten for KAIN
        auto x_new_flat = x_new.flat();
        auto x_old_flat = x.flat();

        // 3. Residual
        auto residual_flat = sub(world, x_new_flat, x_old_flat);
        std::vector<double> rnorms = norm2s(world, residual_flat);
        double res_norm = 0.0;
        for (auto n : rnorms) res_norm += n * n;
        res_norm = std::sqrt(res_norm);

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

        // 8. Compute alpha — collective operation (inner involves MPI)
        double afactor = alpha_factor(type, gs.is_spin_restricted());
        double alpha = afactor * madness::inner(x.flat(), perturbation.flat());

        // 9. Report
        if (print_level >= PrintLevel::Verbose && world.rank() == 0) {
            print("  iter=", iter, " res=", res_norm,
                  " drho=", drho, " alpha=", alpha);
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
