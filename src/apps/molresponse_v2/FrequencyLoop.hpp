#include "ResponseDebugLogger.hpp"
#include "ResponseDebugLoggerMacros.hpp"
#include "ResponseIO.hpp"
#include "ResponseInitializer.hpp"
#include "ResponseManager.hpp"
#include "ResponseRecord.hpp"
#include "ResponseSolver.hpp"
#include "ResponseSolverUtils.hpp"
#include "ResponseState.hpp"

#include <madness/world/world.h>
#define NOT_IMPLEMENTED_THROW                                                  \
  throw std::runtime_error("This solver is not yet implemented.");

template <typename ResponseType>
bool iterate(World &world, const ResponseManager &response_manager,
             const GroundStateData &g_s, const LinearResponseDescriptor &state,
             const LinearResponsePoint &pt, ResponseType &response,
             ResponseDebugLogger &logger, size_t max_iter,
             double conv_thresh) {
  // using Policy = ResponseSolverPolicy<ResponseType>;

  auto &rvec = response;
  auto &all_x = rvec.flat;

  const double dconv =
      std::max(FunctionDefaults<3>::get_thresh() * 10, conv_thresh);
  auto density_target =
      dconv * static_cast<double>(std::max(size_t(5.0), g_s.molecule.natom()));
  const auto x_residual_target = density_target * 10.0;

  auto vp = perturbation_vector(world, g_s, pt);

  auto &phi0 = g_s.orbitals;
  const auto &orbital_energies = g_s.getEnergies();

  // First difference, Make bsh operators is different for each solver
  auto bsh_ops =
      make_bsh_operators(world, response_manager, pt.frequency(),
                         orbital_energies,
                         static_cast<int>(g_s.orbitals.size()), logger, rvec);

  response_solver solver(
      response_vector_allocator(world, static_cast<int>(all_x.size())),
      /*do_printing*/ false);

  auto drho = compute_density(world, rvec, phi0);
  functionT drho_old;

  for (size_t iter = 0; iter < max_iter; ++iter) {
    logger.begin_iteration(iter);
    drho_old = copy(drho);
    // Inner product of response state
    DEBUG_LOG_VALUE(world, &logger, "<x|x>",
                    ResponseSolverUtils::inner(world, rvec.flat, rvec.flat));
    // 1. Coupled-response equations
    vector_real_function_3d x_new;
    DEBUG_TIMED_BLOCK(world, &logger, "compute_rsh", {
      x_new = CoupledResponseEquations(world, g_s, rvec, vp, bsh_ops,
                                       response_manager, logger);
    });
    // 2. Form residual r = x_new - x
    auto residuals = x_new - all_x;
    vector_real_function_3d kain_x;
    DEBUG_TIMED_BLOCK(world, &logger, "KAIN step",
                      { kain_x = solver.update(all_x, residuals); });

    // 3. Compute norm of difference;
    double res_norm = norm2(world, sub(world, all_x, x_new));
    DEBUG_LOG_VALUE(world, &logger, "res_norm", res_norm);
    // 4. Do step restriction
    if (res_norm > response_manager.params().maxrotn() && false) {
      DEBUG_TIMED_BLOCK(world, &logger, "step_restriction", {
        ResponseSolverUtils::do_step_restriction(
            world, all_x, kain_x, res_norm, "a",
            response_manager.params().maxrotn());
      });
    }
    // 5. Update response vector
    rvec.flat = copy(world, kain_x);
    rvec.sync();
    // 6. Compute updated response density
    drho = compute_density(world, rvec, phi0);
    double drho_change = (drho - drho_old).norm2();
    DEBUG_LOG_VALUE(world, &logger, "drho_change", drho_change);
    auto alpha =
        alpha_factor(rvec) * ResponseSolverUtils::inner(world, rvec.flat, vp);
    DEBUG_LOG_VALUE(world, &logger, "alpha", alpha);
    // 7. Convergence check
    logger.end_iteration();
    if (world.rank() == 0) {
      ResponseSolverUtils::print_iteration_line(iter, res_norm, drho_change,
                                                alpha, density_target,
                                                x_residual_target);
    }
    if (drho_change < density_target && res_norm < x_residual_target) {
      rvec.sync();
      return true;
    }
  }
  return false;
};

bool solve_response_vector(World &world, const ResponseManager &response_manager,
                           const GroundStateData &g_s,
                           const LinearResponseDescriptor &state,
                           const LinearResponsePoint &pt,
                           ResponseVector &response_variant,
                           ResponseDebugLogger &logger, size_t max_iter,
                           double conv_thresh);

void promote_response_vector(World &world, const ResponseVector &x_in,
                             ResponseVector &x_out);

void computeFrequencyLoop(World &world, const ResponseManager &response_manager,
                          const LinearResponseDescriptor &state_desc,
                          size_t thresh_index,
                          const GroundStateData &ground_state,
                          ResponseRecord2 &response_record,
                          ResponseDebugLogger &logger,
                          bool at_final_protocol);
