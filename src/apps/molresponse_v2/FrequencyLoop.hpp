#pragma once

#include "ResponseDebugLogger.hpp"
#include "ResponseDebugLoggerMacros.hpp"
#include "ResponseIO.hpp"
#include "ResponseInitializer.hpp"
#include "ResponseManager.hpp"
#include "ResponseSolver.hpp"
#include "ResponseSolverUtils.hpp"
#include "ResponseState.hpp"

#include <cstdlib>
#include <cmath>
#include <limits>
#include <optional>
#include <string>

#include <madness/world/world.h>
#define NOT_IMPLEMENTED_THROW                                                  \
  throw std::runtime_error("This solver is not yet implemented.");

/// Final point-level solver diagnostics persisted into response metadata.
struct ResponseSolveDiagnostics {
  bool converged = false;
  size_t iterations_performed = 0;
  double final_residual_norm = std::numeric_limits<double>::quiet_NaN();
  double final_density_change = std::numeric_limits<double>::quiet_NaN();
  double final_alpha = std::numeric_limits<double>::quiet_NaN();
  size_t max_consecutive_negative_alpha = 0;
  bool reached_iteration_limit = false;
  bool remove_from_frequency_set = false;
  double residual_remove_cutoff = std::numeric_limits<double>::quiet_NaN();
  std::string failure_reason = "not_evaluated";
};

inline double response_point_stall_timeout_seconds() {
  constexpr double k_default_point_stall_timeout_seconds =
      std::numeric_limits<double>::infinity();
  const char *env_value = std::getenv("MADQC_POINT_STALL_TIMEOUT_S");
  if (env_value == nullptr) {
    return k_default_point_stall_timeout_seconds;
  }
  char *end_ptr = nullptr;
  const double parsed = std::strtod(env_value, &end_ptr);
  if (end_ptr == env_value || parsed <= 0.0) {
    return k_default_point_stall_timeout_seconds;
  }
  return parsed;
}

/// Nonlinear iteration driver for one response vector instance.
///
/// This routine performs the KAIN update loop and returns convergence metadata
/// used by `StateSolvePersistence`.
template <typename ResponseType>
ResponseSolveDiagnostics
iterate(World &world, const ResponseManager &response_manager,
        const GroundStateData &g_s, const LinearResponseDescriptor &state,
        const LinearResponsePoint &pt, ResponseType &response,
        ResponseDebugLogger &logger, size_t max_iter, double conv_thresh) {
  // using Policy = ResponseSolverPolicy<ResponseType>;

  auto &rvec = response;
  auto &all_x = rvec.flat;
  const double point_solve_wall_start = madness::wall_time();
  const double point_stall_timeout_s = response_point_stall_timeout_seconds();

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
  ResponseSolveDiagnostics diagnostics;
  size_t consecutive_negative_alpha = 0;
  std::optional<double> previous_residual_norm;
  size_t large_residual_streak = 0;
  size_t explosive_residual_streak = 0;

  const double safety_residual_cutoff =
      std::max(200.0 * x_residual_target, 50.0);
  const double safety_density_cutoff = std::max(200.0 * density_target, 1.0);
  constexpr size_t min_iter_before_safety_check = 4;
  constexpr size_t large_residual_streak_limit = 3;
  constexpr size_t explosive_residual_streak_limit = 2;
  constexpr double explosive_residual_growth_factor = 3.0;
  bool iteration_table_open = false;

  const bool use_verbose_iteration_table = TimedValueLogger::console_enabled();

  auto emit_iteration_table_start = [&]() {
    if (world.rank() != 0 || iteration_table_open ||
        !use_verbose_iteration_table) {
      return;
    }
    madness::print("ITERATION_TABLE_START");
    madness::print("  state = ", pt.perturbationDescription());
    madness::print("  protocol = ", pt.threshold());
    madness::print("  frequency = ", pt.frequency());
    madness::print("  target_residual = ", x_residual_target);
    madness::print("  target_drho = ", density_target);
    iteration_table_open = true;
  };
  auto emit_iteration_table_end = [&](const std::string &reason) {
    if (world.rank() != 0 || !iteration_table_open ||
        !use_verbose_iteration_table) {
      return;
    }
    ResponseSolverUtils::print_iteration_table_border();
    madness::print("ITERATION_TABLE_END");
    madness::print("ITERATION_END");
    madness::print("  state = ", pt.perturbationDescription());
    madness::print("  protocol = ", pt.threshold());
    madness::print("  frequency = ", pt.frequency());
    madness::print("  reason = ", reason);
    madness::print("  iterations = ", diagnostics.iterations_performed);
    madness::print("  converged = ", diagnostics.converged);
    madness::print("------------------------------------------------------------");
    iteration_table_open = false;
  };

  for (size_t iter = 0; iter < max_iter; ++iter) {
    TimedValueLogger::set_iteration_context(iter);
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
    assign_flat_and_sync(rvec, copy(world, kain_x));
    // 6. Compute updated response density
    drho = compute_density(world, rvec, phi0);
    double drho_change = (drho - drho_old).norm2();
    DEBUG_LOG_VALUE(world, &logger, "drho_change", drho_change);
    auto alpha =
        alpha_factor(rvec) * ResponseSolverUtils::inner(world, rvec.flat, vp);
    DEBUG_LOG_VALUE(world, &logger, "alpha", alpha);
    consecutive_negative_alpha = (alpha < 0.0) ? (consecutive_negative_alpha + 1)
                                               : 0;
    diagnostics.max_consecutive_negative_alpha =
        std::max(diagnostics.max_consecutive_negative_alpha,
                 consecutive_negative_alpha);
    // 7. Convergence check
    logger.end_iteration();
    
    if (world.rank() == 0) {
      if (use_verbose_iteration_table) {
        if (iter == 0) {
          emit_iteration_table_start();
        }
        ResponseSolverUtils::print_iteration_line(iter, res_norm, drho_change,
                                                  alpha, density_target,
                                                  x_residual_target);
      } else {
        madness::print("ITERATION_RESIDUAL state=", pt.perturbationDescription(),
                       " protocol=", pt.threshold(), " frequency=",
                       pt.frequency(), " iter=", iter, " residual=", res_norm,
                       " drho=", drho_change, " alpha=", alpha);
      }
    }
    diagnostics.iterations_performed = iter + 1;
    diagnostics.final_residual_norm = res_norm;
    diagnostics.final_density_change = drho_change;
    diagnostics.final_alpha = alpha;
    const double point_elapsed_s = madness::wall_time() - point_solve_wall_start;
    const bool timeout_enabled =
        std::isfinite(point_stall_timeout_s) && point_stall_timeout_s > 0.0;
    if (timeout_enabled && point_elapsed_s > point_stall_timeout_s) {
      diagnostics.reached_iteration_limit = false;
      diagnostics.failure_reason = "point_stall_timeout";
      emit_iteration_table_end("point_stall_timeout");
      TimedValueLogger::clear_iteration_context();
      if (world.rank() == 0) {
        madness::print("WARN POINT_STALL_TIMEOUT state=",
                       pt.perturbationDescription(), " protocol=",
                       pt.threshold(), " frequency=", pt.frequency(),
                       " elapsed_s=", point_elapsed_s,
                       " timeout_s=", point_stall_timeout_s,
                       " iterations=", diagnostics.iterations_performed);
      }
      return diagnostics;
    }

    const bool has_non_finite_metric =
        !std::isfinite(res_norm) || !std::isfinite(drho_change) ||
        !std::isfinite(alpha);
    if (has_non_finite_metric) {
      diagnostics.reached_iteration_limit = false;
      emit_iteration_table_end("non_finite_metric");
      TimedValueLogger::clear_iteration_context();
      if (world.rank() == 0) {
        madness::print(
            "WARN SAFETY_STOP reason=non_finite_metric state=",
            pt.perturbationDescription(), " protocol=", pt.threshold(),
            " frequency=", pt.frequency(),
            " iterations=", diagnostics.iterations_performed);
      }
      return diagnostics;
    }

    const bool residual_is_large = res_norm > safety_residual_cutoff;
    const bool density_change_is_large = drho_change > safety_density_cutoff;
    if (iter + 1 >= min_iter_before_safety_check && residual_is_large &&
        density_change_is_large) {
      ++large_residual_streak;
    } else {
      large_residual_streak = 0;
    }

    const bool residual_is_explosive =
        previous_residual_norm.has_value() && *previous_residual_norm > 0.0 &&
        res_norm >
            (*previous_residual_norm * explosive_residual_growth_factor) &&
        residual_is_large;
    if (iter + 1 >= min_iter_before_safety_check && residual_is_explosive) {
      ++explosive_residual_streak;
    } else {
      explosive_residual_streak = 0;
    }
    previous_residual_norm = res_norm;

    const bool trigger_safety_stop =
        large_residual_streak >= large_residual_streak_limit ||
        explosive_residual_streak >= explosive_residual_streak_limit;
    if (trigger_safety_stop) {
      diagnostics.reached_iteration_limit = false;
      emit_iteration_table_end("divergent_residual_trend");
      TimedValueLogger::clear_iteration_context();
      if (world.rank() == 0) {
        madness::print(
            "WARN SAFETY_STOP reason=divergent_residual_trend state=",
            pt.perturbationDescription(), " protocol=", pt.threshold(),
            " frequency=", pt.frequency(),
            " iterations=", diagnostics.iterations_performed,
            " residual=", res_norm, " drho=", drho_change,
            " residual_cutoff=", safety_residual_cutoff,
            " drho_cutoff=", safety_density_cutoff);
      }
      return diagnostics;
    }

    if (drho_change < density_target && res_norm < x_residual_target) {
      rvec.sync();
      diagnostics.converged = true;
      emit_iteration_table_end("converged");
      TimedValueLogger::clear_iteration_context();
      return diagnostics;
    }
  }
  diagnostics.reached_iteration_limit =
      (diagnostics.iterations_performed >= max_iter) && !diagnostics.converged;
  emit_iteration_table_end("max_iterations_reached");
  TimedValueLogger::clear_iteration_context();
  return diagnostics;
};

/// Variant-dispatch wrapper over `iterate(...)` for `ResponseVector`.
ResponseSolveDiagnostics
solve_response_vector(World &world, const ResponseManager &response_manager,
                      const GroundStateData &g_s,
                      const LinearResponseDescriptor &state,
                      const LinearResponsePoint &pt,
                      ResponseVector &response_variant,
                      ResponseDebugLogger &logger, size_t max_iter,
                      double conv_thresh);

/// Promote/copy static or restricted response vectors into a target variant shape.
void promote_response_vector(World &world, const ResponseVector &x_in,
                             ResponseVector &x_out);

/// Persistence abstraction used by frequency solve routines.
///
/// Implementations provide restart checks and recording backends (JSON metadata,
/// debug logs, or future alternatives) without coupling solve code to storage.
class StateSolvePersistence {
public:
  virtual ~StateSolvePersistence() = default;
  /// True if a point archive exists.
  [[nodiscard]] virtual bool is_saved(const LinearResponsePoint &pt) const = 0;
  /// True if metadata marks the point as converged.
  [[nodiscard]] virtual bool
  is_converged(const LinearResponsePoint &pt) const = 0;
  /// True if failure policy marked this point for frequency removal.
  [[nodiscard]] virtual bool
  is_removed_from_frequency_set(const LinearResponsePoint &pt) const = 0;
  /// True when removed frequencies should be retried instead of skipped.
  [[nodiscard]] virtual bool force_retry_removed_frequencies() const = 0;
  /// Persist `saved`/`converged` status.
  virtual void record_status(const LinearResponsePoint &pt, bool c) = 0;
  /// Persist final residual, density change, and iteration count.
  virtual void record_solver_diagnostics(const LinearResponsePoint &pt,
                                         const ResponseSolveDiagnostics &d,
                                         bool used_fallback_retry) = 0;
  /// Persist wall/cpu timing for one solved point.
  virtual void record_timing(const LinearResponsePoint &pt, double wall_seconds,
                             double cpu_seconds) = 0;
  /// Persist how the initial guess was seeded for this point.
  virtual void record_restart_provenance(
      const LinearResponsePoint &pt, const std::string &source_kind,
      bool loaded_from_disk, bool promoted_from_static,
      const std::optional<double> &source_protocol,
      const std::optional<double> &source_frequency) = 0;
  /// Access iteration logger bound to the persistence backend.
  virtual ResponseDebugLogger &logger() = 0;
  /// Flush buffered debug logs to backing storage.
  virtual void flush_debug_log(World &world) = 0;
  /// Optional progress/metadata polling hook (no-op by default).
  virtual void maybe_poll_progress(World &world, bool force = false) {
    (void)world;
    (void)force;
  }
};

/// Solve an entire frequency series for one descriptor/protocol index.
///
/// Used in channel-series ownership mode where a subgroup/lane owns complete
/// channels and benefits from nearest-frequency continuation.
void computeFrequencyLoop(World &world, const ResponseManager &response_manager,
                          const LinearResponseDescriptor &state_desc,
                          size_t thresh_index,
                          const GroundStateData &ground_state,
                          StateSolvePersistence &persistence,
                          bool at_final_protocol);

/// Solve one independent `(state, protocol, frequency)` point.
///
/// Used in channel-point ownership mode to expose fine-grained parallelism.
void computeFrequencyPoint(World &world,
                           const ResponseManager &response_manager,
                           const LinearResponseDescriptor &state_desc,
                           size_t thresh_index, size_t freq_index,
                           const GroundStateData &ground_state,
                           StateSolvePersistence &persistence,
                           bool at_final_protocol);
