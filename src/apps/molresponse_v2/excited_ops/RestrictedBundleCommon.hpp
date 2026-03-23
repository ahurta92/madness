#pragma once

#include "../ExcitedResponse.hpp"
#include "../ResponseDebugLoggerMacros.hpp"
#include "../ResponseSolverConstants.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace excited_ops {

template <typename ResponseType>
struct BundlePotentials {
    std::vector<ResponseType> lambda;
    std::vector<ResponseType> v0;
    std::vector<ResponseType> gamma;
};

struct StateIterationMetrics {
    double residual       = 0.0;
    double density_change = 0.0;
};

template <typename ResponseType>
struct BundleModel;

inline bool excited_debug_enabled(const ExcitedSolverParams &params) {
    return params.debug_logger != nullptr || params.print_level >= 3;
}

inline std::string state_phase_key(size_t state_index, const std::string &phase) {
    return "state_" + std::to_string(state_index) + "." + phase;
}

inline void log_phase_marker(madness::World &world,
                             ResponseDebugLogger *logger,
                             bool enabled,
                             const std::string &key,
                             double value = 1.0) {
    if (!enabled)
        return;
    TimedValueLogger marker(world, key, logger);
    marker.log(value);
}

template <typename ValueType>
inline void log_phase_value(madness::World &world,
                            ResponseDebugLogger *logger,
                            bool enabled,
                            const std::string &key,
                            const ValueType &value) {
    if (!enabled)
        return;
    TimedValueLogger marker(world, key, logger);
    marker.log(value);
}

template <typename Func>
decltype(auto) timed_phase(madness::World &world,
                           ResponseDebugLogger *logger,
                           bool enabled,
                           const std::string &key,
                           Func &&func) {
    if (!enabled)
        return std::forward<Func>(func)();

    TimedValueLogger timer(world, key, logger);
    if constexpr (std::is_void_v<std::invoke_result_t<Func &>>) {
        std::forward<Func>(func)();
        timer.log();
    } else {
        auto result = std::forward<Func>(func)();
        timer.log();
        return result;
    }
}

template <typename ResponseType>
[[nodiscard]] inline std::vector<ResponseType>
add_kinetic_to_lambda(madness::World &world,
                      const std::vector<ResponseType> &lambdas) {
    std::vector<ResponseType> lambdas_with_kinetic;
    lambdas_with_kinetic.reserve(lambdas.size());
    for (const auto &lambda : lambdas) {
        auto kinetic_flat = apply_kinetic_flat(world, response_all(lambda));
        auto lambda_all  = response_all(lambda);
        madness::gaxpy(world, 1.0, lambda_all, 1.0, kinetic_flat);
        lambdas_with_kinetic.push_back(
            clone_with_all(world, lambda, std::move(lambda_all)));
    }
    world.gop.fence();
    return lambdas_with_kinetic;
}

template <typename ResponseType>
[[nodiscard]] inline bool
rotate_excited_space_common(madness::World &world,
                            std::vector<ResponseType> &states,
                            BundlePotentials<ResponseType> &potentials,
                            std::vector<double> &omega,
                            int print_level,
                            ResponseDebugLogger *logger) {
    const bool debug_enabled = (logger != nullptr) || (print_level >= 3);

    log_phase_marker(world, logger, debug_enabled, "rotate.add_kinetic.begin");
    auto lambdas_with_kinetic = timed_phase(
        world, logger, debug_enabled, "rotate.add_kinetic", [&] {
            return add_kinetic_to_lambda(world, potentials.lambda);
        });

    madness::Tensor<double> overlap;
    madness::Tensor<double> energy;
    log_phase_marker(world, logger, debug_enabled, "rotate.build_matrices.begin");
    timed_phase(world, logger, debug_enabled, "rotate.build_matrices", [&] {
        build_rotation_matrices(world, states, lambdas_with_kinetic, overlap, energy);
    });

    if (print_level >= 2 && world.rank() == 0) {
        madness::print("  S matrix:");
        madness::print(overlap);
        madness::print("  A matrix:");
        madness::print(energy);
    }

    log_phase_marker(world, logger, debug_enabled, "rotate.diagonalize.begin");
    auto diag = timed_phase(world, logger, debug_enabled, "rotate.diagonalize", [&] {
        return diagonalize_bundle(world, overlap, energy, print_level);
    });
    if (!diag.success) {
        log_phase_value(world, logger, debug_enabled, "rotate.success", 0.0);
        if (print_level > 0 && world.rank() == 0)
            madness::print("iterate_excited: diagonalisation failed");
        return false;
    }

    omega = diag.omega;
    log_phase_value(world, logger, debug_enabled, "rotate.success", 1.0);
    log_phase_value(world, logger, debug_enabled, "omega.post_diag", omega);

    log_phase_marker(world, logger, debug_enabled, "rotate.states.begin");
    timed_phase(world, logger, debug_enabled, "rotate.states", [&] {
        states = rotate_bundle(world, states, diag.U);
    });

    log_phase_marker(world, logger, debug_enabled, "rotate.lambda.begin");
    timed_phase(world, logger, debug_enabled, "rotate.lambda", [&] {
        potentials.lambda = rotate_bundle(world, potentials.lambda, diag.U);
    });

    log_phase_marker(world, logger, debug_enabled, "rotate.v0.begin");
    timed_phase(world, logger, debug_enabled, "rotate.v0", [&] {
        potentials.v0 = rotate_bundle(world, potentials.v0, diag.U);
    });

    log_phase_marker(world, logger, debug_enabled, "rotate.gamma.begin");
    timed_phase(world, logger, debug_enabled, "rotate.gamma", [&] {
        potentials.gamma = rotate_bundle(world, potentials.gamma, diag.U);
    });

    if (print_level >= 2 && world.rank() == 0) {
        madness::print("  omega after diag:");
        for (size_t i = 0; i < omega.size(); ++i)
            madness::print("    state=", i,
                           " omega=", omega[i],
                           " eV=", omega[i] * 27.2114);
    }
    return true;
}

template <typename ResponseType>
[[nodiscard]] inline std::vector<madness::vector_real_function_3d>
copy_bundle_flats(madness::World &world,
                  const std::vector<ResponseType> &states) {
    std::vector<madness::vector_real_function_3d> prev_flats;
    prev_flats.reserve(states.size());
    for (const auto &state : states)
        prev_flats.push_back(copy(world, response_all(state), false));
    world.gop.fence();
    return prev_flats;
}

template <typename ResponseType>
[[nodiscard]] inline std::vector<madness::real_function_3d>
compute_bundle_densities(madness::World &world,
                         const std::vector<ResponseType> &states,
                         const GroundStateData &gs) {
    std::vector<madness::real_function_3d> densities;
    densities.reserve(states.size());
    for (const auto &state : states)
        densities.push_back(compute_response_density(world, state, gs));
    return densities;
}

template <typename ResponseType>
[[nodiscard]] inline ResponseBundle<ResponseType>
make_bundle_from_flats(madness::World &world,
                       const std::vector<ResponseType> &prototypes,
                       const std::vector<madness::vector_real_function_3d> &flats) {
    MADNESS_ASSERT(prototypes.size() == flats.size());
    std::vector<ResponseType> states;
    states.reserve(prototypes.size());
    for (size_t i = 0; i < prototypes.size(); ++i)
        states.push_back(
            clone_with_all(world, prototypes[i], copy(world, flats[i], false)));
    world.gop.fence();
    return ResponseBundle<ResponseType>(std::move(states));
}

template <typename ResponseType>
[[nodiscard]] inline StateIterationMetrics
compute_iteration_metrics_common(madness::World &world,
                                 double residual,
                                 const ResponseType &accepted_state,
                                 const madness::real_function_3d &density_before,
                                 const GroundStateData &gs) {
    StateIterationMetrics metrics;
    metrics.residual = residual;

    const auto density_after = compute_response_density(world, accepted_state, gs);
    const auto delta_rho     = density_after - density_before;
    metrics.density_change   = std::sqrt(std::max(0.0, delta_rho.inner(delta_rho)));
    return metrics;
}

template <typename ResponseType>
inline void print_potential_diagnostics(madness::World &world,
                                        const std::vector<ResponseType> &states,
                                        const BundlePotentials<ResponseType> &potentials,
                                        int print_level) {
    if (print_level < 2 || world.rank() != 0)
        return;

    for (size_t i = 0; i < states.size(); ++i) {
        madness::print("  POT_DIAG_BEGIN state=", i, " term=<x|gamma>");
        const auto xgx  = inner(world, response_all(states[i]), response_all(potentials.gamma[i]));
        madness::print("  POT_DIAG_DONE state=", i, " term=<x|gamma> value=", xgx.sum());

        madness::print("  POT_DIAG_BEGIN state=", i, " term=<x|v0>");
        const auto xv0x = inner(world, response_all(states[i]), response_all(potentials.v0[i]));
        madness::print("  POT_DIAG_DONE state=", i, " term=<x|v0> value=", xv0x.sum());

        madness::print("  POT_DIAG_BEGIN state=", i, " term=<x|lambda>");
        const auto xlx  = inner(world, response_all(states[i]), response_all(potentials.lambda[i]));
        madness::print("  POT_DIAG_DONE state=", i, " term=<x|lambda> value=", xlx.sum());

        madness::print("  DIAG state=", i,
                       " <x|gamma>=", xgx.sum(),
                       " <x|v0>=", xv0x.sum(),
                       " <x|lambda>=", xlx.sum());
    }
}

template <typename ResponseType>
ExcitedIterDiagnostics
iterate_trial_common(madness::World &world,
                     ResponseBundle<ResponseType> &bundle,
                     std::vector<double> &omega,
                     const GroundStateData &gs,
                     const ExcitedSolverParams &params) {
    static_assert(!response_is_unrestricted_v<ResponseType>,
                  "iterate_excited: unrestricted not yet implemented");

    using Model = BundleModel<ResponseType>;

    auto &states = bundle.states();
    const size_t num_states = states.size();

    ExcitedIterDiagnostics result;
    if (num_states == 0)
        return result;

    if (omega.size() < num_states)
        omega.resize(num_states, 0.0);

    const bool debug_enabled = excited_debug_enabled(params);
    auto *logger = params.debug_logger;

    const size_t num_orbitals = bundle.num_orbitals();
    bundle_solver<ResponseType> solver(
        ResponseBundleAllocator<ResponseType>{world, num_states, num_orbitals},
        false);
    solver.set_maxsub(static_cast<int>(params.maxsub));

    for (size_t iter = 0; iter < params.maxiter; ++iter) {
        TimedValueLogger::set_iteration_context(iter);
        if (logger != nullptr)
            logger->begin_iteration(iter);

        log_phase_marker(world, logger, debug_enabled, "iteration.begin");
        log_phase_value(world, logger, debug_enabled, "iteration.num_states",
                        static_cast<double>(num_states));
        log_phase_value(world, logger, debug_enabled, "omega.pre_iter", omega);

        log_phase_marker(world, logger, debug_enabled, "compute_potentials.begin");
        auto potentials = timed_phase(world, logger, debug_enabled, "compute_potentials", [&] {
            return Model::compute_bundle_potentials(
                world, states, gs, params.print_level, logger);
        });
        print_potential_diagnostics(world, states, potentials, params.print_level);

        log_phase_marker(world, logger, debug_enabled, "rotate.begin");
        const bool rotate_ok = timed_phase(world, logger, debug_enabled, "rotate", [&] {
            return Model::rotate_excited_space(
                world, states, potentials, omega, params.print_level, logger);
        });
        if (!rotate_ok) {
            if (logger != nullptr)
                logger->end_iteration();
            TimedValueLogger::clear_iteration_context();
            break;
        }

        log_phase_marker(world, logger, debug_enabled, "copy_previous_flats.begin");
        auto previous_flats = timed_phase(world, logger, debug_enabled,
                                          "copy_previous_flats", [&] {
            return copy_bundle_flats(world, states);
        });

        log_phase_marker(world, logger, debug_enabled, "compute_density_before.begin");
        auto density_before = timed_phase(world, logger, debug_enabled,
                                          "compute_density_before", [&] {
            return compute_bundle_densities(world, states, gs);
        });

        std::vector<ResponseType> updated_states;
        std::vector<double> residuals(num_states, 0.0);
        updated_states.reserve(num_states);

        for (size_t i = 0; i < num_states; ++i) {
            const std::string theta_key = state_phase_key(i, "build_theta");
            log_phase_marker(world, logger, debug_enabled, theta_key + ".begin");
            auto theta_state = timed_phase(world, logger, debug_enabled, theta_key, [&] {
                return Model::build_theta(
                    world, states[i], potentials.v0[i], potentials.gamma[i], gs);
            });

            const std::string bsh_key = state_phase_key(i, "bsh_update");
            log_phase_marker(world, logger, debug_enabled, bsh_key + ".begin");
            auto updated_state = timed_phase(world, logger, debug_enabled, bsh_key, [&] {
                return Model::bsh_update_excited(world, states[i], theta_state, omega[i], gs);
            });

            const std::string residual_flat_key = state_phase_key(i, "residual_flat");
            log_phase_marker(world, logger, debug_enabled, residual_flat_key + ".begin");
            auto residual_flat = timed_phase(world, logger, debug_enabled, residual_flat_key, [&] {
                return sub(world, response_all(updated_state), previous_flats[i]);
            });

            residuals[i] = timed_phase(world, logger, debug_enabled,
                                       state_phase_key(i, "residual_norm"), [&] {
                return state_norm(world, residual_flat);
            });
            log_phase_value(world, logger, debug_enabled,
                            state_phase_key(i, "residual_value"), residuals[i]);
            updated_states.push_back(std::move(updated_state));
        }

        timed_phase(world, logger, debug_enabled, "update_loop.fence", [&] {
            world.gop.fence();
        });

        log_phase_marker(world, logger, debug_enabled, "make_prev_bundle.begin");
        auto prev_bundle = timed_phase(world, logger, debug_enabled, "make_prev_bundle", [&] {
            return make_bundle_from_flats(world, states, previous_flats);
        });

        ResponseBundle<ResponseType> updated_bundle(std::move(updated_states));

        log_phase_marker(world, logger, debug_enabled, "bundle_residual.begin");
        auto residual_bundle = timed_phase(world, logger, debug_enabled, "bundle_residual", [&] {
            return updated_bundle - prev_bundle;
        });

        log_phase_marker(world, logger, debug_enabled, "kain_update.begin");
        auto kain_bundle = timed_phase(world, logger, debug_enabled, "kain_update", [&] {
            return solver.update(prev_bundle, residual_bundle);
        });

        const ResponseBundle<ResponseType> &candidate_bundle =
            (iter >= 2) ? kain_bundle : updated_bundle;
        log_phase_value(world, logger, debug_enabled, "candidate_bundle.is_kain",
                        (iter >= 2) ? 1.0 : 0.0);

        std::vector<double> density_changes(num_states, 0.0);
        for (size_t i = 0; i < num_states; ++i) {
            const std::string copy_key = state_phase_key(i, "candidate_copy");
            log_phase_marker(world, logger, debug_enabled, copy_key + ".begin");
            auto candidate_flat = timed_phase(world, logger, debug_enabled, copy_key, [&] {
                return copy(world, response_all(candidate_bundle[i]), false);
            });

            timed_phase(world, logger, debug_enabled, state_phase_key(i, "candidate_copy.fence"), [&] {
                world.gop.fence();
            });

            if (params.max_rotation > 0.0) {
                timed_phase(world, logger, debug_enabled,
                            state_phase_key(i, "step_restriction"), [&] {
                    auto step_flat = sub(world, candidate_flat, response_all(prev_bundle[i]));
                    const double step_norm = state_norm(world, step_flat);
                    log_phase_value(world, logger, debug_enabled,
                                    state_phase_key(i, "step_norm"), step_norm);
                    if (step_norm > params.max_rotation) {
                        const double scale = params.max_rotation / step_norm;
                        madness::gaxpy(world, scale, candidate_flat,
                                       1.0 - scale, response_all(prev_bundle[i]));
                        log_phase_value(world, logger, debug_enabled,
                                        state_phase_key(i, "step_clamped"), 1.0);
                    } else {
                        log_phase_value(world, logger, debug_enabled,
                                        state_phase_key(i, "step_clamped"), 0.0);
                    }
                });
            }

            const std::string assign_key = state_phase_key(i, "assign_project_normalize");
            log_phase_marker(world, logger, debug_enabled, assign_key + ".begin");
            timed_phase(world, logger, debug_enabled, assign_key, [&] {
                assign_all_and_sync(states[i], std::move(candidate_flat));
                project_response_channels(world, states[i], gs);
                normalize_response_metric(world, states[i]);
            });

            const std::string metrics_key = state_phase_key(i, "metrics");
            log_phase_marker(world, logger, debug_enabled, metrics_key + ".begin");
            const auto metrics = timed_phase(world, logger, debug_enabled, metrics_key, [&] {
                return Model::compute_iteration_metrics(
                    world, residuals[i], states[i], density_before[i], gs);
            });
            residuals[i]       = metrics.residual;
            density_changes[i] = metrics.density_change;

            log_phase_value(world, logger, debug_enabled,
                            state_phase_key(i, "residual_post"), residuals[i]);
            log_phase_value(world, logger, debug_enabled,
                            state_phase_key(i, "density_change"), density_changes[i]);
            log_phase_value(world, logger, debug_enabled,
                            state_phase_key(i, "omega"), omega[i]);
        }

        timed_phase(world, logger, debug_enabled, "post_update.fence", [&] {
            world.gop.fence();
        });

        const double max_residual =
            *std::max_element(residuals.begin(), residuals.end());
        const double max_density_change =
            *std::max_element(density_changes.begin(), density_changes.end());
        result.iteration_max_residuals.push_back(max_residual);
        result.iteration_max_density_changes.push_back(max_density_change);
        result.iterations_used = iter + 1;

        log_phase_value(world, logger, debug_enabled, "iteration.max_residual",
                        max_residual);
        log_phase_value(world, logger, debug_enabled,
                        "iteration.max_density_change", max_density_change);

        if (logger != nullptr)
            logger->end_iteration();
        TimedValueLogger::clear_iteration_context();

        if (params.print_level > 0 && world.rank() == 0) {
            madness::print("iterate_excited"
                           " iter=", iter,
                           " max_res=", max_residual,
                           " max_drho=", max_density_change);
            for (size_t i = 0; i < num_states; ++i)
                madness::print("  state=", i,
                               " omega=", omega[i],
                               " eV=", omega[i] * 27.2114,
                               " res=", residuals[i]);
        }

        if (max_residual <= params.dconv) {
            result.converged = true;
            break;
        }
    }

    result.residual_norms.assign(num_states, 0.0);
    result.density_change_norms.assign(num_states, 0.0);
    if (!result.iteration_max_residuals.empty()) {
        const double last_residual = result.iteration_max_residuals.back();
        const double last_density_change =
            result.iteration_max_density_changes.empty()
                ? 0.0
                : result.iteration_max_density_changes.back();
        std::fill(result.residual_norms.begin(),
                  result.residual_norms.end(),
                  last_residual);
        std::fill(result.density_change_norms.begin(),
                  result.density_change_norms.end(),
                  last_density_change);
    }
    return result;
}

} // namespace excited_ops
