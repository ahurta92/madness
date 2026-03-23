#pragma once

#include "RestrictedBundleCommon.hpp"

namespace excited_ops {

template <>
struct BundleModel<TDARestrictedResponse> {
    using response_type = TDARestrictedResponse;

    [[nodiscard]] static BundlePotentials<response_type>
    compute_bundle_potentials(madness::World &world,
                              const std::vector<response_type> &states,
                              const GroundStateData &gs,
                              int print_level,
                              ResponseDebugLogger *logger) {
        const bool debug_enabled = (logger != nullptr) || (print_level >= 3);

        BundlePotentials<response_type> potentials;
        potentials.lambda.reserve(states.size());
        potentials.v0.reserve(states.size());
        potentials.gamma.reserve(states.size());

        for (size_t state_index = 0; state_index < states.size(); ++state_index) {
            const auto phase_key = state_phase_key(state_index, "compute_potentials");
            log_phase_marker(world, logger, debug_enabled, phase_key + ".begin");

            response_type lambda_out;
            response_type v0_out;
            response_type gamma_out;
            timed_phase(world, logger, debug_enabled, phase_key, [&] {
                compute_state_response_potentials(
                    world, states[state_index], gs, lambda_out, v0_out, gamma_out);
            });
            potentials.lambda.push_back(std::move(lambda_out));
            potentials.v0.push_back(std::move(v0_out));
            potentials.gamma.push_back(std::move(gamma_out));
        }
        return potentials;
    }

    [[nodiscard]] static bool
    rotate_excited_space(madness::World &world,
                         std::vector<response_type> &states,
                         BundlePotentials<response_type> &potentials,
                         std::vector<double> &omega,
                         int print_level,
                         ResponseDebugLogger *logger) {
        return rotate_excited_space_common(
            world, states, potentials, omega, print_level, logger);
    }

    [[nodiscard]] static response_type
    build_theta(madness::World &world,
                const response_type &state,
                const response_type &v0_state,
                const response_type &gamma_state,
                const GroundStateData &gs) {
        auto epsilon_flat = apply_hamiltonian_no_diag(world, state, gs);
        auto theta_all = response_all(v0_state) - epsilon_flat + response_all(gamma_state);
        return clone_with_all(world, state, std::move(theta_all));
    }

    [[nodiscard]] static response_type
    bsh_update_excited(madness::World &world,
                       const response_type &state,
                       const response_type &theta_state,
                       double omega,
                       const GroundStateData &gs) {
        std::vector<double> orbital_shifts;
        auto bsh_ops = make_excited_bsh_operators(world, state, omega, gs, orbital_shifts);

        auto theta_all = response_all(theta_state);
        for (size_t p = 0; p < theta_all.size() && p < orbital_shifts.size(); ++p)
            theta_all[p] = ResponseSolverConstants::k_bsh_residual_prefactor *
                           (theta_all[p] + orbital_shifts[p] * response_all(state)[p]);

        auto updated_all = apply(world, bsh_ops, theta_all);
        truncate(world, updated_all);
        world.gop.fence();
        return clone_with_all(world, state, std::move(updated_all));
    }

    [[nodiscard]] static StateIterationMetrics
    compute_iteration_metrics(madness::World &world,
                              double residual,
                              const response_type &accepted_state,
                              const madness::real_function_3d &density_before,
                              const GroundStateData &gs) {
        return compute_iteration_metrics_common(
            world, residual, accepted_state, density_before, gs);
    }

    [[nodiscard]] static ExcitedIterDiagnostics
    iterate_trial(madness::World &world,
                  ResponseBundle<response_type> &bundle,
                  std::vector<double> &omega,
                  const GroundStateData &gs,
                  const ExcitedSolverParams &params) {
        return iterate_trial_common(world, bundle, omega, gs, params);
    }

private:
    static void compute_state_response_potentials(
        madness::World &world,
        const response_type &state,
        const GroundStateData &gs,
        response_type &lambda_out,
        response_type &v0_out,
        response_type &gamma_out) {
        const auto &x_state = response_x(state);
        const auto &phi0    = gs.orbitals;
        const double c_xc   = gs.xcf_.hf_exchange_coefficient();
        const double thresh = madness::FunctionDefaults<3>::get_thresh();

        auto k0x   = K(world, phi0, phi0)(x_state);
        auto xphi  = mul(world, x_state, phi0, true);
        auto rho   = sum(world, xphi, true);
        auto j_rho = apply(madness::CoulombOperator(
                               world,
                               ResponseSolverConstants::k_coulomb_lo,
                               thresh),
                           rho);
        auto gamma_flat = gs.Qhat(
            2.0 * (j_rho * phi0) - c_xc * K(world, phi0, x_state)(phi0));
        auto v0_all      = gs.V_local * response_all(state) - c_xc * k0x;
        auto epsilon_all = apply_hamiltonian_no_diag(world, state, gs);
        auto lambda_all  = v0_all - epsilon_all + gamma_flat;

        v0_out     = clone_with_all(world, state, std::move(v0_all));
        gamma_out  = clone_with_all(world, state, std::move(gamma_flat));
        lambda_out = clone_with_all(world, state, std::move(lambda_all));
    }

    [[nodiscard]] static std::vector<poperatorT>
    make_excited_bsh_operators(madness::World &world,
                               const response_type &state,
                               double omega,
                               const GroundStateData &gs,
                               std::vector<double> &orbital_shifts_out) {
        constexpr double lo           = 1.0e-8;
        constexpr double shift_factor = 0.05;
        const double thresh = madness::FunctionDefaults<3>::get_thresh();
        const size_t n      = gs.orbitals.size();

        std::vector<poperatorT> ops;
        orbital_shifts_out.clear();
        ops.reserve(n);
        orbital_shifts_out.reserve(n);

        for (size_t p = 0; p < n; ++p) {
            const double energy = gs.getEnergies()(long(p));
            double shift = 0.0;
            if ((energy + omega) > 0.0)
                shift = -(energy + omega + shift_factor);
            const double denom      = energy + omega + shift;
            const double stabilized = std::min(denom, -1.0e-8);
            const double mu         = std::sqrt(std::max(1.0e-16, -2.0 * stabilized));
            ops.emplace_back(poperatorT(BSHOperatorPtr3D(world, mu, lo, thresh)));
            orbital_shifts_out.push_back(shift);
        }

        const auto all_size = response_all(state).size();
        if (ops.size() != all_size) {
            orbital_shifts_out.resize(all_size, 0.0);
            if (ops.size() > all_size)
                ops.resize(all_size);
        }
        return ops;
    }
};

inline ExcitedIterDiagnostics
iterate_trial(madness::World &world,
              ResponseBundle<TDARestrictedResponse> &bundle,
              std::vector<double> &omega,
              const GroundStateData &gs,
              const ExcitedSolverParams &params) {
    return BundleModel<TDARestrictedResponse>::iterate_trial(
        world, bundle, omega, gs, params);
}

} // namespace excited_ops
