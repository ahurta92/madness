#pragma once
#include <madness/chem/CalculationParameters.h>

namespace madness {
inline CalculationParameters read_scf_parameters(World& world, const commandlineparser& parser) {
    CalculationParameters result;
    result.enable_scf_inputs();
    result.read_input_and_commandline_options(world, parser, "dft");
    return result;
}

inline iteration::Request scf_iteration_request(const CalculationParameters& parameters) {
    auto p = parameters;
    p.enable_scf_inputs();
    const auto written = p.to_json_if_precedence("defined");
    for (auto i = written.begin(); i != written.end(); ++i) {
        // Retirement checks apply to programmatic requests as well. Numeric
        // validity is checked below against typed, resolved values.
        const auto& key = i.key();
        if (key == "dconv" || key == "econv" || key == "save" || key == "convergence_criteria" ||
            key == "conv_only_dens" || key.starts_with("accuracy.") || key.starts_with("progress.") ||
            key == "convergence.use_default_criteria" || key == "diagnostic.iterations_per_rung" ||
            key == "restart.start_rung" || key == "maxiter_per_rung" ||
            iteration::retired_convergence_input(key))
            p.validate_input_value(key, "0");
    }
    // Only written targets enter the request; the model supplies the default
    // criteria, their targets and the limit rule (frozen occupied: BSH only).
    iteration::Request request;
    for (auto metric : {iteration::Metric::bsh_residual, iteration::Metric::density_change,
                        iteration::Metric::energy_change, iteration::Metric::kinetic_change}) {
        const auto key = std::string("convergence.") + iteration::metric_name(metric);
        if (!p.is_user_defined(key)) continue;
        const double target = p.get<double>(key);
        if (!iteration::finite_nonnegative(target)) throw std::invalid_argument("invalid target: " + key);
        request.targets[metric] = target;
    }
    if (p.is_user_defined("convergence.floor_factor")) {
        request.floor_factor = p.get<double>("convergence.floor_factor");
        if (!iteration::finite_nonnegative(*request.floor_factor))
            throw std::invalid_argument("convergence.floor_factor must be finite and not negative");
    }
    const int interval = p.get<int>("checkpoint.every_iterations");
    if (p.maxiter() < 0 || interval < 0) throw std::invalid_argument("negative SCF update budget or checkpoint cadence");
    request.max_updates = static_cast<std::size_t>(p.maxiter());
    request.checkpoint_interval = static_cast<std::size_t>(interval);
    return request;
}
} // namespace madness
