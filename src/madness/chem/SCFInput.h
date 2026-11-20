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
            key == "restart.start_rung" || key == "maxiter_per_rung")
            p.validate_input_value(key, "0");
    }
    iteration::Request request;
    const auto selection = p.get<std::vector<std::string>>("convergence.criteria");
    for (const auto& name : selection) {
        bool found = false;
        for (auto metric : {iteration::Metric::bsh_residual, iteration::Metric::density_change,
                            iteration::Metric::energy_change, iteration::Metric::kinetic_change})
            if (name == iteration::metric_name(metric)) { request.criteria.push_back(metric); found = true; break; }
        if (!found) throw std::invalid_argument("unknown convergence metric: " + name);
    }
    if (selection.empty()) {
        if (p.is_user_defined("convergence.criteria")) throw std::invalid_argument("empty convergence selection");
        request.criteria = p.freeze_occupied() ? std::vector{iteration::Metric::bsh_residual} :
            std::vector{iteration::Metric::bsh_residual, iteration::Metric::density_change, iteration::Metric::energy_change};
    }
    for (auto metric : {iteration::Metric::bsh_residual, iteration::Metric::density_change,
                        iteration::Metric::energy_change, iteration::Metric::kinetic_change}) {
        const auto key = std::string("convergence.") + iteration::metric_name(metric);
        iteration::Tolerance limit;
        if (p.is_user_defined(key)) limit.absolute = p.get<double>(key);
        limit.relative = p.get<double>(key + ".relative");
        if ((limit.absolute && !iteration::finite_nonnegative(*limit.absolute)) ||
            !iteration::finite_nonnegative(limit.relative)) throw std::invalid_argument("invalid tolerance: " + key);
        if (limit.absolute || limit.relative != 0.0) request.tolerances[metric] = limit;
    }
    const int interval = p.get<int>("checkpoint.every_iterations");
    if (p.maxiter() < 0 || interval < 0) throw std::invalid_argument("negative SCF update budget or checkpoint cadence");
    request.max_updates = static_cast<std::size_t>(p.maxiter());
    request.checkpoint_interval = static_cast<std::size_t>(interval);
    return request;
}
} // namespace madness
