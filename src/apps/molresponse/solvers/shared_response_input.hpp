#pragma once
#include <madness/chem/ResponseParameters.hpp>

namespace molresponse_v3 {
inline std::optional<madness::iteration::Request>
shared_response_request(const ResponseParameters &parameters) {
  namespace it = madness::iteration;
  const auto written = parameters.to_json_if_precedence("defined");
  for (const auto &[key, value] : written.items()) {
    (void)value;
    if (key == "dconv" || key == "excited.maxiter" ||
        key.starts_with("accuracy.") || key.starts_with("progress.") ||
        key == "convergence.use_default_criteria" ||
        key == "diagnostic.iterations_per_rung" || key.starts_with("stall.") ||
        key == "seed.start_rung")
      throw std::invalid_argument("unknown input for shared response: " + key);
    if (const auto retired = it::retired_convergence_input(key))
      throw std::invalid_argument(*retired);
  }
  // Only written targets enter the request. Each solver's model supplies its
  // default criteria and targets; a target it cannot judge is reported.
  it::Request request;
  for (auto metric : {it::Metric::bsh_residual, it::Metric::density_change,
                      it::Metric::energy_change, it::Metric::kinetic_change}) {
    const auto key = std::string("convergence.") + it::metric_name(metric);
    if (!parameters.is_user_defined(key))
      continue;
    const double target = parameters.get<double>(key);
    if (!it::finite_nonnegative(target))
      throw std::invalid_argument("invalid convergence target: " + key);
    request.targets[metric] = target;
  }
  if (parameters.is_user_defined("convergence.floor_factor")) {
    request.floor_factor = parameters.get<double>("convergence.floor_factor");
    if (!it::finite_nonnegative(*request.floor_factor))
      throw std::invalid_argument(
          "convergence.floor_factor must be finite and not negative");
  }
  const int interval = parameters.get<int>("checkpoint.every_iterations");
  if (interval < 0)
    throw std::invalid_argument("negative checkpoint cadence");
  request.max_updates = parameters.maxiter();
  request.checkpoint_interval = static_cast<std::size_t>(interval);
  return request;
}
} // namespace molresponse_v3
