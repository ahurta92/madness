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
  }
  it::Request request;
  for (const auto &name :
       parameters.get<std::vector<std::string>>("convergence.criteria"))
    request.criteria.push_back(it::metric_from_name(name));
  if (request.criteria.empty()) {
    if (parameters.is_user_defined("convergence.criteria"))
      throw std::invalid_argument("empty convergence selection");
    request.criteria = {it::Metric::bsh_residual, it::Metric::density_change,
                        it::Metric::energy_change};
  }
  for (auto metric : {it::Metric::bsh_residual, it::Metric::density_change,
                      it::Metric::energy_change, it::Metric::kinetic_change}) {
    const auto key = std::string("convergence.") + it::metric_name(metric);
    it::Tolerance tolerance;
    const auto absolute = parameters.get<double>(key);
    if (parameters.is_user_defined(key) || absolute >= 0.0)
      tolerance.absolute = absolute;
    tolerance.relative = parameters.get<double>(key + ".relative");
    if ((tolerance.absolute && !it::finite_nonnegative(*tolerance.absolute)) ||
        !it::finite_nonnegative(tolerance.relative))
      throw std::invalid_argument("invalid convergence tolerance: " + key);
    if (tolerance.absolute || tolerance.relative != 0.0)
      request.tolerances[metric] = tolerance;
  }
  const int interval = parameters.get<int>("checkpoint.every_iterations");
  if (interval < 0)
    throw std::invalid_argument("negative checkpoint cadence");
  request.max_updates = parameters.maxiter();
  request.checkpoint_interval = static_cast<std::size_t>(interval);
  return request;
}
} // namespace molresponse_v3
