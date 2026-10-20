#pragma once

#include "../kernels/response_space_ops.hpp"
#include <madness/chem/FunctionResidualReport.h>

namespace molresponse_v3 {
struct FunctionResidualSample {
  int iteration;
  std::vector<madness::FunctionResidualBlocks> channels;
  std::vector<double> aggregates;
};

// Enumerate only stored blocks: static response has no Y; restricted has no beta.
template<class Storage>
madness::FunctionResidualBlocks measure_function_residuals(madness::World& world,
    const Storage& input, const Storage& raw_update) {
  using madness::function_residual_norms;
  madness::FunctionResidualBlocks values{
      {"x_alpha", function_residual_norms(world, input.x_alpha, raw_update.x_alpha)}};
  if constexpr (rs::detail::has_y_alpha<Storage>::value)
    values["y_alpha"] = function_residual_norms(world, input.y_alpha, raw_update.y_alpha);
  if constexpr (rs::detail::has_x_beta<Storage>::value)
    values["x_beta"] = function_residual_norms(world, input.x_beta, raw_update.x_beta);
  if constexpr (rs::detail::has_y_beta<Storage>::value)
    values["y_beta"] = function_residual_norms(world, input.y_beta, raw_update.y_beta);
  return values;
}
} // namespace molresponse_v3
