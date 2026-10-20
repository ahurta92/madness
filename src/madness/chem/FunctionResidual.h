#pragma once

#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <map>
#include <string>
#include <vector>

namespace madness {
// Block names and function indices belong to the application's working basis.
using FunctionResidualBlocks = std::map<std::string, std::vector<double>>;

// Shared definition: raw L2 norm of each fixed-point difference. No weights,
// normalization, or reduction across functions are part of this measurement.
inline std::vector<double> function_residual_norms(
    World& world, const vector_real_function_3d& residuals) {
    return norm2s(world, residuals);
}

inline std::vector<double> function_residual_norms(World& world,
    const vector_real_function_3d& input, const vector_real_function_3d& raw_update) {
    MADNESS_CHECK_THROW(input.size() == raw_update.size(),
                       "per-function residual requires matching function counts");
    // Function handles are shallow; preserve the caller's trees as well as values.
    auto before = copy(world, input);
    auto after = copy(world, raw_update);
    return function_residual_norms(world, sub(world, before, after));
}
} // namespace madness
