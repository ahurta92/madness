#pragma once

#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>

namespace madness {
struct DensityChanges {
    double charge = 0.0;
    double spin = 0.0;
};

// The application supplies its physical spin densities, including occupations.
// These are norms of field differences, not sums/differences of spin norms.
// No normalization, tolerance or convergence decision belongs to this helper.
inline DensityChanges density_change_norms(World& world,
    const real_function_3d& alpha, const real_function_3d& beta,
    const real_function_3d& previous_alpha, const real_function_3d& previous_beta) {
    // Function handles are shallow: preserve caller values and tree states.
    auto current = copy(world, vector_real_function_3d{alpha, beta});
    auto previous = copy(world, vector_real_function_3d{previous_alpha, previous_beta});
    auto delta = sub(world, current, previous);
    return {(delta[0] + delta[1]).norm2(), (delta[0] - delta[1]).norm2()};
}
} // namespace madness
