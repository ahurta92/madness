#ifndef MADNESS_CHEM_KINETICENERGY_H
#define MADNESS_CHEM_KINETICENERGY_H

#include <madness/mra/mra.h>
#include <madness/tensor/distributed_matrix.h>
#include <array>
#include <memory>
#include <vector>

namespace madness {

using KineticGradients3D = std::array<vector_real_function_3d, 3>;

namespace detail {
// The common kinetic definition: T_ij = 1/2 sum_d <D_d f_i | D_d f_j>.
// Only the contraction layout differs (distributed matrix or its diagonal).
// No occupations, spin factors, normalization, or convergence policy belongs here.
template<class Result, class Contract>
Result contract_kinetic_gradients(Result result, const Contract& contract) {
    for (int axis = 0; axis < 3; ++axis) result += contract(axis);
    result *= 0.5;
    return result;
}
} // namespace detail

// Collective. Accept already prepared gradients so SCF retains its configured
// derivatives, halo staging, batching, and timing. No replicated N-by-N matrix.
inline DistributedMatrix<double> kinetic_energy_matrix_from_gradients(
        World& world, const KineticGradients3D& gradients) {
    const auto n = gradients[0].size();
    MADNESS_CHECK(gradients[1].size() == n && gradients[2].size() == n);
    auto matrix = column_distributed_matrix<double>(world, static_cast<long>(n), static_cast<long>(n));
    const auto& distribution = matrix.distribution();
    return detail::contract_kinetic_gradients(matrix,
        [&distribution, &gradients](int axis) {
            return matrix_inner(distribution, gradients[axis], gradients[axis], true);
        });
}

// Collective. Return one unnormalized scalar per function; compute no
// off-diagonal integrals. Copy the handles deeply before derivative operations.
inline std::vector<double> kinetic_energy_diagonal(World& world,
        const vector_real_function_3d& functions,
        const std::vector<std::shared_ptr<real_derivative_3d>>& derivatives) {
    if (functions.empty()) return {};
    MADNESS_CHECK(derivatives.size() == 3);
    auto input = copy(world, functions);
    // Keep only one derivative forest at a time. norm2s works on the
    // reconstructed coefficients without compressing a temporary forest.
    const auto diagonal = detail::contract_kinetic_gradients(
        Tensor<double>(static_cast<long>(functions.size())),
        [&world, &input, &derivatives](int axis) {
            MADNESS_CHECK(derivatives[axis]);
            auto gradient = apply(world, *derivatives[axis], input);
            const auto norms = norm2s(world, gradient);
            Tensor<double> integrals(static_cast<long>(norms.size()));
            for (std::size_t i = 0; i < norms.size(); ++i)
                integrals(static_cast<long>(i)) = norms[i] * norms[i];
            return integrals;
        });
    std::vector<double> values(functions.size());
    for (std::size_t i = 0; i < values.size(); ++i) values[i] = diagonal(static_cast<long>(i));
    return values;
}

// Default derivative discretization used by the driven-response kernels.
inline std::vector<double> kinetic_energy_diagonal(
        World& world, const vector_real_function_3d& functions) {
    return kinetic_energy_diagonal(world, functions, gradient_operator<double, 3>(world));
}

} // namespace madness
#endif
