#include "../kernels/response_space_ops.hpp"
#include <cmath>
#include <iostream>

using namespace madness;

static double gaussian(const coord_3d& r) {
  return std::exp(-(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]));
}
static double diffuse(const coord_3d& r) {
  return std::exp(-0.5*(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]));
}

int main(int argc, char** argv) {
  World& world = initialize(argc, argv);
  int result = 0;
  {
    startup(world, argc, argv, true);
    FunctionDefaults<3>::set_cubic_cell(-12.0, 12.0);
    FunctionDefaults<3>::set_k(8);
    FunctionDefaults<3>::set_thresh(1e-8);
    real_function_3d f = real_factory_3d(world).f(gaussian);
    real_function_3d g = real_factory_3d(world).f(diffuse);
    f.scale(1.0/f.norm2());
    g.scale(1.0/g.norm2());
    auto doubled = copy(f); doubled.scale(2.0);
    vector_real_function_3d functions{f, doubled, g};
    auto before = copy(world, functions);
    const auto values = molresponse_v3::rs::kinetic_diagonal(world, functions);
    // Analytic normalized exp(-a*r^2): <T> = 3*a/2; doubling amplitude
    // multiplies the quadratic contribution by four, without normalization.
    const std::vector<double> expected{1.5, 6.0, 0.75};
    MADNESS_CHECK(values.size() == expected.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
      MADNESS_CHECK(std::abs(values[i] - expected[i]) < 2e-6);
      MADNESS_CHECK((functions[i] - before[i]).norm2() < 1e-12);
    }
    MADNESS_CHECK(molresponse_v3::rs::kinetic_diagonal(world, {}).empty());
    if (world.rank() == 0) std::cout << "per-function kinetic diagonals passed\n";
  }
  finalize();
  return result;
}
