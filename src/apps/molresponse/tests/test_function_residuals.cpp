#include "../solvers/function_residuals.hpp"
#include <cmath>
#include <limits>
#include <iostream>

using namespace madness;

static double gaussian(const coord_3d& r) {
  return std::exp(-(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]));
}
struct Blocks { vector_real_function_3d x_alpha, y_alpha; };

int main(int argc, char** argv) {
  World& world = initialize(argc, argv);
  {
    startup(world, argc, argv, true);
    FunctionDefaults<3>::set_cubic_cell(-12.0, 12.0);
    FunctionDefaults<3>::set_k(8);
    FunctionDefaults<3>::set_thresh(1e-8);
    real_function_3d f = real_factory_3d(world).f(gaussian);
    f.scale(1.0/f.norm2());
    auto four = copy(f); four.scale(4.0);
    auto five = copy(f); five.scale(5.0);
    Blocks input{{f, f}, {f, f}};
    Blocks raw{{four, f}, {f, five}};
    const auto values = molresponse_v3::measure_function_residuals(world, input, raw);
    // A normalized Gaussian gives differences [3,0] in X and [0,4] in Y.
    // These distinguish per-function values, a maximum (4), and quadrature (5).
    const auto& x = values.at("x_alpha");
    const auto& y = values.at("y_alpha");
    MADNESS_CHECK(x.size() == 2 && y.size() == 2 && values.size() == 2);
    MADNESS_CHECK(std::abs(x[0] - 3.0) < 1e-10 && x[1] < 1e-10);
    MADNESS_CHECK(y[0] < 1e-10 && std::abs(y[1] - 4.0) < 1e-10);
    MADNESS_CHECK(std::abs(f.norm2() - 1.0) < 1e-12);
    MADNESS_CHECK(std::abs(four.norm2() - 4.0) < 1e-12);
    MADNESS_CHECK(std::abs(five.norm2() - 5.0) < 1e-12);
    MADNESS_CHECK(function_residual_norms(world, {}).empty());
    bool rejected = false;
    try { function_residual_norms(world, input.x_alpha, vector_real_function_3d{f}); }
    catch (const MadnessException&) { rejected = true; }
    MADNESS_CHECK(rejected);

    auto record = function_residual_record(values);
    MADNESS_CHECK(valid_function_residual_record(record));
    record["blocks"]["x_alpha"][0]["value"] = nullptr;
    MADNESS_CHECK(!valid_function_residual_record(record));
    auto invalid = values;
    invalid["x_alpha"][0] = std::numeric_limits<double>::infinity();
    MADNESS_CHECK(!valid_function_residual_record(function_residual_record(invalid)));
    if (world.rank() == 0) std::cout << "per-function BSH residuals passed\n";
  }
  finalize();
}
