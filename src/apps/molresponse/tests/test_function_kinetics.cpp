#include "../kernels/response_space_ops.hpp"
#include <madness/chem/KineticEnergy.h>
#include <madness/chem/SCF.h>
#include <cmath>
#include <numeric>
#include "../solvers/function_kinetics.hpp"
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
    CalculationParameters parameters;
    parameters.set_user_defined_value("l", 12.0);
    parameters.set_user_defined_value("k", 8);
    parameters.set_user_defined_value("print_level", 0);
    parameters.set_user_defined_value<std::string>("xc", "hf");
    Molecule molecule;
    molecule.add_atom(0.0, 0.0, 0.0, 2.0, 2);
    SCF scf(world, parameters, molecule);
    scf.set_protocol<3>(world, 1e-8);
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
    // A unitary rotation redistributes diagonal entries but preserves trace.
    auto rotated = copy(world, vector_real_function_3d{f, g});
    rotated[0] = (f + g) * (1.0/std::sqrt(2.0));
    rotated[1] = (f - g) * (1.0/std::sqrt(2.0));
    const auto rotated_diagonal = molresponse_v3::rs::kinetic_diagonal(world, rotated);
    MADNESS_CHECK(std::abs(rotated_diagonal[0] + rotated_diagonal[1] - values[0] - values[2]) < 2e-6);
    MADNESS_CHECK(std::abs(rotated_diagonal[0] - values[0]) > 1e-2);
    auto record = molresponse_v3::function_kinetic_record({{"x_alpha", values}}, 0, "rotation-test", 1e-8, 8,
        true, "saved_response_before_archive");
    MADNESS_CHECK(molresponse_v3::valid_function_kinetic_record(record));
    MADNESS_CHECK(std::abs(record.at("trace").at("value").get<double>() - 8.25) < 4e-6);
    record["trace"]["value"] = 4.5;
    MADNESS_CHECK(!molresponse_v3::valid_function_kinetic_record(record));
    // Compare with moldft's real distributed matrix path, including its
    // off-diagonal entries. A wrong factor or an unintended normalization
    // must also fail the independent Gaussian expectations above.
    Tensor<double> matrix(3, 3);
    auto matrix_functions = copy(world, functions);
    scf.kinetic_energy_matrix(world, matrix_functions).copy_to_replicated(matrix);
    for (std::size_t i = 0; i < values.size(); ++i)
      MADNESS_CHECK(std::abs(values[i] - matrix(static_cast<long>(i), static_cast<long>(i))) < 1e-10);
    MADNESS_CHECK(std::abs(matrix(0, 1) - 3.0) < 4e-6);

    // The actual SCF Fock builder performs the separate occupation-weighted
    // sum. Occupations [1, 1/2, 0] give 1.5 + 3.0 = 4.5, per spin.
    Tensor<double> occupations(3);
    occupations(0L) = 1.0; occupations(1L) = 0.5; occupations(2L) = 0.0;
    auto potential = zero_functions<double, 3>(world, 3);
    double total = 0.0;
    scf.make_fock_matrix(world, matrix_functions, potential, occupations, total);
    MADNESS_CHECK(std::abs(total - 4.5) < 4e-6);

    // Definition and derivative discretization are distinct. The shared
    // helper must honor SCF's configured operators, including nondefault ones.
    for (const std::string derivative : {"ble", "bspline"}) {
      scf.param.set_user_defined_value("deriv", derivative);
      scf.set_protocol<3>(world, 1e-8);
      const auto diagonal = kinetic_energy_diagonal(world, functions, scf.gradop);
      scf.kinetic_energy_matrix(world, matrix_functions).copy_to_replicated(matrix);
      for (std::size_t i = 0; i < values.size(); ++i) {
        MADNESS_CHECK(std::abs(diagonal[i] - matrix(static_cast<long>(i), static_cast<long>(i))) < 1e-10);
        MADNESS_CHECK((functions[i] - before[i]).norm2() < 1e-12);
      }
    }
    MADNESS_CHECK(molresponse_v3::rs::kinetic_diagonal(world, {}).empty());
    if (world.rank() == 0) {
      std::cout << "KINETIC_TRACE_EVIDENCE " << nlohmann::json{
          {"units", "hartree"}, {"diagonals", values}, {"unweighted_trace", values[0]+values[1]+values[2]},
          {"occupation_weighted_energy_per_spin", total}, {"rotated_diagonals", rotated_diagonal},
          {"rotation_trace_error", rotated_diagonal[0]+rotated_diagonal[1]-values[0]-values[2]}}
          .dump() << '\n';
      std::cout << "per-function kinetic diagonals passed\n";
    }
  }
  finalize();
  return result;
}
