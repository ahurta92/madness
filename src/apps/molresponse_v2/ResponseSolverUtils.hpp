#pragma once
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <madness/chem/SCF.h>
#include <madness/mra/nonlinsol.h>

#include <madness/external/nlohmann_json/json.hpp>

#include "ResponseVector.hpp"

using json = nlohmann::json;

namespace ResponseSolverUtils {

using namespace madness;

inline void print_iteration_table_border() {
  std::cout
      << "+------+--------------+--------------+--------------+--------------+"
         "--------------+\n";
}

inline void print_iteration_line(int iter, double residual, double deltaE, double xVp, double density_target, double x_residual_target) {
  auto format_sci = [](double value) -> std::string {
    std::ostringstream os;
    os << std::scientific << std::setprecision(3) << std::setw(12) << value;
    return os.str();
  };

  if (iter == 0) {
    std::cout << "ITERATION_TABLE_COLUMNS "
              << "iter,residual,drho,x_vp,target_residual,target_drho\n";
    print_iteration_table_border();
    std::cout << "| " << std::setw(4) << "iter"
              << " | " << std::setw(12) << "residual"
              << " | " << std::setw(12) << "drho"
              << " | " << std::setw(12) << "x_vp"
              << " | " << std::setw(12) << "target_res"
              << " | " << std::setw(12) << "target_drho"
              << " |\n";
    print_iteration_table_border();
  }

  std::cout << "| " << std::setw(4) << iter
            << " | " << format_sci(residual)
            << " | " << format_sci(deltaE)
            << " | " << format_sci(xVp)
            << " | " << format_sci(x_residual_target)
            << " | " << format_sci(density_target)
            << " |\n";

  // Repeat the border periodically so long traces stay readable.
  if ((iter + 1) % 10 == 0) {
    print_iteration_table_border();
  }
}

inline std::vector<poperatorT> make_bsh_operators_response(
    World &world, const double shift, const double omega,
    const Tensor<double> &ground_energies, const double &lo) {
  double tol = FunctionDefaults<3>::get_thresh();
  // Sizes inferred from ground and omega
  size_t num_orbitals = ground_energies.size();  // number of orbitals
  std::vector<poperatorT> ops(num_orbitals);
  // Run over occupied components
  int p = 0;
  std::for_each(ops.begin(), ops.end(), [&](auto &operator_p) {
    double mu = sqrt(-2.0 * (ground_energies(p++) + omega + shift));
    operator_p = poperatorT(BSHOperatorPtr3D(world, mu, lo, tol));
  });
  return ops;
  // End timer
}

inline double inner(World &world, const vector_real_function_3d &x, const vector_real_function_3d &y) {
  double result = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    result += x[i].inner(y[i]);
  }
  return result;
}

inline int infer_function_k(const vector_real_function_3d &functions) {
  for (const auto &f : functions) {
    return f.k();
  }
  return FunctionDefaults<3>::get_k();
}

inline int infer_state_bundle_k(const std::vector<vector_real_function_3d> &states) {
  for (const auto &state : states) {
    if (!state.empty()) {
      return infer_function_k(state);
    }
  }
  return FunctionDefaults<3>::get_k();
}

inline bool align_function_vector_protocol(World &world,
                                           vector_real_function_3d &functions,
                                           int target_k,
                                           double target_thresh) {
  if (functions.empty()) {
    return false;
  }
  const int source_k = infer_function_k(functions);
  if (source_k != target_k) {
    reconstruct(world, functions);
    for (auto &f : functions) {
      f = project(f, target_k, target_thresh, true);
    }
    truncate(world, functions, target_thresh, true);
    return true;
  }
  truncate(world, functions, target_thresh, true);
  return false;
}

inline bool align_state_bundle_protocol(
    World &world, std::vector<vector_real_function_3d> &states, int target_k,
    double target_thresh) {
  bool projected = false;
  for (auto &state : states) {
    projected = align_function_vector_protocol(world, state, target_k,
                                               target_thresh) || projected;
  }
  return projected;
}

inline bool align_response_vector_protocol(World &world, ResponseVector &response,
                                          int target_k, double target_thresh) {
  return std::visit(
      [&](auto &typed_response) {
        const bool projected = align_function_vector_protocol(
            world, typed_response.flat, target_k, target_thresh);
        typed_response.sync();
        return projected;
      },
      response);
}

inline bool align_response_bundle_protocol(World &world,
                                          std::vector<ResponseVector> &responses,
                                          int target_k,
                                          double target_thresh) {
  bool projected = false;
  for (auto &response : responses) {
    projected = align_response_vector_protocol(world, response, target_k,
                                               target_thresh) || projected;
  }
  return projected;
}

inline void do_step_restriction(World &world, const vecfuncT &x, vecfuncT &x_new, const double &anorm, const std::string &spin, const double &maxrotn) {
  // int nres = 0;
  if (anorm > maxrotn) {
    if (world.rank() == 0)
      print("STEP_RESTRICTION_APPLIED norm_change=", anorm,
            " max_rotation=", maxrotn);
    double s = maxrotn / anorm;
    gaxpy(s, x_new, 1.0 - s, x, true);
  }

  world.gop.fence();
  if (world.rank() == 0)
    print("STEP_RESTRICTION_NORM spin=", spin, " norm_change=", anorm);
}
}  // namespace ResponseSolverUtils
