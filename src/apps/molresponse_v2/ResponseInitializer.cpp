#include "ResponseInitializer.hpp"

ResponseVector initialize_guess_vector(World &world, const GroundStateData &gs,
                                       const LinearResponsePoint &pt) {
  size_t num_orbitals = gs.getNumOrbitals();
  bool is_static = pt.is_static();
  bool is_unrestricted = !gs.isSpinRestricted();

  vector_real_function_3d Vp =
      zero_functions_compressed<double, 3>(world, num_orbitals);
  // vector_real_function_3d Vp = perturbation_vector(world, gs, pt);

  if (is_static && !is_unrestricted) {
    StaticRestrictedResponse response(num_orbitals);
    for (size_t i = 0; i < num_orbitals; ++i)
      response.x_alpha[i] = Vp[i];
    response.flatten();
    return response;

  } else if (!is_static && !is_unrestricted) {
    DynamicRestrictedResponse response(num_orbitals);
    for (size_t i = 0; i < num_orbitals; ++i) {
      response.x_alpha[i] = Vp[i];
      response.y_alpha[i] = Vp[i];
    }
    response.flatten();
    return response;

  } else if (is_static && is_unrestricted) {
    StaticUnrestrictedResponse response(num_orbitals);
    for (size_t i = 0; i < num_orbitals; ++i) {
      response.x_alpha[i] = Vp[i];
      response.x_beta[i] = Vp[i];
    }
    response.flatten();
    return response;

  } else if (!is_static && is_unrestricted) {
    DynamicUnrestrictedResponse response(num_orbitals);
    for (size_t i = 0; i < num_orbitals; ++i) {
      response.x_alpha[i] = Vp[i];
      response.x_beta[i] = Vp[i];
    }
    response.flatten();
    return response;
  }

  throw std::runtime_error("Unknown response vector configuration.");
}

