
#include <madness/constants.h>
#include <madness/mra/mra.h>
#include <madness/mra/nonlinsol.h>  // The kain solver
#include <madness/mra/operator.h>
#include <math.h>
#include <stdio.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "../chem/SCFOperators.h"
#include "../chem/molecule.h"
#include "../chem/xcfunctional.h"
#include "molresponse/TDDFT.h"
#include "molresponse/basic_operators.h"
#include "molresponse/ground_parameters.h"
#include "molresponse/load_balance.h"
#include "molresponse/property.h"
#include "molresponse/response_functions.h"
#include "molresponse/response_parameters.h"
#include "molresponse/response_potential.h"
#include "molresponse/timer.h"
#include "molresponse/x_space.h"

X_space TDDFT::Compute_Theta_X(World& world,
                               X_space& Chi,
                               XCOperator<double, 3> xc,
                               std::string calc_type) {
  bool compute_Y = calc_type.compare("full") == 0;
  X_space Theta_X = X_space(world, Chi.num_states(), Chi.num_orbitals());
  // compute
  X_space V0X = compute_V0X(world, Chi, xc, compute_Y);

  V0X.truncate();
  if (r_params.print_level() >= 20) {  // LEGACY_PATCH: collective first, print on rank 0
    auto xv0x = inner(Chi, V0X);  // all ranks must participate (collective)
    if (world.rank() == 0) { print("---------------Theta ----------------"); print("<X|V0|X>"); print(xv0x); }
  }

  X_space E0X(world, Chi.num_states(), Chi.num_orbitals());
  if (r_params.localize().compare("canon") == 0) {
    E0X = Chi.copy();
    E0X.truncate();
    E0X.X = E0X.X * ham_no_diag;
    if (compute_Y) {
      E0X.Y = E0X.Y * ham_no_diag;
    }

    E0X.truncate();
  }

  if (r_params.print_level() >= 20) {
    // LEGACY_PATCH: inner() is collective; compute on all ranks, print on rank 0
    auto xe0x = inner(Chi, E0X);
    if (world.rank() == 0) { print("<X|(E0-diag(E0)|X>"); print(xe0x); }
  }

  X_space gamma;
  // compute
  if (calc_type.compare("full") == 0) {
    gamma = compute_gamma_full(world, Chi, xc);
  } else if (calc_type.compare("static") == 0) {
    gamma = compute_gamma_static(world, Chi, xc);
  } else {
    gamma = compute_gamma_tda(world, Chi, xc);
  }

  Theta_X = (V0X - E0X) + gamma;
  Theta_X.truncate();

  if (r_params.print_level() >= 20) {
    // LEGACY_PATCH: inner() is collective; compute on all ranks, print on rank 0
    auto xthx = inner(Chi, Theta_X);
    if (world.rank() == 0) { print("<X|Theta|X>"); print(xthx); }
  }

  return Theta_X;
}
