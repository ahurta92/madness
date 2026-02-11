#pragma once

#include <Drivers.hpp>
#include <MoldftLib.hpp>
#include <MolresponseLib.hpp>
#include <ParameterManager.hpp>

namespace madness::workflow_builders {

inline void add_response_workflow_drivers(World &world, Params &pm,
                                          qcapp::Workflow &wf) {
  pm.get<CalculationParameters>().set_derived_value("save", true);

  auto reference = std::make_shared<SCFApplication<moldft_lib>>(world, pm);
  wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(reference));
  wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(
      std::make_unique<ResponseApplication<molresponse_lib>>(world, pm,
                                                             reference->calc())));
}

} // namespace madness::workflow_builders
