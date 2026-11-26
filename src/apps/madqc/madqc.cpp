/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
*/

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) &&                   \
    defined(HAVE_UNISTD_H)

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

// static inline int file_exists(const char *inpname) {
//     struct stat buffer;
//     int rc = stat(inpname, &buffer);
//     return (rc == 0);
// }

#endif

#include <madness/misc/info.h>
#include <madness/world/worldmem.h>

#include <fstream>
#include <iostream>

#include "WorkflowDriver.hpp"
#include <madness/chem/RootCoordination.h>
#include <madness/chem/SCFRun.h>
#include <madness_exception.h>

using namespace madness;

/// TODO:
///  - restart options
///  - nemo projections
///  - numerical parameters

void help(const std::string &wf) {
  print("Usage: madqc [options] [input_file]");
  print("  --prefix=<directory>        : workflow output root");
  print("  --continue                  : require an existing shared SCF "
        "checkpoint");
  print("\nOptions:");
  print("  --help=<workflow>           : show this help message");
  print("  --print_parameters=<group>  : print all parameters and exit");
  print("  --restart_info=<directory>  : describe a shared SCF checkpoint and "
        "exit");
  print("  --workflow=<name>           : specify the workflow to run (default: "
        "scf)");
  print("  --optimize                  : optimize the geometry of the "
        "workflow's reference (scf|nemo)");
  print("\nAvailable workflows: " +
        workflow_builders::runnable_workflow_list());
  print("Parameter groups (for --print_parameters): dft, nemo, response, cc2, "
        "cis, oep, pcm, optimization, geometry");
  print("");
  if (wf == "scf") {
    print("madqc --wf=scf");
    print("molecular Hartree-Fock and DFT code");
    print("\nexamples: ");
    print("madqc --wf=scf --geometry=h2o.xyz");
    print("madqc --wf=scf --geometry=he --dft=\"k=8\",xc=lda");
    print("madqc --wf=scf --geometry=h2o --dft=\"k=8; "
          "convergence.energy_change=1.e-4; L=30\"");
  } else if (wf == "nemo") {
    print("madqc --wf=nemo");
    print("molecular Hartree-Fock and DFT code with regularized orbitals");
    print("\nexamples: ");
    print("madqc --wf=nemo --geometry=h2o ");
  } else if (wf == "response") {
    print("madqc --wf=response");
    print("Response theory for DFT and Hartree-Fock");
    print("(runs the ground-state SCF first, then linear/quadratic response)");
    print("\nexamples: ");
    print("madqc --wf=response h2o_response.in");
    print("see  madqc --print_parameters=response  for all response knobs");
  } else if (wf == "mp2" or wf == "cc2") {
    print("madqc --wf=cc2");
    print("MP2/CC2 code for correlated wave functions");
    print("\nexamples: ");
    print("madqc --wf=cc2 --geometry=h2o --dft=\"k=5; econv=1.e-4; L=30\" "
          "--cc2=\"freeze=1\"");
  } else if (wf == "cis") {
    print("madqc --wf=cis");
    print("CIS code for excited states");
  } else if (wf == "oep") {
    print("madqc --wf=oep");
    print("Optimized effective potential code for DFT");
    print("\nexamples: ");
    print("madqc --wf=oep --geometry=h2o");
  } else if (wf == "optimize") {
    print("madqc --optimize --wf=<scf|nemo>");
    print("Geometry optimization: --wf names the reference method, --optimize "
          "asks for its");
    print("geometry to be optimized. (`--wf=optimize` is not a workflow.)");
    print("\nexamples: ");
    print("madqc --optimize --wf=scf --geometry=h2o");
    print("madqc --optimize --wf=nemo --geometry=lih");
    print(
        "madqc --optimize --wf=scf --geometry=h2o --optimization=\"gtol=1.e-4; "
        "maxiter=10\"");
    print("\nsee  madqc --print_parameters=optimization  for all knobs.");
    print(
        "the in-SCF form (--dft=\"gopt=1\") has been removed; its keyvals are");
    print("retired and error with a pointer to the `optimization` group.");
  }
}

void print_parameters(World &world, const commandlineparser &parser,
                      const std::string &group) {
  Params pm;
  if (group.empty()) {
    print("please specify a data group to print parameters for");
    print("\n  --print_parameters=<group>  : print all parameters and exit");
    print("\nAvailable data groups: dft, nemo, response, cc2, cis, oep, pcm, "
          "optimization, geometry");
  } else if (group == "dft" || group == "scf") {
    print("Available parameters for data group: dft");
    auto parameters = pm.get<CalculationParameters>();
    parameters.enable_scf_inputs();
    parameters.print();
  } else if (group == "nemo") {
    print("Available parameters for data group: nemo");
    pm.get<Nemo::NemoCalculationParameters>().print();
  } else if (group == "response") {
    print("Available parameters for data group: response");
    pm.get<ResponseParameters>().print();
  } else if (group == "cc2") {
    print("Available parameters for data group: cc2");
    pm.get<CCParameters>().print();
  } else if (group == "cis") {
    print("Available parameters for data group: cis");
    pm.get<TDHFParameters>().print();
  } else if (group == "oep") {
    print("Available parameters for data group: oep");
    pm.get<OEP_Parameters>().print();
  } else if (group == "optimization") {
    print("Available parameters for data group: optimization");
    pm.get<OptimizationParameters>().print(OptimizationParameters::tag, "end");
  } else if (group == "pcm") {
    print("Available parameters for data group: pcm");
    pm.get<PCMParameters>().print();
  } else if (group == "geometry") {
    Molecule::GeometryParameters geometryparam;
    geometryparam.print("geometry", "end");
  } else {
    std::string msg = "Unknown data group: " + group +
                      "\nAvailable data group are: dft, nemo, response, cc2, "
                      "cis, oep, pcm, optimization, geometry\n";
    print(msg);
  }
}

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  int rc = 0;
  bool solver_started = false;
  std::optional<madqc::RunSummaryTarget> summary;
  if (world.rank() == 0)
    print_header1("MADQC -- Multiresolution Quantum Chemistry Code");
  try {
    commandlineparser parser(argc, argv);
    if (parser.key_exists("help")) {
      help(parser.value("help"));
    } else if (parser.key_exists("print_parameters")) {
      print_parameters(world, parser, parser.value("print_parameters"));
    } else if (parser.key_exists("restart_info")) {
      std::string path = parser.value_raw("restart_info");
      if (path.empty() || path == "restart_info")
        path = ".";
      CalculationParameters inspection;
      inspect_scf(world, std::filesystem::absolute(path), inspection);
    } else {
      startup(world, argc, argv, true);
      solver_started = true;
      if (world.rank() == 0) {
        print(info::print_revision_information());
        print("    git source commit ...     ", info::git_commit());
      }
      madqc::run_workflow(world, parser, &summary);
    }
  } catch (const MadnessException &e) {
    rc = 1;
    if (world.rank() == 0)
      print(e.what(), e.filename, e.msg, e.line);
  } catch (const std::exception &e) {
    rc = 1;
    if (world.rank() == 0)
      print(e.what());
  }
  world.gop.fence();
  if (solver_started)
    print_stats(world);
  // Last in the log: the text of summary.txt, from results.json.
  madqc::print_run_summary(world, summary);
  finalize();
  return rc;
}
