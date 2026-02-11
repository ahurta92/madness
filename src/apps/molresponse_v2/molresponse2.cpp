#include <madness/chem/Drivers.hpp>
#include <madness/chem/ParameterManager.hpp>
#include <madness/chem/WorkflowBuilders.hpp>
#include <madness/misc/info.h>
#include <madness/world/worldmem.h>
#include <madness_exception.h>

using namespace madness;

namespace {

void print_help() {
  print("Usage: molresponse2 [options] [input_file]");
  print("Compatibility wrapper for response workflow.");
  print("");
  print("This executable now runs the same code path as:");
  print("  madqc --wf=response [options] [input_file]");
  print("");
  print("Common options:");
  print("  --help");
  print("  --prefix=<name>");
  print("  --molecule=\"...\"");
  print("  --dft=\"...\"");
  print("  --response=\"...\"");
}

} // namespace

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  commandlineparser parser(argc, argv);

  if (parser.key_exists("help")) {
    if (world.rank() == 0)
      print_help();
    finalize();
    return 0;
  }

  if (world.rank() == 0) {
    print_header1("MOLRESPONSE2 -- compatibility wrapper");
    print("Note: this executable is deprecated; use madqc --wf=response.");
  }

  try {
    startup(world, argc, argv, true);
    if (world.rank() == 0)
      print(info::print_revision_information());

    std::string requested_workflow = "response";
    if (parser.key_exists("workflow"))
      requested_workflow = parser.value("workflow");
    else if (parser.key_exists("wf"))
      requested_workflow = parser.value("wf");
    else if (parser.key_exists("w"))
      requested_workflow = parser.value("w");

    if (requested_workflow != "response") {
      std::string msg =
          "molresponse2 only supports workflow=response. "
          "Use madqc for other workflows.";
      MADNESS_EXCEPTION(msg.c_str(), 1);
    }

    Params pm(world, parser);
    qcapp::Workflow wf;
    workflow_builders::add_response_workflow_drivers(world, pm, wf);
    wf.run(pm.prefix());
  } catch (const MadnessException &e) {
    if (world.rank() == 0) {
      print_header2("caught a MADNESS exception in molresponse2");
      print(e.what(), e.filename, e.msg, e.line);
    }
  } catch (const std::exception &e) {
    if (world.rank() == 0) {
      print_header2("caught an exception in molresponse2");
      print(e.what());
    }
  }

  world.gop.fence();
  world.gop.fence();
  print_stats(world);
  finalize();
  return 0;
}
