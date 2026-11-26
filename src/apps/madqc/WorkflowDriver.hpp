#pragma once

#include "WorkflowInput.h"
#include <Drivers.hpp>
#include <ParameterManager.hpp>
#include <WorkflowBuilders.hpp>
#include <apps/molresponse/madqc_adapter.hpp>
#include <madness/chem/IterationOutput.h>
#include <madness/chem/ResultsSummary.hpp>
#include <madness/chem/RunFiles.h>
#include <madness/chem/SCFRun.h>
#include <madness/chem/VizManifest.hpp>
#include <optional>

namespace madqc {

// Never silently start fresh beside, or read, a run folder an older version
// wrote. Validation finishes collectively before creating any run output.
inline void require_shared_run_layout(madness::World& world,
                                      const std::filesystem::path& root) {
  std::string error;
  if (world.rank() == 0) {
    try {
      if (std::filesystem::exists(root / ".madqc-layout-migration.json") ||
          std::filesystem::exists(root / "task_0" / "moldft"))
        error = "the run folder " + root.string() +
                " was written by an older version of madqc; use a new output prefix";
    } catch (const std::exception& e) { error = e.what(); }
  }
  world.gop.broadcast_serializable(error, 0);
  if (!error.empty()) throw std::runtime_error(error);
  madness::iteration::require_current_run_folder(world, root);
}


inline std::shared_ptr<SCFApplication<moldft_lib>>
make_scf_reference(madness::World &world, const Params &params) {
  auto reference = std::make_shared<SCFApplication<moldft_lib>>(world, params);
  const auto directory = molresponse_v3::effective_dalton_dir(params);
  if (!directory.empty())
    reference->set_pre_run_hook(
        [directory](madness::World &w, const Params &p,
                    const std::filesystem::path &output) {
          molresponse_v3::seed_gs_from_dalton_dir(w, p, output, directory);
        });
  return reference;
}

// One application-level driver. Workflow code owns dependencies; each engine
// retains its numerical state and solver implementation.
/// The run folder whose summary madqc prints at the end of the log, and the
/// verbosity that decides whether it prints. summary.txt itself is written
/// with results.json at every verbosity.
struct RunSummaryTarget {
  std::filesystem::path directory;
  madness::iteration::Verbosity verbosity = madness::iteration::Verbosity::silent;
};

/// Rank zero. Print the summary of the run folder's results.json, the same
/// text as its summary.txt, when the verbosity asks for it.
inline void print_run_summary(madness::World &world,
                              const std::optional<RunSummaryTarget> &target) {
  namespace it = madness::iteration;
  if (world.rank() != 0 || !target ||
      !it::prints(target->verbosity, it::needed::run_summary))
    return;
  const auto path = target->directory / "results.json";
  if (!std::filesystem::exists(path))
    return;
  try {
    std::cout << "\n"
              << it::results_summary(madness::checkpoint::detail::read_json(path))
              << std::flush;
  } catch (const std::exception &e) {
    madness::print("summary not printed:", e.what());
  }
}

/// `summary`, when given, is set to the run folder of a workflow that writes
/// results.json, before any work starts, so a failed run still prints it.
inline void run_workflow(madness::World &world,
                         const madness::commandlineparser &parser,
                         std::optional<RunSummaryTarget> *summary = nullptr) {
  using namespace madness;
  validate_workflow_input(world, parser);
  // Create workflow
  qcapp::Workflow wf;
  std::string user_workflow = "scf";
  if (parser.key_exists("workflow"))
    user_workflow = parser.value("workflow");
  else if (parser.key_exists("wf"))
    user_workflow = parser.value("wf");
  else if (parser.key_exists("w"))
    user_workflow = parser.value("w");

  if (world.rank() == 0) {
    print("input file  :", parser.value("input"));
    print("workflow    :", user_workflow);
    print("");
  }

  // read in all parameters from the input file and the command line
  // logic and interdependent parameter follow later
  Params pm(world, parser);
  const auto &scf_parameters = pm.get<CalculationParameters>();
  const bool shared_scf = user_workflow == "scf" || user_workflow == "response";
  const bool shared_output = shared_scf && !parser.key_exists("optimize");
  const auto run_directory = std::filesystem::absolute(pm.prefix());
  const auto scf_directory = shared_output ? run_directory
      : run_directory / "task_0" / "moldft";
  if (shared_output) require_shared_run_layout(world, run_directory);
  if (shared_output && summary)
    *summary = RunSummaryTarget{
        run_directory,
        madness::iteration::verbosity(
            user_workflow == "response"
                ? pm.get<ResponseParameters>().print_level()
                : scf_parameters.print_level())};
  if (shared_scf && pm.get<IOParameters>().hdf5())
    throw std::runtime_error("shared SCF checkpoints require native archives");
  if (parser.key_exists("continue")) {
    if (!shared_scf || (scf_parameters.restart() != "auto" &&
                        scf_parameters.restart() != "iterate"))
      throw std::runtime_error(
          "--continue requires an SCF workflow with restart auto or iterate");
    require_scf_restart(world, scf_directory);
  }
  if (shared_scf && scf_parameters.restart() == "read_only") {
    if (user_workflow != "scf" || parser.key_exists("optimize") ||
        scf_parameters.derivatives() || scf_parameters.dipole())
      throw std::runtime_error("read_only supports inspection and plotting, "
                               "not workflow computation");
    inspect_scf(world, scf_directory, scf_parameters);
  } else {

    // A relative `dalton.dir` names a directory beside the deck, but every
    // engine runs from its own diagnostics directory, where the
    // same relative path finds nothing. Resolve it once, against the launch
    // directory, before any stage changes directory.
    auto absolute_dalton_dir = [](auto &params) {
      const auto d = params.template get<std::string>("dalton.dir");
      if (!d.empty() && std::filesystem::path(d).is_relative())
        params.template set_user_defined_value<std::string>(
            "dalton.dir",
            std::filesystem::absolute(d).lexically_normal().string());
    };
    absolute_dalton_dir(pm.get<IOParameters>());
    absolute_dalton_dir(pm.get<ResponseParameters>());

    // Deck-level `io` block (roadmap change 5): run-wide restart-backend
    // selection. `io backend hdf5` turns the HDF5 opt-in on for every task
    // that supports it; the response-only `response.hdf5 true` keeps
    // working as an alias (applied here and re-checked in the response
    // adapter). The effective config is stamped into calc_info for ALL
    // tasks via Workflow::set_provenance below.
    // HDF5 restart I/O is a molresponse-only path today, so only honor the
    // deck flag (and stamp the provenance below) for the response workflow —
    // a stray `response.hdf5`/`io backend hdf5` in a cc2/cis/oep/moldft deck
    // must NOT flip the global flag or mis-stamp backend=hdf5 (review LOW).
    // `--optimize` turns the workflow's reference into an optimization task
    // (WorkflowBuilders::add_optimize_workflow_drivers).
    const bool optimize_geometry = parser.key_exists("optimize");
    const bool is_response_wf =
        workflow_builders::workflow_kind_from_name(user_workflow) ==
        workflow_builders::WorkflowKind::Response;
    const bool deck_io_hdf5 =
        pm.get<IOParameters>().hdf5() || pm.get<ResponseParameters>().hdf5();
    if (deck_io_hdf5 && is_response_wf) {
#ifdef MADNESS_HAS_HDF5
      molresponse_v3::set_hdf5_io_enabled(true);
#else
      throw std::runtime_error(
          "io backend hdf5 (or response.hdf5) requested but this build has "
          "no HDF5 support — configure with -DMADNESS_ENABLE_HDF5=ON");
#endif
    } else if (deck_io_hdf5 && !is_response_wf && world.rank() == 0) {
      print("NOTE: `io backend hdf5` / `response.hdf5` only affects the "
            "response workflow's MRA restart I/O — ignored for --wf=",
            user_workflow);
    }

    // The response workflow is built HERE in the app — MADchem's
    // WorkflowBuilders cannot reference v3 (circular lib dependency), and the
    // v2 engine was removed (M1 decoupling Stage 2). molresponse_v3 is THE
    // response engine; `engine v2` in old decks fails loudly below.
    if (is_response_wf && optimize_geometry) {
      // The response step would have to run at the optimized geometry, which
      // needs the downstream consumer work of ARCHITECTURE_ROADMAP change 2.
      throw std::runtime_error(
          "--optimize is not supported for --wf=response; optimize the "
          "geometry first (--optimize --wf=scf) and run the response "
          "calculation on the result");
    }

    if (is_response_wf) {
      const std::string engine = pm.get<ResponseParameters>().engine();
      if (engine == "v2") {
        if (world.rank() == 0)
          print("ERROR: response engine 'v2' (MolresponseLib) was removed.\n"
                "       Delete the `engine v2` line (v3 is the default) — the\n"
                "       input deck is otherwise unchanged; see\n"
                "       molresponse_v3/MIGRATION_FROM_V2.md.");
        throw std::runtime_error("response engine v2 removed — use v3");
      }
      if (world.rank() == 0)
        print("response engine : molresponse_v3");
      auto reference = make_scf_reference(world, pm);
      wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(reference));
      wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(
          std::make_unique<ResponseApplication<molresponse_v3_lib>>(
              world, pm, reference->calc())));
    } else if (user_workflow == "scf" && !optimize_geometry) {
      wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(
          make_scf_reference(world, pm)));
    } else {
      workflow_builders::add_workflow_drivers(world, pm, user_workflow, wf,
                                              optimize_geometry);
    }

// io provenance for ALL tasks (today only response writes MRA restart
// state, but the stamp is run-level: an HDF5 run and a native run must
// be distinguishable from any task's calc_info).
#ifdef MADNESS_HAS_HDF5
    wf.set_provenance(
        "io",
        {{"backend", molresponse_v3::hdf5_io_enabled() ? "hdf5" : "native"},
         {"hdf5_compiled", true}});
    madness::iteration::archive_io_backend() =
        molresponse_v3::hdf5_io_enabled() ? "hdf5" : "native";
#else
    wf.set_provenance("io", {{"backend", "native"}, {"hdf5_compiled", false}});
#endif

    // Echo the effective input of the whole run before any work starts: the
    // molecule, then each driver's own parameter group. This is the only
    // record of what was actually computed (defaults and derived values
    // included) for a run whose deck sets just a handful of keys.
    if (world.rank() == 0) {
      print_header1("Calculation parameters");
      pm.get<Molecule>().print();
      wf.print_parameters(world);
      print("");
      print_header1("Starting calculations");
    }

    std::string prefix = pm.prefix();
    wf.run(prefix, shared_output);

    // Human-readable chemistry report: to stdout and <prefix>.out. The
    // machine-readable <prefix>.calc_info.json remains the source of truth.
    if (world.rank() == 0) {
      qcapp::write_results_summary(std::cout, wf.results());
      // APPEND, never truncate: with the usual `#SBATCH --output=<prefix>.out`
      // this file IS the job's stdout, and opening it for writing wiped the
      // whole run log at the end of every successful run (closeout attempt
      // 12, 2026-09-10: a 97-minute LiH log reduced to the summary).
      std::ofstream report(prefix + ".out", std::ios::app);
      qcapp::write_results_summary(report, wf.results());
      print("Wrote results summary :", prefix + ".out");

      // Index any visualization artifacts (cube/dx/line plots) the run
      // produced into <prefix>.viz_manifest.json for gecko/ParaView/VMD.
      qcapp::write_viz_manifest(prefix, wf.results());
    }
  } // solver workflow
}
} // namespace madqc
