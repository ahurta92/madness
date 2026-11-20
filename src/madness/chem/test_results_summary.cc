#include <madness/chem/ResultsSummary.hpp>
#include <madness/mra/mra.h>
#include <madness/world/test_utilities.h>

#include <sstream>
#include <string>

using namespace madness;

static bool has_text(const std::string &s, const std::string &needle) {
  return s.find(needle) != std::string::npos;
}

int main(int argc, char **argv) {
  World &world = madness::initialize(argc, argv);
  int success = 0;

  if (world.rank() == 0) {
    test_output t("results summary citations");

    const nlohmann::json base_scf = {
        {"model", "scf"},
        {"properties", {{"energy", -1.0}}},
        {"citations", {{"dftd3", false}, {"pcm", false}, {"libxc", false}}}};

    {
      std::ostringstream os;
      qcapp::write_results_summary(os, {{"tasks", nlohmann::json::array({base_scf})}});
      const auto out = os.str();
      const bool ok = !has_text(out, "Citations for external modules");
      t.checkpoint(ok, "omits citations block when no module is requested");
      success += !ok;
    }

    {
      std::ostringstream os;
      nlohmann::json scf = base_scf;
      scf["citations"]["dftd3"] = true;
      scf["citations"]["pcm"] = true;
      scf["citations"]["libxc"] = true;
      qcapp::write_results_summary(os, {{"tasks", nlohmann::json::array({scf})}});
      const auto out = os.str();
      const bool ok = has_text(out, "Citations for external modules") &&
                      has_text(out, "DFT-D3 / simple-dftd3") &&
                      has_text(out, "PCMSolver (PCM)") &&
                      has_text(out, "Libxc");
      t.checkpoint(ok, "prints module citations only when requested");
      success += !ok;
    }

    {
      auto scf = base_scf;
      scf["convergence_info"] = {{"converged_for_thresh", 1e10}, {"converged_for_dconv", 1e10}};
      scf["outcome"] = {{"complete", false}, {"reason", "iteration_limit"}};
      scf["protocol"] = {{"threshold", 1e-6}, {"k", 8}};
      scf["metadata"] = {{"mpi_size", 2}};
      std::ostringstream os;
      qcapp::write_results_summary(os, {{"tasks", nlohmann::json::array({scf})}});
      auto out = os.str();
      bool ok = has_text(out, "NOT converged") && has_text(out, "iteration_limit") &&
          has_text(out, "k = 8") && !has_text(out, "dconv") &&
          !has_text(out, "1e+10") && !has_text(out, "Wall time");
      t.checkpoint(ok, "shared summary uses outcome and never invents elapsed time");
      success += !ok;
      scf["outcome"] = {{"complete", true}, {"reason", "final_criteria_met"}};
      scf["metadata"] = {{"elapsed_time", 12.345}, {"mpi_size", 2}, {"nthreads", 4}};
      os.str(""); os.clear();
      qcapp::write_results_summary(os, {{"tasks", nlohmann::json::array({scf})}});
      out = os.str();
      ok = has_text(out, "Converged") && !has_text(out, "NOT converged") &&
          has_text(out, "final_criteria_met") && has_text(out, "12.345") &&
          has_text(out, "2 MPI x 4 threads") && !has_text(out, "dconv");
      t.checkpoint(ok, "shared summary prints measured invocation time and execution");
      success += !ok;
    }
    {
      nlohmann::json response = {{"type", "response"}, {"metadata", {
        {"elapsed_time", 3.125}, {"v3_timing", {
          {"plan_build", {{"wall_s", 0.125}}}, {"solve", {{"wall_s", 2.5}}},
          {"assemble", {{"wall_s", 0.25}}}}}}}};
      std::ostringstream os;
      qcapp::write_results_summary(os, {{"tasks", nlohmann::json::array({response})}});
      const auto out = os.str();
      const bool ok = has_text(out, "3.125") && has_text(out, "solve 2.500") &&
          has_text(out, "assemble 0.250") && !has_text(out, "MPI x");
      t.checkpoint(ok, "response summary exposes phase times without fabricated thread count");
      success += !ok;
    }

    {
      std::ostringstream os;
      nlohmann::json failed = {{"type", "task_failed"},
                               {"status", "failed"},
                               {"citations",
                                {{"dftd3", false}, {"pcm", true}, {"libxc", false}}}};
      qcapp::write_results_summary(os, {{"tasks", nlohmann::json::array({failed})}});
      const auto out = os.str();
      const bool ok =
          has_text(out, "Citations for external modules") && has_text(out, "PCMSolver (PCM)");
      t.checkpoint(ok, "prints citations for requested modules even when task is failed");
      success += !ok;
    }
  }

  int total_failures = success;
  world.gop.sum(total_failures);
  madness::finalize();
  return total_failures;
}
