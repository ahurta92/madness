// Unit test for response_task_envelope(): the task-record envelope derived from a
// response_metadata document. Pure C++ -- no MPI, no World.
#include "../orchestrator/response_envelope.hpp"

#include <cstdio>

namespace {
int failed = 0;
#define EXPECT(cond, label)                                                \
  do {                                                                     \
    if (cond) { std::printf("  [PASS]  %s\n", label); }                    \
    else      { std::printf("  [FAIL]  %s\n", label); ++failed; }          \
  } while (0)

nlohmann::json two_rung_metadata() {
  using json = nlohmann::json;
  json md;
  md["protocols"] = {{"1e-04_k6", {{"index", 0}, {"thresh", 1e-4}, {"k", 6}}},
                     {"1e-06_k8", {{"index", 1}, {"thresh", 1e-6}, {"k", 8}}}};
  md["fd_states"]["dipole_z"]["1e-04_k6"]["f0.00000"] = {{"converged", true}, {"iter", 5}};
  md["fd_states"]["dipole_z"]["1e-06_k8"]["f0.00000"] = {{"converged", true}, {"iter", 7}};
  md["fd_states"]["dipole_z"]["1e-06_k8"]["f0.02000"] = {{"converged", true}, {"iter", 9}};
  md["excited_states"]["1e-06_k8"] = {{"converged", true}, {"iter", 12}, {"n_roots", 2}};
  md["vbc_states"]["vbc:dipole_z__dipole_z@f0.00000_f0.00000"]["1e-06_k8"] =
      {{"converged", true}, {"metrics", {{"iters", 3}}}};
  md["run_summary"] = {{"stop_reason", "complete"}, {"dropped_work", json::array()}};
  return md;
}
} // namespace

int main() {
  using molresponse_v3::response_task_envelope;

  {
    const auto env = response_task_envelope(two_rung_metadata());
    EXPECT(env.precision["k"] == 8, "precision.k is the finest rung's k");
    EXPECT(env.precision["thresh"] == 1e-6, "precision.thresh is the finest rung's thresh");
    EXPECT(env.precision["protocol_key"] == "1e-06_k8", "precision.protocol_key names the finest rung");
    EXPECT(env.precision["protocol"] == nlohmann::json({1e-4, 1e-6}), "precision.protocol lists rungs in index order");
    EXPECT(env.convergence["status"] == "converged", "all states converged + complete -> converged");
    EXPECT(env.convergence["iterations"] == 12, "iterations is the max over finest-rung states (ES 12)");
    EXPECT(env.convergence["n_states"] == 4, "n_states counts 2 FD + 1 ES + 1 VBC at the finest rung");
    EXPECT(env.convergence["n_unconverged"] == 0, "n_unconverged is 0");
    EXPECT(env.convergence["stop_reason"] == "complete", "stop_reason is copied");
  }
  {
    auto md = two_rung_metadata();
    md["fd_states"]["dipole_z"]["1e-06_k8"]["f0.02000"]["converged"] = false;
    const auto env = response_task_envelope(md);
    EXPECT(env.convergence["status"] == "unconverged", "one unconverged FD state -> unconverged");
    EXPECT(env.convergence["n_unconverged"] == 1, "n_unconverged counts it");
  }
  {
    auto md = two_rung_metadata();
    md["fd_states"]["dipole_z"]["1e-04_k6"]["f0.00000"]["converged"] = false;
    const auto env = response_task_envelope(md);
    EXPECT(env.convergence["status"] == "converged", "a coarse-rung state does not count");
  }
  {
    auto md = two_rung_metadata();
    md["run_summary"]["stop_reason"] = "wall_limit";
    const auto env = response_task_envelope(md);
    EXPECT(env.convergence["status"] == "unconverged", "stop_reason != complete -> unconverged");
  }
  {
    const auto env = response_task_envelope(nlohmann::json::object());
    EXPECT(env.precision.is_null(), "no protocols -> precision null");
    EXPECT(env.convergence["status"] == "unknown", "no protocols -> status unknown");
    EXPECT(env.convergence["iterations"] == -1, "no protocols -> iterations -1");
    EXPECT(env.convergence["stop_reason"] == "unknown", "no run_summary -> stop_reason unknown");
  }

  std::printf("%d failure(s)\n", failed);
  return failed ? 1 : 0;
}
