#include "../solvers/es_save_load.hpp"
#include <madness/world/MADworld.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <unistd.h>

using namespace madness;
using namespace molresponse_v3;

static double gaussian(const coord_3d& r) { return std::exp(-(r[0]*r[0] + r[1]*r[1] + r[2]*r[2])); }

int main(int argc, char** argv) {
  World& world = initialize(argc, argv);
  int result = 0;
  try {
    startup(world, argc, argv, true);
    FunctionDefaults<3>::set_cubic_cell(-10.0, 10.0);
    FunctionDefaults<3>::set_k(8);
    FunctionDefaults<3>::set_thresh(1e-8);
    std::string dir = world.rank() == 0 ? "test-es-kinetic-" + std::to_string(getpid()) : "";
    world.gop.broadcast_serializable(dir, 0);
    real_function_3d f = real_factory_3d(world).f(gaussian);
    f.scale(1.0/f.norm2());
    ESSolver<Full, ClosedShell>::State state;
    ResponseStateXY<ClosedShell> root;
    root.x_alpha = {f};
    root.y_alpha = {copy(f)};
    root.y_alpha[0].scale(0.25);
    state.roots.push_back(root);
    state.stable_index = {7};
    state.omega = Tensor<double>(1);
    state.omega(0L) = 0.5;
    state.iter = 1;
    // Synthetic finite residual metadata required by the existing native loader.
    state.last_bsh_residual = {0.125};
    state.last_density_residual = {0.25};
    const auto functions = function_kinetic_record(measure_function_kinetics(world, root),
        1, "", 1e-8, 8, false, "solver_step_output");
    state.kinetic_history.push_back(es_kinetic_evidence(functions.at("blocks"),
        1, 0, 1, 7, "full", "closed_shell", 1e-8, 8, false));
    save_es_roots<Full, ClosedShell>(world, state, dir, false, 0.0, false);
    auto loaded = load_es_roots<Full, ClosedShell>(world, dir);
    MADNESS_CHECK(loaded.iter == 0 && loaded.stable_index == state.stable_index);
    MADNESS_CHECK(loaded.kinetic_history.size() == 1);
    MADNESS_CHECK(loaded.kinetic_history[0]["kinetic_energy"]["sample"]["attempt_directory"] == dir);
    const auto measured = function_kinetic_record(measure_function_kinetics(world, loaded.roots[0]),
        0, dir, 1e-8, 8, true, "loaded_root");
    MADNESS_CHECK(std::abs(measured.at("trace").at("value").get<double>() -
        functions.at("trace").at("value").get<double>()) < 1e-12);
    nlohmann::json original;
    const auto edit = [&](auto operation) {
      std::string error;
      if (world.rank() == 0) {
        try { operation(); }
        catch (const std::exception& e) { error = e.what(); }
      }
      world.gop.broadcast_serializable(error, 0);
      if (!error.empty()) throw std::runtime_error(error);
      world.gop.fence();
    };
    edit([&] { std::ifstream in(dir + "/roots.json"); in >> original; });
    for (const std::string corruption : {"trace", "model", "shell", "slot"}) {
      edit([&] {
        auto damaged = original;
        auto& report = damaged["kinetic_history"][0];
        if (corruption == "trace") report["kinetic_energy"]["trace"]["value"] = 0.0;
        else {
          auto blocks = report["kinetic_energy"]["blocks"];
          if (corruption == "model") blocks.erase("y_alpha");
          if (corruption == "shell") {
            blocks["x_beta"] = blocks["x_alpha"];
            blocks["y_beta"] = blocks["y_alpha"];
          }
          report = es_kinetic_evidence(blocks, 1, corruption == "slot" ? 1 : 0,
              corruption == "slot" ? 2 : 1, 7, corruption == "model" ? "tda" : "full",
              corruption == "shell" ? "open_shell" : "closed_shell", 1e-8, 8, false);
          report["kinetic_energy"]["sample"]["attempt_directory"] = dir;
        }
        std::ofstream out(dir + "/roots.json");
        out.exceptions(std::ios::failbit | std::ios::badbit);
        out << damaged.dump(2); out.close();
      });
      bool rejected = false;
      try { auto bad = load_es_roots<Full, ClosedShell>(world, dir); }
      catch (const std::runtime_error&) { rejected = true; }
      MADNESS_CHECK(rejected);
    }
    edit([&] { std::ofstream out(dir + "/roots.json"); out << original.dump(2); });
    // Final energy sorting permutes the roots and bookkeeping labels; old
    // sample slots retain their historical meaning across archive roundtrip.
    auto permuted = state;
    ResponseStateXY<ClosedShell> second;
    second.x_alpha = {copy(f)};
    second.y_alpha = {copy(f)};
    second.y_alpha[0].scale(0.5);
    permuted.roots.push_back(second);
    permuted.stable_index = {7, 9};
    permuted.omega = Tensor<double>(2);
    permuted.omega(0L) = 0.7;
    permuted.omega(1L) = 0.5;
    permuted.last_bsh_residual = {0.125, 0.25};
    permuted.last_density_residual = {0.25, 0.5};
    permuted.kinetic_history[0]["identity"]["root_count"] = 2;
    const auto second_functions = function_kinetic_record(measure_function_kinetics(world, second),
        1, "", 1e-8, 8, false, "solver_step_output");
    permuted.kinetic_history.push_back(es_kinetic_evidence(second_functions.at("blocks"),
        1, 1, 2, 9, "full", "closed_shell", 1e-8, 8, false));
    const auto historical = permuted.kinetic_history;
    ESProblem<Full, ClosedShell> problem;
    problem.n_roots = 2;
    ESSolver<Full, ClosedShell> sorter(world, std::move(problem), ConvergencePolicy{}, PrintLevel::Silent);
    sorter.sort_state_by_omega(permuted);
    MADNESS_CHECK(permuted.stable_index[0] == 9 && permuted.kinetic_history == historical);
    save_es_roots<Full, ClosedShell>(world, permuted, dir + "/permuted", false, 0.0, false);
    auto roundtrip = load_es_roots<Full, ClosedShell>(world, dir + "/permuted");
    MADNESS_CHECK(roundtrip.stable_index == permuted.stable_index);
    MADNESS_CHECK(roundtrip.kinetic_history[0]["identity"]["slot"] == 0 &&
        roundtrip.kinetic_history[0]["identity"]["bookkeeping_index"] == 7);
    const auto first_loaded = function_kinetic_record(measure_function_kinetics(world, roundtrip.roots[0]),
        0, dir, 1e-8, 8, true, "loaded_root");
    MADNESS_CHECK(std::abs(first_loaded.at("trace").at("value").get<double>() -
        second_functions.at("trace").at("value").get<double>()) < 1e-12);
    if (world.rank() == 0) print("ES kinetic archive checks passed; artifacts:", dir);
  } catch (const std::exception& error) {
    if (world.rank() == 0) std::cerr << error.what() << '\n';
    result = 1;
  }
  finalize();
  return result;
}
