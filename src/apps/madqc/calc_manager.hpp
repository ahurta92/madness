#ifndef CALC_MANAGER_HPP
#define CALC_MANAGER_HPP

#include <madchem.h>
#include <madness/chem/CalculationParameters.h>
#include <madness/chem/SCF.h>
#include <madness/chem/molecule.h>
#include <apps/molresponse/FrequencyResponse.hpp>
#include <apps/molresponse/ResponseExceptions.hpp>
#include <filesystem>
#include <madness/external/nlohmann_json/json.hpp>
#include <memory>
#include <utility>
#include <vector>
#include "path_manager.hpp"
#include "utils.hpp"

using json = nlohmann::json;
using path = std::filesystem::path;

class CalculationStrategy {
 public:
  virtual void runCalculation(World& world, const json& paths) = 0;
  virtual ~CalculationStrategy() = default;
};

class CompositeCalculationStrategy : public CalculationStrategy {
 private:
  std::vector<std::unique_ptr<CalculationStrategy>> strategies;

 public:
  void addStrategy(std::unique_ptr<CalculationStrategy> strategy) {
    strategies.push_back(std::move(strategy));
  }

  void runCalculation(World& world, const json& paths) override {
    for (const auto& strategy : strategies) {
      strategy->runCalculation(world, paths);
    }
  }
};

class MoldftCalculationStrategy : public CalculationStrategy {

  CalculationParameters parameters;
  json input_json;
  Molecule molecule;

 public:
  void runCalculation(World& world, const json& paths) override {

    // the first step is to look for the moldft paths
    json moldft_paths = paths["moldft"];

    // Get the paths from the json object
    path moldft_path = moldft_paths["calculation"];
    path moldft_restart = moldft_paths["restart"];
    path moldft_calc_info_path = moldft_paths["output"]["calc_info"];

    // print the relevant paths
    if (world.rank() == 0) {
      std::cout << "Calculation Path: " << moldft_paths.dump(4) << std::endl;
    }

    if (!std::filesystem::exists(moldft_path)) {
      std::filesystem::create_directory(moldft_path);
    }

    std::filesystem::current_path(moldft_path);

    json calcInfo;
    auto param1 = parameters;
    world.gop.broadcast_serializable(param1, 0);

    if (world.rank() == 0) {
      ::print("-------------Running moldft------------");
    }

    if (std::filesystem::exists(moldft_restart) &&
        std::filesystem::exists(moldft_calc_info_path)) {
      // if both exist, read the calc_info json
      std::ifstream ifs(moldft_calc_info_path);

      auto moldft_calc_info = json::parse(ifs);
      if (world.rank() == 0) {
        std::cout << "time: " << moldft_calc_info["time"] << std::endl;
        std::cout << "MOLDFT return energy: "
                  << moldft_calc_info["return_energy"] << std::endl;
      }
    } else {
      // if params are different run and if restart exists and if im asking to
      if (std::filesystem::exists(moldft_restart)) {
        param1.set_user_defined_value<bool>("restart", true);
      }
      world.gop.fence();
      if (world.rank() == 0) {

        json moldft_input_json = {};
        moldft_input_json["dft"] = parameters.to_json_if_precedence("defined");
        moldft_input_json["molecule"] = molecule.to_json();

        std::ofstream ofs("moldft.in");
        write_moldft_input(moldft_input_json, ofs);
        ofs.close();
      }
      world.gop.fence();
      commandlineparser parser;
      parser.set_keyval("input", "moldft.in");

      if (world.rank() == 0)
        ::print("input filename: ", parser.value("input"));

      print_meminfo(world.rank(), "startup");
      FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3>>(world)));

      std::cout.precision(6);
      SCF calc(world, parser);
      if (world.rank() == 0) {
        ::print("\n\n");
        ::print(
            " MADNESS Hartree-Fock and Density Functional Theory "
            "Program");
        ::print(
            " ----------------------------------------------------------"
            "\n");
        calc.param.print("dft");
      }
      if (world.size() > 1) {
        calc.set_protocol<3>(world, 1e-4);
        calc.make_nuclear_potential(world);
        calc.initial_load_bal(world);
      }
      // vama
      calc.set_protocol<3>(world, calc.param.protocol()[0]);
      // calc.set_protocol<3>(world, 1e-4);
      world.gop.fence();
      MolecularEnergy ME(world, calc);
      world.gop.fence();
      // double energy=ME.value(calc.molecule.get_all_coords().flat()); // ugh!
      ME.value(calc.molecule.get_all_coords().flat());  // ugh!
      world.gop.fence();
      const real_function_3d rho =
          2.0 * calc.make_density(world, calc.aocc, calc.amo);
      auto dipole_t = calc.dipole(world, rho);
      std::map<std::string, double> results;
      results["scf_energy"] = calc.current_energy;
      world.gop.fence();
      if (world.rank() == 0) {
        calc.output_scf_info_schema(results, dipole_t);
        ME.output_calc_info_schema();
      }
    }

    // Add actual calculation logic here
  }

  MoldftCalculationStrategy(const CalculationParameters& params, Molecule mol)
      : parameters(params), molecule(std::move(mol)){};
};

class MP2CalculationStrategy : public CalculationStrategy {
 public:
  void runCalculation(World& world, const json& paths) override {
    std::cout << "Running MP2 Calculation\n";
    std::cout << "Calc Directory: " << paths["mp2_calc_dir"] << "\n";
    std::cout << "Restart Path: " << paths["mp2_restart"] << "\n";
    std::cout << "Calc Info Path: " << paths["mp2_calc_info"] << "\n";
    // Add actual calculation logic here
  }
};

class ResponseCalculationStrategy : public CalculationStrategy {

  ResponseParameters parameters;
  std::string op;
  std::string xc;
  std::vector<double> freqs;

 public:
  explicit ResponseCalculationStrategy(const ResponseParameters& params)
      : parameters(params) {}

  static void append_to_alpha_json(const double& omega,
                                   const std::vector<std::string>& ij,
                                   const Tensor<double>& alpha,
                                   nlohmann::ordered_json& alpha_json) {
    auto num_unique_elements = ij.size();
    for (int i = 0; i < num_unique_elements; i++) {
      alpha_json["omega"].push_back(omega);
      alpha_json["ij"].push_back(ij[i]);
      alpha_json["alpha"].push_back(alpha[i]);
    }
  }
  static void add_alpha_i_to_json(nlohmann::ordered_json& alpha_i,
                                  nlohmann::ordered_json& alpha_json) {
    print(alpha_json.dump(4));

    alpha_json["omega"].insert(alpha_json["omega"].end(),
                               alpha_i["omega"].begin(),
                               alpha_i["omega"].end());
    alpha_json["ij"].insert(alpha_json["ij"].end(), alpha_i["ij"].begin(),
                            alpha_i["ij"].end());
    alpha_json["alpha"].insert(alpha_json["alpha"].end(),
                               alpha_i["alpha"].begin(),
                               alpha_i["alpha"].end());
  }

  /**
     *
     * @param world
     * @param filename
     * @param frequency
     * @param property
     * @param xc
     * @param moldft_path
     * @param restart_path
     * @return
     */
  bool runFrequency(World& world, ResponseParameters& r_params,
                    double frequency) {

    op = r_params.perturbation();

    // Set the response parameters

    if (world.rank() == 0) {

      json input_json = {};
      input_json["response"] = r_params.to_json_if_precedence("defined");
      std::ofstream out("response.in");
      write_json_to_input_file(input_json, {"response"}, out);
    }
    bool converged = false;
    // if rbase exists and converged I just return save path and true
    if (world.rank() == 0) {
      ::print("Checking if response has converged for frequency: ", frequency);

      if (std::filesystem::exists("response_base.json")) {
        {
          std::ifstream ifs("response_base.json");
          json response_base;
          ifs >> response_base;
          if (response_base["converged"] &&
              response_base["precision"]["dconv"] == r_params.dconv())
            converged = true;
          ::print("Response has converged: ", converged);
        }
      }
    }
    world.gop.broadcast(converged, 0);
    // add logic to compute alpha.json if needed

    if (converged) {
      return true;
    } else {

      world.gop.fence();
      if (world.rank() == 0) {
        ::print("Running response calculation for frequency: ", frequency);
      }

      GroundStateCalculation ground_calculation{world, r_params.archive()};
      Molecule molecule = ground_calculation.molecule();
      r_params.set_ground_state_calculation_data(ground_calculation);
      r_params.set_derived_values(world, molecule);
      CalcParams calc_params = {ground_calculation, molecule, r_params};

      RHS_Generator rhs_generator;
      if (op == "dipole") {
        rhs_generator = dipole_generator;
      } else {
        rhs_generator = nuclear_generator;
      }
      FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3>>(world)));
      FrequencyResponse calc(world, calc_params, frequency, rhs_generator);
      if (world.rank() == 0) {
        ::print("\n\n");
        ::print(
            " MADNESS Time-Dependent Density Functional Theory Response "
            "Program");
        ::print(
            " ----------------------------------------------------------\n");
        ::print("\n");
        calc_params.molecule.print();
        ::print("\n");
        calc_params.response_parameters.print("response");
        // put the response parameters in a j_molrespone json object
        calc_params.response_parameters.to_json(calc.j_molresponse);
      }
      calc.solve(world);
      auto [omega, polar_omega] = calc.get_response_data();
      // flatten polar_omega

      auto alpha = polar_omega.flat();

      std::vector<std::string> ij{"XX", "XY", "XZ", "YX", "YY",
                                  "YZ", "ZX", "ZY", "ZZ"};
      nlohmann::ordered_json alpha_json;
      alpha_json["omega"] = {};
      alpha_json["ij"] = {};
      alpha_json["alpha"] = {};

      append_to_alpha_json(omega, ij, alpha, alpha_json);
      calc.j_molresponse["properties"] = {};
      calc.j_molresponse["properties"]["alpha"] = alpha_json;

      // set protocol to the first
      if (world.rank() == 0) {
        // calc.time_data.to_json(calc.j_molresponse);
        calc.output_json();
      }
      // calc.time_data.print_data();
      return calc.j_molresponse["converged"];
    }
  }

  void runCalculation(World& world, const json& paths) override {

    auto& moldft_paths = paths["moldft"];
    auto moldft_restart = moldft_paths["restart"].get<std::string>();
    moldft_restart = std::string(path(moldft_restart).replace_extension(""));

    auto& response_paths = paths["response"];
    auto& calc_paths = response_paths["calculation"];

    auto restart_paths =
        response_paths["restart"].get<std::vector<std::string>>();
    auto output_paths =
        response_paths["output"].get<std::vector<std::string>>();
    auto alpha_path = response_paths["properties"]["alpha"].get<std::string>();

    size_t num_freqs = calc_paths.size();
    // I need to analyze alpha.json to see which frequencies are missing.
    // Or I just write them to seperate files and then combine them later?

    auto freqs = parameters.freq_range();

    bool last_converged = false;
    nlohmann::ordered_json alpha_json;
    alpha_json["omega"] = json::array();
    alpha_json["ij"] = json::array();
    alpha_json["alpha"] = json::array();

    for (size_t i = 0; i < num_freqs; i++) {

      auto freq_i = freqs[i];

      if (!std::filesystem::exists(calc_paths[i])) {
        std::filesystem::create_directory(calc_paths[i]);
      }
      std::filesystem::current_path(calc_paths[i]);

      print("current path: ", std::filesystem::current_path());
      print("calc path: ", calc_paths[i]);
      print("freq: ", freq_i);
      // if the last converged is true, then we can restart from the last save path

      bool restart = true;
      path save_path = restart_paths[i];  // current restart path aka save path
      path restart_path = (i > 0) ? restart_paths[i - 1] : restart_paths[i];
      std::string save_string = save_path.filename().stem();
      if (world.rank() == 0) {
        ::print("-------------Running response------------ at frequency: ",
                freq_i);
        ::print("moldft restart path", moldft_restart);
        ::print("restart path", restart_path);
        ::print("save path", save_path);
        ::print("save string", save_string);
      }

      ResponseParameters r_params = parameters;

      r_params.set_user_defined_value("omega", freq_i);
      r_params.set_user_defined_value("archive", moldft_restart);
      if (last_converged || i == 0) {

        if (world.rank() == 0) {
          r_params.set_user_defined_value("save", true);
          r_params.set_user_defined_value("save_file", save_string);
          if (restart) {  // if we are trying a restart calculation
            if (std::filesystem::exists(save_path)) {
              // if the save path exists then we know we can
              //  restart from the previous save
              r_params.set_user_defined_value("restart", true);
              r_params.set_user_defined_value("restart_file", save_string);
            } else if (std::filesystem::exists(restart_path)) {
              ::print("restart path exists", restart_path);
              r_params.set_user_defined_value("restart", true);
              ::print(restart_path.parent_path().stem());
              ::print(restart_path.filename().stem());

              // get the directory of the restart path
              auto new_restart_path = path("../") /
                                      restart_path.parent_path().stem() /
                                      restart_path.filename().stem();

              // format restart path to be ../restart_path/restart_path

              //
              ::print("new restart file: ", restart_path);
              r_params.set_user_defined_value("restart_file",
                                              new_restart_path.string());
              // Then we restart from the previous file instead
            } else {
              r_params.set_user_defined_value("restart", false);
            }
            // neither file exists therefore you need to start from fresh
          }
        }
      } else {
        throw Response_Convergence_Error{};
      }
      world.gop.broadcast_serializable(r_params, 0);
      last_converged = runFrequency(world, r_params, freq_i);

      nlohmann::ordered_json alpha_i;

      if (last_converged) {
        // read output_paths[i]
        std::ifstream ifs(output_paths[i]);
        json response_base_i;
        ifs >> response_base_i;
        alpha_i = response_base_i["properties"]["alpha"];
      }
      if (world.rank() == 0) {
        print("last converged: ", last_converged);
        print(alpha_i.dump(4));
      }

      // combine alpha_json_i with alpha_json
      add_alpha_i_to_json(alpha_i, alpha_json);
    }

    std::ofstream out_file(alpha_path);

    if (world.rank() == 0) {
      out_file << alpha_json.dump(4);
    }
  }
};

class CalcManager {
 private:
  std::unique_ptr<CalculationStrategy> strategy;
  PathManager pathManager;

 public:
  void setCalculationStrategy(
      std::unique_ptr<CalculationStrategy> newStrategy) {
    strategy = std::move(newStrategy);
  }

  void runCalculations(World& world, const path& root) {
    json paths = pathManager.generateCalcPaths(root);
    // output the paths to disk
    if (world.rank() == 0) {
      std::ofstream ofs("paths.json");
      ofs << paths.dump(4);
      ofs.close();
    }

    if (world.rank() == 0) {
      print(paths.dump(4));
    }

    strategy->runCalculation(world, paths);
  }

  void addPathStrategy(std::unique_ptr<PathStrategy> strategy) {
    pathManager.addStrategy(std::move(strategy));
  }

  CalcManager() = default;
};

#endif
