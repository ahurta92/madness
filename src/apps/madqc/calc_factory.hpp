#ifndef CALC_STRATS_HPP
#define CALC_STRATS_HPP

#include "calc_manager.hpp"
#include "parameter_manager.hpp"

// The driver is the calculation type (e.g. energy, gradient, hessian, etc.)
// The model is the combination of theory and basis set
// For us, the basis set is adaptive, so we don't need to worry about that
// we simply compute to a certain precision

// Let's start of defining the energy driver.  I need to create a calc_manager
// The models are going to reflect which application we are running with
// - Moldft
// - Response (moldft + molresponse)
// - MP2 (moldft + mp2)
// - CIS (moldft + cis)
// - OEP (moldft + oep)
// - MP3 (moldft + mp3)
//
//
//
//
//

// define the properties that can be calculated for each model
// This example is for the response model
struct response_property_map {
  bool alpha;
  bool beta;
  bool shg;
};

void to_json(nlohmann::json& j, const response_property_map& r) {
  j = nlohmann::json{{"alpha", r.alpha}, {"beta", r.beta}, {"shg", r.shg}};
}
void from_json(const nlohmann::json& j, response_property_map& r) {
  j.at("alpha").get_to(r.alpha);
  j.at("beta").get_to(r.beta);
  j.at("shg").get_to(r.shg);
}

// A property map is a map of properties to be calculated for each model
using model = std::string;
using model_properties = std::map<std::string, bool>;
//maybe just use json for this?
using property_map = std::map<model, model_properties>;

std::unique_ptr<CalcManager> createCalcManager(const std::string& model_name, ParameterManager& pm,
                                               property_map properties) {
  // Create a new CalcManager
  auto calc_manager = std::make_unique<CalcManager>();
  auto molecule = pm.get_molecule();
  // All calculations start with a reference
  auto params = pm.get_moldft_params();
  auto moldft_properties = properties["moldft"];
  auto moldft = std::make_unique<MoldftCalculationStrategy>(params, molecule, "moldft", moldft_properties);
  calc_manager->addStrategy(std::move(moldft));

  if (model_name == "moldft") {

    return calc_manager;

  } else if (model_name == "response") {

    auto& response_params = pm.get_molresponse_params();
    auto response_properties = properties.at("response");  // map of properties to be calculated for response

    auto perturbations = response_params.perturbations();
    // say you want to do moldft hf and response lda... well it's possible here TODO: vecotrize this so you can do mixed?
    auto xc = response_params.xc();
    auto freq_range = response_params.freq_range();

    vector<std::string> input_names = {"moldft"};
    for (auto const& perturbation : perturbations) {
      ResponseInput r_input = std::make_tuple(perturbation, xc, freq_range);
      auto response_calc =
          std::make_unique<LinearResponseStrategy>(response_params, r_input, perturbation, input_names);
      calc_manager->addStrategy(std::move(response_calc));
    }

    if (response_properties["beta"]) {

      auto set_freqs = [&]() {
        vector<double> freqs_copy = freq_range;
        auto num_freqs = freq_range.size();
        auto compare_freqs = [](double x, double y) {
          return std::abs(x - y) < 1e-3;
        };

        for (int i = 0; i < num_freqs; i++) {    // for i=0:n-1
          for (int j = i; j < num_freqs; j++) {  // for j = i  omega_3=-(omega_1+omega_2)
            auto omega_1 = freq_range[i];
            auto omega_2 = freq_range[j];
            auto omega_3 = omega_1 + omega_2;

            // if you can find omega_3 in freq_copy skip
            if (omega_2 == 0.0)
              continue;
            if (std::find_if(freqs_copy.begin(), freqs_copy.end(),
                             [&](double x) { return compare_freqs(x, omega_3); }) != freqs_copy.end()) {
              continue;
            } else {
              freqs_copy.push_back(omega_3);
            }
          }
        }
        return freqs_copy;
      };

      freq_range = set_freqs();
      print("Frequency Range: ", freq_range);
      // this is where I we create our calculation
      ResponseInput hyp_input = std::make_tuple("dipole", xc, freq_range);
      std::vector<std::string> r_input_names = {"moldft"};
      std::vector<std::string> h_input_names = {"moldft", "dipole"};
      auto response_hyper =
          std::make_unique<LinearResponseStrategy>(response_params, hyp_input, "dipole", r_input_names);
      auto hyper_calc = std::make_unique<ResponseHyper>(response_params, hyp_input, "hyper", h_input_names);
      calc_manager->addStrategy(std::move(response_hyper));
      calc_manager->addStrategy(std::move(hyper_calc));
    }
    if (response_properties["shg"]) {
      throw std::invalid_argument("SHG not implemented yet");
    }
  } else if (model_name == "MP2") {
    throw std::invalid_argument("MP2 not implemented yet");
    /*auto moldft_strategy = std::make_unique<MoldftCalculationStrategy>(params, molecule);*/
    /*auto mp2_strategy = std::make_unique<MP2CalculationStrategy>(params, molecule);*/
    /*calc_manager->addStrategy(std::move(moldft_strategy));*/
    /*calc_manager->addStrategy(std::move(mp2_strategy));*/
  } else if (model_name == "CIS") {
    throw std::invalid_argument("CIS not implemented yet");
    /*auto moldft_strategy = std::make_unique<MoldftCalculationStrategy>(params, molecule);*/
    /*auto cis_strategy = std::make_unique<CISCalculationStrategy>(params, molecule);*/
    /*calc_manager->addStrategy(std::move(moldft_strategy));*/
    /*calc_manager->addStrategy(std::move(cis_strategy));*/
  } else if (model_name == "OEP") {
    throw std::invalid_argument("OEP not implemented yet");
    /*auto moldft_strategy = std::make_unique<MoldftCalculationStrategy>(params, molecule);*/
    /*auto oep_strategy = std::make_unique<OEPCalculationStrategy>(params, molecule);*/
    /*calc_manager->addStrategy(std::move(moldft_strategy));*/
    /*calc_manager->addStrategy(std::move(oep_strategy));*/
  } else if (model_name == "MP3") {
    throw std::invalid_argument("MP3 not implemented yet");
    /*auto moldft_strategy = std::make_unique<MoldftCalculationStrategy>(params, molecule);*/
    /*auto mp3_strategy = std::make_unique<MP3CalculationStrategy>(params, molecule);*/
    /*calc_manager->addStrategy(std::move(moldft_strategy));*/
    /*calc_manager->addStrategy(std::move(mp3_strategy));*/
  } else {
    throw std::invalid_argument("Unknown model name: " + model_name);
  }
  // Return the configured CalcManager
  return calc_manager;
}

// Below here we define more complex calculation strategies
DynamicInput example_calculation(Molecule& molecule) {

  // Example setup calculation function
  auto setupCalculation = [&](CalcManager& calc_manager, const std::string& iter_name) {
    // Create strategies and set them up for the current iteration
    auto params = CalculationParameters();
    auto moldft_strategy = std::make_unique<MoldftCalculationStrategy>(params, molecule);
    calc_manager.addStrategy(std::move(moldft_strategy));
  };

  return DynamicInput{};
}

#endif
