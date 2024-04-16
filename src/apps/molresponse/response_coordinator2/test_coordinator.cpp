//
// Created by adrianhurtado on 1/1/22.
//
#include "CalculationParameters.h"

#define CATCH_CONFIG_RUNNER

#include "FrequencyResponse.hpp"
#include "coordinator.hpp"
#include "madness/external/catch/catch.hpp"
#include "string"
#include <filesystem>


using path = std::filesystem::path;

int main(int argc, char *argv[]) {
    World &world = madness::initialize(argc, argv);
    int result = 0;
    world.gop.fence();
    startup(world, argc, argv);
    { result = Catch::Session().run(argc, argv); }

    return result;

}


std::pair<CalculationParameters, ResponseParameters> get_params(const std::string &json_input_path) {

    std::ifstream input_stream(json_input_path);
    json input_json = json::parse(input_stream);

    json moldft_json = input_json["dft"];
    json molresponse_json = input_json["response"];
    CalculationParameters moldft_params;


    moldft_params.from_json(moldft_json);
    ResponseParameters molresponse_params;
    molresponse_params.from_json(molresponse_json);
    return std::make_pair(moldft_params, molresponse_params);

}

TEST_CASE("Basic Response Manager Test") {
    // Set up the run directories
    using namespace madness;

    World &world = World::get_default();
    std::cout.precision(6);

    path json_input_path("resources/inputs/input.json");
    path molecule_path("resources/molecules/H2O.mol");
    // print full path
    // Read in json

    auto [moldft_params, molresponse_params] = get_params(json_input_path);
    auto response_manager = ResponseCalcManager(world, molecule_path, moldft_params, molresponse_params);
    response_manager.print();
    // The json is converted into a temporary getKW file which is then read by the parser.
    // Now we need to write a function
}

TEST_CASE("Param Hash Generation Test") {
    // Set up the run directories
    using namespace madness;

    World &world = World::get_default();
    std::cout.precision(6);

    path mol1_path("resources/molecules/H2O.mol");
    path mol2_path("resources/molecules/Be.mol");
    path json_input_path("resources/inputs/input.json");

    auto [moldft_params, molresponse_params] = get_params(json_input_path);


    auto response_manager = ResponseCalcManager(world, mol1_path, moldft_params, molresponse_params);
    auto response_manager2 = ResponseCalcManager(world, mol2_path, moldft_params, molresponse_params);

    // Test that the param directories are the same
    CHECK(response_manager.get_param_path() == response_manager2.get_param_path());
    // Test that the molecule directories are different
    CHECK(response_manager.get_moldft_path() != response_manager2.get_moldft_path());

}

TEST_CASE("Test Manager Run MOLDFT") {
    using namespace madness;

    World &world = World::get_default();
    std::cout.precision(6);

    path mol1_path("resources/molecules/He.mol");
    path json_input_path("resources/inputs/input.json");

    auto [moldft_params, molresponse_params] = get_params(json_input_path);


    auto response_manager = ResponseCalcManager(world, mol1_path, moldft_params, molresponse_params);

    response_manager.run_moldft(world, true);
}

// Create function which converts json input to and outputs to a output stream

void write_molecule_json_to_input_file(const json &molecule_json, std::ostream &output_stream) {

    output_stream << "geometry" << std::endl;
    output_stream << "    " << "units" << " " << std::string(molecule_json["units"]) << std::endl;
    output_stream << "    " << "eprec" << " " << molecule_json["eprec"] << std::endl;
    // for each key within the block write the key value pair to the file line by line
    // grab the geometry and symbols as two separate arrays
    std::vector<std::string> symbols;
    std::vector<std::vector<double>> geometry;

    auto coords = molecule_json["geometry"];
    auto symbols_json = molecule_json["symbols"];

    auto n = symbols_json.size();

    std::vector<std::pair<std::string, std::vector<double>>> lines(n);

    for (int i = 0; i < n; i++) {
        symbols.push_back(symbols_json[i]);
        geometry.push_back(coords[i]);
    }

    for (int i = 0; i < n; i++) {
        output_stream << "    " << symbols[i] << " " << geometry[i][0] << " " << geometry[i][1] << " " << geometry[i][2]
                      << std::endl;
    }
    output_stream << "end" << std::endl << std::endl;

}

void
write_json_to_input_file(const json &input_json, const std::vector<std::string> &keys, std::ostream &output_stream) {

    for (auto key: keys) {

        if (input_json.find(key) == input_json.end()) {
            throw std::runtime_error("Key not found in input json");
        }

        output_stream << key << std::endl;
        // for each key within the block write the key value pair to the file line by line
        for (auto &[key, value]: input_json[key].items()) {
            output_stream << "    " << key << " " << value << std::endl;
        }
        output_stream << "end" << std::endl << std::endl;

    }

}

void write_moldft_input(const json &input_json, std::ostream &out) {
    write_json_to_input_file(input_json, {"dft"}, out);
    write_molecule_json_to_input_file(input_json["molecule"], out);
}
void write_response_input(const json &input_json, std::ostream &out) {
    write_json_to_input_file(input_json, {"response"}, out);
}

TEST_CASE("Testing Function which writes input file from json input") {

    using namespace madness;

    World &world = World::get_default();
    std::cout.precision(6);
    path mol1_path("resources/molecules/H2.mol");
    path json_input_path("resources/inputs/input.json");

    auto [moldft_params, molresponse_params] = get_params(json_input_path);

    std::ifstream molecule_stream(mol1_path);
    Molecule molecule;
    if (!molecule_stream.is_open()) {
        throw std::runtime_error("Could not open molecule file");
    } else {
        molecule.read(molecule_stream);
        molecule_stream.close();
    }
    ::print("Molecule read from file: ");
    molecule.print();


    json all_input_blocks = {};

    all_input_blocks["dft"] = moldft_params.to_json_if_precedence("defined");
    all_input_blocks["response"] = molresponse_params.to_json_if_precedence("defined");
    all_input_blocks["molecule"] = molecule.to_json();



    print(all_input_blocks.dump(4));
    write_json_to_input_file(all_input_blocks, {"dft", "response"}, std::cout);
    write_molecule_json_to_input_file(all_input_blocks["molecule"], std::cout);

    std::ofstream moldft_input_file("example_moldft.in");
    write_moldft_input(all_input_blocks, moldft_input_file);
    moldft_input_file.close();

    std::ofstream response_input_file("example_response.in");
    write_response_input(all_input_blocks, response_input_file);
    response_input_file.close();

}
