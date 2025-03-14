#ifndef UTILS_HPP
#define UTILS_HPP
#include <apps/molresponse/response_parameters.h>
#include <madchem.h>
#include <madness/chem/CalculationParameters.h>
#include <madness/chem/SCF.h>
#include <madness/chem/molecule.h>
#include <algorithm>
#include <filesystem>
#include <madness/external/nlohmann_json/json.hpp>
#include <utility>

using path = std::filesystem::path;
using json = nlohmann::json;
using commandlineparser = madness::commandlineparser;

auto split(const std::string& s, char delim) -> vector<std::string> {
  vector<std::string> result;
  std::stringstream ss(s);
  std::string item;

  while (getline(ss, item, delim)) {
    result.push_back(item);
  }

  return result;
}

// Add a path to a root path
auto addPath(const path& root, const std::string& branch) -> path {
  path p_branch = root;
  p_branch += branch;
  return p_branch;
}

void write_molecule_json_to_input_file(const json& molecule_json,
                                       std::ostream& output_stream) {
  output_stream << "geometry" << std::endl;
  // first write the parameters
  auto parameters = molecule_json["parameters"];
  for (auto& [key, value] : parameters.items()) {
    output_stream << "    " << key << " " << value << std::endl;
  }

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
    output_stream << "    " << symbols[i] << " " << geometry[i][0] << " "
                  << geometry[i][1] << " " << geometry[i][2] << std::endl;
  }
  output_stream << "end" << std::endl << std::endl;
}

void write_json_to_input_file(const json& input_json,
                              const std::vector<std::string>& keys,
                              std::ostream& output_stream) {
  for (const auto& key : keys) {
    if (input_json.find(key) == input_json.end()) {
      throw std::runtime_error("Key not found in input json");
    }

    output_stream << key << std::endl;
    // for each key within the block write the key value pair to the file line
    // by line
    for (auto& [key, value] : input_json[key].items()) {
      output_stream << "    " << key << " " << value << std::endl;
    }
    output_stream << "end" << std::endl << std::endl;
  }
}

void write_moldft_input(const json& input_json, std::ostream& out) {
  write_json_to_input_file(input_json, {"dft"}, out);
  write_molecule_json_to_input_file(input_json["molecule"], out);
}

void write_response_input(const json& input_json, std::ostream& out) {
  write_json_to_input_file(input_json, {"response"}, out);
}
#endif
