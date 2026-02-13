#pragma once

#include "ResponseState.hpp"
#include "StateGenerator.hpp"
#include "../../madness/chem/ResponseParameters.hpp"

#include <madness/external/nlohmann_json/json.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <functional>
#include <iomanip>
#include <optional>
#include <set>
#include <stdexcept>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

using json = nlohmann::json;

struct DerivedStateDependency {
  std::string perturbation_description;
  double frequency = 0.0;

  [[nodiscard]] json to_json() const {
    return {{"perturbation", perturbation_description},
            {"frequency", frequency}};
  }
};

struct DerivedStateRequest {
  std::string derived_state_id;
  std::string target_perturbation;
  std::string vbc_source_id;
  double target_frequency = 0.0;
  double threshold = 0.0;
  std::string source_b_perturbation;
  std::string source_c_perturbation;
  double source_b_frequency = 0.0;
  double source_c_frequency = 0.0;
  std::vector<DerivedStateDependency> dependencies;

  [[nodiscard]] json to_json() const {
    json deps = json::array();
    for (const auto &d : dependencies)
      deps.push_back(d.to_json());

    return {{"id", derived_state_id},
            {"target_perturbation", target_perturbation},
            {"vbc_source", vbc_source_id},
            {"target_frequency", target_frequency},
            {"threshold", threshold},
            {"source_b_perturbation", source_b_perturbation},
            {"source_c_perturbation", source_c_perturbation},
            {"source_b_frequency", source_b_frequency},
            {"source_c_frequency", source_c_frequency},
            {"dependencies", deps}};
  }
};

struct DerivedStatePlan {
  std::vector<DerivedStateRequest> requests;

  [[nodiscard]] json to_json() const {
    json reqs = json::array();
    for (const auto &req : requests)
      reqs.push_back(req.to_json());
    return {{"num_requests", requests.size()}, {"requests", reqs}};
  }
};

struct DerivedStateGateEntry {
  std::string derived_state_id;
  bool ready = false;
  std::vector<DerivedStateDependency> missing_dependencies;

  [[nodiscard]] json to_json() const {
    json missing = json::array();
    for (const auto &d : missing_dependencies)
      missing.push_back(d.to_json());
    return {{"id", derived_state_id},
            {"ready", ready},
            {"missing_dependencies", missing}};
  }
};

struct DerivedStateGateReport {
  size_t total_requests = 0;
  size_t ready_requests = 0;
  size_t blocked_requests = 0;
  size_t total_missing_dependencies = 0;
  std::vector<DerivedStateGateEntry> entries;

  [[nodiscard]] bool all_ready() const { return blocked_requests == 0; }

  [[nodiscard]] json to_json() const {
    json rows = json::array();
    for (const auto &entry : entries)
      rows.push_back(entry.to_json());
    return {{"total_requests", total_requests},
            {"ready_requests", ready_requests},
            {"blocked_requests", blocked_requests},
            {"total_missing_dependencies", total_missing_dependencies},
            {"all_ready", all_ready()},
            {"entries", rows}};
  }
};

class DerivedStatePlanner {
public:
  static DerivedStatePlan
  build_vbc_driven_quadratic_plan(const ResponseParameters &response_params,
                                  const Molecule &molecule,
                                  bool spin_restricted,
                                  const std::vector<double> &protocol) {
    DerivedStatePlan plan;
    std::unordered_set<std::string> seen_ids;

    const double threshold =
        protocol.empty() ? 1.0e-6 : *std::min_element(protocol.begin(), protocol.end());

    const auto props = response_params.requested_properties();
    bool need_beta = false;
    bool need_raman = false;
    for (const auto &raw : props) {
      auto p = canonical_property(raw);
      need_beta = need_beta || (p == "hyperpolarizability");
      need_raman = need_raman || (p == "raman");
    }

    const auto dipole_freqs = response_params.dipole_frequencies();
    const auto dipole_dirs = response_params.dipole_directions();

    if (need_beta) {
      for (char a : dipole_dirs) {
        for (char b : dipole_dirs) {
          for (char c : dipole_dirs) {
            for (double freq_b : dipole_freqs) {
              for (double freq_c : dipole_freqs) {
                add_request(plan, seen_ids, DipolePerturbation{a},
                            DipolePerturbation{b}, DipolePerturbation{c},
                            freq_b, freq_c, threshold, spin_restricted);
              }
            }
          }
        }
      }
    }

    if (need_raman) {
      const auto natoms = molecule.natom();
      const auto nuc_dirs = response_params.nuclear_directions();
      for (char a : dipole_dirs) {
        for (char b : dipole_dirs) {
          for (long atom = 0; atom < natoms; ++atom) {
            for (char n : nuc_dirs) {
              for (double freq_b : dipole_freqs) {
                const double freq_c = 0.0;
                add_request(plan, seen_ids, DipolePerturbation{a},
                            DipolePerturbation{b},
                            NuclearDisplacementPerturbation{
                                static_cast<int>(atom), n},
                            freq_b,
                            freq_c, threshold, spin_restricted);
              }
            }
          }
        }
      }
    }

    return plan;
  }

  template <typename DependencyReadyFn>
  static DerivedStateGateReport
  evaluate_dependency_gate(const DerivedStatePlan &plan,
                           const GeneratedStateData &generated_state_data,
                           size_t threshold_index,
                           DependencyReadyFn &&dependency_ready) {
    DerivedStateGateReport report;
    report.total_requests = plan.requests.size();
    report.entries.reserve(plan.requests.size());

    for (const auto &req : plan.requests) {
      DerivedStateGateEntry entry;
      entry.derived_state_id = req.derived_state_id;
      entry.ready = true;

      for (const auto &dep : req.dependencies) {
        auto st_it =
            generated_state_data.state_map.find(dep.perturbation_description);
        if (st_it == generated_state_data.state_map.end()) {
          entry.missing_dependencies.push_back(dep);
          entry.ready = false;
          continue;
        }

        auto freq_index = find_frequency_index(st_it->second, dep.frequency);
        if (!freq_index.has_value()) {
          entry.missing_dependencies.push_back(dep);
          entry.ready = false;
          continue;
        }

        LinearResponsePoint pt{st_it->second, threshold_index, *freq_index};
        if (!dependency_ready(pt)) {
          entry.missing_dependencies.push_back(dep);
          entry.ready = false;
        }
      }

      if (entry.ready) {
        ++report.ready_requests;
      } else {
        ++report.blocked_requests;
        report.total_missing_dependencies += entry.missing_dependencies.size();
      }

      report.entries.push_back(std::move(entry));
    }

    return report;
  }

  static VBCResponseState make_vbc_state(const DerivedStateRequest &req,
                                         bool spin_restricted) {
    const auto b_desc =
        !req.source_b_perturbation.empty()
            ? req.source_b_perturbation
            : (req.dependencies.empty()
                   ? std::string{}
                   : req.dependencies.front().perturbation_description);
    const auto c_desc =
        !req.source_c_perturbation.empty()
            ? req.source_c_perturbation
            : (req.dependencies.size() < 2
                   ? std::string{}
                   : req.dependencies[1].perturbation_description);

    const double freq_b = !req.source_b_perturbation.empty()
                              ? req.source_b_frequency
                              : (req.dependencies.empty()
                                     ? req.source_b_frequency
                                     : req.dependencies.front().frequency);
    const double freq_c = !req.source_c_perturbation.empty()
                              ? req.source_c_frequency
                              : (req.dependencies.size() < 2
                                     ? req.source_c_frequency
                                     : req.dependencies[1].frequency);

    if (b_desc.empty() || c_desc.empty()) {
      throw std::runtime_error("Derived request " + req.derived_state_id +
                               " is missing VBC source dependencies");
    }

    return VBCResponseState(parse_perturbation(b_desc),
                            parse_perturbation(c_desc), freq_b, freq_c,
                            req.threshold, spin_restricted);
  }

private:
  static char canonical_direction(char raw) {
    return static_cast<char>(
        std::tolower(static_cast<unsigned char>(raw)));
  }

  static Perturbation parse_perturbation(const std::string &description) {
    if (description.rfind("Dipole_", 0) == 0) {
      if (description.size() < std::string("Dipole_x").size()) {
        throw std::runtime_error("Invalid dipole perturbation description: " +
                                 description);
      }
      return DipolePerturbation{canonical_direction(description.back())};
    }

    if (description.rfind("Nuc_", 0) == 0) {
      if (description.size() < std::string("Nuc_0x").size()) {
        throw std::runtime_error("Invalid nuclear perturbation description: " +
                                 description);
      }
      const std::string atom_token =
          description.substr(4, description.size() - 5);
      if (atom_token.empty()) {
        throw std::runtime_error("Missing atom index in perturbation: " +
                                 description);
      }
      const int atom_index = std::stoi(atom_token);
      return NuclearDisplacementPerturbation{
          atom_index, canonical_direction(description.back())};
    }

    if (description.rfind("Mag_", 0) == 0) {
      if (description.size() < std::string("Mag_x").size()) {
        throw std::runtime_error("Invalid magnetic perturbation description: " +
                                 description);
      }
      return MagneticPerturbation{canonical_direction(description.back())};
    }

    throw std::runtime_error("Unsupported perturbation description: " +
                             description);
  }

  static std::string canonical_property(const std::string &raw) {
    if (raw.size() >= 2 &&
        ((raw.front() == '"' && raw.back() == '"') ||
         (raw.front() == '\'' && raw.back() == '\''))) {
      return raw.substr(1, raw.size() - 2);
    }
    return raw;
  }

  static std::string make_request_id(const Perturbation &a, const Perturbation &b,
                                     const Perturbation &c, double freq_b,
                                     double freq_c) {
    std::ostringstream oss;
    oss << "Derived_" << describe_perturbation(a) << "__"
        << describe_perturbation(b) << "_" << describe_perturbation(c)
        << "_f" << std::fixed << std::setprecision(3) << freq_b << "_"
        << std::fixed << std::setprecision(3) << freq_c;
    return oss.str();
  }

  static void add_request(DerivedStatePlan &plan,
                          std::unordered_set<std::string> &seen_ids,
                          const Perturbation &a, const Perturbation &b,
                          const Perturbation &c, double freq_b, double freq_c,
                          double threshold, bool spin_restricted) {
    const std::string request_id = make_request_id(a, b, c, freq_b, freq_c);
    if (!seen_ids.insert(request_id).second) {
      return;
    }

    VBCResponseState vbc_state(b, c, freq_b, freq_c, threshold,
                               spin_restricted);
    DerivedStateRequest req;
    req.derived_state_id = request_id;
    req.target_perturbation = describe_perturbation(a);
    req.vbc_source_id = vbc_state.response_filename();
    req.target_frequency = freq_b + freq_c;
    req.threshold = threshold;
    req.source_b_perturbation = describe_perturbation(b);
    req.source_c_perturbation = describe_perturbation(c);
    req.source_b_frequency = freq_b;
    req.source_c_frequency = freq_c;
    req.dependencies = {{describe_perturbation(b), freq_b},
                        {describe_perturbation(c), freq_c}};
    plan.requests.push_back(std::move(req));
  }

  static std::optional<size_t>
  find_frequency_index(const LinearResponseDescriptor &desc,
                       double target_frequency) {
    auto it = desc.frequency_map.find(target_frequency);
    if (it != desc.frequency_map.end())
      return it->second;

    constexpr double tol = 1.0e-9;
    for (size_t i = 0; i < desc.frequencies.size(); ++i) {
      if (std::abs(desc.frequencies[i] - target_frequency) < tol) {
        return i;
      }
    }
    return std::nullopt;
  }
};
