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

/*
Developer Overview
- Purpose: Build and gate derived quadratic-state requests (currently VBC-driven)
  from response inputs before property assembly.
- Planning algorithm: Expand requested properties into (A,B,C,wb,wc) tuples,
  deduplicate by stable derived_state_id, and store dependencies on linear
  response points for B(wb) and C(wc).
- Dependency gate: For each request, verify that required linear points exist
  and are marked ready by the caller-provided predicate.
- Execution bridge: Reconstruct VBCResponseState from request metadata so solve
  code can run derived requests without re-deriving tuple logic.
*/

using json = nlohmann::json;

/// Dependency on a linear response point required before a derived request can run.
struct DerivedStateDependency {
  // Perturbation label of the linear state dependency (e.g. Dipole_x).
  std::string perturbation_description;
  // Frequency for the dependency point that must be available.
  double frequency = 0.0;

  [[nodiscard]] json to_json() const {
    return {{"perturbation", perturbation_description},
            {"frequency", frequency}};
  }
};

/// One executable derived-state task reconstructed from input/property expansion.
///
/// Big-picture role:
/// - Stage 1: produced by `build_vbc_driven_quadratic_plan(...)`.
/// - Stage 2d: dependency-gated and then turned back into `VBCResponseState`.
struct DerivedStateRequest {
  // Stable ID used in logs/metadata for this derived-state request.
  std::string derived_state_id;
  // Target perturbation for the derived solve (the "A" in A(B,C)).
  std::string target_perturbation;
  // Backing VBC response filename key produced from (B,C,wb,wc).
  std::string vbc_source_id;
  // Frequency for the target derived state (wb + wc).
  double target_frequency = 0.0;
  // Protocol threshold used when creating the derived state object.
  double threshold = 0.0;
  // Explicit source perturbations/frequencies used to reconstruct VBC state.
  std::string source_b_perturbation;
  std::string source_c_perturbation;
  double source_b_frequency = 0.0;
  double source_c_frequency = 0.0;
  // Dependency gate inputs that must be converged before this request runs.
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

/// Container for all derived requests emitted during planning.
struct DerivedStatePlan {
  // All derived requests discovered from input/property expansion.
  std::vector<DerivedStateRequest> requests;

  [[nodiscard]] json to_json() const {
    json reqs = json::array();
    for (const auto &req : requests)
      reqs.push_back(req.to_json());
    return {{"num_requests", requests.size()}, {"requests", reqs}};
  }
};

/// One row of dependency-gate output for a request.
struct DerivedStateGateEntry {
  // Mirrors DerivedStateRequest::derived_state_id.
  std::string derived_state_id;
  // True only if every dependency was found and marked ready.
  bool ready = false;
  // Missing points explaining why this request is blocked.
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

/// Aggregate dependency-gate report used by orchestration metadata and dispatch.
struct DerivedStateGateReport {
  // Summary counters for dependency-gate diagnostics.
  size_t total_requests = 0;
  size_t ready_requests = 0;
  size_t blocked_requests = 0;
  size_t total_missing_dependencies = 0;
  // Per-request gate status aligned with plan.requests order.
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

/// Planner for VBC-driven derived response requests (beta/Raman pathways).
///
/// This class keeps tuple expansion logic (A,B,C,wb,wc), deduplication, and
/// dependency gating out of the stage-2 solver loop so execution can be
/// restart-aware and inspectable through metadata.
class DerivedStatePlanner {
public:
  /// Expand response settings into derived requests.
  ///
  /// Generates one request per unique `(A, B, C, wb, wc)` tuple required by
  /// current property settings, after applying beta triplet policy and
  /// frequency canonicalization.
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

    auto dipole_freqs = response_params.dipole_frequencies();
    for (double &freq : dipole_freqs) {
      freq = canonicalize_response_frequency(freq);
    }
    std::sort(dipole_freqs.begin(), dipole_freqs.end());
    dipole_freqs.erase(std::unique(dipole_freqs.begin(), dipole_freqs.end()),
                       dipole_freqs.end());
    const auto dipole_dirs = response_params.dipole_directions();

    if (need_beta) {
      const auto beta_frequency_pairs = select_beta_frequency_pairs(
          dipole_freqs, response_params.beta_shg(), response_params.beta_or(),
          response_params.beta_all_triplets());
      for (char a : dipole_dirs) {
        for (char b : dipole_dirs) {
          for (char c : dipole_dirs) {
            for (const auto &[freq_b, freq_c] : beta_frequency_pairs) {
              add_request(plan, seen_ids, DipolePerturbation{a},
                          DipolePerturbation{b}, DipolePerturbation{c}, freq_b,
                          freq_c, threshold, spin_restricted);
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

  /// Evaluate whether each request is ready at a specific protocol index.
  ///
  /// `dependency_ready(pt)` is provided by the caller (usually metadata-backed)
  /// and decides if each linear prerequisite point is currently usable.
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

  /// Reconstruct a runtime `VBCResponseState` from persisted request metadata.
  ///
  /// This is the bridge between planning data and stage-2 derived execution.
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
  /// Normalize axis labels (`X/Y/Z`) to canonical lower-case.
  static char canonical_direction(char raw) {
    return static_cast<char>(
        std::tolower(static_cast<unsigned char>(raw)));
  }

  /// Parse string perturbation ids (e.g. `Dipole_x`, `Nuc_0z`) into typed variants.
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

  /// Remove optional quote wrappers from property names.
  static std::string canonical_property(const std::string &raw) {
    if (raw.size() >= 2 &&
        ((raw.front() == '"' && raw.back() == '"') ||
         (raw.front() == '\'' && raw.back() == '\''))) {
      return raw.substr(1, raw.size() - 2);
    }
    return raw;
  }

  /// Select beta frequency pairs from policy flags.
  ///
  /// Precedence is:
  /// 1. `beta_all_triplets`: full cartesian `(wb,wc)`
  /// 2. `beta_or`: `(0,wc)`
  /// 3. `beta_shg`: `(w,w)`
  static std::vector<std::pair<double, double>> select_beta_frequency_pairs(
      const std::vector<double> &frequencies, bool beta_shg, bool beta_or,
      bool beta_all_triplets) {
    std::vector<std::pair<double, double>> pairs;
    if (beta_all_triplets) {
      pairs.reserve(frequencies.size() * frequencies.size());
      for (const double freq_b : frequencies) {
        for (const double freq_c : frequencies) {
          pairs.emplace_back(freq_b, freq_c);
        }
      }
      return pairs;
    }

    if (beta_or) {
      pairs.reserve(frequencies.size());
      for (const double freq_c : frequencies) {
        pairs.emplace_back(0.0, freq_c);
      }
      return pairs;
    }

    if (beta_shg) {
      pairs.reserve(frequencies.size());
      for (const double freq : frequencies) {
        pairs.emplace_back(freq, freq);
      }
      return pairs;
    }

    throw std::runtime_error(
        "No beta frequency-triplet mode enabled. Set one of: beta.shg, "
        "beta.or, beta.all_triplets.");
  }

  /// Build a stable identifier used for dedupe and metadata keys.
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

  /// Canonicalize and append one derived request if it has not been seen before.
  static void add_request(DerivedStatePlan &plan,
                          std::unordered_set<std::string> &seen_ids,
                          const Perturbation &a, const Perturbation &b,
                          const Perturbation &c, double freq_b, double freq_c,
                          double threshold, bool spin_restricted) {
    freq_b = canonicalize_response_frequency(freq_b);
    freq_c = canonicalize_response_frequency(freq_c);
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

  /// Find a frequency index in a descriptor with map-first and tolerance fallback.
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
