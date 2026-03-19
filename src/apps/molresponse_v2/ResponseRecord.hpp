#pragma once
#include "ResponseState.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <madness/external/nlohmann_json/json.hpp>
#include <map>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

// ============================================================================
// Excited-state identity and I/O types
//
// These live here because ResponseRecord2 uses them directly for persisting
// per-root metadata. ExcitedResponse.hpp includes this header for the same
// types.
//
// Terminology: a collection of M excited-state response vectors is called a
// "response space" (not "bundle") because we diagonalise to find eigenvectors
// of that space.
// ============================================================================

struct ExcitedRootDescriptor {
    std::string root_id;       // stable string key, e.g. "es_root_0001"
    size_t      stable_index = 0;  // permanent integer identity assigned at first appearance
    size_t      slot_index   = 0;  // current position in the response space vector
    double      energy       = 0.0;
    std::string display_name;
};

inline void to_json(nlohmann::json &j, const ExcitedRootDescriptor &root) {
    j = nlohmann::json{{"root_id",      root.root_id},
                       {"stable_index", root.stable_index},
                       {"slot_index",   root.slot_index},
                       {"energy",       root.energy},
                       {"display_name", root.display_name}};
}

inline void from_json(const nlohmann::json &j, ExcitedRootDescriptor &root) {
    root.root_id      = j.value("root_id", std::string{});
    root.stable_index = j.value("stable_index", static_cast<size_t>(0));
    root.slot_index   = j.value("slot_index",
                                j.value("root_index", root.stable_index));
    root.energy       = j.value("energy", 0.0);
    root.display_name = j.value("display_name",
                                j.value("name", std::string{}));
}

// Input to ExcitedResponse::solve_protocol for one protocol level.
struct ExcitedProtocolInput {
    double threshold         = 0.0;
    double dconv             = 0.0;
    size_t protocol_index    = 0;
    bool   restart_saved     = false;
    bool   restart_converged = false;
    size_t owner_group       = 0;
    bool   tda               = false;
    size_t num_states        = 1;
    size_t guess_max_iter    = 5;
    size_t maxiter           = 20;
    size_t maxsub            = 8;
};

// Result from ExcitedResponse::solve_protocol for one protocol level.
struct ExcitedProtocolResult {
    bool        attempted                   = false;
    bool        saved                       = false;
    bool        converged                   = false;
    bool        failed                      = false;
    bool        skipped                     = false;
    bool        restart_reused              = false;
    std::string stage_status                = "pending";
    std::string response_variant            = "unknown";
    std::string restart_support_mode        = "guess_only";
    std::string restart_source              = "none";
    std::string snapshot_kind               = "none";
    bool        bundle_state_present        = false;  // TODO: rename response_space_present
    bool        restart_capable             = false;
    double      restart_source_threshold    = 0.0;
    size_t      iterations                  = 0;
    std::vector<double>                  energies;
    std::vector<std::string>             state_names;
    std::vector<ExcitedRootDescriptor>   roots;
    std::vector<size_t>                  slot_permutation;
    std::vector<double>                  residual_norms;
    std::vector<double>                  density_change_norms;
    std::vector<double>                  relative_residual_norms;
    std::vector<double>                  iteration_max_residuals;
    std::vector<double>                  iteration_max_density_changes;
    std::vector<double>                  iteration_max_relative_residuals;
    std::string convergence_mode            = "max_residual";
    std::string accelerator_mode            = "none";
    size_t      accelerator_subspace        = 0;
    double      density_convergence_target  = 0.0;
    double      relative_convergence_target = 0.0;
    double      max_rotation                = 0.0;
};

inline void to_json(nlohmann::json &j, const ExcitedProtocolResult &r) {
    j = nlohmann::json{
        {"attempted",                        r.attempted},
        {"saved",                            r.saved},
        {"converged",                        r.converged},
        {"failed",                           r.failed},
        {"skipped",                          r.skipped},
        {"restart_reused",                   r.restart_reused},
        {"stage_status",                     r.stage_status},
        {"response_variant",                 r.response_variant},
        {"restart_support_mode",             r.restart_support_mode},
        {"restart_source",                   r.restart_source},
        {"snapshot_kind",                    r.snapshot_kind},
        {"bundle_state_present",             r.bundle_state_present},
        {"restart_capable",                  r.restart_capable},
        {"restart_source_threshold",         r.restart_source_threshold},
        {"iterations",                       r.iterations},
        {"energies",                         r.energies},
        {"state_names",                      r.state_names},
        {"roots",                            r.roots},
        {"slot_permutation",                 r.slot_permutation},
        {"residual_norms",                   r.residual_norms},
        {"density_change_norms",             r.density_change_norms},
        {"relative_residual_norms",          r.relative_residual_norms},
        {"iteration_max_residuals",          r.iteration_max_residuals},
        {"iteration_max_density_changes",    r.iteration_max_density_changes},
        {"iteration_max_relative_residuals", r.iteration_max_relative_residuals},
        {"convergence_mode",                 r.convergence_mode},
        {"accelerator_mode",                 r.accelerator_mode},
        {"accelerator_subspace",             r.accelerator_subspace},
        {"density_convergence_target",       r.density_convergence_target},
        {"relative_convergence_target",      r.relative_convergence_target},
        {"max_rotation",                     r.max_rotation}};
}

inline void from_json(const nlohmann::json &j, ExcitedProtocolResult &r) {
    r.attempted              = j.value("attempted",              false);
    r.saved                  = j.value("saved",                  false);
    r.converged              = j.value("converged",              false);
    r.failed                 = j.value("failed",                 false);
    r.skipped                = j.value("skipped",                false);
    r.restart_reused         = j.value("restart_reused",         false);
    r.stage_status           = j.value("stage_status",           std::string("pending"));
    r.response_variant       = j.value("response_variant",       std::string("unknown"));
    r.restart_support_mode   = j.value("restart_support_mode",   std::string("guess_only"));
    r.restart_source         = j.value("restart_source",         std::string("none"));
    r.snapshot_kind          = j.value("snapshot_kind",          std::string("none"));
    r.bundle_state_present   = j.value("bundle_state_present",   false);
    r.restart_capable        = j.value("restart_capable",        false);
    r.restart_source_threshold = j.value("restart_source_threshold", 0.0);
    r.iterations             = j.value("iterations",             static_cast<size_t>(0));
    if (j.contains("energies") && j["energies"].is_array())
        r.energies = j["energies"].get<std::vector<double>>();
    if (j.contains("state_names") && j["state_names"].is_array())
        r.state_names = j["state_names"].get<std::vector<std::string>>();
    if (j.contains("roots") && j["roots"].is_array())
        r.roots = j["roots"].get<std::vector<ExcitedRootDescriptor>>();
    if (j.contains("slot_permutation") && j["slot_permutation"].is_array())
        r.slot_permutation = j["slot_permutation"].get<std::vector<size_t>>();
    if (j.contains("residual_norms") && j["residual_norms"].is_array())
        r.residual_norms = j["residual_norms"].get<std::vector<double>>();
    if (j.contains("density_change_norms") && j["density_change_norms"].is_array())
        r.density_change_norms = j["density_change_norms"].get<std::vector<double>>();
    if (j.contains("relative_residual_norms") && j["relative_residual_norms"].is_array())
        r.relative_residual_norms = j["relative_residual_norms"].get<std::vector<double>>();
    if (j.contains("iteration_max_residuals") && j["iteration_max_residuals"].is_array())
        r.iteration_max_residuals = j["iteration_max_residuals"].get<std::vector<double>>();
    if (j.contains("iteration_max_density_changes") && j["iteration_max_density_changes"].is_array())
        r.iteration_max_density_changes = j["iteration_max_density_changes"].get<std::vector<double>>();
    if (j.contains("iteration_max_relative_residuals") && j["iteration_max_relative_residuals"].is_array())
        r.iteration_max_relative_residuals = j["iteration_max_relative_residuals"].get<std::vector<double>>();
    r.convergence_mode           = j.value("convergence_mode",           std::string("max_residual"));
    r.accelerator_mode           = j.value("accelerator_mode",           std::string("none"));
    r.accelerator_subspace       = j.value("accelerator_subspace",       static_cast<size_t>(0));
    r.density_convergence_target = j.value("density_convergence_target", 0.0);
    r.relative_convergence_target= j.value("relative_convergence_target",0.0);
    r.max_rotation               = j.value("max_rotation",               0.0);
}

using json = nlohmann::json;
namespace fs = std::filesystem;

class ResponseRecord2 {
public:
  using json = nlohmann::json;

  struct MissingItem {
    std::string state;    // perturbationDescription()
    std::string freq;     // canonical "0.500"
    std::string protocol; // canonical "1e-06"
    bool saved_found = false;
    bool converged_found = false;
  };

  ResponseRecord2(World &world, const std::string &filepath)
      : world_(world), path_(filepath) {
    std::string json_string;
    if (world_.rank() == 0) {
      if (fs::exists(path_)) {
        std::ifstream in(path_);
        if (in) {
          std::stringstream buf;
          buf << in.rdbuf();
          json_string = buf.str();
        } else {
          std::cerr << "Error opening file: " << path_ << std::endl;
          json_string = "{}";
        }
      } else {
        json_string = "{}";
      }
    }
    world_.gop.fence();
    world_.gop.broadcast_serializable(json_string, 0);
    world_.gop.fence();

    if (json_string.empty())
      data_ = json::object();
    else
      data_ = json::parse(json_string, nullptr, /*allow_exceptions=*/true);

    ensure_root(data_);
  }

  // --- Initialization for a generated set of states ---
  void initialize_states(const std::vector<LinearResponseDescriptor> &states) {
    for (const auto &st : states) {
      const std::string state_id = st.perturbationDescription();
      ensure_state(data_, state_id);

      // NB: canonicalize keys up front
      for (double thr : st.thresholds) {
        const std::string pkey = protocol_key(thr);
        ensure_protocol(data_, state_id, pkey);

        for (double f : st.frequencies) {
          const std::string fkey = freq_key(f);
          auto &node = data_["states"][state_id]["protocols"][pkey];
          if (!node["saved"].contains(fkey))
            node["saved"][fkey] = false;
          if (!node["converged"].contains(fkey))
            node["converged"][fkey] = false;
        }
      }
    }
    write(); // sync to disk + broadcast
  }

  void initialize_excited_bundle(bool enabled, size_t num_states, bool tda,
                                 size_t guess_max_iter, size_t maxiter,
                                 size_t maxsub, size_t owner_group,
                                 const std::vector<double> &protocols) {
    ensure_excited_root(data_);
    auto &plan = data_["excited_states"]["plan"];
    if (!plan.is_object()) {
      plan = json::object();
    }
    plan["enabled"] = enabled;
    plan["num_states"] = num_states;
    plan["tda"] = tda;
    plan["guess_max_iter"] = guess_max_iter;
    plan["maxiter"] = maxiter;
    plan["maxsub"] = maxsub;
    plan["owner_group"] = owner_group;
    plan["protocols"] = json::array();

    for (double protocol : protocols) {
      const std::string pkey = protocol_key(protocol);
      plan["protocols"].push_back(pkey);
      ensure_excited_protocol(data_, pkey);
    }
    write();
  }

  void record_excited_protocol_status(double protocol, bool saved,
                                      bool converged) {
    const std::string pkey = protocol_key(protocol);
    ensure_excited_protocol(data_, pkey);
    auto &node = data_["excited_states"]["protocols"][pkey];
    node["saved"] = saved;
    node["converged"] = converged;
    write();
  }

  void record_excited_protocol_timing(double protocol, double wall_seconds,
                                      double cpu_seconds) {
    const std::string pkey = protocol_key(protocol);
    ensure_excited_protocol(data_, pkey);
    auto &node = data_["excited_states"]["protocols"][pkey];
    node["timings"] = {{"wall_seconds", wall_seconds},
                       {"cpu_seconds", cpu_seconds}};
    write();
  }

  void record_excited_protocol_energies(double protocol,
                                        const std::vector<double> &energies) {
    const std::string pkey = protocol_key(protocol);
    ensure_excited_protocol(data_, pkey);
    auto &node = data_["excited_states"]["protocols"][pkey];
    node["energies"] = energies;
    write();
  }

  [[nodiscard]] static ExcitedProtocolResult
  normalize_excited_protocol_result(double protocol,
                                    ExcitedProtocolResult result) {
    const size_t root_count = std::max(
        result.energies.size(),
        std::max(result.state_names.size(), result.roots.size()));
    if (root_count == 0) {
      result.state_names.clear();
      result.roots.clear();
      result.slot_permutation.clear();
      return result;
    }

    if (result.roots.size() < root_count) {
      size_t next_stable_index = 0;
      for (const auto &root : result.roots) {
        next_stable_index = std::max(next_stable_index, root.stable_index + 1);
      }
      result.roots.reserve(root_count);
      for (size_t i = result.roots.size(); i < root_count; ++i) {
        ExcitedRootDescriptor root;
        root.stable_index = next_stable_index++;
        root.root_id = make_excited_root_id(root.stable_index);
        result.roots.push_back(std::move(root));
      }
    } else if (result.roots.size() > root_count) {
      result.roots.resize(root_count);
    }

    const size_t original_energy_count = result.energies.size();
    if (result.energies.size() < root_count) {
      result.energies.resize(root_count, 0.0);
    }

    size_t next_stable_index = 0;
    for (const auto &root : result.roots) {
      next_stable_index = std::max(next_stable_index, root.stable_index + 1);
    }

    std::vector<size_t> used_stable_indices;
    used_stable_indices.reserve(result.roots.size());
    for (auto &root : result.roots) {
      const bool duplicate_stable_index =
          std::find(used_stable_indices.begin(), used_stable_indices.end(),
                    root.stable_index) != used_stable_indices.end();
      if (duplicate_stable_index) {
        root.stable_index = next_stable_index++;
      }
      if (root.root_id.empty()) {
        root.root_id = make_excited_root_id(root.stable_index);
      }
      used_stable_indices.push_back(root.stable_index);
    }

    std::sort(result.roots.begin(), result.roots.end(),
              [](const ExcitedRootDescriptor &lhs,
                 const ExcitedRootDescriptor &rhs) {
                if (lhs.stable_index != rhs.stable_index) {
                  return lhs.stable_index < rhs.stable_index;
                }
                return lhs.root_id < rhs.root_id;
              });

    std::vector<bool> slot_used(root_count, false);
    for (auto &root : result.roots) {
      if (root.slot_index >= root_count || slot_used[root.slot_index]) {
        root.slot_index = root_count;
        continue;
      }
      slot_used[root.slot_index] = true;
    }
    size_t next_slot = 0;
    for (auto &root : result.roots) {
      if (root.slot_index < root_count) {
        continue;
      }
      while (next_slot < root_count && slot_used[next_slot]) {
        ++next_slot;
      }
      root.slot_index = (next_slot < root_count) ? next_slot : 0;
      if (next_slot < root_count) {
        slot_used[next_slot] = true;
      }
    }

    bool state_names_complete = result.state_names.size() == root_count;
    if (state_names_complete) {
      for (const auto &name : result.state_names) {
        if (name.empty()) {
          state_names_complete = false;
          break;
        }
      }
    }
    std::vector<std::string> slot_names =
        state_names_complete ? result.state_names
                             : assign_excited_state_names(result.energies, protocol);

    for (const auto &root : result.roots) {
      if (root.slot_index >= root_count) {
        continue;
      }
      if (root.slot_index >= original_energy_count) {
        result.energies[root.slot_index] = root.energy;
      }
      if (!root.display_name.empty()) {
        slot_names[root.slot_index] = root.display_name;
      }
    }

    result.slot_permutation.assign(root_count, 0);
    for (auto &root : result.roots) {
      if (root.slot_index >= root_count) {
        continue;
      }
      root.energy = result.energies[root.slot_index];
      if (root.display_name.empty()) {
        root.display_name = slot_names[root.slot_index];
      }
      slot_names[root.slot_index] = root.display_name;
      result.slot_permutation[root.slot_index] = root.stable_index;
    }

    result.state_names = std::move(slot_names);
    return result;
  }

  void record_excited_protocol_result(double protocol, size_t owner_group,
                                      const ExcitedProtocolResult &result,
                                      double wall_seconds,
                                      double cpu_seconds) {
    const auto normalized =
        normalize_excited_protocol_result(protocol, result);
    const std::string pkey = protocol_key(protocol);
    ensure_excited_protocol(data_, pkey);
    auto &node = data_["excited_states"]["protocols"][pkey];
    node["saved"] = normalized.saved;
    node["converged"] = normalized.converged;
    node["attempted"] = normalized.attempted;
    node["failed"] = normalized.failed;
    node["skipped"] = normalized.skipped;
    node["restart_reused"] = normalized.restart_reused;
    node["stage_status"] = normalized.stage_status;
    node["response_variant"] = normalized.response_variant;
    node["restart_support_mode"] = normalized.restart_support_mode;
    node["restart_source"] = normalized.restart_source;
    node["snapshot_kind"] = normalized.snapshot_kind;
    node["bundle_state_present"] = normalized.bundle_state_present;
    node["restart_capable"] = normalized.restart_capable;
    node["restart_source_threshold"] = normalized.restart_source_threshold;
    node["owner_group"] = owner_group;
    node["iterations"] = normalized.iterations;
    node["timings"] = {{"wall_seconds", wall_seconds},
                       {"cpu_seconds", cpu_seconds}};
    node["energies"] = normalized.energies;
    node["state_names"] = normalized.state_names;
    node["roots"] = normalized.roots;
    node["slot_permutation"] = normalized.slot_permutation;
    node["residual_norms"] = normalized.residual_norms;
    node["density_change_norms"] = normalized.density_change_norms;
    node["relative_residual_norms"] = normalized.relative_residual_norms;
    node["iteration_max_residuals"] = normalized.iteration_max_residuals;
    node["iteration_max_density_changes"] =
        normalized.iteration_max_density_changes;
    node["iteration_max_relative_residuals"] =
        normalized.iteration_max_relative_residuals;
    node["convergence_mode"] = normalized.convergence_mode;
    node["accelerator_mode"] = normalized.accelerator_mode;
    node["accelerator_subspace"] = normalized.accelerator_subspace;
    node["density_convergence_target"] =
        normalized.density_convergence_target;
    node["relative_convergence_target"] =
        normalized.relative_convergence_target;
    node["max_rotation"] = normalized.max_rotation;
    write();
  }

  // --- Pretty table ---
  void print_summary(int proto_digits = 0, int freq_decimals = 3) const {
    using std::cout;
    using std::left;
    using std::setw;
    constexpr int W_ROW = 5, W_STATE = 32, W_PROTO = 12, W_FREQ = 10,
                  W_SAVED = 7, W_CONV = 10;

    cout << " Response State Summary\n";
    cout << setw(W_ROW) << "#" << "  " << setw(W_STATE) << left << "State"
         << setw(W_PROTO) << "Protocol" << setw(W_FREQ) << "Freq"
         << setw(W_SAVED) << "Saved" << setw(W_CONV) << "Converged" << "\n";
    cout << std::string(
                W_ROW + 2 + W_STATE + W_PROTO + W_FREQ + W_SAVED + W_CONV, '-')
         << "\n";

    size_t row = 0;
    if (!data_.contains("states") || !data_["states"].is_object()) {
      cout << "(no states)\n";
      return;
    }

    for (const auto &kv : data_["states"].items()) {
      const std::string &state_id = kv.key();
      const auto &entry = kv.value();
      const auto &protos =
          (entry.contains("protocols") && entry["protocols"].is_object())
              ? entry["protocols"]
              : json::object();
      if (protos.empty()) {
        cout << setw(W_ROW) << row++ << "  " << setw(W_STATE) << left
             << state_id << setw(W_PROTO) << "-" << setw(W_FREQ) << "-"
             << setw(W_SAVED) << "-" << setw(W_CONV) << "-" << "\n";
        continue;
      }

      auto proto_keys = numeric_keys(protos); // (double, string) sorted asc
      for (const auto &[pnum, pkey] : proto_keys) {
        const auto &node = protos.at(pkey);
        const auto &saved =
            (node.contains("saved") && node["saved"].is_object())
                ? node["saved"]
                : json::object();
        const auto &conv =
            (node.contains("converged") && node["converged"].is_object())
                ? node["converged"]
                : json::object();

        json union_obj = json::object();
        for (const auto &it : saved.items())
          union_obj[it.key()] = true;
        for (const auto &it : conv.items())
          union_obj[it.key()] = true;

        auto freq_keys = numeric_keys(union_obj); // numeric sort
        if (freq_keys.empty()) {
          cout << setw(W_ROW) << row++ << "  " << setw(W_STATE) << left
               << state_id << setw(W_PROTO) << fmt_sci(pnum, proto_digits)
               << setw(W_FREQ) << "-" << setw(W_SAVED) << "-" << setw(W_CONV)
               << "-" << "\n";
          continue;
        }

        for (const auto &[fnum, fkey] : freq_keys) {
          const bool s =
              saved.contains(fkey) ? saved.at(fkey).get<bool>() : false;
          const bool c =
              conv.contains(fkey) ? conv.at(fkey).get<bool>() : false;
          cout << setw(W_ROW) << row++ << "  " << setw(W_STATE) << left
               << state_id << setw(W_PROTO) << fmt_sci(pnum, proto_digits)
               << setw(W_FREQ) << fmt_fixed(fnum, freq_decimals)
               << setw(W_SAVED) << (s ? "yes" : "no") << setw(W_CONV)
               << (c ? "yes" : "no") << "\n";
        }
      }
    }
  }

  // --- Queries (overloads take doubles; we canonicalize to keys) ---
  [[nodiscard]] bool is_saved(const std::string &state_id, double protocol,
                double freq) const {
    return get_flag(state_id, protocol_key(protocol), freq_key(freq), "saved");
  }
  [[nodiscard]] bool is_converged(const std::string &state_id, double protocol,
                    double freq) const {
    return get_flag(state_id, protocol_key(protocol), freq_key(freq),
                    "converged");
  }
  [[nodiscard]] bool
  is_marked_for_frequency_removal(const std::string &state_id, double protocol,
                                  double freq) const {
    const std::string pkey = protocol_key(protocol);
    const std::string fkey = freq_key(freq);
    if (!data_.contains("states") || !data_["states"].is_object()) {
      return false;
    }
    const auto sit = data_["states"].find(state_id);
    if (sit == data_["states"].end() || !(*sit).is_object() ||
        !(*sit).contains("protocols") || !(*sit)["protocols"].is_object()) {
      return false;
    }
    const auto pit = (*sit)["protocols"].find(pkey);
    if (pit == (*sit)["protocols"].end() || !pit->is_object() ||
        !pit->contains("solver_diagnostics") ||
        !(*pit)["solver_diagnostics"].is_object()) {
      return false;
    }
    const auto dit = (*pit)["solver_diagnostics"].find(fkey);
    if (dit == (*pit)["solver_diagnostics"].end() || !dit->is_object() ||
        !dit->contains("remove_from_frequency_set")) {
      return false;
    }
    return (*dit)["remove_from_frequency_set"].is_boolean() &&
           (*dit)["remove_from_frequency_set"].get<bool>();
  }
  // --- Mutations ---
  void mark_saved(const std::string &state_id, double protocol, double freq) {
    set_flag(state_id, protocol_key(protocol), freq_key(freq), "saved", true);
  }
  void mark_converged(const std::string &state_id, double protocol,
                      double freq, bool c) {
    set_flag(state_id, protocol_key(protocol), freq_key(freq), "converged", c);
  }
  void record_status(const LinearResponsePoint &pt, bool c) {
    const std::string sid = pt.perturbationDescription();
    const std::string p = protocol_key(pt.threshold());
    const std::string f = freq_key(pt.frequency());
    ensure_protocol(data_, sid, p);
    data_["states"][sid]["protocols"][p]["saved"][f] = true;
    data_["states"][sid]["protocols"][p]["converged"][f] = c;
    write();
  }

  void record_timing(const LinearResponsePoint &pt, double wall_seconds,
                     double cpu_seconds) {
    const std::string sid = pt.perturbationDescription();
    const std::string p = protocol_key(pt.threshold());
    const std::string f = freq_key(pt.frequency());
    ensure_protocol(data_, sid, p);
    auto &node = data_["states"][sid]["protocols"][p];
    if (!node.contains("timings") || !node["timings"].is_object()) {
      node["timings"] = json::object();
    }
    node["timings"][f] = {{"wall_seconds", wall_seconds},
                          {"cpu_seconds", cpu_seconds}};
    write();
  }

  void record_solver_diagnostics(const LinearResponsePoint &pt, bool converged,
                                 size_t iterations_performed,
                                 double final_residual_norm,
                                 double final_density_change,
                                 bool used_fallback_retry,
                                 double final_alpha =
                                     std::numeric_limits<double>::quiet_NaN(),
                                 size_t max_consecutive_negative_alpha = 0,
                                 bool reached_iteration_limit = false,
                                 bool remove_from_frequency_set = false,
                                 double residual_remove_cutoff =
                                     std::numeric_limits<double>::quiet_NaN(),
                                 const std::string &failure_reason =
                                     "not_evaluated") {
    const std::string sid = pt.perturbationDescription();
    const std::string p = protocol_key(pt.threshold());
    const std::string f = freq_key(pt.frequency());
    ensure_protocol(data_, sid, p);
    auto &node = data_["states"][sid]["protocols"][p];
    if (!node.contains("solver_diagnostics") ||
        !node["solver_diagnostics"].is_object()) {
      node["solver_diagnostics"] = json::object();
    }
    node["solver_diagnostics"][f] = {
        {"converged", converged},
        {"iterations_performed", iterations_performed},
        {"final_residual_norm", final_residual_norm},
        {"final_density_change", final_density_change},
        {"used_fallback_retry", used_fallback_retry},
        {"final_alpha", final_alpha},
        {"max_consecutive_negative_alpha", max_consecutive_negative_alpha},
        {"reached_iteration_limit", reached_iteration_limit},
        {"remove_from_frequency_set", remove_from_frequency_set},
        {"residual_remove_cutoff", residual_remove_cutoff},
        {"failure_reason", failure_reason}};
    write();
  }

  void record_restart_provenance(
      const LinearResponsePoint &pt, const std::string &source_kind,
      bool loaded_from_disk, bool promoted_from_static,
      const std::optional<double> &source_protocol,
      const std::optional<double> &source_frequency) {
    const std::string sid = pt.perturbationDescription();
    const std::string p = protocol_key(pt.threshold());
    const std::string f = freq_key(pt.frequency());
    ensure_protocol(data_, sid, p);
    auto &node = data_["states"][sid]["protocols"][p];
    if (!node.contains("restart_provenance") ||
        !node["restart_provenance"].is_object()) {
      node["restart_provenance"] = json::object();
    }
    json provenance = {{"kind", source_kind},
                       {"loaded_from_disk", loaded_from_disk},
                       {"promoted_from_static", promoted_from_static},
                       {"source_protocol", nullptr},
                       {"source_frequency", nullptr}};
    if (source_protocol.has_value()) {
      provenance["source_protocol"] = protocol_key(*source_protocol);
    }
    if (source_frequency.has_value()) {
      provenance["source_frequency"] = freq_key(*source_frequency);
    }
    node["restart_provenance"][f] = std::move(provenance);
  }

  // --- Property gating helpers ---
  static std::string
  final_protocol_key_from(const std::vector<double> &protos) {
    if (protos.empty())
      return "inf"; // sentinel; will never match
    return protocol_key(*std::min_element(protos.begin(), protos.end()));
  }

  // If converged at multiple protocols, return the most accurate (smallest
  // numeric) protocol key.
  [[nodiscard]] std::optional<std::string>
  best_converged_protocol(const std::string &state_id, double freq) const {
    if (!data_.contains("states"))
      return std::nullopt;
    const auto sit = data_["states"].find(state_id);
    if (sit == data_["states"].end())
      return std::nullopt;

    const auto &protos = (*sit)["protocols"];
    std::vector<std::string> pkeys;
    pkeys.reserve(protos.size());
    for (auto it = protos.begin(); it != protos.end(); ++it)
      pkeys.push_back(it.key());
    std::sort(pkeys.begin(), pkeys.end(),
              [](const std::string &a, const std::string &b) {
                return protocol_numeric(a) < protocol_numeric(b);
              });

    const std::string fk = freq_key(freq);
    for (const auto &p : pkeys) {
      const auto &node = protos.at(p);
      if (node.contains("converged")) {
        const auto &conv = node["converged"];
        if (conv.contains(fk) && conv.at(fk).get<bool>())
          return p;
      }
    }
    return std::nullopt;
  }

  // Check saved && converged at *final* protocol for each (state, freq)
  [[nodiscard]] std::vector<MissingItem>
  missing_at_final_protocol(const std::vector<std::string> &state_ids,
                            const std::vector<double> &freqs,
                            const std::string &final_proto_key) const {
    std::vector<MissingItem> out;
    for (const auto &sid : state_ids) {
      for (double f : freqs) {
        const std::string fk = freq_key(f);
        bool s = false, c = false;
        if (data_.contains("states")) {
          auto sit = data_["states"].find(sid);
          if (sit != data_["states"].end()) {
            auto pit = (*sit)["protocols"].find(final_proto_key);
            if (pit != (*sit)["protocols"].end()) {
              const auto &node = (*sit)["protocols"].at(final_proto_key);
              if (node.contains("saved") && node["saved"].contains(fk))
                s = node["saved"].at(fk).get<bool>();
              if (node.contains("converged") && node["converged"].contains(fk))
                c = node["converged"].at(fk).get<bool>();
            }
          }
        }
        if (!(s && c))
          out.push_back({sid, fk, final_proto_key, s, c});
      }
    }
    return out;
  }

  // Throw if not ready; call this right before property computations.
  void enforce_ready_for_properties(const std::vector<std::string> &state_ids,
                                    const std::vector<double> &freqs,
                                    const std::string &final_proto_key) const {
    auto missing = missing_at_final_protocol(state_ids, freqs, final_proto_key);
    if (!missing.empty()) {
      if (world_.rank() == 0) {
        std::ostringstream msg;
        msg << "Property gate failed at final protocol " << final_proto_key
            << ":\n";
        for (const auto &m : missing) {
          msg << "  - " << m.state << " @ " << m.freq
              << " (saved=" << (m.saved_found ? "true" : "false")
              << ", converged=" << (m.converged_found ? "true" : "false")
              << ")\n";
        }
        std::cerr << msg.str();
      }
      throw std::runtime_error(
          "Required states/frequencies not ready at final protocol.");
    }
  }

  // --- File I/O (rank 0 writes; all ranks sync) ---
  void write() {
    std::string json_string;
    if (world_.rank() == 0) {
      std::ofstream out(path_);
      if (!out)
        throw std::runtime_error("Cannot open " + path_ + " for writing");
      out << std::setw(2) << data_ << "\n";
      json_string = data_.dump();
    }
    world_.gop.fence();
    world_.gop.broadcast_serializable(json_string, 0);
    world_.gop.fence();
    if (world_.rank() != 0) {
      data_ = json::parse(json_string);
    }
  }

  [[nodiscard]] json to_json() const { return data_; }

  // --- Optional: final_saved compatibility flag you had ---
  void mark_final_saved(const std::string &state_id, bool flag = true) {
    ensure_state(data_, state_id);
    data_["states"][state_id]["final_saved"] = flag;
    write();
  }

  // --- Flat rows (useful for CSV/logs) ---
  struct Row {
    std::string state, freq, protocol;
    bool saved = false, converged = false;
  };
  [[nodiscard]] std::vector<Row> to_rows() const {
    std::vector<Row> rows;
    if (!data_.contains("states"))
      return rows;
    for (auto sit = data_["states"].begin(); sit != data_["states"].end();
         ++sit) {
      const std::string& sid = sit.key();
      if (!(*sit).contains("protocols"))
        continue;
      const auto &protos = (*sit)["protocols"];
      for (auto pit = protos.begin(); pit != protos.end(); ++pit) {
        const std::string& pkey = pit.key();
        const auto &node = pit.value();
        const auto &saved = node.value("saved", json::object());
        const auto &conv = node.value("converged", json::object());
        std::map<std::string, bool> k;
        for (auto it = saved.begin(); it != saved.end(); ++it)
          k[it.key()] = true;
        for (auto it = conv.begin(); it != conv.end(); ++it)
          k[it.key()] = true;
        for (const auto &kv : k) {
          const auto &fk = kv.first;
          bool s = saved.contains(fk) ? saved.at(fk).get<bool>() : false;
          bool c = conv.contains(fk) ? conv.at(fk).get<bool>() : false;
          rows.push_back({sid, fk, pkey, s, c});
        }
      }
    }
    std::sort(rows.begin(), rows.end(), [](const Row &a, const Row &b) {
      if (a.state != b.state)
        return a.state < b.state;
      double af = std::stod(a.freq), bf = std::stod(b.freq);
      if (af != bf)
        return af < bf;
      return protocol_numeric(a.protocol) < protocol_numeric(b.protocol);
    });
    return rows;
  }
  // -------- Canonicalization & shape --------
  static std::string freq_key(double f) {
    std::ostringstream os;
    os << std::fixed << std::setprecision(3) << f;
    return os.str();
  }
  static std::string protocol_key(double thr) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.0e", thr); // "1e-06"
    std::string s(buf);
    // normalize "+": "1e+00" -> "1e+00" (keep plus; stod handles it)
    return s;
  }
  static double protocol_numeric(const std::string &p) {
    try {
      return std::stod(p);
    } catch (...) {
      return std::numeric_limits<double>::infinity();
    }
  }

private:
  World &world_;
  std::string path_;
  json data_;

  static void ensure_root(json &root) {
    if (!root.is_object())
      root = json::object();
    if (!root.contains("states"))
      root["states"] = json::object();
    ensure_excited_root(root);
  }

  static void ensure_excited_root(json &root) {
    if (!root.contains("excited_states") || !root["excited_states"].is_object()) {
      root["excited_states"] = json::object();
    }
    auto &excited = root["excited_states"];
    if (!excited.contains("plan") || !excited["plan"].is_object()) {
      excited["plan"] = json::object();
    }
    if (!excited.contains("protocols") || !excited["protocols"].is_object()) {
      excited["protocols"] = json::object();
    }
  }

  static void ensure_excited_protocol(json &root, const std::string &pkey) {
    ensure_excited_root(root);
    auto &protocols = root["excited_states"]["protocols"];
    if (!protocols.contains(pkey) || !protocols[pkey].is_object()) {
      protocols[pkey] = json::object();
    }
    auto &node = protocols[pkey];
    if (!node.contains("saved") || !node["saved"].is_boolean()) {
      node["saved"] = false;
    }
    if (!node.contains("converged") || !node["converged"].is_boolean()) {
      node["converged"] = false;
    }
    if (!node.contains("timings") || !node["timings"].is_object()) {
      node["timings"] = {{"wall_seconds", 0.0}, {"cpu_seconds", 0.0}};
    }
    if (!node.contains("energies") || !node["energies"].is_array()) {
      node["energies"] = json::array();
    }
    if (!node.contains("state_names") || !node["state_names"].is_array()) {
      node["state_names"] = json::array();
    }
    if (!node.contains("roots") || !node["roots"].is_array()) {
      node["roots"] = json::array();
    }
    if (!node.contains("slot_permutation") ||
        !node["slot_permutation"].is_array()) {
      node["slot_permutation"] = json::array();
    }
    if (!node.contains("iterations") || !node["iterations"].is_number_unsigned()) {
      node["iterations"] = 0;
    }
    if (!node.contains("stage_status") || !node["stage_status"].is_string()) {
      node["stage_status"] = "placeholder_pending_solver";
    }
    if (!node.contains("response_variant") ||
        !node["response_variant"].is_string()) {
      node["response_variant"] = "unknown";
    }
    if (!node.contains("restart_support_mode") ||
        !node["restart_support_mode"].is_string()) {
      node["restart_support_mode"] = "guess_only";
    }
    if (!node.contains("restart_source") || !node["restart_source"].is_string()) {
      node["restart_source"] = "none";
    }
    if (!node.contains("snapshot_kind") || !node["snapshot_kind"].is_string()) {
      node["snapshot_kind"] = "none";
    }
    if (!node.contains("bundle_state_present") ||
        !node["bundle_state_present"].is_boolean()) {
      node["bundle_state_present"] = false;
    }
    if (!node.contains("restart_capable") ||
        !node["restart_capable"].is_boolean()) {
      node["restart_capable"] = false;
    }
    if (!node.contains("restart_source_threshold") ||
        !node["restart_source_threshold"].is_number()) {
      node["restart_source_threshold"] = 0.0;
    }
    if (!node.contains("owner_group") || !node["owner_group"].is_number_unsigned()) {
      node["owner_group"] = 0;
    }
    if (!node.contains("attempted") || !node["attempted"].is_boolean()) {
      node["attempted"] = false;
    }
    if (!node.contains("failed") || !node["failed"].is_boolean()) {
      node["failed"] = false;
    }
    if (!node.contains("skipped") || !node["skipped"].is_boolean()) {
      node["skipped"] = false;
    }
    if (!node.contains("restart_reused") || !node["restart_reused"].is_boolean()) {
      node["restart_reused"] = false;
    }
    if (!node.contains("residual_norms") || !node["residual_norms"].is_array()) {
      node["residual_norms"] = json::array();
    }
    if (!node.contains("density_change_norms") ||
        !node["density_change_norms"].is_array()) {
      node["density_change_norms"] = json::array();
    }
    if (!node.contains("relative_residual_norms") ||
        !node["relative_residual_norms"].is_array()) {
      node["relative_residual_norms"] = json::array();
    }
    if (!node.contains("iteration_max_residuals") ||
        !node["iteration_max_residuals"].is_array()) {
      node["iteration_max_residuals"] = json::array();
    }
    if (!node.contains("iteration_max_density_changes") ||
        !node["iteration_max_density_changes"].is_array()) {
      node["iteration_max_density_changes"] = json::array();
    }
    if (!node.contains("iteration_max_relative_residuals") ||
        !node["iteration_max_relative_residuals"].is_array()) {
      node["iteration_max_relative_residuals"] = json::array();
    }
    if (!node.contains("convergence_mode") ||
        !node["convergence_mode"].is_string()) {
      node["convergence_mode"] = "max_residual";
    }
    if (!node.contains("accelerator_mode") ||
        !node["accelerator_mode"].is_string()) {
      node["accelerator_mode"] = "none";
    }
    if (!node.contains("accelerator_subspace") ||
        !node["accelerator_subspace"].is_number_unsigned()) {
      node["accelerator_subspace"] = 0;
    }
    if (!node.contains("density_convergence_target") ||
        !node["density_convergence_target"].is_number()) {
      node["density_convergence_target"] = 0.0;
    }
    if (!node.contains("relative_convergence_target") ||
        !node["relative_convergence_target"].is_number()) {
      node["relative_convergence_target"] = 0.0;
    }
    if (!node.contains("max_rotation") ||
        !node["max_rotation"].is_number()) {
      node["max_rotation"] = 0.0;
    }
  }
  static void ensure_state(json &root, const std::string &sid) {
    ensure_root(root);
    auto &states = root["states"];
    if (!states.contains(sid)) {
      states[sid] = json::object();
      states[sid]["protocols"] = json::object();
    } else if (!states[sid].contains("protocols")) {
      states[sid]["protocols"] = json::object();
    }
  }
  static void ensure_protocol(json &root, const std::string &sid,
                              const std::string &pkey) {
    ensure_state(root, sid);
    auto &protos = root["states"][sid]["protocols"];
    if (!protos.contains(pkey)) {
      protos[pkey] = json::object(
          {{"saved", json::object()},
           {"converged", json::object()},
           {"solver_diagnostics", json::object()},
           {"restart_provenance", json::object()}});
    }
  }

  // -------- Printing helpers / sorting --------
  static std::string fmt_sci(double x, int digits_after_decimal = 0) {
    std::ostringstream os;
    os << std::scientific << std::setprecision(digits_after_decimal) << x;
    return os.str();
  }
  static std::string fmt_fixed(double x, int decimals = 3) {
    std::ostringstream os;
    os << std::fixed << std::setprecision(decimals) << x;
    return os.str();
  }

  static std::string make_excited_root_id(size_t stable_index) {
    std::ostringstream os;
    os << "es_root_" << std::setw(4) << std::setfill('0') << stable_index;
    return os.str();
  }

  static std::string excited_alpha_suffix(size_t index) {
    std::string suffix;
    size_t value = index;
    do {
      const char c = static_cast<char>('a' + (value % 26));
      suffix.push_back(c);
      if (value < 26) {
        break;
      }
      value = (value / 26) - 1;
    } while (true);
    std::reverse(suffix.begin(), suffix.end());
    return suffix;
  }

  static std::vector<std::string>
  assign_excited_state_names(const std::vector<double> &energies,
                             double protocol_threshold) {
    const size_t n = energies.size();
    std::vector<std::string> names(n);
    if (n == 0) {
      return names;
    }

    const double grouping_tol =
        std::max(1.0e-12, 10.0 * std::abs(protocol_threshold));
    size_t group_root_index = 1;
    size_t i = 0;
    while (i < n) {
      size_t j = i + 1;
      const double anchor = energies[i];
      while (j < n && std::abs(energies[j] - anchor) <= grouping_tol) {
        ++j;
      }
      const size_t group_size = j - i;
      if (group_size == 1) {
        names[i] = "es" + std::to_string(group_root_index);
      } else {
        for (size_t k = 0; k < group_size; ++k) {
          names[i + k] = "es" + std::to_string(group_root_index) +
                         excited_alpha_suffix(k);
        }
      }
      ++group_root_index;
      i = j;
    }
    return names;
  }

  static std::vector<std::pair<double, std::string>>
  numeric_keys(const json &obj) {
    std::vector<std::pair<double, std::string>> out;
    out.reserve(obj.size());
    for (const auto &kv : obj.items())
      out.emplace_back(std::stod(kv.key()), kv.key());
    std::sort(out.begin(), out.end(),
              [](const auto &a, const auto &b) { return a.first < b.first; });
    return out;
  }

  // -------- core flag ops (with canonical keys precomputed) --------
 [[nodiscard]] bool get_flag(const std::string &sid, const std::string &pkey,
                const std::string &fkey, const std::string &which) const {
    if (!data_.contains("states"))
      return false;
    const auto sit = data_["states"].find(sid);
    if (sit == data_["states"].end())
      return false;
    const auto pit = (*sit)["protocols"].find(pkey);
    if (pit == (*sit)["protocols"].end())
      return false;
    const auto &node = (*pit);
    if (!node.contains(which) || !node[which].is_object())
      return false;
    const auto &sub = node[which];
    if (!sub.contains(fkey))
      return false;
    return sub.at(fkey).get<bool>();
  }

  void set_flag(const std::string &sid, const std::string &pkey,
                const std::string &fkey, const std::string &which, bool value) {
    ensure_protocol(data_, sid, pkey);
    data_["states"][sid]["protocols"][pkey][which][fkey] = value;
    write();
  }
};
