#include "ExcitedStateBundleSolver.hpp"

#include "GroundStateData.hpp"
#include "ResponseSolver.hpp"
#include "ResponseSolverUtils.hpp"
#include "ResponseVector.hpp"
#include "../molresponse/response_macrotask.hpp"

#include <madness/external/nlohmann_json/json.hpp>
#include <madness/chem/molecular_functors.h>
#include <madness/world/worldgop.h>
#include <madness/mra/funcdefaults.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

namespace fs = std::filesystem;
using json = nlohmann::json;

std::string protocol_key(double threshold) {
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%.0e", threshold);
  return std::string(buf);
}

std::string format_vector_head(const std::vector<double> &values,
                               size_t max_items = 6,
                               int precision = 4) {
  std::ostringstream os;
  os << "[";
  const size_t n = std::min(max_items, values.size());
  for (size_t i = 0; i < n; ++i) {
    if (i > 0) {
      os << ", ";
    }
    os << std::scientific << std::setprecision(precision) << values[i];
  }
  if (values.size() > n) {
    os << ", ...";
  }
  os << "]";
  return os.str();
}

std::string format_string_head(const std::vector<std::string> &values,
                               size_t max_items = 6) {
  std::ostringstream os;
  os << "[";
  const size_t n = std::min(max_items, values.size());
  for (size_t i = 0; i < n; ++i) {
    if (i > 0) {
      os << ", ";
    }
    os << values[i];
  }
  if (values.size() > n) {
    os << ", ...";
  }
  os << "]";
  return os.str();
}

std::pair<size_t, double>
max_abs_value_with_index(const std::vector<double> &values) {
  if (values.empty()) {
    return {0, 0.0};
  }
  size_t max_idx = 0;
  double max_val = std::fabs(values[0]);
  for (size_t i = 1; i < values.size(); ++i) {
    const double candidate = std::fabs(values[i]);
    if (candidate > max_val) {
      max_val = candidate;
      max_idx = i;
    }
  }
  return {max_idx, max_val};
}

std::string alpha_suffix(size_t index) {
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

std::vector<std::string> assign_excited_state_names(
    const std::vector<double> &energies, double protocol_threshold) {
  const size_t n = energies.size();
  std::vector<std::string> names(n);
  if (n == 0) {
    return names;
  }

  const double grouping_tol = std::max(1.0e-12, 10.0 * protocol_threshold);
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
                       alpha_suffix(k);
      }
    }
    ++group_root_index;
    i = j;
  }
  return names;
}

std::string make_excited_root_id(size_t stable_index) {
  std::ostringstream os;
  os << "es_root_" << std::setw(4) << std::setfill('0') << stable_index;
  return os.str();
}

template <std::size_t NDIM>
class LocalizedGaussianGuess : public FunctionFunctorInterface<double, NDIM> {
  using coordT = Vector<double, NDIM>;

public:
  LocalizedGaussianGuess(const coordT &origin, double exponent,
                         std::vector<int> powers)
      : origin_(origin), exponent_(exponent), powers_(std::move(powers)) {
    if (powers_.size() != NDIM) {
      powers_.assign(NDIM, 0);
    }
  }

  double operator()(const coordT &xyz) const override {
    double radial2 = 0.0;
    double prefactor = 1.0;
    for (std::size_t i = 0; i < NDIM; ++i) {
      const double x = xyz[i] - origin_[i];
      radial2 += x * x;
      if (powers_[i] > 0) {
        prefactor *= std::pow(x, powers_[i]);
      }
    }
    return prefactor * std::exp(-exponent_ * radial2);
  }

  std::vector<coordT> special_points() const override { return {origin_}; }

private:
  coordT origin_;
  double exponent_ = 1.0;
  std::vector<int> powers_;
};

bool file_exists_world(madness::World &world, const std::string &path) {
  bool exists = false;
  if (world.rank() == 0) {
    exists = fs::exists(path);
  }
  world.gop.broadcast_serializable(exists, 0);
  return exists;
}

struct RestartSnapshot {
  // Stage-2c restart contract.
  //
  // Full restart-capable snapshots carry the active typed ResponseVector bundle
  // plus stable root metadata. Trial-space x states remain as a weaker
  // guess-only fallback, but they are not treated as sufficient for exact
  // restart of dynamic/full-response variants.
  bool has_data = false;
  size_t schema_version = 2;
  bool converged = false;
  bool bundle_state_present = false;
  bool restart_capable = false;
  std::string snapshot_kind = "unknown";
  std::string response_variant = "unknown";
  std::string restart_support_mode = "guess_only";
  size_t protocol_index = 0;
  double protocol_threshold = 0.0;
  int protocol_k = -1;
  size_t owner_group = 0;
  size_t archive_world_size = 0;
  size_t iterations = 0;
  size_t state_count = 0;
  size_t num_orbitals = 0;
  std::vector<double> energies;
  std::vector<std::string> state_names;
  std::vector<ExcitedRootDescriptor> roots;
  std::vector<size_t> slot_permutation;
  std::vector<double> residual_norms;
  std::vector<double> density_change_norms;
  std::vector<double> relative_residual_norms;
  std::vector<double> iteration_max_residuals;
  std::vector<double> iteration_max_density_changes;
  std::vector<double> iteration_max_relative_residuals;
  std::string convergence_mode = "max_residual";
  std::string accelerator_mode = "none";
  size_t accelerator_subspace = 0;
  double density_convergence_target = 0.0;
  double relative_convergence_target = 0.0;
  double max_rotation = 0.0;
  bool has_trial_states = false;
  int trial_space_k = -1;
  size_t trial_space_num_orbitals = 0;
  std::vector<vector_real_function_3d> trial_states;
  int bundle_k = -1;
  size_t bundle_num_orbitals = 0;
  std::vector<ResponseVector> response_bundle;
};

std::string trial_space_archive_path(const std::string &snapshot_path) {
  return snapshot_path + ".trial_states";
}

std::string response_bundle_archive_path(const std::string &snapshot_path) {
  return snapshot_path + ".response_bundle";
}

bool archive_exists_world(madness::World &world, const std::string &path) {
  bool exists = false;
  if (world.rank() == 0) {
    exists = fs::exists(path) || fs::exists(path + ".00000");
  }
  world.gop.broadcast_serializable(exists, 0);
  return exists;
}

void read_restart_trial_states(madness::World &world, const std::string &path,
                               RestartSnapshot &snapshot) {
  const std::string trial_path = trial_space_archive_path(path);
  if (!archive_exists_world(world, trial_path)) {
    return;
  }

  archive::ParallelInputArchive ar(world, trial_path.c_str());
  int loaded_k = -1;
  size_t num_states = 0;
  ar & loaded_k;
  ar & num_states;

  std::vector<vector_real_function_3d> trial_states;
  trial_states.reserve(num_states);
  size_t max_orbitals = 0;
  for (size_t i = 0; i < num_states; ++i) {
    size_t state_orbitals = 0;
    ar & state_orbitals;
    max_orbitals = std::max(max_orbitals, state_orbitals);
    auto state = zero_functions_compressed<double, 3>(
        world, static_cast<int>(state_orbitals));
    for (size_t j = 0; j < state_orbitals; ++j) {
      ar & state[j];
    }
    trial_states.push_back(std::move(state));
  }
  world.gop.fence();

  snapshot.has_trial_states = !trial_states.empty();
  snapshot.trial_space_k = loaded_k;
  snapshot.trial_space_num_orbitals = max_orbitals;
  snapshot.trial_states = std::move(trial_states);
}

void write_restart_trial_states(madness::World &world, const std::string &path,
                                const RestartSnapshot &snapshot) {
  if (!snapshot.has_trial_states || snapshot.trial_states.empty()) {
    return;
  }
  const std::string trial_path = trial_space_archive_path(path);
  archive::ParallelOutputArchive ar(world, trial_path.c_str(), 1);
  const size_t num_states = snapshot.trial_states.size();
  ar & snapshot.trial_space_k;
  ar & num_states;
  for (const auto &state : snapshot.trial_states) {
    const size_t state_orbitals = state.size();
    ar & state_orbitals;
    for (size_t j = 0; j < state_orbitals; ++j) {
      ar & state[j];
    }
  }
  world.gop.fence();
}

void read_restart_response_bundle(madness::World &world, const std::string &path,
                                  RestartSnapshot &snapshot) {
  const std::string bundle_path = response_bundle_archive_path(path);
  if (!snapshot.bundle_state_present || !archive_exists_world(world, bundle_path)) {
    snapshot.bundle_state_present = false;
    snapshot.response_bundle.clear();
    return;
  }

  archive::ParallelInputArchive ar(world, bundle_path.c_str());
  std::string variant;
  size_t num_states = 0;
  size_t num_orbitals = 0;
  int loaded_k = -1;
  ar & variant;
  ar & num_states;
  ar & num_orbitals;
  ar & loaded_k;

  bool compatible = response_variant_is_valid(variant);
  std::vector<ResponseVector> response_bundle;
  response_bundle.reserve(num_states);
  for (size_t i = 0; i < num_states; ++i) {
    size_t flat_size = 0;
    ar & flat_size;
    auto loaded_flat = zero_functions_compressed<double, 3>(
        world, static_cast<int>(flat_size));
    for (size_t j = 0; j < flat_size; ++j) {
      ar & loaded_flat[j];
    }
    if (!compatible) {
      continue;
    }
    auto response = make_response_vector_from_variant(num_orbitals, variant);
    auto &flat = get_flat(response);
    if (flat.size() != flat_size) {
      compatible = false;
      response_bundle.clear();
      continue;
    }
    flat = std::move(loaded_flat);
    sync_response(response);
    response_bundle.push_back(std::move(response));
  }
  world.gop.fence();

  if (!compatible) {
    snapshot.bundle_state_present = false;
    snapshot.restart_capable = false;
    snapshot.response_bundle.clear();
    return;
  }

  snapshot.response_variant = variant;
  snapshot.bundle_k = loaded_k;
  snapshot.bundle_num_orbitals = num_orbitals;
  snapshot.response_bundle = std::move(response_bundle);
  snapshot.bundle_state_present = !snapshot.response_bundle.empty();
  if (!snapshot.bundle_state_present) {
    snapshot.restart_capable = false;
  }
  if (snapshot.state_count == 0) {
    snapshot.state_count = snapshot.response_bundle.size();
  }
  if (snapshot.num_orbitals == 0) {
    snapshot.num_orbitals = snapshot.bundle_num_orbitals;
  }
}

void write_restart_response_bundle(madness::World &world, const std::string &path,
                                   const RestartSnapshot &snapshot) {
  if (!snapshot.bundle_state_present || snapshot.response_bundle.empty()) {
    return;
  }

  const std::string bundle_path = response_bundle_archive_path(path);
  archive::ParallelOutputArchive ar(world, bundle_path.c_str(), 1);
  const std::string variant = response_variant_is_valid(snapshot.response_variant)
                                  ? snapshot.response_variant
                                  : response_variant_name(snapshot.response_bundle.front());
  const size_t num_states = snapshot.response_bundle.size();
  const size_t num_orbitals =
      snapshot.bundle_num_orbitals > 0 ? snapshot.bundle_num_orbitals
                                       : response_num_orbitals(snapshot.response_bundle.front());
  const int bundle_k = (snapshot.bundle_k >= 0)
                           ? snapshot.bundle_k
                           : madness::FunctionDefaults<3>::get_k();
  std::string archived_variant = variant;
  size_t archived_num_states = num_states;
  size_t archived_num_orbitals = num_orbitals;
  int archived_bundle_k = bundle_k;
  ar & archived_variant;
  ar & archived_num_states;
  ar & archived_num_orbitals;
  ar & archived_bundle_k;
  for (const auto &response : snapshot.response_bundle) {
    const auto &flat = get_flat(response);
    const size_t flat_size = flat.size();
    size_t archived_flat_size = flat_size;
    ar & archived_flat_size;
    for (size_t j = 0; j < flat_size; ++j) {
      auto archived_function = flat[j];
      ar & archived_function;
    }
  }
  world.gop.fence();
}

RestartSnapshot read_restart_snapshot(madness::World &world,
                                      const std::string &path) {
  RestartSnapshot snapshot;
  if (world.rank() == 0 && fs::exists(path)) {
    std::ifstream in(path);
    if (in.good()) {
      try {
        json j;
        in >> j;
        snapshot.has_data = true;
        snapshot.schema_version = j.value("schema_version", static_cast<size_t>(1));
        snapshot.converged = j.value("converged", false);
        snapshot.bundle_state_present = j.value("bundle_state_present", false);
        snapshot.restart_capable = j.value("restart_capable", false);
        snapshot.snapshot_kind =
            j.value("snapshot_kind", std::string("legacy_trial_seed"));
        snapshot.response_variant =
            j.value("response_variant", std::string("unknown"));
        snapshot.restart_support_mode =
            j.value("restart_support_mode", std::string("guess_only"));
        snapshot.protocol_index = j.value("protocol_index", static_cast<size_t>(0));
        snapshot.protocol_threshold = j.value("protocol_threshold", 0.0);
        snapshot.protocol_k = j.value("protocol_k", -1);
        snapshot.owner_group = j.value("owner_group", static_cast<size_t>(0));
        snapshot.archive_world_size =
            j.value("archive_world_size", static_cast<size_t>(0));
        snapshot.iterations = j.value("iterations", static_cast<size_t>(0));
        snapshot.state_count = j.value("state_count", static_cast<size_t>(0));
        snapshot.num_orbitals = j.value("num_orbitals", static_cast<size_t>(0));
        snapshot.bundle_k = j.value("bundle_k", -1);
        snapshot.bundle_num_orbitals =
            j.value("bundle_num_orbitals", static_cast<size_t>(0));
        if (j.contains("energies") && j["energies"].is_array()) {
          snapshot.energies = j["energies"].get<std::vector<double>>();
        }
        if (j.contains("state_names") && j["state_names"].is_array()) {
          snapshot.state_names = j["state_names"].get<std::vector<std::string>>();
        }
        if (j.contains("roots") && j["roots"].is_array()) {
          snapshot.roots = j["roots"].get<std::vector<ExcitedRootDescriptor>>();
        }
        if (j.contains("slot_permutation") && j["slot_permutation"].is_array()) {
          snapshot.slot_permutation =
              j["slot_permutation"].get<std::vector<size_t>>();
        }
        if (j.contains("residual_norms") && j["residual_norms"].is_array()) {
          snapshot.residual_norms =
              j["residual_norms"].get<std::vector<double>>();
        }
        if (j.contains("density_change_norms") &&
            j["density_change_norms"].is_array()) {
          snapshot.density_change_norms =
              j["density_change_norms"].get<std::vector<double>>();
        }
        if (j.contains("relative_residual_norms") &&
            j["relative_residual_norms"].is_array()) {
          snapshot.relative_residual_norms =
              j["relative_residual_norms"].get<std::vector<double>>();
        }
        if (j.contains("iteration_max_residuals") &&
            j["iteration_max_residuals"].is_array()) {
          snapshot.iteration_max_residuals =
              j["iteration_max_residuals"].get<std::vector<double>>();
        }
        if (j.contains("iteration_max_density_changes") &&
            j["iteration_max_density_changes"].is_array()) {
          snapshot.iteration_max_density_changes =
              j["iteration_max_density_changes"].get<std::vector<double>>();
        }
        if (j.contains("iteration_max_relative_residuals") &&
            j["iteration_max_relative_residuals"].is_array()) {
          snapshot.iteration_max_relative_residuals =
              j["iteration_max_relative_residuals"].get<std::vector<double>>();
        }
        snapshot.convergence_mode =
            j.value("convergence_mode", std::string("max_residual"));
        snapshot.accelerator_mode =
            j.value("accelerator_mode", std::string("none"));
        snapshot.accelerator_subspace =
            j.value("accelerator_subspace", static_cast<size_t>(0));
        snapshot.density_convergence_target =
            j.value("density_convergence_target", 0.0);
        snapshot.relative_convergence_target =
            j.value("relative_convergence_target", 0.0);
        snapshot.max_rotation = j.value("max_rotation", 0.0);
      } catch (...) {
        snapshot = RestartSnapshot{};
      }
    }
  }
  world.gop.broadcast_serializable(snapshot.has_data, 0);
  world.gop.broadcast_serializable(snapshot.schema_version, 0);
  world.gop.broadcast_serializable(snapshot.converged, 0);
  world.gop.broadcast_serializable(snapshot.bundle_state_present, 0);
  world.gop.broadcast_serializable(snapshot.restart_capable, 0);
  world.gop.broadcast_serializable(snapshot.snapshot_kind, 0);
  world.gop.broadcast_serializable(snapshot.response_variant, 0);
  world.gop.broadcast_serializable(snapshot.restart_support_mode, 0);
  world.gop.broadcast_serializable(snapshot.protocol_index, 0);
  world.gop.broadcast_serializable(snapshot.protocol_threshold, 0);
  world.gop.broadcast_serializable(snapshot.protocol_k, 0);
  world.gop.broadcast_serializable(snapshot.owner_group, 0);
  world.gop.broadcast_serializable(snapshot.archive_world_size, 0);
  world.gop.broadcast_serializable(snapshot.iterations, 0);
  world.gop.broadcast_serializable(snapshot.state_count, 0);
  world.gop.broadcast_serializable(snapshot.num_orbitals, 0);
  world.gop.broadcast_serializable(snapshot.energies, 0);
  world.gop.broadcast_serializable(snapshot.state_names, 0);
  world.gop.broadcast_serializable(snapshot.slot_permutation, 0);
  world.gop.broadcast_serializable(snapshot.residual_norms, 0);
  world.gop.broadcast_serializable(snapshot.density_change_norms, 0);
  world.gop.broadcast_serializable(snapshot.relative_residual_norms, 0);
  world.gop.broadcast_serializable(snapshot.iteration_max_residuals, 0);
  world.gop.broadcast_serializable(snapshot.iteration_max_density_changes, 0);
  world.gop.broadcast_serializable(snapshot.iteration_max_relative_residuals,
                                   0);
  world.gop.broadcast_serializable(snapshot.convergence_mode, 0);
  world.gop.broadcast_serializable(snapshot.accelerator_mode, 0);
  world.gop.broadcast_serializable(snapshot.accelerator_subspace, 0);
  world.gop.broadcast_serializable(snapshot.density_convergence_target, 0);
  world.gop.broadcast_serializable(snapshot.relative_convergence_target, 0);
  world.gop.broadcast_serializable(snapshot.max_rotation, 0);
  world.gop.broadcast_serializable(snapshot.bundle_k, 0);
  world.gop.broadcast_serializable(snapshot.bundle_num_orbitals, 0);
  std::string roots_json = "[]";
  if (world.rank() == 0) {
    roots_json = json(snapshot.roots).dump();
  }
  world.gop.broadcast_serializable(roots_json, 0);
  if (world.rank() != 0) {
    snapshot.roots =
        json::parse(roots_json).get<std::vector<ExcitedRootDescriptor>>();
  }
  read_restart_trial_states(world, path, snapshot);
  read_restart_response_bundle(world, path, snapshot);
  return snapshot;
}

void write_restart_snapshot(madness::World &world, const std::string &path,
                            const RestartSnapshot &snapshot) {
  if (world.rank() == 0) {
    json j = {{"schema_version", snapshot.schema_version},
              {"converged", snapshot.converged},
              {"bundle_state_present", snapshot.bundle_state_present},
              {"restart_capable", snapshot.restart_capable},
              {"snapshot_kind", snapshot.snapshot_kind},
              {"response_variant", snapshot.response_variant},
              {"restart_support_mode", snapshot.restart_support_mode},
              {"protocol_index", snapshot.protocol_index},
              {"protocol_threshold", snapshot.protocol_threshold},
              {"protocol_k", snapshot.protocol_k},
              {"owner_group", snapshot.owner_group},
              {"archive_world_size", snapshot.archive_world_size},
              {"iterations", snapshot.iterations},
              {"state_count", snapshot.state_count},
              {"num_orbitals", snapshot.num_orbitals},
              {"energies", snapshot.energies},
              {"state_names", snapshot.state_names},
              {"roots", snapshot.roots},
              {"slot_permutation", snapshot.slot_permutation},
              {"bundle_k", snapshot.bundle_k},
              {"bundle_num_orbitals", snapshot.bundle_num_orbitals},
              {"residual_norms", snapshot.residual_norms},
              {"density_change_norms", snapshot.density_change_norms},
              {"relative_residual_norms", snapshot.relative_residual_norms},
              {"iteration_max_residuals", snapshot.iteration_max_residuals},
              {"iteration_max_density_changes",
               snapshot.iteration_max_density_changes},
              {"iteration_max_relative_residuals",
               snapshot.iteration_max_relative_residuals},
              {"convergence_mode", snapshot.convergence_mode},
              {"accelerator_mode", snapshot.accelerator_mode},
              {"accelerator_subspace", snapshot.accelerator_subspace},
              {"density_convergence_target",
               snapshot.density_convergence_target},
              {"relative_convergence_target",
               snapshot.relative_convergence_target},
              {"max_rotation", snapshot.max_rotation}};
    std::ofstream out(path);
    out << j.dump(2) << "\n";
  }
  world.gop.fence();
  write_restart_trial_states(world, path, snapshot);
  write_restart_response_bundle(world, path, snapshot);
}

struct ExcitedTrialSpace {
  size_t num_states = 0;
  size_t num_orbitals = 0;
  std::vector<vector_real_function_3d> x_states;
  std::vector<double> omega;
};

// Default adapter used until the legacy excited-state solver is fully
// reintegrated into the molresponse_v2 build.
class RestartAwareExcitedScaffoldSolver final : public ExcitedStateBundleSolver {
public:
  explicit RestartAwareExcitedScaffoldSolver(ExcitedBundleSolverConfig config)
      : config_(std::move(config)), workflow_(config_) {}

  [[nodiscard]] std::string name() const override {
    return "restart_aware_scaffold_adapter";
  }

  [[nodiscard]] ExcitedBundleProtocolResult
  solve_protocol(madness::World &world,
                 const ExcitedBundleProtocolInput &input) const override {
    return workflow_.solve_protocol(world, input);
  }

private:
  class ExcitedProtocolWorkflow {
  public:
    explicit ExcitedProtocolWorkflow(ExcitedBundleSolverConfig config)
        : config_(std::move(config)) {}

    [[nodiscard]] ExcitedBundleProtocolResult
    solve_protocol(madness::World &world,
                   const ExcitedBundleProtocolInput &input) {
      ExcitedBundleProtocolResult result;
      const std::string protocol_file = protocol_restart_file(input.threshold);
      const auto restart_snapshot = read_restart_snapshot(world, protocol_file);

      prepare_protocol(world, input.threshold);
      ensure_ground_data(world);
      const bool full_protocol_restart_ready =
          snapshot_supports_full_restart(restart_snapshot, input);

      // Same-protocol skip requires a full restart-capable bundle snapshot, not
      // just legacy x-state metadata.
      if (input.restart_saved && input.restart_converged &&
          full_protocol_restart_ready) {
        load_restart_seed(restart_snapshot);
        ensure_trial_space_matches_guess(world, input);
        result.attempted = false;
        result.saved = true;
        result.converged = true;
        result.failed = false;
        result.skipped = true;
        result.restart_reused = true;
        result.stage_status = "restart_ready_skip";
        result.response_variant = restart_snapshot.response_variant;
        result.restart_support_mode = restart_snapshot.restart_support_mode;
        result.restart_source = "current_protocol_snapshot";
        result.snapshot_kind = restart_snapshot.snapshot_kind;
        result.bundle_state_present = restart_snapshot.bundle_state_present;
        result.restart_capable = restart_snapshot.restart_capable;
        result.restart_source_threshold = restart_snapshot.protocol_threshold;
        result.iterations = restart_snapshot.iterations;
        result.energies = restart_snapshot.energies;
        result.state_names = restart_snapshot.state_names;
        result.roots = restart_snapshot.roots;
        result.slot_permutation = restart_snapshot.slot_permutation;
        if (result.state_names.size() != result.energies.size()) {
          result.state_names =
              assign_excited_state_names(result.energies, input.threshold);
        }
        if (result.roots.empty() && !result.energies.empty()) {
          result.roots.reserve(result.energies.size());
          for (size_t i = 0; i < result.energies.size(); ++i) {
            const std::string display_name =
                (i < result.state_names.size()) ? result.state_names[i]
                                                : std::string{};
            result.roots.push_back(
                {make_excited_root_id(i), i, i, result.energies[i],
                 display_name});
          }
        }
        if (result.slot_permutation.empty() && !result.roots.empty()) {
          result.slot_permutation.reserve(result.roots.size());
          for (const auto &root : result.roots) {
            result.slot_permutation.push_back(root.stable_index);
          }
        }
        result.residual_norms = restart_snapshot.residual_norms;
        result.density_change_norms = restart_snapshot.density_change_norms;
        result.relative_residual_norms =
            restart_snapshot.relative_residual_norms;
        result.iteration_max_residuals = restart_snapshot.iteration_max_residuals;
        result.iteration_max_density_changes =
            restart_snapshot.iteration_max_density_changes;
        result.iteration_max_relative_residuals =
            restart_snapshot.iteration_max_relative_residuals;
        result.convergence_mode = restart_snapshot.convergence_mode;
        result.accelerator_mode = restart_snapshot.accelerator_mode;
        result.accelerator_subspace = restart_snapshot.accelerator_subspace;
        result.density_convergence_target =
            restart_snapshot.density_convergence_target;
        result.relative_convergence_target =
            restart_snapshot.relative_convergence_target;
        result.max_rotation = restart_snapshot.max_rotation;
        return result;
      }

      if (input.restart_saved && input.restart_converged && world.rank() == 0 &&
          !full_protocol_restart_ready) {
        madness::print(
            "WARN EXCITED_RESTART_SKIP_BLOCKED threshold=", input.threshold,
            " reason=missing_full_bundle_snapshot variant=",
            restart_snapshot.response_variant,
            " snapshot_kind=", restart_snapshot.snapshot_kind,
            " restart_capable=", restart_snapshot.restart_capable);
      }

      const auto initialize_result =
          initialize_protocol_guess(world, input, restart_snapshot);

      result = iterate(world, input, initialize_result.restart_reused);
      if (!initialize_result.stage_status.empty()) {
        result.stage_status =
            initialize_result.stage_status + "__" + result.stage_status;
      }
      if (result.response_variant.empty() || result.response_variant == "unknown") {
        result.response_variant = initialize_result.response_variant;
      }
      if (result.restart_support_mode.empty()) {
        result.restart_support_mode = initialize_result.restart_support_mode;
      }
      result.restart_source = initialize_result.restart_source;
      if (result.restart_source_threshold <= 0.0) {
        result.restart_source_threshold =
            initialize_result.restart_source_threshold;
      }
      return result;
    }

  private:
    struct ProtocolInitResult {
      std::string stage_status;
      bool restart_reused = false;
      std::string restart_source = "none";
      std::string snapshot_kind = "none";
      std::string response_variant = "unknown";
      std::string restart_support_mode = "guess_only";
      bool bundle_state_present = false;
      bool restart_capable = false;
      double restart_source_threshold = 0.0;
    };

    struct RestartSeedChoice {
      RestartSnapshot snapshot;
      bool found = false;
      bool from_current_protocol = false;
      bool from_lower_protocol = false;
      bool from_guess_archive = false;
      bool uses_full_bundle = false;
      double source_threshold = 0.0;
    };

    ExcitedBundleSolverConfig config_;
    bool has_active_guess_ = false;
    ExcitedTrialSpace trial_space_;
    std::vector<ResponseVector> active_response_bundle_;
    std::vector<double> omega_;
    std::vector<std::string> state_names_;
    std::vector<ExcitedRootDescriptor> root_descriptors_;
    std::vector<size_t> slot_permutation_;
    size_t next_root_stable_index_ = 0;
    std::vector<double> residual_norms_;
    std::vector<double> density_change_norms_;
    std::vector<double> relative_residual_norms_;
    size_t last_iteration_count_ = 0;
    std::vector<double> iteration_max_residuals_;
    std::vector<double> iteration_max_density_changes_;
    std::vector<double> iteration_max_relative_residuals_;
    std::string convergence_mode_ = "max_residual";
    std::string accelerator_mode_ = "none";
    size_t accelerator_subspace_ = 0;
    double density_convergence_target_ = 0.0;
    double relative_convergence_target_ = 0.0;
    double max_rotation_ = 0.0;
    std::unique_ptr<GroundStateData> ground_data_;
    int ground_data_k_ = -1;
    double ground_data_thresh_ = std::numeric_limits<double>::infinity();
    int preliminaries_k_ = -1;
    double preliminaries_thresh_ = std::numeric_limits<double>::infinity();
    size_t ground_num_orbitals_ = 0;
    Tensor<double> ground_energies_;
    vector_real_function_3d ground_orbitals_;
    int trial_space_k_ = -1;
    double trial_space_thresh_ = std::numeric_limits<double>::infinity();
    int active_bundle_k_ = -1;
    double active_bundle_thresh_ = std::numeric_limits<double>::infinity();

    [[nodiscard]] std::string guess_archive_file() const {
      std::ostringstream os;
      os << config_.output_prefix << ".excited_bundle.guess.restartdata";
      return os.str();
    }

    [[nodiscard]] std::string expected_response_variant(bool tda) const {
      const bool unrestricted =
          (ground_data_ != nullptr) && !ground_data_->isSpinRestricted();
      return response_variant_name(tda, unrestricted);
    }

    [[nodiscard]] static std::string restart_support_mode_for_variant(
        const std::string &variant) {
      // Phase 2 boundary: exact full-bundle restart is implemented for the
      // restricted-shell variants. Unrestricted snapshots still persist all
      // channels, but they remain guess-only until that path is validated.
      return response_variant_is_unrestricted(variant)
                 ? "guess_only_unrestricted_variant"
                 : "full_bundle_resume";
    }

    [[nodiscard]] static bool variant_supports_full_restart(
        const std::string &variant) {
      return restart_support_mode_for_variant(variant) == "full_bundle_resume";
    }

    [[nodiscard]] std::string active_bundle_variant() const {
      if (!active_response_bundle_.empty()) {
        return response_variant_name(active_response_bundle_.front());
      }
      return "unknown";
    }

    [[nodiscard]] bool active_bundle_matches_input(
        const ExcitedBundleProtocolInput &input) const {
      return !active_response_bundle_.empty() &&
             active_response_bundle_.size() == input.num_states &&
             active_bundle_variant() == expected_response_variant(input.tda);
    }

    struct IterationContract {
      double density_target = 0.0;
      double relative_target = 0.0;
      double max_rotation = 0.0;
      size_t accelerator_subspace = 0;
      bool use_dual_gate = false;
      bool accelerator_active = false;
      bool step_restriction_active = false;
      std::string convergence_mode = "max_residual";
      std::string accelerator_mode = "none";
    };

    [[nodiscard]] static double legacy_max_rotation(double threshold) {
      double max_rotation = 0.5;
      if (threshold >= 1.0e-2) {
        max_rotation = 2.0;
      } else if (threshold >= 1.0e-4) {
        max_rotation = 0.25;
      } else if (threshold >= 1.0e-6) {
        max_rotation = 0.10;
      } else if (threshold >= 1.0e-8) {
        max_rotation = 0.05;
      }
      return max_rotation;
    }

    [[nodiscard]] static IterationContract
    make_iteration_contract(const ExcitedBundleProtocolInput &input,
                            const std::string &variant) {
      IterationContract contract;
      contract.density_target =
          std::max(100.0 * input.threshold, input.dconv);
      contract.relative_target =
          std::max(50.0 * input.threshold, 0.5 * input.dconv);
      contract.max_rotation = legacy_max_rotation(input.threshold);
      contract.accelerator_subspace = input.maxsub;

      if (response_variant_is_unrestricted(variant)) {
        contract.use_dual_gate = false;
        contract.convergence_mode = "max_residual_unrestricted_scaffold";
        contract.accelerator_active = false;
        contract.step_restriction_active = false;
        contract.accelerator_mode =
            (input.maxsub > 1) ? "unsupported_unrestricted_variant"
                               : "disabled_maxsub_le_1";
      } else {
        contract.use_dual_gate = true;
        contract.convergence_mode = "density_relative_dual_gate";
        contract.accelerator_active = input.maxsub > 1;
        // Legacy excited-state update_response computes max_rotation but leaves
        // step restriction disabled in the main bundle path.
        contract.step_restriction_active = false;
        contract.accelerator_mode =
            contract.accelerator_active ? "kain_per_root"
                                        : "disabled_maxsub_le_1";
      }
      return contract;
    }

    void configure_iteration_contract(const ExcitedBundleProtocolInput &input,
                                      const std::string &variant) {
      const auto contract = make_iteration_contract(input, variant);
      convergence_mode_ = contract.convergence_mode;
      accelerator_mode_ = contract.accelerator_mode;
      accelerator_subspace_ = contract.accelerator_subspace;
      density_convergence_target_ = contract.density_target;
      relative_convergence_target_ = contract.relative_target;
      max_rotation_ = contract.max_rotation;
    }

    void populate_snapshot_header(madness::World &world,
                                  RestartSnapshot &snapshot,
                                  const ExcitedBundleProtocolInput &input,
                                  const std::string &snapshot_kind) const {
      snapshot.schema_version = 2;
      snapshot.snapshot_kind = snapshot_kind;
      snapshot.protocol_index = input.protocol_index;
      snapshot.protocol_threshold = input.threshold;
      snapshot.protocol_k = current_protocol_k();
      snapshot.owner_group = input.owner_group;
      snapshot.archive_world_size = static_cast<size_t>(world.size());
      snapshot.state_count = omega_.size();
      snapshot.num_orbitals = ground_num_orbitals_;
      snapshot.response_variant = active_bundle_variant();
      if (snapshot.response_variant == "unknown") {
        snapshot.response_variant = expected_response_variant(input.tda);
      }
      snapshot.restart_support_mode =
          restart_support_mode_for_variant(snapshot.response_variant);
      snapshot.bundle_state_present = !active_response_bundle_.empty();
      snapshot.restart_capable =
          snapshot.bundle_state_present &&
          variant_supports_full_restart(snapshot.response_variant);
      snapshot.bundle_k = active_bundle_k_;
      snapshot.bundle_num_orbitals =
          (snapshot.bundle_state_present && !active_response_bundle_.empty())
              ? response_num_orbitals(active_response_bundle_.front())
              : ground_num_orbitals_;
      snapshot.response_bundle = active_response_bundle_;
    }

    void write_guess_archive(madness::World &world, bool converged,
                             size_t iterations,
                             const ExcitedBundleProtocolInput &input) const {
      RestartSnapshot snapshot;
      snapshot.has_data = true;
      snapshot.converged = converged;
      snapshot.iterations = iterations;
      populate_snapshot_header(
          world, snapshot, input,
          active_response_bundle_.empty() ? "guess_trial_seed" : "guess_bundle");
      snapshot.energies = omega_;
      snapshot.state_names = state_names_;
      snapshot.roots = stable_sorted_root_descriptors();
      snapshot.slot_permutation = slot_permutation_;
      snapshot.residual_norms = residual_norms_;
      snapshot.density_change_norms = density_change_norms_;
      snapshot.relative_residual_norms = relative_residual_norms_;
      snapshot.iteration_max_residuals = iteration_max_residuals_;
      snapshot.iteration_max_density_changes = iteration_max_density_changes_;
      snapshot.iteration_max_relative_residuals =
          iteration_max_relative_residuals_;
      snapshot.convergence_mode = convergence_mode_;
      snapshot.accelerator_mode = accelerator_mode_;
      snapshot.accelerator_subspace = accelerator_subspace_;
      snapshot.density_convergence_target = density_convergence_target_;
      snapshot.relative_convergence_target = relative_convergence_target_;
      snapshot.max_rotation = max_rotation_;
      snapshot.has_trial_states = !trial_space_.x_states.empty();
      snapshot.trial_space_k = trial_space_k_;
      snapshot.trial_space_num_orbitals = ground_num_orbitals_;
      snapshot.trial_states = trial_space_.x_states;
      write_restart_snapshot(world, guess_archive_file(), snapshot);
    }

    static int protocol_k(double thresh) {
      if (thresh >= 0.9e-2)
        return 4;
      if (thresh >= 0.9e-4)
        return 6;
      if (thresh >= 0.9e-6)
        return 8;
      if (thresh >= 0.9e-8)
        return 10;
      return 12;
    }

    void prepare_protocol(madness::World &world, double thresh) const {
      // Basic protocol prep that mirrors the old solve-path intent
      // (set thresh + k before initialize/iterate). Full operator/ground setup
      // will move here during legacy solver translation.
      madness::FunctionDefaults<3>::set_k(protocol_k(thresh));
      madness::FunctionDefaults<3>::set_thresh(thresh);
      if (world.rank() == 0) {
        madness::print("EXCITED_PROTOCOL_PREP threshold=", thresh,
                       " k=", madness::FunctionDefaults<3>::get_k());
      }
      world.gop.fence();
    }

    [[nodiscard]] static int current_protocol_k() {
      return madness::FunctionDefaults<3>::get_k();
    }

    [[nodiscard]] static double current_protocol_thresh() {
      return madness::FunctionDefaults<3>::get_thresh();
    }

    [[nodiscard]] bool snapshot_variant_matches_input(
        const RestartSnapshot &snapshot,
        const ExcitedBundleProtocolInput &input) const {
      return response_variant_is_valid(snapshot.response_variant) &&
             snapshot.response_variant == expected_response_variant(input.tda);
    }

    [[nodiscard]] bool snapshot_supports_full_restart(
        const RestartSnapshot &snapshot,
        const ExcitedBundleProtocolInput &input) const {
      return snapshot.has_data && snapshot.bundle_state_present &&
             snapshot.restart_capable &&
             snapshot_variant_matches_input(snapshot, input) &&
             snapshot.response_bundle.size() == input.num_states &&
             snapshot.bundle_num_orbitals == ground_num_orbitals_;
    }

    [[nodiscard]] bool snapshot_can_seed_guess(
        const RestartSnapshot &snapshot,
        const ExcitedBundleProtocolInput &input) const {
      if (!snapshot.has_data) {
        return false;
      }
      if (snapshot.bundle_state_present && snapshot_variant_matches_input(snapshot, input) &&
          snapshot.response_bundle.size() == input.num_states &&
          snapshot.bundle_num_orbitals == ground_num_orbitals_) {
        return true;
      }
      if (snapshot.has_trial_states && snapshot.trial_states.size() >= input.num_states) {
        return true;
      }
      return !snapshot.energies.empty();
    }

    void load_restart_seed(const RestartSnapshot &restart_snapshot) {
      omega_ = restart_snapshot.energies;
      if (omega_.size() != restart_snapshot.state_count &&
          restart_snapshot.state_count > 0) {
        omega_.resize(restart_snapshot.state_count);
      }
      state_names_ = restart_snapshot.state_names;
      if (state_names_.size() != omega_.size()) {
        state_names_.clear();
      }
      root_descriptors_ = restart_snapshot.roots;
      std::sort(root_descriptors_.begin(), root_descriptors_.end(),
                [](const ExcitedRootDescriptor &lhs,
                   const ExcitedRootDescriptor &rhs) {
                  if (lhs.slot_index != rhs.slot_index) {
                    return lhs.slot_index < rhs.slot_index;
                  }
                  return lhs.stable_index < rhs.stable_index;
                });
      slot_permutation_ = restart_snapshot.slot_permutation;
      next_root_stable_index_ = 0;
      for (const auto &root : root_descriptors_) {
        next_root_stable_index_ =
            std::max(next_root_stable_index_, root.stable_index + 1);
      }
      if (root_descriptors_.empty() && !omega_.empty()) {
        root_descriptors_.reserve(omega_.size());
        for (size_t i = 0; i < omega_.size(); ++i) {
          const std::string display_name =
              (i < state_names_.size()) ? state_names_[i] : std::string{};
          root_descriptors_.push_back(
              {make_excited_root_id(i), i, i, omega_[i], display_name});
        }
        next_root_stable_index_ = omega_.size();
      }
      residual_norms_ = restart_snapshot.residual_norms;
      if (residual_norms_.size() != omega_.size()) {
        residual_norms_.assign(omega_.size(), 1.0e-2);
      }
      density_change_norms_ = restart_snapshot.density_change_norms;
      if (density_change_norms_.size() != omega_.size()) {
        density_change_norms_.assign(omega_.size(), 0.0);
      }
      relative_residual_norms_ = restart_snapshot.relative_residual_norms;
      if (relative_residual_norms_.size() != omega_.size()) {
        relative_residual_norms_.assign(omega_.size(), 0.0);
      }
      last_iteration_count_ = restart_snapshot.iterations;
      iteration_max_residuals_ = restart_snapshot.iteration_max_residuals;
      iteration_max_density_changes_ =
          restart_snapshot.iteration_max_density_changes;
      iteration_max_relative_residuals_ =
          restart_snapshot.iteration_max_relative_residuals;
      convergence_mode_ = restart_snapshot.convergence_mode;
      accelerator_mode_ = restart_snapshot.accelerator_mode;
      accelerator_subspace_ = restart_snapshot.accelerator_subspace;
      density_convergence_target_ =
          restart_snapshot.density_convergence_target;
      relative_convergence_target_ =
          restart_snapshot.relative_convergence_target;
      max_rotation_ = restart_snapshot.max_rotation;

      active_response_bundle_.clear();
      if (restart_snapshot.bundle_state_present &&
          !restart_snapshot.response_bundle.empty()) {
        active_response_bundle_ = restart_snapshot.response_bundle;
        active_bundle_k_ = restart_snapshot.bundle_k;
        active_bundle_thresh_ = restart_snapshot.protocol_threshold;
      } else {
        active_bundle_k_ = -1;
        active_bundle_thresh_ = std::numeric_limits<double>::infinity();
      }

      if (restart_snapshot.has_trial_states &&
          !restart_snapshot.trial_states.empty()) {
        trial_space_.x_states = restart_snapshot.trial_states;
        trial_space_.num_states = trial_space_.x_states.size();
        trial_space_.num_orbitals = restart_snapshot.trial_space_num_orbitals;
        trial_space_.omega = omega_;
        trial_space_k_ = restart_snapshot.trial_space_k;
        trial_space_thresh_ = restart_snapshot.protocol_threshold;
      } else if (!active_response_bundle_.empty()) {
        trial_space_.x_states.clear();
        trial_space_.num_states = 0;
        trial_space_.num_orbitals = ground_num_orbitals_;
        trial_space_.omega = omega_;
      }
      has_active_guess_ = !omega_.empty();
    }

    void normalize_root_descriptors() {
      std::set<size_t> used_stable_indices;
      size_t next_index = next_root_stable_index_;
      for (const auto &root : root_descriptors_) {
        next_index = std::max(next_index, root.stable_index + 1);
      }
      for (auto &root : root_descriptors_) {
        bool reassigned_stable_index = false;
        if (used_stable_indices.count(root.stable_index) > 0) {
          root.stable_index = next_index++;
          reassigned_stable_index = true;
        }
        used_stable_indices.insert(root.stable_index);
        if (root.root_id.empty() || reassigned_stable_index) {
          root.root_id = make_excited_root_id(root.stable_index);
        }
      }
      next_root_stable_index_ = next_index;
    }

    void sync_root_descriptors_for_protocol(double threshold) {
      if (omega_.empty()) {
        root_descriptors_.clear();
        state_names_.clear();
        slot_permutation_.clear();
        return;
      }

      if (root_descriptors_.size() != omega_.size()) {
        std::vector<ExcitedRootDescriptor> resized;
        const size_t preserve = std::min(root_descriptors_.size(), omega_.size());
        resized.reserve(omega_.size());
        for (size_t i = 0; i < preserve; ++i) {
          resized.push_back(root_descriptors_[i]);
        }
        for (size_t i = preserve; i < omega_.size(); ++i) {
          ExcitedRootDescriptor root;
          root.stable_index = next_root_stable_index_++;
          root.root_id = make_excited_root_id(root.stable_index);
          resized.push_back(std::move(root));
        }
        root_descriptors_ = std::move(resized);
      }

      normalize_root_descriptors();

      state_names_ = assign_excited_state_names(omega_, threshold);
      slot_permutation_.clear();
      slot_permutation_.reserve(root_descriptors_.size());
      for (size_t slot = 0; slot < root_descriptors_.size(); ++slot) {
        auto &root = root_descriptors_[slot];
        root.slot_index = slot;
        root.energy = (slot < omega_.size()) ? omega_[slot] : 0.0;
        if (root.display_name.empty()) {
          root.display_name = (slot < state_names_.size()) ? state_names_[slot]
                                                           : std::string{};
        }
        slot_permutation_.push_back(root.stable_index);
      }
    }

    [[nodiscard]] std::vector<ExcitedRootDescriptor>
    stable_sorted_root_descriptors() const {
      auto roots = root_descriptors_;
      std::sort(roots.begin(), roots.end(),
                [](const ExcitedRootDescriptor &lhs,
                   const ExcitedRootDescriptor &rhs) {
                  if (lhs.stable_index != rhs.stable_index) {
                    return lhs.stable_index < rhs.stable_index;
                  }
                  return lhs.root_id < rhs.root_id;
                });
      return roots;
    }

    void apply_slot_reordering(const std::vector<size_t> &new_slot_to_old_slot,
                               double threshold) {
      if (omega_.empty()) {
        root_descriptors_.clear();
        state_names_.clear();
        slot_permutation_.clear();
        return;
      }
      if (root_descriptors_.empty() ||
          new_slot_to_old_slot.size() != root_descriptors_.size()) {
        sync_root_descriptors_for_protocol(threshold);
        return;
      }

      const auto previous = root_descriptors_;
      std::vector<ExcitedRootDescriptor> reordered(root_descriptors_.size());
      std::vector<bool> consumed(previous.size(), false);
      size_t next_unused = 0;
      auto take_next_unused = [&]() -> ExcitedRootDescriptor {
        while (next_unused < previous.size() && consumed[next_unused]) {
          ++next_unused;
        }
        if (next_unused < previous.size()) {
          consumed[next_unused] = true;
          return previous[next_unused++];
        }
        ExcitedRootDescriptor root;
        root.stable_index = next_root_stable_index_++;
        root.root_id = make_excited_root_id(root.stable_index);
        return root;
      };

      for (size_t slot = 0; slot < new_slot_to_old_slot.size(); ++slot) {
        const size_t old_slot = new_slot_to_old_slot[slot];
        if (old_slot < previous.size() && !consumed[old_slot]) {
          reordered[slot] = previous[old_slot];
          consumed[old_slot] = true;
        } else {
          reordered[slot] = take_next_unused();
        }
      }

      root_descriptors_ = std::move(reordered);
      sync_root_descriptors_for_protocol(threshold);
    }

    void ensure_ground_data(madness::World &world) {
      const int target_k = current_protocol_k();
      const double target_thresh = current_protocol_thresh();
      if (ground_data_ != nullptr && ground_data_k_ == target_k &&
          std::abs(ground_data_thresh_ - target_thresh) <= 1.0e-14) {
        const bool prelim_ready =
            (preliminaries_k_ == target_k &&
             std::abs(preliminaries_thresh_ - target_thresh) <= 1.0e-14);
        if (prelim_ready) {
          return;
        }
      }
      if (ground_data_ == nullptr) {
        ground_data_ =
            std::make_unique<GroundStateData>(world, config_.archive_file, Molecule());
      }
      ground_data_->prepareOrbitals(world, target_k, target_thresh);
      const double vtol = std::max(1.0e-12, 0.1 * target_thresh);
      const auto coulop = CoulombOperator(world, 1.0e-8, target_thresh);
      ground_data_->computePreliminaries(world, coulop, vtol, "moldft.fock.json");
      ground_num_orbitals_ = static_cast<size_t>(ground_data_->getNumOrbitals());
      ground_energies_ = ground_data_->getEnergies();
      ground_orbitals_ = copy(world, ground_data_->getOrbitals(), true);
      ground_data_k_ = target_k;
      ground_data_thresh_ = target_thresh;
      preliminaries_k_ = target_k;
      preliminaries_thresh_ = target_thresh;
      world.gop.fence();
    }

    [[nodiscard]] static double state_inner(madness::World &world,
                                            const vector_real_function_3d &a,
                                            const vector_real_function_3d &b) {
      return ResponseSolverUtils::inner(world, a, b);
    }

    [[nodiscard]] static double state_norm(madness::World &world,
                                           const vector_real_function_3d &a) {
      return std::sqrt(std::max(0.0, state_inner(world, a, a)));
    }

    void project_and_orthonormalize(madness::World &world,
                                    std::vector<vector_real_function_3d> &states) {
      if (states.empty()) {
        return;
      }
      QProjector<double, 3> projector(ground_orbitals_);
      for (auto &state : states) {
        state = projector(state);
      }

      constexpr double tiny_norm = 1.0e-11;
      std::vector<vector_real_function_3d> filtered;
      filtered.reserve(states.size());
      for (auto &state : states) {
        if (state_norm(world, state) > tiny_norm) {
          filtered.push_back(state);
        }
      }
      states.swap(filtered);
      if (states.empty()) {
        return;
      }

      for (size_t i = 0; i < states.size(); ++i) {
        const double ni = state_norm(world, states[i]);
        if (ni <= tiny_norm) {
          continue;
        }
        scale(world, states[i], 1.0 / ni, false);
        for (size_t j = i + 1; j < states.size(); ++j) {
          const double proj = state_inner(world, states[i], states[j]);
          gaxpy(world, 1.0, states[j], -proj, states[i], false);
        }
      }
      for (auto &state : states) {
        const double ni = state_norm(world, state);
        if (ni > tiny_norm) {
          scale(world, state, 1.0 / ni, false);
        }
      }
      for (auto &state : states) {
        truncate(world, state, madness::FunctionDefaults<3>::get_thresh());
      }
      world.gop.fence();
    }

    [[nodiscard]] std::vector<double>
    estimate_state_energies(madness::World &world,
                            const std::vector<vector_real_function_3d> &states) const {
      std::vector<double> energies(states.size(), 0.0);
      for (size_t i = 0; i < states.size(); ++i) {
        double weighted = 0.0;
        double denom = 0.0;
        for (size_t orb = 0; orb < states[i].size(); ++orb) {
          const double coeff_norm2 =
              std::max(0.0, states[i][orb].inner(states[i][orb]));
          const double e_occ = std::abs(ground_energies_(long(orb)));
          weighted += e_occ * coeff_norm2;
          denom += coeff_norm2;
        }
        energies[i] = (denom > 1.0e-14)
                          ? (weighted / denom + 0.02 * static_cast<double>(i + 1))
                          : (0.10 + 0.05 * static_cast<double>(i + 1));
      }
      return energies;
    }

    struct StateUpdateMetrics {
      double residual = 0.0;
      double density_change = 0.0;
      double relative_residual = 0.0;
      double step_norm = 0.0;
      bool accelerator_applied = false;
      bool step_restricted = false;
    };

    static void align_component_to_ground(
        madness::World &world, vector_real_function_3d &component,
        size_t n_orbitals) {
      if (component.size() < n_orbitals) {
        auto zeros = zero_functions_compressed<double, 3>(
            world, static_cast<int>(n_orbitals - component.size()));
        component.insert(component.end(), zeros.begin(), zeros.end());
      } else if (component.size() > n_orbitals) {
        component.resize(n_orbitals);
      }
      auto zero_replacements = zero_functions_compressed<double, 3>(
          world, static_cast<int>(component.size()));
      for (size_t i = 0; i < component.size(); ++i) {
        if (!component[i].is_initialized()) {
          component[i] = zero_replacements[i];
        }
      }
    }

    static void ensure_initialized_flat(madness::World &world,
                                        vector_real_function_3d &functions) {
      if (functions.empty()) {
        return;
      }
      auto zeros = zero_functions_compressed<double, 3>(
          world, static_cast<int>(functions.size()));
      for (size_t i = 0; i < functions.size(); ++i) {
        if (!functions[i].is_initialized()) {
          functions[i] = zeros[i];
        }
      }
    }

    template <typename ResponseType>
    void align_response_to_ground(madness::World &world,
                                  ResponseType &response) const {
      align_component_to_ground(world, response.x_alpha, ground_num_orbitals_);
      if constexpr (std::is_same_v<ResponseType, DynamicRestrictedResponse> ||
                    std::is_same_v<ResponseType, DynamicUnrestrictedResponse>) {
        align_component_to_ground(world, response.y_alpha, ground_num_orbitals_);
      }
      if constexpr (std::is_same_v<ResponseType, StaticUnrestrictedResponse> ||
                    std::is_same_v<ResponseType, DynamicUnrestrictedResponse>) {
        align_component_to_ground(world, response.x_beta, ground_num_orbitals_);
      }
      if constexpr (std::is_same_v<ResponseType, DynamicUnrestrictedResponse>) {
        align_component_to_ground(world, response.y_beta, ground_num_orbitals_);
      }
      response.flatten();
    }

    template <typename ResponseType>
    [[nodiscard]] static size_t
    count_uninitialized_flat(const ResponseType &response) {
      size_t count = 0;
      for (const auto &f : response.flat) {
        if (!f.is_initialized()) {
          ++count;
        }
      }
      return count;
    }

    template <typename ResponseType>
    [[nodiscard]] vector_real_function_3d
    apply_hamiltonian_no_diag(madness::World &world,
                              const ResponseType &response) const {
      if constexpr (std::is_same_v<ResponseType, StaticRestrictedResponse>) {
        return transform(world, response.x_alpha, ground_data_->Hamiltonian_no_diag,
                         true);
      } else if constexpr (std::is_same_v<ResponseType,
                                          DynamicRestrictedResponse>) {
        auto epsilon = transform(world, response.x_alpha,
                                 ground_data_->Hamiltonian_no_diag, true);
        auto eps_y = transform(world, response.y_alpha,
                               ground_data_->Hamiltonian_no_diag, true);
        epsilon.insert(epsilon.end(), eps_y.begin(), eps_y.end());
        return epsilon;
      } else if constexpr (std::is_same_v<ResponseType,
                                          StaticUnrestrictedResponse>) {
        auto epsilon = transform(world, response.x_alpha,
                                 ground_data_->Hamiltonian_no_diag, true);
        auto eps_b = transform(world, response.x_beta,
                               ground_data_->Hamiltonian_no_diag, true);
        epsilon.insert(epsilon.end(), eps_b.begin(), eps_b.end());
        return epsilon;
      } else {
        auto epsilon = transform(world, response.x_alpha,
                                 ground_data_->Hamiltonian_no_diag, true);
        auto eps_ya = transform(world, response.y_alpha,
                                ground_data_->Hamiltonian_no_diag, true);
        auto eps_xb = transform(world, response.x_beta,
                                ground_data_->Hamiltonian_no_diag, true);
        auto eps_yb = transform(world, response.y_beta,
                                ground_data_->Hamiltonian_no_diag, true);
        epsilon.insert(epsilon.end(), eps_ya.begin(), eps_ya.end());
        epsilon.insert(epsilon.end(), eps_xb.begin(), eps_xb.end());
        epsilon.insert(epsilon.end(), eps_yb.begin(), eps_yb.end());
        return epsilon;
      }
    }

    template <typename ResponseType>
    [[nodiscard]] real_function_3d
    compute_response_density(madness::World &world,
                             const ResponseType &response) const {
      if constexpr (std::is_same_v<ResponseType, StaticRestrictedResponse> ||
                    std::is_same_v<ResponseType, DynamicRestrictedResponse>) {
        return compute_density(world, response, ground_orbitals_);
      } else if constexpr (std::is_same_v<ResponseType,
                                          StaticUnrestrictedResponse>) {
        auto xphi_a = mul(world, response.x_alpha, ground_orbitals_, true);
        auto xphi_b = mul(world, response.x_beta, ground_orbitals_, true);
        return sum(world, xphi_a, true) + sum(world, xphi_b, true);
      } else {
        auto alpha_xy =
            gaxpy_oop(1.0, response.x_alpha, 1.0, response.y_alpha, true);
        auto beta_xy =
            gaxpy_oop(1.0, response.x_beta, 1.0, response.y_beta, true);
        auto xphi_a = mul(world, alpha_xy, ground_orbitals_, true);
        auto xphi_b = mul(world, beta_xy, ground_orbitals_, true);
        return sum(world, xphi_a, true) + sum(world, xphi_b, true);
      }
    }

    template <typename ResponseType>
    [[nodiscard]] double response_metric_norm2(
        madness::World &world, const ResponseType &response) const {
      double metric = ResponseSolverUtils::inner(world, response.x_alpha,
                                                 response.x_alpha);
      if constexpr (std::is_same_v<ResponseType, StaticUnrestrictedResponse> ||
                    std::is_same_v<ResponseType, DynamicUnrestrictedResponse>) {
        metric += ResponseSolverUtils::inner(world, response.x_beta,
                                             response.x_beta);
      }
      if constexpr (std::is_same_v<ResponseType, DynamicRestrictedResponse>) {
        metric -= ResponseSolverUtils::inner(world, response.y_alpha,
                                             response.y_alpha);
      }
      if constexpr (std::is_same_v<ResponseType, DynamicUnrestrictedResponse>) {
        metric -= ResponseSolverUtils::inner(world, response.y_alpha,
                                             response.y_alpha);
        metric -= ResponseSolverUtils::inner(world, response.y_beta,
                                             response.y_beta);
      }
      return metric;
    }

    template <typename ResponseType>
    [[nodiscard]] double response_metric_inner(
        madness::World &world, const ResponseType &lhs,
        const ResponseType &rhs) const {
      double metric =
          ResponseSolverUtils::inner(world, lhs.x_alpha, rhs.x_alpha);
      if constexpr (std::is_same_v<ResponseType, StaticUnrestrictedResponse> ||
                    std::is_same_v<ResponseType, DynamicUnrestrictedResponse>) {
        metric += ResponseSolverUtils::inner(world, lhs.x_beta, rhs.x_beta);
      }
      if constexpr (std::is_same_v<ResponseType, DynamicRestrictedResponse>) {
        metric -= ResponseSolverUtils::inner(world, lhs.y_alpha, rhs.y_alpha);
      }
      if constexpr (std::is_same_v<ResponseType, DynamicUnrestrictedResponse>) {
        metric -= ResponseSolverUtils::inner(world, lhs.y_alpha, rhs.y_alpha);
        metric -= ResponseSolverUtils::inner(world, lhs.y_beta, rhs.y_beta);
      }
      return metric;
    }

    template <typename ResponseType>
    [[nodiscard]] double response_plain_inner(madness::World &world,
                                              const ResponseType &lhs,
                                              const ResponseType &rhs) const {
      double value = ResponseSolverUtils::inner(world, lhs.x_alpha, rhs.x_alpha);
      if constexpr (std::is_same_v<ResponseType, DynamicRestrictedResponse> ||
                    std::is_same_v<ResponseType, DynamicUnrestrictedResponse>) {
        value += ResponseSolverUtils::inner(world, lhs.y_alpha, rhs.y_alpha);
      }
      if constexpr (std::is_same_v<ResponseType, StaticUnrestrictedResponse> ||
                    std::is_same_v<ResponseType, DynamicUnrestrictedResponse>) {
        value += ResponseSolverUtils::inner(world, lhs.x_beta, rhs.x_beta);
      }
      if constexpr (std::is_same_v<ResponseType, DynamicUnrestrictedResponse>) {
        value += ResponseSolverUtils::inner(world, lhs.y_beta, rhs.y_beta);
      }
      return value;
    }

    template <typename ResponseType>
    void normalize_response_metric(madness::World &world,
                                   ResponseType &response) const {
      const double metric = response_metric_norm2(world, response);
      if (!std::isfinite(metric) || metric <= 1.0e-14) {
        return;
      }
      const double inv_norm = 1.0 / std::sqrt(metric);
      scale(world, response.x_alpha, inv_norm, false);
      if constexpr (std::is_same_v<ResponseType, DynamicRestrictedResponse> ||
                    std::is_same_v<ResponseType, DynamicUnrestrictedResponse>) {
        scale(world, response.y_alpha, inv_norm, false);
      }
      if constexpr (std::is_same_v<ResponseType, StaticUnrestrictedResponse> ||
                    std::is_same_v<ResponseType, DynamicUnrestrictedResponse>) {
        scale(world, response.x_beta, inv_norm, false);
      }
      if constexpr (std::is_same_v<ResponseType, DynamicUnrestrictedResponse>) {
        scale(world, response.y_beta, inv_norm, false);
      }
      response.flatten();
    }

    template <typename ResponseType>
    void project_response_channels(madness::World &world,
                                   ResponseType &response) const {
      QProjector<double, 3> projector(ground_orbitals_);
      response.x_alpha = projector(response.x_alpha);
      if constexpr (std::is_same_v<ResponseType, DynamicRestrictedResponse> ||
                    std::is_same_v<ResponseType, DynamicUnrestrictedResponse>) {
        response.y_alpha = projector(response.y_alpha);
      }
      if constexpr (std::is_same_v<ResponseType, StaticUnrestrictedResponse> ||
                    std::is_same_v<ResponseType, DynamicUnrestrictedResponse>) {
        response.x_beta = projector(response.x_beta);
      }
      if constexpr (std::is_same_v<ResponseType, DynamicUnrestrictedResponse>) {
        response.y_beta = projector(response.y_beta);
      }
      response.flatten();
    }

    template <typename ResponseType> struct BundleResponsePotentials {
      std::vector<ResponseType> lambda;
      std::vector<ResponseType> v0;
      std::vector<ResponseType> gamma;
    };

    template <typename ResponseType>
    [[nodiscard]] ResponseType
    make_response_from_flat(madness::World &world, const ResponseType &prototype,
                            vector_real_function_3d flat) const {
      auto response = prototype;
      response.flat = std::move(flat);
      ensure_initialized_flat(world, response.flat);
      response.sync();
      align_response_to_ground(world, response);
      return response;
    }

    template <typename ResponseType>
    void compute_state_response_potentials(
        madness::World &world, const ResponseType &state, ResponseType &lambda_out,
        ResponseType &v0_out, ResponseType &gamma_out) const {
      auto current = state;
      align_response_to_ground(world, current);
      if (current.flat.empty()) {
        lambda_out = current;
        v0_out = current;
        gamma_out = current;
        return;
      }

      std::vector<int> orbital_index(current.flat.size());
      std::iota(orbital_index.begin(), orbital_index.end(), 0);
      std::vector<int> state_index(current.flat.size(), 0);

      ResponseComputeGroundExchange ground_exchange_task_impl;
      MacroTask ground_exchange_task(world, ground_exchange_task_impl);
      ResponseComputeGammaX response_exchange_task_impl;
      MacroTask response_exchange_task(world, response_exchange_task_impl);

      const bool use_tda_exchange = !response_has_y_channel_v<ResponseType>;
      auto k0 = ground_exchange_task(orbital_index, state_index, current.flat,
                                     ground_orbitals_, use_tda_exchange);
      auto gx = response_exchange_task(orbital_index, state_index, current.flat,
                                       ground_orbitals_, use_tda_exchange);

      const double c_xc = ground_data_->xcf_.hf_exchange_coefficient();
      auto v_local = ground_data_->V_local * current.flat;
      auto v0_flat = v_local - c_xc * k0;
      auto e0_flat = apply_hamiltonian_no_diag(world, current);
      auto lambda_flat = v0_flat - e0_flat + gx;

      v0_out = make_response_from_flat(world, current, std::move(v0_flat));
      gamma_out = make_response_from_flat(world, current, std::move(gx));
      lambda_out =
          make_response_from_flat(world, current, std::move(lambda_flat));
    }

    template <typename ResponseType>
    [[nodiscard]] BundleResponsePotentials<ResponseType>
    compute_bundle_response_potentials(
        madness::World &world, const std::vector<ResponseType> &states,
        size_t iter) const {
      BundleResponsePotentials<ResponseType> potentials;
      potentials.lambda.resize(states.size());
      potentials.v0.resize(states.size());
      potentials.gamma.resize(states.size());
      for (size_t i = 0; i < states.size(); ++i) {
        compute_state_response_potentials(world, states[i], potentials.lambda[i],
                                          potentials.v0[i], potentials.gamma[i]);
      }
      if (world.rank() == 0 && config_.print_level > 1) {
        madness::print("EXCITED_POTENTIALS iter=", iter,
                       " states=", states.size());
      }
      return potentials;
    }

    template <typename ResponseType>
    void build_rotation_matrices(
        madness::World &world, const std::vector<ResponseType> &states,
        const std::vector<ResponseType> &lambda_states, Tensor<double> &S,
        Tensor<double> &A) const {
      const size_t n_states = states.size();
      S = Tensor<double>(n_states, n_states);
      A = Tensor<double>(n_states, n_states);
      for (size_t i = 0; i < n_states; ++i) {
        for (size_t j = i; j < n_states; ++j) {
          const double sij = response_metric_inner(world, states[i], states[j]);
          S(i, j) = sij;
          S(j, i) = sij;
          const double aij =
              response_plain_inner(world, states[i], lambda_states[j]);
          const double aji =
              response_plain_inner(world, states[j], lambda_states[i]);
          A(i, j) = aij;
          A(j, i) = aji;
        }
      }
      S = 0.5 * (S + transpose(S));
      A = 0.5 * (A + transpose(A));
    }

    struct OverlapConditioningResult {
      Tensor<double> overlap;
      size_t floored_singular_values = 0;
      double floor = 0.0;
      double min_singular_value = 0.0;
      bool applied = false;
      bool success = false;
    };

    [[nodiscard]] static OverlapConditioningResult
    condition_excited_overlap_matrix(const Tensor<double> &overlap,
                                     double thresh_degenerate) {
      OverlapConditioningResult result;
      const size_t n_states = static_cast<size_t>(overlap.dim(0));
      if (n_states == 0) {
        result.success = true;
        result.overlap = overlap;
        return result;
      }

      Tensor<double> left_vectors;
      Tensor<double> singular_values;
      Tensor<double> right_vectors;
      auto overlap_copy = copy(overlap);
      try {
        svd(overlap_copy, left_vectors, singular_values, right_vectors);
      } catch (...) {
        return result;
      }

      const double floor =
          std::max(1.0e-12, 10.0 * std::max(thresh_degenerate, 1.0e-14));
      result.floor = floor;
      result.min_singular_value =
          singular_values.dim(0) > 0 ? singular_values(long(0)) : 0.0;
      if (singular_values.dim(0) > 0) {
        result.min_singular_value = singular_values(long(0));
        for (int64_t i = 1; i < singular_values.dim(0); ++i) {
          result.min_singular_value =
              std::min(result.min_singular_value, singular_values(i));
        }
      }

      std::vector<double> clipped_singulars(n_states, floor);
      for (size_t i = 0; i < n_states && i < static_cast<size_t>(singular_values.dim(0));
           ++i) {
        const double sigma = singular_values(long(i));
        if (!std::isfinite(sigma) || sigma < floor) {
          clipped_singulars[i] = floor;
          ++result.floored_singular_values;
        } else {
          clipped_singulars[i] = sigma;
        }
      }

      Tensor<double> conditioned(n_states, n_states);
      for (size_t row = 0; row < n_states; ++row) {
        for (size_t col = 0; col < n_states; ++col) {
          double value = 0.0;
          for (size_t k = 0; k < n_states; ++k) {
            value += left_vectors(long(row), long(k)) * clipped_singulars[k] *
                     left_vectors(long(col), long(k));
          }
          conditioned(long(row), long(col)) = value;
        }
      }

      result.overlap = 0.5 * (conditioned + transpose(conditioned));
      result.applied = result.floored_singular_values > 0;
      result.success = true;
      return result;
    }

    template <typename ResponseType>
    void rotate_bundle_component(
        madness::World &world, const std::vector<ResponseType> &source_states,
        std::vector<ResponseType> &rotated_states, const Tensor<double> &U,
        vector_real_function_3d ResponseType::*component_ptr) const {
      const size_t n_states = source_states.size();
      if (n_states == 0) {
        return;
      }
      const size_t n_orbitals = (source_states[0].*component_ptr).size();
      for (size_t orb = 0; orb < n_orbitals; ++orb) {
        vector_real_function_3d by_state(n_states);
        for (size_t s = 0; s < n_states; ++s) {
          by_state[s] = (source_states[s].*component_ptr)[orb];
        }
        const auto rotated_orbital = transform(world, by_state, U, false);
        for (size_t s = 0; s < n_states; ++s) {
          (rotated_states[s].*component_ptr)[orb] = rotated_orbital[s];
        }
      }
    }

    template <typename ResponseType>
    [[nodiscard]] std::vector<ResponseType>
    rotate_bundle_states(madness::World &world,
                         const std::vector<ResponseType> &states,
                         const Tensor<double> &U) const {
      std::vector<ResponseType> rotated = states;
      for (auto &state : rotated) {
        align_response_to_ground(world, state);
      }
      rotate_bundle_component(world, states, rotated, U, &ResponseType::x_alpha);
      if constexpr (std::is_same_v<ResponseType, DynamicRestrictedResponse> ||
                    std::is_same_v<ResponseType, DynamicUnrestrictedResponse>) {
        rotate_bundle_component(world, states, rotated, U,
                                &ResponseType::y_alpha);
      }
      if constexpr (std::is_same_v<ResponseType, StaticUnrestrictedResponse> ||
                    std::is_same_v<ResponseType, DynamicUnrestrictedResponse>) {
        rotate_bundle_component(world, states, rotated, U, &ResponseType::x_beta);
      }
      if constexpr (std::is_same_v<ResponseType, DynamicUnrestrictedResponse>) {
        rotate_bundle_component(world, states, rotated, U, &ResponseType::y_beta);
      }
      for (auto &state : rotated) {
        state.flatten();
        align_response_to_ground(world, state);
      }
      world.gop.fence();
      return rotated;
    }

    [[nodiscard]] bool diagonalize_excited_bundle(
        madness::World &world, Tensor<double> S, Tensor<double> A,
        std::vector<double> &omega, Tensor<double> &U_out,
        std::vector<size_t> *slot_order_out = nullptr) const {
      const size_t n_states = static_cast<size_t>(S.dim(0));
      if (n_states == 0) {
        return false;
      }

      const double diag_floor =
          std::max(1.0e-12, 10.0 * madness::FunctionDefaults<3>::get_thresh());
      for (size_t i = 0; i < n_states; ++i) {
        if (!std::isfinite(S(i, i)) || std::abs(S(i, i)) < diag_floor) {
          S(i, i) = diag_floor;
        }
      }
      S = 0.5 * (S + transpose(S));
      A = 0.5 * (A + transpose(A));
      const auto response_input = copy(A);

      Tensor<double> eigenvectors(n_states, n_states);
      Tensor<double> eigenvalues(n_states);
      bool overlap_conditioned_retry = false;
      try {
        sygvp(world, A, S, 1, eigenvectors, eigenvalues);
      } catch (...) {
        const auto conditioning = condition_excited_overlap_matrix(
            S, madness::FunctionDefaults<3>::get_thresh());
        if (!conditioning.success) {
          if (world.rank() == 0 && config_.print_level > 0) {
            madness::print("EXCITED_ROTATE_FAIL reason=sygvp_exception");
          }
          return false;
        }
        if (world.rank() == 0 && config_.print_level > 0) {
          madness::print("EXCITED_ROTATE_RETRY reason=condition_overlap_svd",
                         " floored_singular_values=",
                         conditioning.floored_singular_values,
                         " overlap_floor=", conditioning.floor,
                         " min_singular_value=",
                         conditioning.min_singular_value);
        }
        try {
          auto conditioned_overlap = conditioning.overlap;
          auto conditioned_response = copy(response_input);
          sygvp(world, conditioned_response, conditioned_overlap, 1,
                eigenvectors, eigenvalues);
          overlap_conditioned_retry = conditioning.applied;
        } catch (...) {
          if (world.rank() == 0 && config_.print_level > 0) {
            madness::print("EXCITED_ROTATE_FAIL reason=sygvp_exception",
                           " retry=condition_overlap_svd");
          }
          return false;
        }
      }

      auto U = eigenvectors;
      omega.resize(n_states);
      for (size_t i = 0; i < n_states; ++i) {
        omega[i] = eigenvalues(long(i));
      }

      // Prefer an eigenvector ordering that keeps large coefficients near the
      // diagonal before final sorting, as in legacy excited_eig.
      bool switched = true;
      while (switched) {
        switched = false;
        for (size_t i = 0; i < n_states; ++i) {
          for (size_t j = i + 1; j < n_states; ++j) {
            const double sold = U(long(i), long(i)) * U(long(i), long(i)) +
                                U(long(j), long(j)) * U(long(j), long(j));
            const double snew = U(long(i), long(j)) * U(long(i), long(j)) +
                                U(long(j), long(i)) * U(long(j), long(i));
            if (snew <= sold) {
              continue;
            }
            for (size_t row = 0; row < n_states; ++row) {
              std::swap(U(long(row), long(i)), U(long(row), long(j)));
            }
            std::swap(omega[i], omega[j]);
            switched = true;
          }
        }
      }

      // Fix phases so diagonal components are positive when possible.
      for (size_t i = 0; i < n_states; ++i) {
        if (U(long(i), long(i)) < 0.0) {
          for (size_t row = 0; row < n_states; ++row) {
            U(long(row), long(i)) *= -1.0;
          }
        }
      }

      // Undo arbitrary rotations inside near-degenerate clusters.
      const double thresh_degenerate =
          madness::FunctionDefaults<3>::get_thresh();
      size_t ilo = 0;
      while (ilo + 1 < n_states) {
        size_t ihi = ilo;
        while ((ihi + 1) < n_states) {
          const double cluster_tol =
              thresh_degenerate * 10.0 * std::max(std::abs(omega[ilo]), 1.0);
          if (std::abs(omega[ilo] - omega[ihi + 1]) >= cluster_tol) {
            break;
          }
          ++ihi;
        }
        const size_t nclus = ihi - ilo + 1;
        if (nclus > 1) {
          Tensor<double> q(nclus, nclus);
          for (size_t r = 0; r < nclus; ++r) {
            for (size_t c = 0; c < nclus; ++c) {
              q(long(r), long(c)) = U(long(ilo + r), long(ilo + c));
            }
          }
          Tensor<double> W(nclus, nclus);
          Tensor<double> sigma(nclus);
          Tensor<double> VH(nclus, nclus);
          try {
            svd(q, W, sigma, VH);
            Tensor<double> q_rot(nclus, nclus);
            for (size_t r = 0; r < nclus; ++r) {
              for (size_t c = 0; c < nclus; ++c) {
                double value = 0.0;
                for (size_t k = 0; k < nclus; ++k) {
                  value += W(long(r), long(k)) * VH(long(k), long(c));
                }
                q_rot(long(c), long(r)) = value; // transpose(W * VH)
              }
            }

            Tensor<double> block(n_states, nclus);
            for (size_t row = 0; row < n_states; ++row) {
              for (size_t c = 0; c < nclus; ++c) {
                block(long(row), long(c)) = U(long(row), long(ilo + c));
              }
            }
            for (size_t row = 0; row < n_states; ++row) {
              for (size_t c = 0; c < nclus; ++c) {
                double value = 0.0;
                for (size_t k = 0; k < nclus; ++k) {
                  value += block(long(row), long(k)) * q_rot(long(k), long(c));
                }
                U(long(row), long(ilo + c)) = value;
              }
            }
          } catch (...) {
            if (world.rank() == 0 && config_.print_level > 0) {
              madness::print("EXCITED_ROTATE_WARN reason=degenerate_cluster_svd_exception",
                             " ilo=", ilo, " ihi=", ihi);
            }
          }
        }
        ilo = ihi + 1;
      }

      std::vector<size_t> order(n_states);
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(), [&](size_t lhs, size_t rhs) {
        return omega[lhs] < omega[rhs];
      });

      Tensor<double> sorted_U(n_states, n_states);
      std::vector<double> sorted_omega(n_states, 0.0);
      for (size_t col = 0; col < n_states; ++col) {
        const size_t old_col = order[col];
        sorted_omega[col] = omega[old_col];
        for (size_t row = 0; row < n_states; ++row) {
          sorted_U(long(row), long(col)) = U(long(row), long(old_col));
        }
        if (sorted_U(long(col), long(col)) < 0.0) {
          for (size_t row = 0; row < n_states; ++row) {
            sorted_U(long(row), long(col)) *= -1.0;
          }
        }
      }

      omega = std::move(sorted_omega);
      if (slot_order_out != nullptr) {
        *slot_order_out = order;
      }
      U_out = std::move(sorted_U);
      if (overlap_conditioned_retry && world.rank() == 0 &&
          config_.print_level > 1) {
        madness::print("EXCITED_ROTATE_RECOVERED reason=condition_overlap_svd",
                       " nstates=", n_states);
      }
      return true;
    }

    template <typename ResponseType>
    [[nodiscard]] bool rotate_excited_bundle_states(
        madness::World &world, std::vector<ResponseType> &states,
        BundleResponsePotentials<ResponseType> &potentials,
        std::vector<double> &omega, size_t iter,
        std::vector<size_t> *slot_order_out = nullptr) const {
      if (states.empty()) {
        return false;
      }

      for (auto &state : states) {
        align_response_to_ground(world, state);
      }

      Tensor<double> S;
      Tensor<double> A;
      build_rotation_matrices(world, states, potentials.lambda, S, A);

      Tensor<double> U;
      std::vector<double> rotated_omega;
      if (!diagonalize_excited_bundle(world, S, A, rotated_omega, U,
                                      slot_order_out)) {
        return false;
      }
      states = rotate_bundle_states(world, states, U);
      potentials.lambda = rotate_bundle_states(world, potentials.lambda, U);
      potentials.v0 = rotate_bundle_states(world, potentials.v0, U);
      potentials.gamma = rotate_bundle_states(world, potentials.gamma, U);
      omega = std::move(rotated_omega);
      if (world.rank() == 0 && config_.print_level > 1) {
        madness::print("EXCITED_ROTATE iter=", iter,
                       " omega_head=", format_vector_head(omega),
                       " nstates=", states.size());
      }
      return true;
    }

    [[nodiscard]] bool rotate_response_bundle(
        madness::World &world, std::vector<ResponseVector> &responses,
        std::vector<double> &omega, size_t iter) const {
      if (responses.empty()) {
        return false;
      }
      const size_t variant_index = responses.front().index();
      for (const auto &response : responses) {
        if (response.index() != variant_index) {
          if (world.rank() == 0 && config_.print_level > 0) {
            madness::print(
                "EXCITED_ROTATE_SKIP reason=variant_mismatch variant0=",
                variant_index);
          }
          return false;
        }
      }

      bool rotated = false;
      std::visit(
          [&](auto &typed_first) {
            using ResponseType = std::decay_t<decltype(typed_first)>;
            if constexpr (std::is_same_v<ResponseType, StaticUnrestrictedResponse> ||
                          std::is_same_v<ResponseType,
                                         DynamicUnrestrictedResponse>) {
              if (world.rank() == 0 && config_.print_level > 1) {
                madness::print(
                    "EXCITED_ROTATE_SKIP reason=unrestricted_wip iter=", iter);
              }
              rotated = false;
            } else {
              std::vector<ResponseType> typed_states;
              typed_states.reserve(responses.size());
              for (const auto &response : responses) {
                typed_states.push_back(std::get<ResponseType>(response));
              }
              auto potentials = compute_bundle_response_potentials(
                  world, typed_states, iter);
              rotated = rotate_excited_bundle_states(
                  world, typed_states, potentials, omega, iter);
              if (rotated) {
                for (size_t i = 0; i < typed_states.size(); ++i) {
                  responses[i] = std::move(typed_states[i]);
                }
              }
            }
          },
          responses.front());
      return rotated;
    }

    template <typename ResponseType>
    [[nodiscard]] std::vector<poperatorT>
    make_excited_bsh_operators(madness::World &world, const ResponseType &response,
                               double omega,
                               std::vector<double> &orbital_shifts) const {
      constexpr double lo = 1.0e-8;
      constexpr double shift_factor = 0.05;
      const double thresh = madness::FunctionDefaults<3>::get_thresh();
      const size_t n_orbitals = ground_num_orbitals_;

      std::vector<poperatorT> ops;
      orbital_shifts.clear();
      auto append_block = [&](double freq, bool shifted) {
        for (size_t p = 0; p < n_orbitals; ++p) {
          const double e = ground_energies_(long(p));
          double shift = 0.0;
          if (shifted && (e + freq) > 0.0) {
            shift = -(e + freq + shift_factor);
          }
          const double denom = e + freq + shift;
          const double stabilized = std::min(denom, -1.0e-8);
          const double mu = std::sqrt(std::max(1.0e-16, -2.0 * stabilized));
          ops.emplace_back(
              poperatorT(BSHOperatorPtr3D(world, mu, lo, thresh)));
          orbital_shifts.push_back(shift);
        }
      };

      if constexpr (std::is_same_v<ResponseType, StaticRestrictedResponse>) {
        append_block(omega, true);
      } else if constexpr (std::is_same_v<ResponseType,
                                          DynamicRestrictedResponse>) {
        append_block(omega, true);
        append_block(-omega, false);
      } else if constexpr (std::is_same_v<ResponseType,
                                          StaticUnrestrictedResponse>) {
        append_block(omega, true);
        append_block(omega, true);
      } else {
        append_block(omega, true);
        append_block(-omega, false);
        append_block(omega, true);
        append_block(-omega, false);
      }

      if (ops.size() != response.flat.size()) {
        // Keep shifts/indexing valid even if a caller passes an unexpected shape.
        orbital_shifts.resize(response.flat.size(), 0.0);
        if (ops.size() > response.flat.size()) {
          ops.resize(response.flat.size());
        }
      }
      return ops;
    }

    template <typename ResponseType>
    [[nodiscard]] ResponseType bsh_update_from_theta(
        madness::World &world, const ResponseType &state,
        const ResponseType &theta_state, double omega, double protocol_threshold,
        size_t state_id, size_t iter, double &state_residual_out) const {
      state_residual_out = 0.0;
      auto current = state;
      auto theta_response = theta_state;
      const bool stage_trace = (world.rank() == 0 && config_.print_level > 1);
      auto log_stage = [&](const char *stage) {
        if (!stage_trace) {
          return;
        }
        madness::print("EXCITED_STAGE iter=", iter, " state=", state_id,
                       " fn=bsh_update_from_theta stage=", stage,
                       " omega=", omega, " nflat=", current.flat.size());
      };

      log_stage("enter");
      align_response_to_ground(world, current);
      align_response_to_ground(world, theta_response);
      log_stage("after_align_response");
      if (current.flat.empty()) {
        log_stage("return_empty");
        return current;
      }

      std::vector<double> orbital_shifts;
      auto bsh_ops =
          make_excited_bsh_operators(world, current, omega, orbital_shifts);

      auto theta = theta_response.flat;
      const double update_trunc_threshold =
          std::max(1.0e-12, 1.0e-4 * protocol_threshold);
      for (size_t p = 0;
           p < theta.size() && p < current.flat.size() && p < orbital_shifts.size();
           ++p) {
        theta[p] = -2.0 * (theta[p] + orbital_shifts[p] * current.flat[p]);
      }
      truncate(world, theta, update_trunc_threshold, true);
      log_stage("after_theta_shift_truncate");

      auto updated_flat = apply(world, bsh_ops, theta);
      ensure_initialized_flat(world, updated_flat);
      log_stage("after_bsh_apply");
      auto updated = current;
      assign_flat_and_sync(updated, updated_flat);
      project_response_channels(world, updated);
      // truncate() modifies function data in-place through shared_ptrs that
      // project_response_channels already rebuilt via flatten().  The typed
      // channels (x_alpha, y_alpha) share those same ptrs, so no sync needed.
      truncate(world, updated.flat, update_trunc_threshold, true);
      normalize_response_metric(world, updated);
      align_response_to_ground(world, updated);
      log_stage("after_project_truncate_normalize");

      if (world.rank() == 0 && config_.print_level > 1) {
        const auto max_shift_pair = max_abs_value_with_index(orbital_shifts);
        const size_t max_shift_orb = max_shift_pair.first;
        const double max_shift_signed =
            (max_shift_orb < orbital_shifts.size()) ? orbital_shifts[max_shift_orb]
                                                    : 0.0;
        madness::print(
            "EXCITED_BSH_UPDATE iter=", iter, " state=", state_id,
            " omega=", omega,
            " shift_head=", format_vector_head(orbital_shifts),
            " max_abs_shift=", max_shift_pair.second,
            " max_shift_signed=", max_shift_signed,
            " max_shift_orb=", max_shift_orb,
            " bsh_ops_size=", bsh_ops.size(),
            " bsh_update_size=", updated.flat.size(),
            " theta_size=", theta.size(),
            " state_flat_size=", current.flat.size());
      }

      for (size_t p = 0; p < updated.flat.size() && p < current.flat.size(); ++p) {
        auto diff = updated.flat[p] - current.flat[p];
        state_residual_out = std::max(state_residual_out, diff.norm2());
      }
      log_stage("after_residual_eval");
      world.gop.fence();
      log_stage("exit");
      return updated;
    }

    template <typename ResponseType>
    [[nodiscard]] StateUpdateMetrics iterate_state_from_potentials(
        madness::World &world, ResponseType &state, double omega,
        const ResponseType &v0_state, const ResponseType &gamma_state,
        double protocol_threshold, double damping, size_t state_id,
        size_t iter, double max_rotation,
        response_solver *accelerator = nullptr) const {
      StateUpdateMetrics metrics;
      const bool stage_trace = (world.rank() == 0 && config_.print_level > 1);
      auto log_stage = [&](const char *stage) {
        if (!stage_trace) {
          return;
        }
        madness::print("EXCITED_STAGE iter=", iter, " state=", state_id,
                       " fn=iterate_state_from_potentials stage=", stage,
                       " omega=", omega, " nflat=", state.flat.size());
      };
      log_stage("enter");
      align_response_to_ground(world, state);
      log_stage("after_align_response");
      const auto density_before = compute_response_density(world, state);
      log_stage("after_density_before");
      const auto previous_flat = copy(world, state.flat, true);
      const double previous_state_norm =
          std::max(state_norm(world, previous_flat), 1.0e-12);
      log_stage("after_copy_previous");

      auto rotated_e0_flat = apply_hamiltonian_no_diag(world, state);
      auto theta = state;
      theta.flat = v0_state.flat - rotated_e0_flat + gamma_state.flat;
      ensure_initialized_flat(world, theta.flat);
      theta.sync();
      align_response_to_ground(world, theta);
      log_stage("after_theta_build");

      double update_residual = 0.0;
      auto updated =
          bsh_update_from_theta(world, state, theta, omega, protocol_threshold,
                                state_id, iter, update_residual);
      log_stage("after_bsh_update");

      auto candidate_flat = copy(world, updated.flat, true);
      auto residual_flat = sub(world, previous_flat, updated.flat);
      ensure_initialized_flat(world, residual_flat);
      if (accelerator != nullptr && iter > 1) {
        candidate_flat = accelerator->update(previous_flat, residual_flat);
        metrics.accelerator_applied = true;
      } else if (damping < 0.999999) {
        candidate_flat =
            gaxpy_oop(damping, updated.flat, 1.0 - damping, previous_flat, true);
      }
      ensure_initialized_flat(world, candidate_flat);
      metrics.step_norm =
          state_norm(world, sub(world, previous_flat, candidate_flat));
      if (max_rotation > 0.0 && metrics.step_norm > max_rotation) {
        ResponseSolverUtils::do_step_restriction(world, previous_flat,
                                                 candidate_flat,
                                                 metrics.step_norm, "a",
                                                 max_rotation);
        metrics.step_restricted = true;
        metrics.step_norm =
            state_norm(world, sub(world, previous_flat, candidate_flat));
      }
      assign_flat_and_sync(state, std::move(candidate_flat));
      log_stage("after_damped_update");
      project_response_channels(world, state);
      normalize_response_metric(world, state);
      align_response_to_ground(world, state);
      log_stage("after_project_normalize");

      const auto density_after = compute_response_density(world, state);
      log_stage("after_density_after");
      metrics.density_change = (density_after - density_before).norm2();
      metrics.residual = update_residual;
      metrics.relative_residual = metrics.residual / previous_state_norm;
      log_stage("exit");
      return metrics;
    }

    template <typename ResponseType>
    [[nodiscard]] ResponseType bsh_update_state(madness::World &world,
                                                const ResponseType &state,
                                                double omega,
                                                double protocol_threshold,
                                                size_t state_id,
                                                size_t iter,
                                                double &state_residual_out) const {
      state_residual_out = 0.0;
      auto current = state;
      const bool stage_trace = (world.rank() == 0 && config_.print_level > 1);
      auto log_stage = [&](const char *stage) {
        if (!stage_trace) {
          return;
        }
        madness::print("EXCITED_STAGE iter=", iter, " state=", state_id,
                       " fn=bsh_update_state stage=", stage, " omega=", omega,
                       " nflat=", current.flat.size());
      };
      log_stage("enter");
      align_response_to_ground(world, current);
      log_stage("after_align_response");
      if (current.flat.empty()) {
        log_stage("return_empty");
        return current;
      }

      std::vector<int> orbital_index(current.flat.size());
      std::iota(orbital_index.begin(), orbital_index.end(), 0);
      std::vector<int> state_index(current.flat.size(), 0);

      ResponseComputeGroundExchange ground_exchange_task_impl;
      MacroTask ground_exchange_task(world, ground_exchange_task_impl);
      ResponseComputeGammaX response_exchange_task_impl;
      MacroTask response_exchange_task(world, response_exchange_task_impl);

      const bool use_tda_exchange = !response_has_y_channel_v<ResponseType>;
      auto k0 = ground_exchange_task(orbital_index, state_index, current.flat,
                                     ground_orbitals_, use_tda_exchange);
      auto gx = response_exchange_task(orbital_index, state_index, current.flat,
                                       ground_orbitals_, use_tda_exchange);
      log_stage("after_exchange_tasks");

      const double c_xc = ground_data_->xcf_.hf_exchange_coefficient();
      auto v_local = ground_data_->V_local * current.flat;
      auto v0x = v_local - c_xc * k0;
      auto epsilon = apply_hamiltonian_no_diag(world, current);
      auto theta = v0x - epsilon + gx;
      log_stage("after_theta_build");

      std::vector<double> orbital_shifts;
      auto bsh_ops =
          make_excited_bsh_operators(world, current, omega, orbital_shifts);

      const double update_trunc_threshold =
          std::max(1.0e-12, 1.0e-4 * protocol_threshold);
      for (size_t p = 0;
           p < theta.size() && p < current.flat.size() && p < orbital_shifts.size();
           ++p) {
        theta[p] = -2.0 * (theta[p] + orbital_shifts[p] * current.flat[p]);
      }
      truncate(world, theta, update_trunc_threshold, true);
      log_stage("after_theta_shift_truncate");

      auto updated_flat = apply(world, bsh_ops, theta);
      ensure_initialized_flat(world, updated_flat);
      log_stage("after_bsh_apply");
      auto updated = current;
      assign_flat_and_sync(updated, updated_flat);
      project_response_channels(world, updated);
      // truncate() modifies function data in-place through shared_ptrs that
      // project_response_channels already rebuilt via flatten().  The typed
      // channels (x_alpha, y_alpha) share those same ptrs, so no sync needed.
      truncate(world, updated.flat, update_trunc_threshold, true);
      normalize_response_metric(world, updated);
      // Projection/normalization can leave null Function impls in near-zero
      // channels. Re-align to ensure downstream density kernels always see
      // initialized functions with the expected channel shape.
      align_response_to_ground(world, updated);
      log_stage("after_project_truncate_normalize");

      const auto max_shift_pair = max_abs_value_with_index(orbital_shifts);
      const size_t max_shift_orb = max_shift_pair.first;
      const double shift = (max_shift_orb < orbital_shifts.size())
                               ? orbital_shifts[max_shift_orb]
                               : 0.0;
      const double max_abs_shift = max_shift_pair.second;
      const double state_norm_in = state_norm(world, current.flat);
      const double theta_norm = state_norm(world, theta);
      const double updated_norm_post = state_norm(world, updated.flat);
      if (world.rank() == 0 && config_.print_level > 1) {
        madness::print(
            "EXCITED_BSH_UPDATE iter=", iter, " state=", state_id,
            " omega=", omega,
            " shift_head=", format_vector_head(orbital_shifts),
            " max_abs_shift=", max_abs_shift,
            " max_shift_signed=", shift, " max_shift_orb=", max_shift_orb,
            " bsh_ops_size=", bsh_ops.size(),
            " bsh_update_size=", updated.flat.size(),
            " theta_size=", theta.size(), " state_flat_size=", current.flat.size());
      }
      if (world.rank() == 0 && config_.print_level > 1 &&
          updated_norm_post < 1.0e-12) {
        static std::atomic<int> zero_update_reports{0};
        const int report_idx =
            zero_update_reports.fetch_add(1, std::memory_order_relaxed);
        if (report_idx < 48) {
          madness::print(
              "EXCITED_BSH_ZERO_UPDATE report=", report_idx,
              " omega=", omega, " shift=", shift, " state_norm=",
              state_norm_in, " theta_norm=", theta_norm,
              " updated_norm_post=", updated_norm_post,
              " update_trunc_threshold=", update_trunc_threshold);
        }
      }

      for (size_t p = 0; p < updated.flat.size() && p < current.flat.size(); ++p) {
        auto diff = updated.flat[p] - current.flat[p];
        state_residual_out = std::max(state_residual_out, diff.norm2());
      }
      log_stage("after_residual_eval");
      // Release large temporaries, but do not clear `updated_flat` since
      // `updated.flat` may still share underlying Function implementations.
      k0.clear();
      gx.clear();
      v_local.clear();
      v0x.clear();
      epsilon.clear();
      theta.clear();
      world.gop.fence();
      log_stage("exit");
      return updated;
    }

    template <typename ResponseType>
    [[nodiscard]] StateUpdateMetrics iterate_state(
        madness::World &world, ResponseType &state, double omega,
        double protocol_threshold, double damping, size_t state_id,
        size_t iter) const {
      StateUpdateMetrics metrics;
      const bool stage_trace = (world.rank() == 0 && config_.print_level > 1);
      auto log_stage = [&](const char *stage) {
        if (!stage_trace) {
          return;
        }
        madness::print("EXCITED_STAGE iter=", iter, " state=", state_id,
                       " fn=iterate_state stage=", stage, " omega=", omega,
                       " nflat=", state.flat.size());
      };
      log_stage("enter");
      align_response_to_ground(world, state);
      log_stage("after_align_response");
      const auto density_before = compute_response_density(world, state);
      log_stage("after_density_before");
      const auto previous_flat = copy(world, state.flat, true);
      const double previous_state_norm =
          std::max(state_norm(world, previous_flat), 1.0e-12);
      log_stage("after_copy_previous");

      double update_residual = 0.0;
      auto updated = bsh_update_state(world, state, omega, protocol_threshold,
                                      state_id, iter, update_residual);
      log_stage("after_bsh_update");
      state.flat = gaxpy_oop(damping, updated.flat, 1.0 - damping, state.flat, true);
      ensure_initialized_flat(world, state.flat);
      log_stage("after_damped_update");
      state.sync();
      project_response_channels(world, state);
      normalize_response_metric(world, state);
      const size_t uninit_post_normalize = count_uninitialized_flat(state);
      if (config_.print_level > 1 && uninit_post_normalize > 0) {
        madness::print("EXCITED_WARN rank=", world.rank(), " iter=", iter,
                       " state=", state_id,
                       " stage=post_normalize uninitialized_flat=",
                       uninit_post_normalize, " nflat=", state.flat.size());
      }
      // Keep the state in a safe, canonical shape before density evaluation.
      align_response_to_ground(world, state);
      const size_t uninit_post_align = count_uninitialized_flat(state);
      if (config_.print_level > 1 && uninit_post_align > 0) {
        madness::print("EXCITED_WARN rank=", world.rank(), " iter=", iter,
                       " state=", state_id, " stage=post_align uninitialized_flat=",
                       uninit_post_align, " nflat=", state.flat.size());
      }
      log_stage("after_project_normalize");

      const auto density_after = compute_response_density(world, state);
      log_stage("after_density_after");
      metrics.density_change = (density_after - density_before).norm2();
      metrics.residual = update_residual;
      metrics.relative_residual = metrics.residual / previous_state_norm;
      metrics.step_norm =
          state_norm(world, sub(world, previous_flat, state.flat));
      log_stage("exit");
      return metrics;
    }

    void sort_by_energy(std::vector<vector_real_function_3d> &states,
                        std::vector<double> &energies) const {
      std::vector<size_t> order(energies.size());
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(),
                [&](size_t a, size_t b) { return energies[a] < energies[b]; });

      std::vector<vector_real_function_3d> sorted_states;
      std::vector<double> sorted_energies;
      sorted_states.reserve(order.size());
      sorted_energies.reserve(order.size());
      for (const auto idx : order) {
        sorted_states.push_back(states[idx]);
        sorted_energies.push_back(energies[idx]);
      }
      states.swap(sorted_states);
      energies.swap(sorted_energies);
    }

    [[nodiscard]] ExcitedTrialSpace
    make_derivative_trial(madness::World &world, size_t m, size_t direction_offset) {
      ensure_ground_data(world);
      ExcitedTrialSpace trial;
      trial.num_states = m;
      trial.num_orbitals = ground_num_orbitals_;
      trial.x_states.reserve(m);
      trial.omega.assign(m, 0.0);

      real_derivative_3d Dx(world, 0);
      real_derivative_3d Dy(world, 1);
      real_derivative_3d Dz(world, 2);
      const auto dphi_x = apply(world, Dx, ground_orbitals_);
      const auto dphi_y = apply(world, Dy, ground_orbitals_);
      const auto dphi_z = apply(world, Dz, ground_orbitals_);

      for (size_t i = 0; i < m; ++i) {
        const size_t orb = i % std::max<size_t>(ground_num_orbitals_, 1);
        const size_t dir = (i + direction_offset) % 3;
        auto state = zero_functions_compressed<double, 3>(
            world, static_cast<int>(ground_num_orbitals_));
        if (ground_num_orbitals_ > 0) {
          if (dir == 0) {
            state[orb] = copy(dphi_x[orb]);
          } else if (dir == 1) {
            state[orb] = copy(dphi_y[orb]);
          } else {
            state[orb] = copy(dphi_z[orb]);
          }
        }
        trial.x_states.push_back(std::move(state));
      }
      project_and_orthonormalize(world, trial.x_states);
      trial.omega = estimate_state_energies(world, trial.x_states);
      sort_by_energy(trial.x_states, trial.omega);
      trial.num_states = trial.x_states.size();
      return trial;
    }

    [[nodiscard]] ExcitedTrialSpace
    make_projected_ao_trial(madness::World &world, size_t m) {
      ensure_ground_data(world);
      ExcitedTrialSpace trial;
      trial.num_states = 0;
      trial.num_orbitals = ground_num_orbitals_;
      if (ground_num_orbitals_ == 0) {
        return trial;
      }

      const Molecule &molecule = ground_data_->getMolecule();
      if (molecule.natom() == 0) {
        return trial;
      }

      const std::vector<std::string> candidate_bases = {
          "d-aug-cc-pVQZ", "aug-cc-pVQZ", "aug-cc-pVTZ", "aug-cc-pVDZ",
          "6-31g"};
      const double thresh = current_protocol_thresh();
      const size_t generation_target = std::max<size_t>(m, 2 * m);

      for (const auto &basis_name : candidate_bases) {
        try {
          AtomicBasisSet ao_basis(basis_name);
          const int nbf = ao_basis.nbf(molecule);
          if (nbf <= 0) {
            continue;
          }
          if (world.rank() == 0 && config_.print_level > 0) {
            madness::print("EXCITED_TRIAL_AO basis=", basis_name,
                           " nbf=", nbf);
          }

          vector_real_function_3d ao_functions(static_cast<size_t>(nbf));
          for (int i = 0; i < nbf; ++i) {
            auto aofunc = real_functor_3d(new madchem::AtomicBasisFunctor(
                ao_basis.get_atomic_basis_function(molecule, i)));
            ao_functions[static_cast<size_t>(i)] =
                real_factory_3d(world)
                    .functor(aofunc)
                    .truncate_on_project()
                    .nofence()
                    .truncate_mode(1);
          }
          world.gop.fence();
          truncate(world, ao_functions, thresh, true);

          QProjector<double, 3> projector(ground_orbitals_);
          for (auto &ao_fn : ao_functions) {
            ao_fn = projector(ao_fn);
          }
          truncate(world, ao_functions, thresh, true);

          std::vector<vector_real_function_3d> states;
          states.reserve(generation_target);
          for (size_t i = 0;
               i < ao_functions.size() && states.size() < generation_target;
               ++i) {
            for (size_t orb = 0;
                 orb < ground_num_orbitals_ && states.size() < generation_target;
                 ++orb) {
              auto state = zero_functions_compressed<double, 3>(
                  world, static_cast<int>(ground_num_orbitals_));
              state[orb] = copy(ao_functions[i]);
              states.push_back(std::move(state));
            }
          }
          project_and_orthonormalize(world, states);
          if (!states.empty()) {
            trial.x_states = std::move(states);
            trial.omega = estimate_state_energies(world, trial.x_states);
            sort_by_energy(trial.x_states, trial.omega);
            trial.num_states = trial.x_states.size();
            return trial;
          }
        } catch (...) {
          continue;
        }
      }
      return trial;
    }

    [[nodiscard]] ExcitedTrialSpace
    make_localized_gaussian_trial(madness::World &world, size_t m,
                                  unsigned seed) {
      ensure_ground_data(world);
      ExcitedTrialSpace trial;
      trial.num_states = 0;
      trial.num_orbitals = ground_num_orbitals_;
      if (ground_num_orbitals_ == 0) {
        return trial;
      }

      const auto atoms = ground_data_->getMolecule().get_atoms();
      if (atoms.empty()) {
        return trial;
      }

      std::mt19937 rng(seed);
      std::uniform_real_distribution<double> exponent_dist(0.04, 1.20);
      std::uniform_real_distribution<double> shift_dist(-0.35, 0.35);
      std::uniform_int_distribution<int> axis_dist(0, 2);
      std::uniform_int_distribution<int> power_dist(0, 2);

      const size_t generation_target = std::max<size_t>(m, 2 * m);
      trial.x_states.reserve(generation_target);
      for (size_t i = 0; i < generation_target; ++i) {
        const auto &atom = atoms[i % atoms.size()];
        coord_3d center = atom.get_coords();
        center[0] += shift_dist(rng);
        center[1] += shift_dist(rng);
        center[2] += shift_dist(rng);

        std::vector<int> powers(3, 0);
        powers[axis_dist(rng)] = power_dist(rng);
        const double exponent = exponent_dist(rng);

        auto gauss = real_factory_3d(world).functor(real_functor_3d(
            new LocalizedGaussianGuess<3>(center, exponent, powers)));

        auto state = zero_functions_compressed<double, 3>(
            world, static_cast<int>(ground_num_orbitals_));
        const size_t orb = i % ground_num_orbitals_;
        state[orb] = std::move(gauss);
        trial.x_states.push_back(std::move(state));
      }

      project_and_orthonormalize(world, trial.x_states);
      if (!trial.x_states.empty()) {
        trial.omega = estimate_state_energies(world, trial.x_states);
        sort_by_energy(trial.x_states, trial.omega);
      }
      trial.num_states = trial.x_states.size();
      return trial;
    }

    [[nodiscard]] ExcitedTrialSpace
    make_random_trial(madness::World &world,
                      const ExcitedBundleProtocolInput &input) {
      const size_t m = std::max<size_t>(2 * input.num_states, 1);
      auto trial = make_projected_ao_trial(world, m);
      if (trial.x_states.size() < input.num_states) {
        auto derivative_trial = make_derivative_trial(world, m, 0);
        for (auto &state : derivative_trial.x_states) {
          trial.x_states.push_back(std::move(state));
        }
        project_and_orthonormalize(world, trial.x_states);
        trial.omega = estimate_state_energies(world, trial.x_states);
        sort_by_energy(trial.x_states, trial.omega);
        trial.num_states = trial.x_states.size();
      }
      if (trial.x_states.size() < input.num_states) {
        auto gaussian_trial = make_localized_gaussian_trial(
            world, m, 37u + static_cast<unsigned>(input.protocol_index));
        for (auto &state : gaussian_trial.x_states) {
          trial.x_states.push_back(std::move(state));
        }
        project_and_orthonormalize(world, trial.x_states);
        trial.omega = estimate_state_energies(world, trial.x_states);
        sort_by_energy(trial.x_states, trial.omega);
        trial.num_states = trial.x_states.size();
      }
      std::mt19937 rng(17u + static_cast<unsigned>(input.protocol_index));
      std::uniform_real_distribution<double> dist(0.95, 1.05);
      for (auto &w : trial.omega) {
        w *= dist(rng);
      }
      return trial;
    }

    [[nodiscard]] ExcitedTrialSpace
    make_nwchem_trial(madness::World &world,
                      const ExcitedBundleProtocolInput &input) {
      return make_derivative_trial(world, std::max<size_t>(2 * input.num_states, 1), 1);
    }

    [[nodiscard]] ExcitedTrialSpace
    create_trial_functions(madness::World &world,
                           const ExcitedBundleProtocolInput &input) {
      return make_derivative_trial(world, std::max<size_t>(2 * input.num_states, 1), 2);
    }

    [[nodiscard]] ExcitedTrialSpace
    create_trial_functions2(madness::World &world,
                            const ExcitedBundleProtocolInput &input) {
      return make_derivative_trial(world, std::max<size_t>(2 * input.num_states, 1), 0);
    }

    void iterate_trial(madness::World &world,
                       const ExcitedBundleProtocolInput &input,
                       ExcitedTrialSpace &trial) {
      if (trial.x_states.empty()) {
        trial.omega.clear();
        trial.num_states = 0;
        trial.num_orbitals = ground_num_orbitals_;
        return;
      }
      const size_t iters = std::max<size_t>(input.guess_max_iter, 1);
      for (size_t iter = 0; iter < iters; ++iter) {
        project_and_orthonormalize(world, trial.x_states);
        auto estimate = estimate_state_energies(world, trial.x_states);
        double max_delta = 0.0;
        if (trial.omega.size() != estimate.size()) {
          trial.omega = estimate;
          const auto [_, max_abs_estimate] = max_abs_value_with_index(trial.omega);
          max_delta = max_abs_estimate;
        } else {
          for (size_t i = 0; i < trial.omega.size(); ++i) {
            const double delta = estimate[i] - trial.omega[i];
            max_delta = std::max(max_delta, std::fabs(delta));
            trial.omega[i] = 0.6 * trial.omega[i] + 0.4 * estimate[i];
          }
        }
        sort_by_energy(trial.x_states, trial.omega);
        if (world.rank() == 0 && config_.print_level > 1) {
          madness::print("EXCITED_GUESS_ITER iter=", (iter + 1), "/", iters,
                         " max_delta=", max_delta,
                         " omega_head=", format_vector_head(trial.omega));
        }
      }
      trial.num_states = trial.x_states.size();
      trial.num_orbitals = ground_num_orbitals_;
    }

    void sort_trial_space(madness::World &world, ExcitedTrialSpace &trial) const {
      if (trial.x_states.empty()) {
        trial.omega.clear();
        trial.num_states = 0;
        trial.num_orbitals = ground_num_orbitals_;
        return;
      }
      if (trial.omega.size() != trial.x_states.size()) {
        trial.omega = estimate_state_energies(world, trial.x_states);
      }
      sort_by_energy(trial.x_states, trial.omega);
      trial.num_states = trial.x_states.size();
      trial.num_orbitals = ground_num_orbitals_;
    }

    void select_lowest_trial_roots(madness::World &world, ExcitedTrialSpace &trial,
                                   size_t requested_states) const {
      sort_trial_space(world, trial);
      const size_t keep = std::min(requested_states, trial.x_states.size());
      if (trial.x_states.size() > keep) {
        trial.x_states.resize(keep);
      }
      if (trial.omega.size() > keep) {
        trial.omega.resize(keep);
      } else if (trial.omega.size() < trial.x_states.size()) {
        trial.omega = estimate_state_energies(world, trial.x_states);
      }
      trial.num_states = trial.x_states.size();
      trial.num_orbitals = ground_num_orbitals_;
    }

    void seed_active_bundle_from_trial_space(
        madness::World &world, const ExcitedBundleProtocolInput &input) {
      active_response_bundle_ = build_response_bundle_from_trial_space(world, input);
      active_bundle_k_ = current_protocol_k();
      active_bundle_thresh_ = current_protocol_thresh();
    }

    [[nodiscard]] ResponseVector
    make_state_response_vector(madness::World &world,
                               const vector_real_function_3d &x_state,
                               bool tda) const {
      auto x_aligned = copy(world, x_state, true);
      align_component_to_ground(world, x_aligned, ground_num_orbitals_);
      auto zeros = zero_functions_compressed<double, 3>(
          world, static_cast<int>(ground_num_orbitals_));

      const bool unrestricted = !ground_data_->isSpinRestricted();
      if (!unrestricted && tda) {
        StaticRestrictedResponse response;
        response.x_alpha = copy(world, x_aligned, true);
        response.flatten();
        return response;
      }
      if (!unrestricted && !tda) {
        DynamicRestrictedResponse response;
        response.x_alpha = copy(world, x_aligned, true);
        response.y_alpha = copy(world, zeros, true);
        response.flatten();
        return response;
      }
      if (unrestricted && tda) {
        StaticUnrestrictedResponse response;
        response.x_alpha = copy(world, x_aligned, true);
        response.x_beta = copy(world, zeros, true);
        response.flatten();
        return response;
      }

      DynamicUnrestrictedResponse response;
      response.x_alpha = copy(world, x_aligned, true);
      response.y_alpha = copy(world, zeros, true);
      response.x_beta = copy(world, zeros, true);
      response.y_beta = copy(world, zeros, true);
      response.flatten();
      return response;
    }

    [[nodiscard]] vector_real_function_3d
    extract_x_state(madness::World &world, const ResponseVector &response) const {
      return std::visit(
          [&](const auto &typed_response) {
            return copy(world, typed_response.x_alpha, true);
          },
          response);
    }

    [[nodiscard]] std::vector<ResponseVector>
    build_response_bundle_from_trial_space(madness::World &world,
                                           const ExcitedBundleProtocolInput &input) const {
      std::vector<ResponseVector> responses;
      responses.reserve(trial_space_.x_states.size());
      for (const auto &x_state : trial_space_.x_states) {
        responses.push_back(make_state_response_vector(world, x_state, input.tda));
      }
      return responses;
    }

    [[nodiscard]] std::vector<ResponseVector>
    build_response_bundle_seed(madness::World &world,
                               const ExcitedBundleProtocolInput &input) const {
      if (active_bundle_matches_input(input)) {
        return active_response_bundle_;
      }
      return build_response_bundle_from_trial_space(world, input);
    }

    void sync_trial_space_from_response_bundle(
        madness::World &world, const std::vector<ResponseVector> &responses) {
      trial_space_.x_states.clear();
      trial_space_.x_states.reserve(responses.size());
      for (const auto &response : responses) {
        trial_space_.x_states.push_back(extract_x_state(world, response));
      }
      trial_space_.num_states = trial_space_.x_states.size();
      trial_space_.num_orbitals = ground_num_orbitals_;
    }

    [[nodiscard]] StateUpdateMetrics
    iterate_response_state(madness::World &world, ResponseVector &response,
                           double omega, double protocol_threshold,
                           double damping, size_t state_id,
                           size_t iter) const {
      return std::visit(
          [&](auto &typed_response) {
            return iterate_state(world, typed_response, omega, protocol_threshold,
                                 damping, state_id, iter);
          },
          response);
    }

    template <typename ResponseType>
    [[nodiscard]] bool iterate_typed_bundle_legacy_sequence(
        madness::World &world, std::vector<ResponseType> &states,
        std::vector<double> &omega, double protocol_threshold, double damping,
        size_t iter, std::vector<double> &state_residuals,
        std::vector<double> &state_density_changes,
        std::vector<double> &state_relative_residuals, double max_rotation,
        std::vector<response_solver> *state_accelerators = nullptr,
        std::vector<size_t> *slot_order_out = nullptr) const {
      if (states.empty() || omega.empty()) {
        return false;
      }
      auto potentials = compute_bundle_response_potentials(world, states, iter);
      bool bundle_rotated = rotate_excited_bundle_states(
          world, states, potentials, omega, iter, slot_order_out);

      size_t n_states = std::min(states.size(), omega.size());
      n_states = std::min(n_states, potentials.v0.size());
      n_states = std::min(n_states, potentials.gamma.size());
      if (states.size() != n_states) {
        states.resize(n_states);
      }
      if (omega.size() != n_states) {
        omega.resize(n_states);
      }
      if (potentials.v0.size() != n_states) {
        potentials.v0.resize(n_states);
      }
      if (potentials.gamma.size() != n_states) {
        potentials.gamma.resize(n_states);
      }
      if (state_residuals.size() != n_states) {
        state_residuals.assign(n_states, 0.0);
      }
      if (state_density_changes.size() != n_states) {
        state_density_changes.assign(n_states, 0.0);
      }
      if (state_relative_residuals.size() != n_states) {
        state_relative_residuals.assign(n_states, 0.0);
      }

      for (size_t i = 0; i < n_states; ++i) {
        if (world.rank() == 0 && config_.print_level > 1) {
          madness::print("EXCITED_STATE_BEGIN iter=", iter, " state=", i,
                         " omega=", omega[i]);
        }
        response_solver *state_accelerator = nullptr;
        if (state_accelerators != nullptr && i < state_accelerators->size()) {
          state_accelerator = &(*state_accelerators)[i];
        }
        const auto metrics = iterate_state_from_potentials(
            world, states[i], omega[i], potentials.v0[i], potentials.gamma[i],
            protocol_threshold, damping, i, iter, max_rotation,
            state_accelerator);
        state_residuals[i] = metrics.residual;
        state_density_changes[i] = metrics.density_change;
        state_relative_residuals[i] = metrics.relative_residual;
        if (world.rank() == 0 && config_.print_level > 1) {
          madness::print("EXCITED_STATE_DONE iter=", iter, " state=", i,
                         " residual=", metrics.residual,
                         " relative_residual=", metrics.relative_residual,
                         " drho=", metrics.density_change,
                         " step_norm=", metrics.step_norm,
                         " accelerator=", metrics.accelerator_applied,
                         " step_restricted=", metrics.step_restricted);
        }
        // MADNESS function destruction is deferred to a global fence.
        // In excited bundles we update many states per iteration, so fence
        // each state to avoid large transient retention.
        world.gop.fence();
      }
      return bundle_rotated;
    }

    void build_fresh_guess(madness::World &world,
                           const ExcitedBundleProtocolInput &input) {
      const size_t requested_states = std::max<size_t>(input.num_states, 1);
      if (input.tda && input.num_states <= 6) {
        trial_space_ = create_trial_functions2(world, input);
      } else if (!input.tda && input.num_states <= 2) {
        trial_space_ = make_nwchem_trial(world, input);
      } else if (!input.tda && input.num_states <= 6) {
        trial_space_ = create_trial_functions(world, input);
      } else {
        trial_space_ = make_random_trial(world, input);
      }

      if (trial_space_.x_states.size() < requested_states) {
        auto expanded = make_random_trial(world, input);
        if (expanded.x_states.size() > trial_space_.x_states.size()) {
          trial_space_ = std::move(expanded);
        }
      }

      iterate_trial(world, input, trial_space_);
      select_lowest_trial_roots(world, trial_space_, requested_states);

      if (trial_space_.x_states.size() < requested_states) {
        auto expanded = make_random_trial(world, input);
        if (expanded.x_states.size() > trial_space_.x_states.size()) {
          trial_space_ = std::move(expanded);
          iterate_trial(world, input, trial_space_);
          select_lowest_trial_roots(world, trial_space_, requested_states);
        }
      }
      omega_ = trial_space_.omega;
      state_names_.clear();
      root_descriptors_.clear();
      slot_permutation_.clear();
      next_root_stable_index_ = 0;
      active_response_bundle_.clear();
      active_bundle_k_ = -1;
      active_bundle_thresh_ = std::numeric_limits<double>::infinity();
      if (!trial_space_.x_states.empty()) {
        seed_active_bundle_from_trial_space(world, input);
      }
      residual_norms_ = std::vector<double>(omega_.size(), 1.0e-2);
      density_change_norms_.assign(omega_.size(), 0.0);
      relative_residual_norms_.assign(omega_.size(), 0.0);
      last_iteration_count_ = std::max<size_t>(input.guess_max_iter, 1);
      iteration_max_residuals_.clear();
      iteration_max_density_changes_.clear();
      iteration_max_relative_residuals_.clear();
      has_active_guess_ = !omega_.empty();
      trial_space_k_ = current_protocol_k();
      trial_space_thresh_ = current_protocol_thresh();
      if (world.rank() == 0 && config_.print_level > 0) {
        madness::print("EXCITED_FRESH_GUESS_SELECT requested_states=",
                       requested_states, " selected_states=", omega_.size(),
                       " snapshot_kind=",
                       (active_response_bundle_.empty() ? "guess_trial_seed"
                                                        : "guess_bundle"),
                       " omega_head=", format_vector_head(omega_));
      }
    }

    void align_trial_space_protocol(madness::World &world) {
      const int target_k = current_protocol_k();
      const double target_thresh = current_protocol_thresh();
      if (trial_space_.x_states.empty()) {
        trial_space_k_ = target_k;
        trial_space_thresh_ = target_thresh;
        return;
      }

      const int source_k = ResponseSolverUtils::infer_state_bundle_k(trial_space_.x_states);
      const bool projected = ResponseSolverUtils::align_state_bundle_protocol(
          world, trial_space_.x_states, target_k, target_thresh);
      if (projected) {
        if (world.rank() == 0) {
          madness::print("EXCITED_TRIAL_REPROJECT source_k=", source_k,
                         " target_k=", target_k, " threshold=", target_thresh,
                         " states=", trial_space_.x_states.size());
        }
        project_and_orthonormalize(world, trial_space_.x_states);
      }
      trial_space_k_ = target_k;
      trial_space_thresh_ = target_thresh;
    }

    bool align_active_bundle_protocol(madness::World &world,
                                      const ExcitedBundleProtocolInput &input) {
      const int target_k = current_protocol_k();
      const double target_thresh = current_protocol_thresh();
      if (active_response_bundle_.empty()) {
        active_bundle_k_ = target_k;
        active_bundle_thresh_ = target_thresh;
        return false;
      }

      const bool variant_matches = active_bundle_matches_input(input);
      if (!variant_matches) {
        if (world.rank() == 0 && config_.print_level > 0) {
          madness::print("EXCITED_BUNDLE_RESTART_REJECT reason=variant_mismatch",
                         " have=", active_bundle_variant(),
                         " want=", expected_response_variant(input.tda));
        }
        active_response_bundle_.clear();
        active_bundle_k_ = -1;
        active_bundle_thresh_ = std::numeric_limits<double>::infinity();
        return false;
      }

      const int source_k = active_bundle_k_;
      const double source_thresh = active_bundle_thresh_;
      const bool projected = ResponseSolverUtils::align_response_bundle_protocol(
          world, active_response_bundle_, target_k, target_thresh);
      for (auto &response : active_response_bundle_) {
        std::visit(
            [&](auto &typed_response) {
              align_response_to_ground(world, typed_response);
              project_response_channels(world, typed_response);
              normalize_response_metric(world, typed_response);
              align_response_to_ground(world, typed_response);
            },
            response);
      }
      if ((projected ||
           std::abs(source_thresh - target_thresh) > 1.0e-14) &&
          world.rank() == 0) {
        madness::print("EXCITED_BUNDLE_REPROJECT source_k=", source_k,
                       " target_k=", target_k,
                       " source_threshold=", source_thresh,
                       " target_threshold=", target_thresh,
                       " states=", active_response_bundle_.size(),
                       " variant=", active_bundle_variant());
      }
      active_bundle_k_ = target_k;
      active_bundle_thresh_ = target_thresh;
      sync_trial_space_from_response_bundle(world, active_response_bundle_);
      trial_space_.num_states = trial_space_.x_states.size();
      trial_space_.num_orbitals = ground_num_orbitals_;
      trial_space_.omega = omega_;
      trial_space_k_ = target_k;
      trial_space_thresh_ = target_thresh;
      return true;
    }

    void ensure_trial_space_matches_guess(
        madness::World &world, const ExcitedBundleProtocolInput &input) {
      ensure_ground_data(world);
      if (trial_space_.x_states.empty() && !active_response_bundle_.empty()) {
        sync_trial_space_from_response_bundle(world, active_response_bundle_);
        trial_space_.num_states = trial_space_.x_states.size();
        trial_space_.num_orbitals = ground_num_orbitals_;
        trial_space_.omega = omega_;
      }
      if (trial_space_.x_states.empty()) {
        trial_space_ = make_random_trial(world, input);
        iterate_trial(world, input, trial_space_);
        trial_space_k_ = current_protocol_k();
        trial_space_thresh_ = current_protocol_thresh();
      }

      if (trial_space_.x_states.size() < input.num_states) {
        auto expanded_trial = make_random_trial(world, input);
        if (expanded_trial.x_states.size() > trial_space_.x_states.size()) {
          trial_space_ = std::move(expanded_trial);
          iterate_trial(world, input, trial_space_);
          trial_space_k_ = current_protocol_k();
          trial_space_thresh_ = current_protocol_thresh();
        }
      }

      if (trial_space_.x_states.size() > input.num_states) {
        trial_space_.x_states.resize(input.num_states);
      }
      if (!omega_.empty()) {
        if (trial_space_.x_states.size() > omega_.size()) {
          trial_space_.x_states.resize(omega_.size());
        } else if (trial_space_.x_states.size() < omega_.size()) {
          omega_.resize(trial_space_.x_states.size());
        }
        trial_space_.omega = omega_;
      } else {
        omega_ = estimate_state_energies(world, trial_space_.x_states);
        trial_space_.omega = omega_;
      }
      if (residual_norms_.size() != omega_.size()) {
        residual_norms_.assign(omega_.size(), 1.0e-2);
      }
      has_active_guess_ = !omega_.empty();
      align_trial_space_protocol(world);
    }

    void ensure_state_names_for_protocol(double threshold) {
      sync_root_descriptors_for_protocol(threshold);
    }

    [[nodiscard]] RestartSeedChoice
    select_restart_seed(madness::World &world,
                        const ExcitedBundleProtocolInput &input,
                        const RestartSnapshot &current_protocol_snapshot) const {
      RestartSeedChoice selected;
      auto is_stalled_snapshot = [](const RestartSnapshot &snapshot) -> bool {
        if (snapshot.converged || snapshot.iteration_max_residuals.empty()) {
          return false;
        }
        const auto [min_it, max_it] = std::minmax_element(
            snapshot.iteration_max_residuals.begin(),
            snapshot.iteration_max_residuals.end());
        const double best = *min_it;
        const double worst = *max_it;
        const double final_residual = snapshot.iteration_max_residuals.back();
        if (!std::isfinite(best) || !std::isfinite(worst) ||
            !std::isfinite(final_residual)) {
          return true;
        }
        return (best >= 0.95 && worst >= 0.95 && final_residual >= 0.95);
      };
      auto print_skip_stalled = [&](double source_threshold,
                                    const RestartSnapshot &snapshot) {
        if (world.rank() != 0 || config_.print_level <= 0 ||
            snapshot.iteration_max_residuals.empty()) {
          return;
        }
        const auto [min_it, max_it] = std::minmax_element(
            snapshot.iteration_max_residuals.begin(),
            snapshot.iteration_max_residuals.end());
        madness::print("EXCITED_RESTART_SKIP threshold=", source_threshold,
                       " reason=stalled_snapshot min_iter_residual=", *min_it,
                       " max_iter_residual=", *max_it,
                       " final_iter_residual=",
                       snapshot.iteration_max_residuals.back(),
                       " iterations=", snapshot.iterations);
      };
      auto print_skip_incompatible = [&](double source_threshold,
                                         const RestartSnapshot &snapshot,
                                         const char *reason) {
        if (world.rank() != 0 || config_.print_level <= 0) {
          return;
        }
        madness::print("EXCITED_RESTART_SKIP threshold=", source_threshold,
                       " reason=", reason,
                       " snapshot_kind=", snapshot.snapshot_kind,
                       " variant=", snapshot.response_variant,
                       " restart_support=", snapshot.restart_support_mode,
                       " restart_capable=", snapshot.restart_capable);
      };
      if (current_protocol_snapshot.has_data) {
        if (is_stalled_snapshot(current_protocol_snapshot)) {
          print_skip_stalled(input.threshold, current_protocol_snapshot);
        } else if (snapshot_supports_full_restart(current_protocol_snapshot, input)) {
          selected.snapshot = current_protocol_snapshot;
          selected.found = true;
          selected.from_current_protocol = true;
          selected.uses_full_bundle = true;
          selected.source_threshold = input.threshold;
          return selected;
        } else {
          print_skip_incompatible(input.threshold, current_protocol_snapshot,
                                  "current_protocol_not_restart_capable");
        }
      }

      const size_t protocol_count = config_.protocols.size();
      const size_t restart_scan_end = std::min(input.protocol_index, protocol_count);
      for (size_t prev = restart_scan_end; prev > 0; --prev) {
        const double threshold = config_.protocols[prev - 1];
        const auto snapshot =
            read_restart_snapshot(world, protocol_restart_file(threshold));
        if (!snapshot.has_data) {
          continue;
        }
        if (is_stalled_snapshot(snapshot)) {
          print_skip_stalled(threshold, snapshot);
          continue;
        }
        if (!snapshot_supports_full_restart(snapshot, input)) {
          print_skip_incompatible(threshold, snapshot,
                                  "lower_protocol_not_restart_capable");
          continue;
        }
        selected.snapshot = snapshot;
        selected.found = true;
        selected.from_lower_protocol = true;
        selected.uses_full_bundle = true;
        selected.source_threshold = threshold;
        break;
      }
      if (!selected.found) {
        const auto guess_snapshot =
            read_restart_snapshot(world, guess_archive_file());
        if (snapshot_can_seed_guess(guess_snapshot, input)) {
          selected.snapshot = guess_snapshot;
          selected.found = true;
          selected.from_guess_archive = true;
          selected.uses_full_bundle = guess_snapshot.bundle_state_present &&
                                      snapshot_variant_matches_input(guess_snapshot,
                                                                    input);
          selected.source_threshold = guess_snapshot.protocol_threshold;
        }
      }
      return selected;
    }

    [[nodiscard]] ProtocolInitResult
    initialize_protocol_guess(madness::World &world,
                              const ExcitedBundleProtocolInput &input,
                              const RestartSnapshot &current_protocol_snapshot) {
      ensure_ground_data(world);
      const auto restart_choice =
          select_restart_seed(world, input, current_protocol_snapshot);
      if (restart_choice.found) {
        load_restart_seed(restart_choice.snapshot);
        if (restart_choice.uses_full_bundle) {
          align_active_bundle_protocol(world, input);
        }
        ensure_trial_space_matches_guess(world, input);
        ensure_state_names_for_protocol(input.threshold);
        configure_iteration_contract(input, expected_response_variant(input.tda));
        if (world.rank() == 0) {
          if (restart_choice.from_current_protocol) {
            madness::print("EXCITED_INIT strategy=protocol_restart_guess states=",
                           omega_.size(),
                           " variant=", restart_choice.snapshot.response_variant,
                           " snapshot_kind=",
                           restart_choice.snapshot.snapshot_kind,
                           " bundle_seed=", restart_choice.uses_full_bundle);
          } else if (restart_choice.from_lower_protocol) {
            madness::print(
                "EXCITED_INIT strategy=lower_protocol_restart_guess "
                "source_threshold=",
                restart_choice.source_threshold, " states=", omega_.size(),
                " variant=", restart_choice.snapshot.response_variant,
                " snapshot_kind=", restart_choice.snapshot.snapshot_kind,
                " bundle_seed=", restart_choice.uses_full_bundle);
          } else if (restart_choice.from_guess_archive) {
            madness::print("EXCITED_INIT strategy=guess_archive states=",
                           omega_.size(),
                           " variant=", restart_choice.snapshot.response_variant,
                           " snapshot_kind=", restart_choice.snapshot.snapshot_kind,
                           " bundle_seed=", restart_choice.uses_full_bundle);
          }
        }
        ProtocolInitResult init;
        init.stage_status = restart_choice.from_current_protocol
                                ? "protocol_restart_guess"
                                : (restart_choice.from_lower_protocol
                                       ? "lower_protocol_restart_guess"
                                       : "guess_archive");
        init.restart_reused = true;
        init.restart_source = restart_choice.from_current_protocol
                                  ? "current_protocol_snapshot"
                                  : (restart_choice.from_lower_protocol
                                         ? "lower_protocol_snapshot"
                                         : "guess_archive");
        init.snapshot_kind = restart_choice.snapshot.snapshot_kind;
        init.response_variant = restart_choice.snapshot.response_variant;
        init.restart_support_mode = restart_choice.snapshot.restart_support_mode;
        init.bundle_state_present = restart_choice.snapshot.bundle_state_present;
        init.restart_capable = restart_choice.snapshot.restart_capable;
        init.restart_source_threshold = restart_choice.source_threshold;
        return init;
      }

      if (has_active_guess_ && !omega_.empty()) {
        if (!active_response_bundle_.empty()) {
          align_active_bundle_protocol(world, input);
        }
        ensure_trial_space_matches_guess(world, input);
        ensure_state_names_for_protocol(input.threshold);
        configure_iteration_contract(input, expected_response_variant(input.tda));
        if (world.rank() == 0) {
          madness::print("EXCITED_INIT strategy=carryover_guess states=",
                         omega_.size());
        }
        ProtocolInitResult init;
        init.stage_status = "carryover_guess";
        init.restart_reused = true;
        init.restart_source = "carryover";
        init.response_variant = expected_response_variant(input.tda);
        init.restart_support_mode =
            restart_support_mode_for_variant(init.response_variant);
        init.snapshot_kind = active_response_bundle_.empty() ? "none"
                                                             : "carryover_bundle";
        init.bundle_state_present = !active_response_bundle_.empty();
        init.restart_capable = init.bundle_state_present &&
                               variant_supports_full_restart(init.response_variant);
        init.restart_source_threshold = input.threshold;
        return init;
      }

      build_fresh_guess(world, input);
      ensure_state_names_for_protocol(input.threshold);
      configure_iteration_contract(input, expected_response_variant(input.tda));
      write_guess_archive(world, false, last_iteration_count_, input);
      if (world.rank() == 0) {
        madness::print("EXCITED_INIT strategy=fresh_guess states=",
                       omega_.size());
      }
      ProtocolInitResult init;
      init.stage_status = "protocol_fresh_guess";
      init.restart_reused = false;
      init.restart_source = "fresh_guess";
      init.response_variant = expected_response_variant(input.tda);
      init.restart_support_mode =
          restart_support_mode_for_variant(init.response_variant);
      init.snapshot_kind = active_response_bundle_.empty() ? "guess_trial_seed"
                                                           : "guess_bundle";
      init.bundle_state_present = !active_response_bundle_.empty();
      init.restart_capable =
          init.bundle_state_present &&
          variant_supports_full_restart(init.response_variant);
      init.restart_source_threshold = input.threshold;
      return init;
    }

    [[nodiscard]] ExcitedBundleProtocolResult
    iterate(madness::World &world, const ExcitedBundleProtocolInput &input,
            bool restart_seed_reused) {
      ensure_ground_data(world);
      ensure_trial_space_matches_guess(world, input);
      if (omega_.empty()) {
        build_fresh_guess(world, input);
        ensure_trial_space_matches_guess(world, input);
      }
      if (trial_space_.x_states.size() > omega_.size()) {
        trial_space_.x_states.resize(omega_.size());
      } else if (trial_space_.x_states.size() < omega_.size()) {
        omega_.resize(trial_space_.x_states.size());
      }
      trial_space_.omega = omega_;
      ensure_state_names_for_protocol(input.threshold);
      auto response_states = build_response_bundle_seed(world, input);
      if (response_states.size() > omega_.size()) {
        response_states.resize(omega_.size());
      } else if (response_states.size() < omega_.size()) {
        omega_.resize(response_states.size());
      }
      sync_trial_space_from_response_bundle(world, response_states);

      const std::string response_variant =
          response_states.empty() ? expected_response_variant(input.tda)
                                  : response_variant_name(response_states.front());
      const auto contract = make_iteration_contract(input, response_variant);
      configure_iteration_contract(input, response_variant);
      const double effective_max_rotation =
          contract.step_restriction_active ? contract.max_rotation : 0.0;

      const size_t max_iterations = std::max<size_t>(input.maxiter, 1);
      const double convergence_target = std::max(10.0 * input.threshold, 1.0e-7);
      bool converged = false;
      size_t used_iterations = 0;
      residual_norms_.assign(omega_.size(), 0.0);
      density_change_norms_.assign(omega_.size(), 0.0);
      relative_residual_norms_.assign(omega_.size(), 0.0);
      iteration_max_residuals_.clear();
      iteration_max_density_changes_.clear();
      iteration_max_relative_residuals_.clear();
      std::vector<response_solver> state_accelerators;
      if (!response_states.empty() && contract.accelerator_active) {
        state_accelerators.reserve(response_states.size());
        for (const auto &response : response_states) {
          const size_t flat_size =
              std::visit([](const auto &typed_response) {
                return typed_response.flat.size();
              },
                         response);
          state_accelerators.emplace_back(
              response_vector_allocator(world, flat_size), false);
          state_accelerators.back().set_maxsub(
              static_cast<int>(std::max<size_t>(contract.accelerator_subspace,
                                                static_cast<size_t>(1))));
        }
      }
      if (world.rank() == 0 && config_.print_level > 0) {
        madness::print(
            "EXCITED_ITER_START protocol=", input.threshold, " states=",
            omega_.size(), " maxiter=", max_iterations,
            " convergence_mode=", contract.convergence_mode,
            " density_target=", contract.density_target,
            " relative_target=", contract.relative_target,
            " fallback_max_residual_target=", convergence_target,
            " accelerator_mode=", contract.accelerator_mode,
            " accelerator_subspace=", contract.accelerator_subspace,
            " max_rotation=", contract.max_rotation,
            " step_restriction_active=", contract.step_restriction_active,
            " restart_seed_reused=", restart_seed_reused,
            " trial_space_k=", trial_space_k_,
            " solve_variant=", response_variant);
      }

      for (size_t iter = 0; iter < max_iterations; ++iter) {
        used_iterations = iter + 1;
        std::vector<double> state_residuals(omega_.size(), 0.0);
        std::vector<double> state_density_changes(omega_.size(), 0.0);
        std::vector<double> state_relative_residuals(omega_.size(), 0.0);
        if (world.rank() == 0 && config_.print_level > 1) {
          madness::print("EXCITED_ITER_OMEGA_PRE iter=", used_iterations,
                         " omega_head=", format_vector_head(omega_));
        }
        bool bundle_rotated = false;
        std::vector<size_t> bundle_slot_order;
        if (!response_states.empty()) {
          const double damping = (iter == 0 ? 0.85 : 0.65);
          std::visit(
              [&](auto &typed_first) {
                using ResponseType = std::decay_t<decltype(typed_first)>;
                if constexpr (std::is_same_v<ResponseType,
                                             StaticUnrestrictedResponse> ||
                              std::is_same_v<ResponseType,
                                             DynamicUnrestrictedResponse>) {
                  for (size_t i = 0; i < response_states.size() && i < omega_.size();
                       ++i) {
                    if (world.rank() == 0 && config_.print_level > 1) {
                      madness::print("EXCITED_STATE_BEGIN iter=", used_iterations,
                                     " state=", i, " omega=", omega_[i]);
                    }
                    const auto metrics = iterate_response_state(
                        world, response_states[i], omega_[i], input.threshold,
                        damping, i, used_iterations);
                    state_residuals[i] = metrics.residual;
                    state_density_changes[i] = metrics.density_change;
                    state_relative_residuals[i] = metrics.relative_residual;
                    if (world.rank() == 0 && config_.print_level > 1) {
                      madness::print("EXCITED_STATE_DONE iter=", used_iterations,
                                     " state=", i,
                                     " residual=", metrics.residual,
                                     " relative_residual=",
                                     metrics.relative_residual,
                                     " drho=", metrics.density_change);
                    }
                    world.gop.fence();
                  }
                } else {
                  std::vector<ResponseType> typed_states;
                  typed_states.reserve(response_states.size());
                  for (const auto &response : response_states) {
                    typed_states.push_back(std::get<ResponseType>(response));
                  }
                  bundle_rotated = iterate_typed_bundle_legacy_sequence(
                      world, typed_states, omega_, input.threshold, damping,
                      used_iterations, state_residuals, state_density_changes,
                      state_relative_residuals, effective_max_rotation,
                      state_accelerators.empty() ? nullptr
                                                 : &state_accelerators,
                      &bundle_slot_order);
                  if (response_states.size() > typed_states.size()) {
                    response_states.resize(typed_states.size());
                  }
                  for (size_t i = 0; i < typed_states.size(); ++i) {
                    response_states[i] = std::move(typed_states[i]);
                  }
                }
              },
              response_states.front());
          if (bundle_rotated) {
            apply_slot_reordering(bundle_slot_order, input.threshold);
          }
          sync_trial_space_from_response_bundle(world, response_states);
        }

        double max_residual = 0.0;
        double max_density_change = 0.0;
        double max_relative_residual = 0.0;
        if (!bundle_rotated) {
          auto target_omega = trial_space_.x_states.empty()
                                  ? omega_
                                  : estimate_state_energies(world, trial_space_.x_states);
          if (target_omega.size() < omega_.size()) {
            target_omega.resize(omega_.size(), 0.0);
          }
          for (size_t i = 0; i < omega_.size(); ++i) {
            const double delta = target_omega[i] - omega_[i];
            omega_[i] += 0.45 * delta;
            state_residuals[i] = std::max(state_residuals[i], std::fabs(delta));
            max_residual = std::max(max_residual, state_residuals[i]);
            max_density_change =
                std::max(max_density_change, state_density_changes[i]);
            max_relative_residual =
                std::max(max_relative_residual, state_relative_residuals[i]);
          }
        } else {
          for (size_t i = 0; i < omega_.size(); ++i) {
            max_residual = std::max(max_residual, state_residuals[i]);
            max_density_change =
                std::max(max_density_change, state_density_changes[i]);
            max_relative_residual =
                std::max(max_relative_residual, state_relative_residuals[i]);
          }
        }
        sync_root_descriptors_for_protocol(input.threshold);
        residual_norms_ = state_residuals;
        density_change_norms_ = state_density_changes;
        relative_residual_norms_ = state_relative_residuals;
        iteration_max_residuals_.push_back(max_residual);
        iteration_max_density_changes_.push_back(max_density_change);
        iteration_max_relative_residuals_.push_back(max_relative_residual);
        if (world.rank() == 0 && config_.print_level > 1) {
          const auto [worst_state, worst_residual] =
              max_abs_value_with_index(state_residuals);
          const auto [worst_density_state, worst_density] =
              max_abs_value_with_index(state_density_changes);
          const auto [worst_relative_state, worst_relative] =
              max_abs_value_with_index(state_relative_residuals);
          const double worst_omega =
              (worst_state < omega_.size()) ? omega_[worst_state] : 0.0;
          const std::string worst_name =
              (worst_state < state_names_.size()) ? state_names_[worst_state]
                                                  : std::string("es?");
          madness::print(
              "EXCITED_ITER iter=", used_iterations, "/", max_iterations,
              " max_residual=", max_residual,
              " max_density_change=", max_density_change,
              " max_relative_residual=", max_relative_residual,
              " worst_state=", worst_state,
              " worst_state_residual=", worst_residual,
              " worst_state_name=", worst_name,
              " worst_state_omega=", worst_omega,
              " worst_density_state=", worst_density_state,
              " worst_density_change=", worst_density,
              " worst_relative_state=", worst_relative_state,
              " worst_relative_residual=", worst_relative,
              " omega_head=", format_vector_head(omega_),
              " residual_head=", format_vector_head(state_residuals),
              " relative_head=", format_vector_head(state_relative_residuals),
              " drho_head=", format_vector_head(state_density_changes));
        }
        if (contract.use_dual_gate) {
          converged = max_density_change <= contract.density_target &&
                      max_relative_residual <= contract.relative_target;
        } else {
          converged = max_residual <= convergence_target;
        }
        if (converged) {
          converged = true;
          break;
        }
      }
      sync_trial_space_from_response_bundle(world, response_states);
      trial_space_.num_states = trial_space_.x_states.size();
      last_iteration_count_ = used_iterations;
      trial_space_.omega = omega_;
      active_response_bundle_ = response_states;
      active_bundle_k_ = current_protocol_k();
      active_bundle_thresh_ = current_protocol_thresh();
      sync_root_descriptors_for_protocol(input.threshold);
      if (world.rank() == 0 && config_.print_level > 0) {
        const double final_max_residual =
            iteration_max_residuals_.empty() ? 0.0 : iteration_max_residuals_.back();
        const double final_max_density =
            iteration_max_density_changes_.empty()
                ? 0.0
                : iteration_max_density_changes_.back();
        const double final_max_relative =
            iteration_max_relative_residuals_.empty()
                ? 0.0
                : iteration_max_relative_residuals_.back();
        madness::print("EXCITED_ITER_DONE protocol=", input.threshold,
                       " converged=", converged,
                       " iterations=", used_iterations,
                       " final_max_residual=", final_max_residual,
                       " final_max_density_change=", final_max_density,
                       " final_max_relative_residual=", final_max_relative,
                       " convergence_mode=", convergence_mode_,
                       " accelerator_mode=", accelerator_mode_,
                       " state_names_head=", format_string_head(state_names_),
                       " omega_head=", format_vector_head(omega_));
      }

      RestartSnapshot snapshot;
      snapshot.has_data = true;
      snapshot.converged = converged;
      snapshot.iterations = used_iterations;
      populate_snapshot_header(
          world, snapshot, input,
          active_response_bundle_.empty() ? "protocol_trial_seed"
                                          : "protocol_bundle");
      snapshot.energies = omega_;
      snapshot.state_names = state_names_;
      snapshot.roots = stable_sorted_root_descriptors();
      snapshot.slot_permutation = slot_permutation_;
      snapshot.residual_norms = residual_norms_;
      snapshot.density_change_norms = density_change_norms_;
      snapshot.relative_residual_norms = relative_residual_norms_;
      snapshot.iteration_max_residuals = iteration_max_residuals_;
      snapshot.iteration_max_density_changes = iteration_max_density_changes_;
      snapshot.iteration_max_relative_residuals =
          iteration_max_relative_residuals_;
      snapshot.convergence_mode = convergence_mode_;
      snapshot.accelerator_mode = accelerator_mode_;
      snapshot.accelerator_subspace = accelerator_subspace_;
      snapshot.density_convergence_target = density_convergence_target_;
      snapshot.relative_convergence_target =
          relative_convergence_target_;
      snapshot.max_rotation = max_rotation_;
      snapshot.has_trial_states = !trial_space_.x_states.empty();
      snapshot.trial_space_k = trial_space_k_;
      snapshot.trial_space_num_orbitals = ground_num_orbitals_;
      snapshot.trial_states = trial_space_.x_states;
      write_restart_snapshot(world, protocol_restart_file(input.threshold),
                             snapshot);
      write_guess_archive(world, converged, used_iterations, input);

      ExcitedBundleProtocolResult result;
      result.attempted = true;
      result.saved = true;
      result.converged = converged;
      result.failed = !converged;
      result.skipped = false;
      result.restart_reused = restart_seed_reused;
      result.stage_status =
          converged ? "iterate_converged" : "iterate_reached_iter_budget";
      if (result.restart_reused) {
        result.stage_status = "restart_reused__" + result.stage_status;
      }
      result.iterations = used_iterations;
      result.energies = omega_;
      result.state_names = state_names_;
      result.roots = stable_sorted_root_descriptors();
      result.slot_permutation = slot_permutation_;
      result.residual_norms = residual_norms_;
      result.density_change_norms = density_change_norms_;
      result.relative_residual_norms = relative_residual_norms_;
      result.iteration_max_residuals = iteration_max_residuals_;
      result.iteration_max_density_changes = iteration_max_density_changes_;
      result.iteration_max_relative_residuals =
          iteration_max_relative_residuals_;
      result.convergence_mode = convergence_mode_;
      result.accelerator_mode = accelerator_mode_;
      result.accelerator_subspace = accelerator_subspace_;
      result.density_convergence_target = density_convergence_target_;
      result.relative_convergence_target =
          relative_convergence_target_;
      result.max_rotation = max_rotation_;
      result.response_variant = snapshot.response_variant;
      result.restart_support_mode = snapshot.restart_support_mode;
      result.snapshot_kind = snapshot.snapshot_kind;
      result.bundle_state_present = snapshot.bundle_state_present;
      result.restart_capable = snapshot.restart_capable;
      return result;
    }

    [[nodiscard]] std::string protocol_restart_file(double threshold) const {
      std::ostringstream os;
      os << config_.output_prefix << ".excited_bundle."
         << protocol_key(threshold) << ".restartdata";
      return os.str();
    }
  };

  ExcitedBundleSolverConfig config_;
  mutable ExcitedProtocolWorkflow workflow_;
};

} // namespace

std::unique_ptr<ExcitedStateBundleSolver>
make_excited_state_bundle_solver_adapter(const ExcitedBundleSolverConfig &config) {
  if (config.archive_file.empty()) {
    return std::make_unique<ExcitedStateBundleNoopSolver>();
  }
  return std::make_unique<RestartAwareExcitedScaffoldSolver>(config);
}
