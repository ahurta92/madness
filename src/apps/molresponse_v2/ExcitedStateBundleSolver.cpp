#include "ExcitedStateBundleSolver.hpp"

#include "GroundStateData.hpp"
#include "ResponseSolverUtils.hpp"

#include <madness/external/nlohmann_json/json.hpp>
#include <madness/world/worldgop.h>
#include <madness/mra/funcdefaults.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
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

bool file_exists_world(madness::World &world, const std::string &path) {
  bool exists = false;
  if (world.rank() == 0) {
    exists = fs::exists(path);
  }
  world.gop.broadcast_serializable(exists, 0);
  return exists;
}

struct RestartSnapshot {
  bool has_data = false;
  bool converged = false;
  size_t iterations = 0;
  std::vector<double> energies;
  std::vector<double> residual_norms;
  std::vector<double> iteration_max_residuals;
};

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
        snapshot.converged = j.value("converged", false);
        snapshot.iterations = j.value("iterations", static_cast<size_t>(0));
        if (j.contains("energies") && j["energies"].is_array()) {
          snapshot.energies = j["energies"].get<std::vector<double>>();
        }
        if (j.contains("residual_norms") && j["residual_norms"].is_array()) {
          snapshot.residual_norms =
              j["residual_norms"].get<std::vector<double>>();
        }
        if (j.contains("iteration_max_residuals") &&
            j["iteration_max_residuals"].is_array()) {
          snapshot.iteration_max_residuals =
              j["iteration_max_residuals"].get<std::vector<double>>();
        }
      } catch (...) {
        snapshot = RestartSnapshot{};
      }
    }
  }
  world.gop.broadcast_serializable(snapshot.has_data, 0);
  world.gop.broadcast_serializable(snapshot.converged, 0);
  world.gop.broadcast_serializable(snapshot.iterations, 0);
  world.gop.broadcast_serializable(snapshot.energies, 0);
  world.gop.broadcast_serializable(snapshot.residual_norms, 0);
  world.gop.broadcast_serializable(snapshot.iteration_max_residuals, 0);
  return snapshot;
}

void write_restart_snapshot(madness::World &world, const std::string &path,
                            const RestartSnapshot &snapshot) {
  if (world.rank() == 0) {
    json j = {{"converged", snapshot.converged},
              {"iterations", snapshot.iterations},
              {"energies", snapshot.energies},
              {"residual_norms", snapshot.residual_norms},
              {"iteration_max_residuals", snapshot.iteration_max_residuals}};
    std::ofstream out(path);
    out << j.dump(2) << "\n";
  }
  world.gop.fence();
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

      // Matches old response flow: if this protocol is already converged, skip.
      if (input.restart_saved && input.restart_converged) {
        if (restart_snapshot.has_data && !restart_snapshot.energies.empty()) {
          load_restart_guess(restart_snapshot);
          ensure_trial_space_matches_guess(world, input);
        }
        result.attempted = false;
        result.saved = true;
        result.converged = true;
        result.failed = false;
        result.skipped = true;
        result.restart_reused = true;
        result.stage_status = "restart_ready_skip";
        result.iterations = restart_snapshot.iterations;
        result.energies = restart_snapshot.energies;
        result.residual_norms = restart_snapshot.residual_norms;
        result.iteration_max_residuals = restart_snapshot.iteration_max_residuals;
        return result;
      }

      prepare_protocol(world, input.threshold);

      const auto initialize_result =
          initialize_protocol_guess(world, input, restart_snapshot);

      result = iterate(world, input, initialize_result.restart_reused);
      if (!initialize_result.stage_status.empty()) {
        result.stage_status =
            initialize_result.stage_status + "__" + result.stage_status;
      }
      return result;
    }

  private:
    struct ProtocolInitResult {
      std::string stage_status;
      bool restart_reused = false;
    };

    struct RestartSeedChoice {
      RestartSnapshot snapshot;
      bool found = false;
      bool from_current_protocol = false;
      bool from_lower_protocol = false;
      double source_threshold = 0.0;
    };

    ExcitedBundleSolverConfig config_;
    bool has_active_guess_ = false;
    ExcitedTrialSpace trial_space_;
    std::vector<double> omega_;
    std::vector<double> residual_norms_;
    size_t last_iteration_count_ = 0;
    std::vector<double> iteration_max_residuals_;
    bool ground_loaded_ = false;
    size_t ground_num_orbitals_ = 0;
    Tensor<double> ground_energies_;
    vector_real_function_3d ground_orbitals_;

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
        madness::print("Excited bundle protocol prep: thresh=", thresh,
                       " k=", madness::FunctionDefaults<3>::get_k());
      }
      world.gop.fence();
    }

    void ensure_ground_data(madness::World &world) {
      if (ground_loaded_) {
        return;
      }
      GroundStateData ground(world, config_.archive_file, Molecule());
      ground.prepareOrbitals(world, madness::FunctionDefaults<3>::get_k(),
                             madness::FunctionDefaults<3>::get_thresh());
      ground_num_orbitals_ = static_cast<size_t>(ground.getNumOrbitals());
      ground_energies_ = ground.getEnergies();
      ground_orbitals_ = copy(world, ground.getOrbitals(), true);
      ground_loaded_ = true;
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

    [[nodiscard]] double shift_for_state_energy(double omega) const {
      // Keep all denominators negative for the BSH operator.
      double shift = 0.0;
      for (size_t p = 0; p < ground_num_orbitals_; ++p) {
        const double e = ground_energies_(long(p));
        if (e + omega + shift >= 0.0) {
          shift = std::max(shift, -0.05 - (e + omega));
        }
      }
      return shift;
    }

    [[nodiscard]] vector_real_function_3d
    bsh_update_state(madness::World &world, const vector_real_function_3d &state,
                     double omega, double protocol_threshold,
                     double &state_residual_out) const {
      state_residual_out = 0.0;
      if (state.empty()) {
        return state;
      }

      const double shift = shift_for_state_energy(omega);
      const double lo = 1.0e-8;
      auto bsh_ops = ResponseSolverUtils::make_bsh_operators_response(
          world, shift, omega, ground_energies_, lo);

      auto theta = copy(world, state, true);
      for (size_t p = 0; p < theta.size(); ++p) {
        const double coeff =
            -2.0 * (ground_energies_(long(p)) + omega + shift);
        theta[p] = coeff * theta[p];
      }
      truncate(world, theta, protocol_threshold, true);

      auto updated = apply(world, bsh_ops, theta);
      QProjector<double, 3> projector(ground_orbitals_);
      updated = projector(updated);
      truncate(world, updated, protocol_threshold, true);

      for (size_t p = 0; p < state.size(); ++p) {
        auto diff = updated[p] - state[p];
        state_residual_out = std::max(state_residual_out, diff.norm2());
      }
      return updated;
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
    make_random_trial(madness::World &world,
                      const ExcitedBundleProtocolInput &input) {
      const size_t m = std::max<size_t>(2 * input.num_states, 1);
      auto trial = make_derivative_trial(world, m, 0);
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
        return;
      }
      const size_t iters = std::max<size_t>(input.guess_max_iter, 1);
      for (size_t iter = 0; iter < iters; ++iter) {
        project_and_orthonormalize(world, trial.x_states);
        auto estimate = estimate_state_energies(world, trial.x_states);
        if (trial.omega.size() != estimate.size()) {
          trial.omega = estimate;
        } else {
          for (size_t i = 0; i < trial.omega.size(); ++i) {
            trial.omega[i] = 0.6 * trial.omega[i] + 0.4 * estimate[i];
          }
        }
        sort_by_energy(trial.x_states, trial.omega);
      }
    }

    void load_restart_guess(const RestartSnapshot &restart_snapshot) {
      omega_ = restart_snapshot.energies;
      residual_norms_ = restart_snapshot.residual_norms;
      if (residual_norms_.size() != omega_.size()) {
        residual_norms_.assign(omega_.size(), 1.0e-2);
      }
      last_iteration_count_ = restart_snapshot.iterations;
      iteration_max_residuals_ = restart_snapshot.iteration_max_residuals;
      has_active_guess_ = !omega_.empty();
    }

    void build_fresh_guess(madness::World &world,
                           const ExcitedBundleProtocolInput &input) {
      if (input.tda) {
        trial_space_ = create_trial_functions2(world, input);
      } else if (input.num_states <= 2) {
        trial_space_ = make_nwchem_trial(world, input);
      } else if (input.num_states <= 6) {
        trial_space_ = create_trial_functions(world, input);
      } else {
        trial_space_ = make_random_trial(world, input);
      }
      iterate_trial(world, input, trial_space_);
      if (trial_space_.x_states.size() > input.num_states) {
        trial_space_.x_states.resize(input.num_states);
        trial_space_.omega.resize(input.num_states);
      }
      omega_ = trial_space_.omega;
      residual_norms_ = std::vector<double>(omega_.size(), 1.0e-2);
      last_iteration_count_ = std::max<size_t>(input.guess_max_iter, 1);
      iteration_max_residuals_.clear();
      has_active_guess_ = !omega_.empty();
    }

    void ensure_trial_space_matches_guess(
        madness::World &world, const ExcitedBundleProtocolInput &input) {
      if (trial_space_.x_states.empty()) {
        trial_space_ = create_trial_functions(world, input);
        iterate_trial(world, input, trial_space_);
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
        omega_ = trial_space_.omega;
      }
      if (residual_norms_.size() != omega_.size()) {
        residual_norms_.assign(omega_.size(), 1.0e-2);
      }
      has_active_guess_ = !omega_.empty();
    }

    [[nodiscard]] RestartSeedChoice
    select_restart_seed(madness::World &world,
                        const ExcitedBundleProtocolInput &input,
                        const RestartSnapshot &current_protocol_snapshot) const {
      RestartSeedChoice selected;
      if (current_protocol_snapshot.has_data &&
          !current_protocol_snapshot.energies.empty()) {
        selected.snapshot = current_protocol_snapshot;
        selected.found = true;
        selected.from_current_protocol = true;
        selected.source_threshold = input.threshold;
        return selected;
      }

      double best_threshold = -1.0;
      for (const auto threshold : config_.protocols) {
        if (threshold >= input.threshold - 1.0e-14) {
          continue;
        }
        const auto snapshot =
            read_restart_snapshot(world, protocol_restart_file(threshold));
        if (!snapshot.has_data || snapshot.energies.empty()) {
          continue;
        }
        if (!selected.found || threshold > best_threshold) {
          selected.snapshot = snapshot;
          selected.found = true;
          selected.from_lower_protocol = true;
          selected.source_threshold = threshold;
          best_threshold = threshold;
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
        load_restart_guess(restart_choice.snapshot);
        ensure_trial_space_matches_guess(world, input);
        if (world.rank() == 0) {
          if (restart_choice.from_current_protocol) {
            madness::print("Excited bundle init: using protocol restart guess (",
                           omega_.size(), " states)");
          } else if (restart_choice.from_lower_protocol) {
            madness::print(
                "Excited bundle init: using lower-protocol restart guess (",
                "source_thresh=", restart_choice.source_threshold,
                ", states=", omega_.size(), ")");
          }
        }
        return {restart_choice.from_current_protocol
                    ? "protocol_restart_guess"
                    : "lower_protocol_restart_guess",
                true};
      }

      if (has_active_guess_ && !omega_.empty()) {
        ensure_trial_space_matches_guess(world, input);
        if (world.rank() == 0) {
          madness::print("Excited bundle init: using in-memory carryover guess (",
                         omega_.size(), " states)");
        }
        return {"carryover_guess", true};
      }

      build_fresh_guess(world, input);
      if (world.rank() == 0) {
        madness::print("Excited bundle init: fresh trial/guess path (",
                       omega_.size(), " states)");
      }
      return {"protocol_fresh_guess", false};
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

      const size_t max_iterations = std::max<size_t>(input.maxiter, 1);
      const double convergence_target = std::max(10.0 * input.threshold, 1.0e-7);
      bool converged = false;
      size_t used_iterations = 0;
      iteration_max_residuals_.clear();

      for (size_t iter = 0; iter < max_iterations; ++iter) {
        used_iterations = iter + 1;
        std::vector<double> state_residuals(omega_.size(), 0.0);
        if (!trial_space_.x_states.empty()) {
          const double damping = (iter == 0 ? 0.85 : 0.65);
          for (size_t i = 0; i < trial_space_.x_states.size(); ++i) {
            double state_residual = 0.0;
            auto updated = bsh_update_state(
                world, trial_space_.x_states[i], omega_[i], input.threshold,
                state_residual);
            trial_space_.x_states[i] = gaxpy_oop(
                damping, updated, 1.0 - damping, trial_space_.x_states[i], true);
            state_residuals[i] = state_residual;
          }
          project_and_orthonormalize(world, trial_space_.x_states);
        }

        auto target_omega = trial_space_.x_states.empty()
                                ? omega_
                                : estimate_state_energies(world, trial_space_.x_states);
        if (target_omega.size() < omega_.size()) {
          target_omega.resize(omega_.size(), 0.0);
        }

        double max_residual = 0.0;
        for (size_t i = 0; i < omega_.size(); ++i) {
          const double delta = target_omega[i] - omega_[i];
          omega_[i] += 0.45 * delta;
          state_residuals[i] = std::max(state_residuals[i], std::fabs(delta));
          max_residual = std::max(max_residual, state_residuals[i]);
        }
        residual_norms_ = state_residuals;
        iteration_max_residuals_.push_back(max_residual);
        if (max_residual <= convergence_target) {
          converged = true;
          break;
        }
      }
      last_iteration_count_ = used_iterations;

      RestartSnapshot snapshot;
      snapshot.has_data = true;
      snapshot.converged = converged;
      snapshot.iterations = used_iterations;
      snapshot.energies = omega_;
      snapshot.residual_norms = residual_norms_;
      snapshot.iteration_max_residuals = iteration_max_residuals_;
      write_restart_snapshot(world, protocol_restart_file(input.threshold),
                             snapshot);

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
      result.residual_norms = residual_norms_;
      result.iteration_max_residuals = iteration_max_residuals_;
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
