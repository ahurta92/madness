#include "FrequencyLoop.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>

namespace {
struct RestartGuessProvenance {
  std::string source_kind = "initial_guess";
  std::optional<double> source_protocol;
  std::optional<double> source_frequency;
  bool loaded_from_disk = false;
  bool promoted_from_static = false;
};

struct GuessSeedResult {
  ResponseVector guess;
  RestartGuessProvenance provenance;
};

struct PointSolveResult {
  ResponseVector response;
  ResponseSolveDiagnostics diagnostics;
  bool completed_without_exception = true;
};

struct SolverAttemptResult {
  ResponseSolveDiagnostics diagnostics;
  bool completed_without_exception = true;
  std::string exception_message;
};

static bool should_solve_point(StateSolvePersistence &persistence,
                               const LinearResponsePoint &pt,
                               bool at_final_protocol) {
  if (!persistence.force_retry_removed_frequencies()) {
    const LinearResponsePoint protocol0_pt{pt.desc, /*thresh_index=*/0,
                                           pt.freq_index};
    if (persistence.is_removed_from_frequency_set(protocol0_pt)) {
      return false;
    }
  }
  const bool is_saved = persistence.is_saved(pt);
  return !is_saved || (at_final_protocol && !persistence.is_converged(pt));
}

static constexpr double k_remove_frequency_residual_cutoff_final_protocol =
    3.0e-2;
static constexpr double k_remove_frequency_residual_cutoff_protocol0 = 1.0e2;
static constexpr size_t k_negative_alpha_run_cutoff_final_protocol = 2;
static constexpr size_t k_negative_alpha_run_cutoff_protocol0 = 6;

static void apply_failure_policy(ResponseSolveDiagnostics &diagnostics,
                                 bool at_first_protocol,
                                 bool at_final_protocol) {
  diagnostics.remove_from_frequency_set = false;
  diagnostics.residual_remove_cutoff =
      std::numeric_limits<double>::quiet_NaN();

  const bool apply_removal_policy = at_first_protocol || at_final_protocol;
  if (!apply_removal_policy) {
    diagnostics.failure_reason = "policy_not_applied_this_protocol";
    return;
  }

  const bool using_protocol0_policy = at_first_protocol && !at_final_protocol;
  const double residual_cutoff = using_protocol0_policy
                                     ? k_remove_frequency_residual_cutoff_protocol0
                                     : k_remove_frequency_residual_cutoff_final_protocol;
  const size_t negative_alpha_cutoff =
      using_protocol0_policy ? k_negative_alpha_run_cutoff_protocol0
                             : k_negative_alpha_run_cutoff_final_protocol;
  diagnostics.residual_remove_cutoff = residual_cutoff;

  if (diagnostics.converged) {
    diagnostics.failure_reason = "converged";
    return;
  }

  const bool residual_is_high =
      std::isfinite(diagnostics.final_residual_norm) &&
      diagnostics.final_residual_norm >= residual_cutoff;
  if (residual_is_high) {
    diagnostics.remove_from_frequency_set = true;
    diagnostics.failure_reason = using_protocol0_policy
                                     ? "protocol0_residual_above_cutoff"
                                     : "final_residual_above_cutoff";
    return;
  }

  if (diagnostics.max_consecutive_negative_alpha >= negative_alpha_cutoff) {
    diagnostics.remove_from_frequency_set = true;
    diagnostics.failure_reason = using_protocol0_policy
                                     ? "protocol0_consecutive_negative_alpha"
                                     : "consecutive_negative_alpha";
    return;
  }

  diagnostics.failure_reason = using_protocol0_policy
                                   ? "protocol0_unconverged_below_cutoff"
                                   : "unconverged_below_cutoff";
}

static void print_point_solve_report(
    const LinearResponsePoint &pt,
    const ResponseSolveDiagnostics &solve_diagnostics,
    const RestartGuessProvenance &provenance, bool used_fallback_retry,
    bool completed_without_exception,
    const std::string &solver_exception_message, double state_wall_seconds,
    double state_cpu_seconds) {
  auto maybe_double_or_dash = [](const std::optional<double> &value) -> std::string {
    if (!value.has_value()) {
      return "-";
    }
    std::ostringstream os;
    os << std::scientific << std::setprecision(6) << *value;
    return os.str();
  };

  madness::print("SOLVE_REPORT");
  madness::print("  state = ", pt.perturbationDescription());
  madness::print("  protocol = ", pt.threshold());
  madness::print("  frequency = ", pt.frequency());

  madness::print("  [timing]");
  madness::print("    wall_s = ", state_wall_seconds);
  madness::print("    cpu_s = ", state_cpu_seconds);

  madness::print("  [restart]");
  madness::print("    kind = ", provenance.source_kind);
  madness::print("    loaded_from_disk = ", provenance.loaded_from_disk);
  madness::print("    promoted_from_static = ", provenance.promoted_from_static);
  madness::print("    source_protocol = ",
                 maybe_double_or_dash(provenance.source_protocol));
  madness::print("    source_frequency = ",
                 maybe_double_or_dash(provenance.source_frequency));

  madness::print("  [result]");
  madness::print("    converged = ", solve_diagnostics.converged);
  madness::print("    iterations = ", solve_diagnostics.iterations_performed);
  madness::print("    final_residual = ", solve_diagnostics.final_residual_norm);
  madness::print("    final_drho = ", solve_diagnostics.final_density_change);
  madness::print("    final_alpha = ", solve_diagnostics.final_alpha);
  madness::print("    max_negative_alpha_run = ",
                 solve_diagnostics.max_consecutive_negative_alpha);
  madness::print("    reached_iter_limit = ",
                 solve_diagnostics.reached_iteration_limit);
  madness::print("    remove_from_frequency_set = ",
                 solve_diagnostics.remove_from_frequency_set);
  madness::print("    residual_remove_cutoff = ",
                 solve_diagnostics.residual_remove_cutoff);
  madness::print("    failure_reason = ", solve_diagnostics.failure_reason);
  madness::print("    fallback_retry = ", used_fallback_retry);

  if (!completed_without_exception) {
    madness::print("  [exception]");
    madness::print("    message = ", solver_exception_message);
    madness::print("    action = mark_unconverged_continue");
  }
  madness::print("------------------------------------------------------------");
}

static void print_solve_start(World &world, const LinearResponsePoint &pt,
                              const char *execution_path) {
  if (world.rank() != 0) {
    return;
  }
  const double point_stall_timeout_seconds =
      response_point_stall_timeout_seconds();
  const bool timeout_enabled = std::isfinite(point_stall_timeout_seconds) &&
                               point_stall_timeout_seconds > 0.0;
  madness::print("SOLVE_START");
  madness::print("  execution_path = ", execution_path);
  madness::print("  state = ", pt.perturbationDescription());
  madness::print("  protocol_index = ", pt.thresh_index);
  madness::print("  protocol_threshold = ", pt.threshold());
  madness::print("  frequency_index = ", pt.freq_index);
  madness::print("  frequency = ", pt.frequency());
  madness::print("  static_point = ", pt.is_static());
  madness::print("  point_stall_timeout_s = ",
                 timeout_enabled ? std::to_string(point_stall_timeout_seconds)
                                 : std::string("disabled"));
  madness::print("  point_key = ", pt.response_filename());
  madness::print("------------------------------------------------------------");
}

static GuessSeedResult seed_response_guess_for_point(
    World &world, const GroundStateData &ground_state,
    const LinearResponseDescriptor &state_desc, const LinearResponsePoint &pt,
    size_t thresh_index, int num_orbitals,
    StateSolvePersistence &persistence,
    const ResponseVector *previous_response, bool have_previous_freq_response,
    ResponseVector initial_guess) {
  RestartGuessProvenance provenance;
  bool has_loaded_guess = false;
  const bool is_saved = persistence.is_saved(pt);
  if (is_saved) {
    has_loaded_guess = load_response_vector(world, num_orbitals, pt, initial_guess);
    if (has_loaded_guess) {
      provenance.source_kind = "same_protocol_archive";
      provenance.source_protocol = pt.threshold();
      provenance.source_frequency = pt.frequency();
      provenance.loaded_from_disk = true;
    }
  }
  if (!has_loaded_guess && thresh_index > 0) {
    LinearResponsePoint coarser_pt{state_desc, thresh_index - 1, pt.freq_index};
    has_loaded_guess =
        load_response_vector(world, num_orbitals, coarser_pt, initial_guess);
    if (has_loaded_guess) {
      provenance.source_kind = "coarser_protocol_archive";
      provenance.source_protocol = coarser_pt.threshold();
      provenance.source_frequency = coarser_pt.frequency();
      provenance.loaded_from_disk = true;
    }
  }

  if (!has_loaded_guess) {
    if (!pt.is_static()) {
      if (have_previous_freq_response && previous_response != nullptr) {
        initial_guess = *previous_response;
        provenance.source_kind = "previous_frequency_memory";
        provenance.source_protocol = pt.threshold();
        if (pt.freq_index > 0) {
          provenance.source_frequency = state_desc.frequency(pt.freq_index - 1);
        }
        if (pt.freq_index > 0) {
          LinearResponsePoint prev_freq_pt{state_desc, thresh_index,
                                           pt.freq_index - 1};
          if (prev_freq_pt.is_static()) {
            promote_response_vector(world, *previous_response, initial_guess);
            provenance.promoted_from_static = true;
          }
        }
      } else if (pt.freq_index > 0) {
        LinearResponsePoint prev_freq_pt{state_desc, thresh_index,
                                         pt.freq_index - 1};
        if (load_response_vector(world, num_orbitals, prev_freq_pt,
                                 initial_guess)) {
          provenance.source_kind = "previous_frequency_archive";
          provenance.source_protocol = prev_freq_pt.threshold();
          provenance.source_frequency = prev_freq_pt.frequency();
          provenance.loaded_from_disk = true;
          world.gop.fence();
          if (prev_freq_pt.is_static()) {
            promote_response_vector(world, initial_guess, initial_guess);
            provenance.promoted_from_static = true;
          }
        } else {
          initial_guess = initialize_guess_vector(world, ground_state, pt);
        }
      } else {
        initial_guess = initialize_guess_vector(world, ground_state, pt);
      }
    } else {
      initial_guess = initialize_guess_vector(world, ground_state, pt);
    }
  }

  return GuessSeedResult{std::move(initial_guess), std::move(provenance)};
}

static std::optional<GuessSeedResult>
seed_retry_guess_from_previous_frequency_archive(
    World &world, const LinearResponseDescriptor &state_desc,
    const LinearResponsePoint &pt, size_t thresh_index, int num_orbitals,
    bool is_unrestricted, StateSolvePersistence &persistence) {
  if (pt.is_static() || pt.freq_index == 0) {
    return std::nullopt;
  }

  const LinearResponsePoint prev_freq_pt{state_desc, thresh_index,
                                         pt.freq_index - 1};
  if (!persistence.is_converged(prev_freq_pt)) {
    return std::nullopt;
  }

  ResponseVector retry_guess =
      make_response_vector(num_orbitals, /*is_static=*/state_desc.is_static(pt.freq_index),
                           is_unrestricted);
  if (!load_response_vector(world, num_orbitals, prev_freq_pt, retry_guess)) {
    return std::nullopt;
  }
  world.gop.fence();

  RestartGuessProvenance provenance;
  provenance.source_kind = "previous_frequency_archive_retry";
  provenance.source_protocol = prev_freq_pt.threshold();
  provenance.source_frequency = prev_freq_pt.frequency();
  provenance.loaded_from_disk = true;
  if (prev_freq_pt.is_static()) {
    promote_response_vector(world, retry_guess, retry_guess);
    provenance.promoted_from_static = true;
  }

  return GuessSeedResult{std::move(retry_guess), std::move(provenance)};
}

static bool should_retry_from_previous_frequency_seed(
    const LinearResponsePoint &pt, const RestartGuessProvenance &provenance,
    const ResponseSolveDiagnostics &diagnostics, bool at_first_protocol) {
  if (!at_first_protocol || diagnostics.converged) {
    return false;
  }
  if (pt.freq_index == 0 || pt.is_static()) {
    return false;
  }
  // Retry only when the first attempt was truly an initializer run.
  return provenance.source_kind == "initial_guess";
}

static SolverAttemptResult run_solver_attempt(
    World &world, const ResponseManager &response_manager,
    const GroundStateData &ground_state,
    const LinearResponseDescriptor &state_desc, const LinearResponsePoint &pt,
    ResponseVector &response_guess, ResponseDebugLogger &logger,
    bool at_first_protocol, bool at_final_protocol) {
  logger.start_state(pt);
  SolverAttemptResult attempt;
  const auto max_iter = response_manager.params().maxiter();
  const auto conv_thresh = response_manager.params().dconv();
  try {
    attempt.diagnostics = solve_response_vector(
        world, response_manager, ground_state, state_desc, pt, response_guess,
        logger, max_iter, conv_thresh);
  } catch (const std::exception &ex) {
    attempt.completed_without_exception = false;
    attempt.diagnostics.converged = false;
    attempt.diagnostics.iterations_performed = 0;
    attempt.diagnostics.final_residual_norm =
        std::numeric_limits<double>::infinity();
    attempt.diagnostics.final_density_change =
        std::numeric_limits<double>::infinity();
    attempt.diagnostics.final_alpha = std::numeric_limits<double>::quiet_NaN();
    attempt.diagnostics.reached_iteration_limit = false;
    attempt.exception_message = ex.what();
  } catch (...) {
    attempt.completed_without_exception = false;
    attempt.diagnostics.converged = false;
    attempt.diagnostics.iterations_performed = 0;
    attempt.diagnostics.final_residual_norm =
        std::numeric_limits<double>::infinity();
    attempt.diagnostics.final_density_change =
        std::numeric_limits<double>::infinity();
    attempt.diagnostics.final_alpha = std::numeric_limits<double>::quiet_NaN();
    attempt.diagnostics.reached_iteration_limit = false;
    attempt.exception_message = "unknown exception";
  }

  apply_failure_policy(attempt.diagnostics, at_first_protocol, at_final_protocol);
  if (!attempt.completed_without_exception &&
      (at_first_protocol || at_final_protocol) &&
      !attempt.diagnostics.converged) {
    const bool is_await_timeout_exception =
        attempt.exception_message.find("ThreadPool::await() timed out") !=
        std::string::npos;
    attempt.diagnostics.remove_from_frequency_set = true;
    if (is_await_timeout_exception) {
      attempt.diagnostics.failure_reason =
          at_first_protocol && !at_final_protocol
              ? "protocol0_threadpool_await_timeout"
              : "threadpool_await_timeout";
    } else {
      attempt.diagnostics.failure_reason =
          at_first_protocol && !at_final_protocol ? "protocol0_solver_exception"
                                                  : "solver_exception";
    }
  }

  return attempt;
}

static PointSolveResult solve_and_record_point(
    World &world, const ResponseManager &response_manager,
    const GroundStateData &ground_state,
    const LinearResponseDescriptor &state_desc, const LinearResponsePoint &pt,
    StateSolvePersistence &persistence, ResponseVector response_guess,
    const RestartGuessProvenance &provenance, bool at_final_protocol,
    bool at_first_protocol) {
  auto &logger = persistence.logger();

  const double state_wall_start = madness::wall_time();
  const double state_cpu_start = madness::cpu_time();
  const bool is_unrestricted = !ground_state.isSpinRestricted();
  const auto num_orbitals = static_cast<int>(ground_state.getNumOrbitals());

  RestartGuessProvenance final_provenance = provenance;
  bool used_fallback_retry = false;

  auto solver_attempt =
      run_solver_attempt(world, response_manager, ground_state, state_desc, pt,
                         response_guess, logger, at_first_protocol,
                         at_final_protocol);

  const bool fallback_retry_eligible = should_retry_from_previous_frequency_seed(
      pt, provenance, solver_attempt.diagnostics, at_first_protocol);
  if (fallback_retry_eligible) {
    auto retry_seed = seed_retry_guess_from_previous_frequency_archive(
        world, state_desc, pt, pt.thresh_index, num_orbitals, is_unrestricted,
        persistence);
    if (retry_seed.has_value()) {
      used_fallback_retry = true;
      final_provenance = retry_seed->provenance;
      if (world.rank() == 0) {
        madness::print(
            "SOLVE_RETRY state=", pt.perturbationDescription(),
            " protocol=", pt.threshold(), " frequency=", pt.frequency(),
            " seed_kind=previous_frequency_archive_retry");
      }
      response_guess = std::move(retry_seed->guess);
      solver_attempt =
          run_solver_attempt(world, response_manager, ground_state, state_desc, pt,
                             response_guess, logger, at_first_protocol,
                             at_final_protocol);
      if (at_first_protocol && !solver_attempt.diagnostics.converged) {
        solver_attempt.diagnostics.remove_from_frequency_set = true;
        if (solver_attempt.diagnostics.failure_reason ==
            "protocol0_unconverged_below_cutoff") {
          solver_attempt.diagnostics.failure_reason = "protocol0_retry_failed";
        }
      }
    } else if (at_first_protocol && !solver_attempt.diagnostics.converged) {
      solver_attempt.diagnostics.remove_from_frequency_set = true;
      solver_attempt.diagnostics.failure_reason =
          "protocol0_unconverged_no_retry_seed";
    }
  }

  ResponseSolveDiagnostics solve_diagnostics = solver_attempt.diagnostics;
  const bool completed_without_exception =
      solver_attempt.completed_without_exception;
  const std::string &solver_exception_message = solver_attempt.exception_message;

  const bool is_converged = solve_diagnostics.converged;

  const double state_wall_seconds = madness::wall_time() - state_wall_start;
  const double state_cpu_seconds = madness::cpu_time() - state_cpu_start;
  persistence.record_timing(pt, state_wall_seconds, state_cpu_seconds);
  persistence.record_restart_provenance(
      pt, final_provenance.source_kind, final_provenance.loaded_from_disk,
      final_provenance.promoted_from_static, final_provenance.source_protocol,
      final_provenance.source_frequency);
  persistence.record_solver_diagnostics(pt, solve_diagnostics,
                                        used_fallback_retry);
  if (world.rank() == 0) {
    if (TimedValueLogger::console_enabled()) {
      logger.print_timing_table(pt);
      logger.print_values_table(pt);
    }
    print_point_solve_report(
        pt, solve_diagnostics, final_provenance, used_fallback_retry,
        completed_without_exception, solver_exception_message,
        state_wall_seconds, state_cpu_seconds);
    if (solve_diagnostics.remove_from_frequency_set) {
      madness::print("WARN FREQUENCY_REMOVAL_SUMMARY state=",
                     pt.perturbationDescription());
      madness::print("WARN FREQUENCY_REMOVAL_SUMMARY protocol=", pt.threshold());
      madness::print("WARN FREQUENCY_REMOVAL_SUMMARY frequency=", pt.frequency());
      madness::print("WARN FREQUENCY_REMOVAL_SUMMARY reason=",
                     solve_diagnostics.failure_reason);
    }
    if (!completed_without_exception) {
      madness::print("WARN SOLVER_EXCEPTION_SUMMARY state=",
                     pt.perturbationDescription());
      madness::print("WARN SOLVER_EXCEPTION_SUMMARY protocol=", pt.threshold());
      madness::print("WARN SOLVER_EXCEPTION_SUMMARY frequency=", pt.frequency());
      madness::print("WARN SOLVER_EXCEPTION_SUMMARY message=",
                     solver_exception_message);
    }
  }

  world.gop.fence();
  save_response_vector(world, pt, response_guess);
  world.gop.fence();
  persistence.record_status(pt, is_converged);
  persistence.maybe_poll_progress(world, /*force=*/false);
  return PointSolveResult{std::move(response_guess),
                          std::move(solve_diagnostics),
                          completed_without_exception};
}

} // namespace

ResponseSolveDiagnostics
solve_response_vector(World &world, const ResponseManager &response_manager,
                      const GroundStateData &g_s,
                      const LinearResponseDescriptor &state,
                      const LinearResponsePoint &pt,
                      ResponseVector &response_variant,
                      ResponseDebugLogger &logger, size_t max_iter,
                      double conv_thresh) {
  return std::visit(
      overloaded{// Only wire types that are READY:
                 [&](StaticRestrictedResponse &vector) {
                   return iterate(world, response_manager, g_s, state, pt,
                                  vector, logger, max_iter, conv_thresh);
                 },
                 [&](DynamicRestrictedResponse &vector) {
                   return iterate(world, response_manager, g_s, state, pt,
                                  vector, logger, max_iter, conv_thresh);
                 },
                 [&](auto &) -> ResponseSolveDiagnostics {
                   throw std::logic_error(
                       "This response type isnt implemented yet");
                 }},
      response_variant);
}

void promote_response_vector(World &world, const ResponseVector &x_in,
                             ResponseVector &x_out) {
  if (std::holds_alternative<StaticRestrictedResponse>(x_in)) {
    if (world.rank() == 0) {
      madness::print(
          "RESPONSE_PROMOTE from=static_restricted to=dynamic_restricted");
    }
    const auto &prev_resp = std::get<StaticRestrictedResponse>(x_in);
    DynamicRestrictedResponse current_resp;
    current_resp.x_alpha = copy(world, prev_resp.x_alpha);
    current_resp.y_alpha = copy(world, prev_resp.x_alpha);
    current_resp.flatten();
    x_out = current_resp;
  } else if (std::holds_alternative<StaticUnrestrictedResponse>(x_in)) {
    if (world.rank() == 0) {
      madness::print("RESPONSE_PROMOTE from=static_unrestricted "
                     "to=dynamic_unrestricted");
    }
    const auto &prev_resp = std::get<StaticUnrestrictedResponse>(x_in);

    DynamicUnrestrictedResponse current_resp;
    current_resp.x_alpha = copy(world, prev_resp.x_alpha);
    current_resp.x_beta = copy(world, prev_resp.x_beta);
    current_resp.y_alpha = copy(world, prev_resp.x_alpha);
    current_resp.y_beta = copy(world, prev_resp.x_beta);
    current_resp.flatten();
    x_out = current_resp;
  } else if (std::holds_alternative<DynamicRestrictedResponse>(x_in)) {
    if (world.rank() == 0) {
      madness::print("RESPONSE_COPY type=dynamic_restricted");
    }
    const auto &prev_resp = std::get<DynamicRestrictedResponse>(x_in);

    DynamicRestrictedResponse current_resp;
    current_resp.x_alpha = copy(world, prev_resp.x_alpha);
    current_resp.y_alpha = copy(world, prev_resp.y_alpha);
    current_resp.flatten();
    x_out = current_resp;
  } else if (std::holds_alternative<DynamicUnrestrictedResponse>(x_in)) {
    if (world.rank() == 0) {
      madness::print("RESPONSE_COPY type=dynamic_unrestricted");
    }
    const auto &prev_resp = std::get<DynamicUnrestrictedResponse>(x_in);

    DynamicUnrestrictedResponse current_resp;
    current_resp.x_alpha = copy(world, prev_resp.x_alpha);
    current_resp.x_beta = copy(world, prev_resp.x_beta);
    current_resp.y_alpha = copy(world, prev_resp.y_alpha);
    current_resp.y_beta = copy(world, prev_resp.y_beta);
    current_resp.flatten();
    x_out = current_resp;
  } else {
    throw std::runtime_error(
        "Unknown response variant in promote_response_vector");
  }
}

void computeFrequencyLoop(World &world,
                          const ResponseManager &response_manager,
                          const LinearResponseDescriptor &state_desc,
                          size_t thresh_index,
                          const GroundStateData &ground_state,
                          StateSolvePersistence &persistence,
                          bool at_final_protocol) {

  const bool is_unrestricted = !ground_state.isSpinRestricted();
  const auto num_orbitals = static_cast<int>(ground_state.getNumOrbitals());
  const bool at_first_protocol = (thresh_index == 0);

  ResponseVector previous_response = make_response_vector(
      num_orbitals, /*is_static=*/state_desc.is_static(0), is_unrestricted);
  bool has_previous_freq_response = false;
  size_t previous_response_freq_index = std::numeric_limits<size_t>::max();

  for (size_t freq_index = 0; freq_index < state_desc.num_frequencies();
       ++freq_index) {
    LinearResponsePoint pt{state_desc, thresh_index, freq_index};
    print_solve_start(world, pt, "frequency_loop");

    if (!should_solve_point(persistence, pt, at_final_protocol)) {
      continue;
    }

    world.gop.fence();

    ResponseVector initial_guess = make_response_vector(
        num_orbitals, /*is_static=*/state_desc.is_static(freq_index),
        is_unrestricted);
    const bool should_use_previous_frequency_memory =
        has_previous_freq_response &&
        previous_response_freq_index != std::numeric_limits<size_t>::max() &&
        previous_response_freq_index + 1 == freq_index;
    GuessSeedResult seeded_guess = seed_response_guess_for_point(
        world, ground_state, state_desc, pt, thresh_index, num_orbitals,
        persistence, should_use_previous_frequency_memory ? &previous_response
                                                          : nullptr,
        should_use_previous_frequency_memory, std::move(initial_guess));

    PointSolveResult point_result = solve_and_record_point(
        world, response_manager, ground_state, state_desc, pt, persistence,
        std::move(seeded_guess.guess), seeded_guess.provenance,
        at_final_protocol, at_first_protocol);
    if (point_result.completed_without_exception &&
        point_result.diagnostics.converged) {
      previous_response = std::move(point_result.response);
      has_previous_freq_response = true;
      previous_response_freq_index = freq_index;
    }
  }
}

void computeFrequencyPoint(World &world,
                           const ResponseManager &response_manager,
                           const LinearResponseDescriptor &state_desc,
                           size_t thresh_index, size_t freq_index,
                           const GroundStateData &ground_state,
                           StateSolvePersistence &persistence,
                           bool at_final_protocol) {
  const bool is_unrestricted = !ground_state.isSpinRestricted();
  const auto num_orbitals = static_cast<int>(ground_state.getNumOrbitals());
  const bool at_first_protocol = (thresh_index == 0);

  LinearResponsePoint pt{state_desc, thresh_index, freq_index};
  print_solve_start(world, pt, "frequency_point");

  if (!should_solve_point(persistence, pt, at_final_protocol)) {
    return;
  }

  world.gop.fence();
  ResponseVector initial_guess = make_response_vector(
      num_orbitals, /*is_static=*/state_desc.is_static(freq_index),
      is_unrestricted);
  GuessSeedResult seeded_guess = seed_response_guess_for_point(
      world, ground_state, state_desc, pt, thresh_index, num_orbitals,
      persistence, nullptr, false, std::move(initial_guess));

  (void)solve_and_record_point(world, response_manager, ground_state,
                               state_desc, pt, persistence,
                               std::move(seeded_guess.guess),
                               seeded_guess.provenance, at_final_protocol,
                               at_first_protocol);
}
