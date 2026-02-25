#include "FrequencyLoop.hpp"

bool solve_response_vector(World &world,
                           const ResponseManager &response_manager,
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
                 [&](auto &) -> bool {
                   throw std::logic_error(
                       "This response type isn’t implemented yet");
                 }},
      response_variant);
}

void promote_response_vector(World &world, const ResponseVector &x_in,
                             ResponseVector &x_out) {
  if (std::holds_alternative<StaticRestrictedResponse>(x_in)) {
    if (world.rank() == 0) {
      madness::print("🔁 Promoting static restricted → dynamic restricted");
    }
    const auto &prev_resp = std::get<StaticRestrictedResponse>(x_in);
    DynamicRestrictedResponse current_resp;
    current_resp.x_alpha = copy(world, prev_resp.x_alpha);
    current_resp.y_alpha = copy(world, prev_resp.x_alpha);
    current_resp.flatten();
    x_out = current_resp;
  } else if (std::holds_alternative<StaticUnrestrictedResponse>(x_in)) {
    if (world.rank() == 0) {
      madness::print("🔁 Promoting static unrestricted → dynamic unrestricted");
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
      madness::print("📥 Copying dynamic restricted response");
    }
    const auto &prev_resp = std::get<DynamicRestrictedResponse>(x_in);

    DynamicRestrictedResponse current_resp;
    current_resp.x_alpha = copy(world, prev_resp.x_alpha);
    current_resp.y_alpha = copy(world, prev_resp.y_alpha);
    current_resp.flatten();
    x_out = current_resp;
  } else if (std::holds_alternative<DynamicUnrestrictedResponse>(x_in)) {
    if (world.rank() == 0) {
      madness::print("📥 Copying dynamic unrestricted response");
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

  bool is_unrestricted = !ground_state.isSpinRestricted();
  auto num_orbitals = static_cast<int>(ground_state.getNumOrbitals());
  auto &logger = persistence.logger();

  ResponseVector previous_response = make_response_vector(
      num_orbitals, /*is_static=*/state_desc.is_static(0), is_unrestricted);
  bool have_previous_freq_response = false;


  for (size_t freq_index = 0; freq_index < state_desc.num_frequencies();
       ++freq_index) {
  ResponseVector x_0 = make_response_vector(
      num_orbitals, /*is_static=*/state_desc.is_static(freq_index), is_unrestricted);
    LinearResponsePoint pt{state_desc, thresh_index, freq_index};
    if(world.rank() == 0) {
      madness::print("🔔 Starting response solve for ", pt.perturbationDescription(),
                     " at thresh ", pt.threshold(),
                     " freq ", pt.frequency());
    }

    bool is_saved = persistence.is_saved(pt);
    bool should_solve =
        !is_saved ||
        (at_final_protocol && !persistence.is_converged(pt));
    if (!should_solve)
      continue;

    world.gop.fence();

    std::string restart_source_kind = "initial_guess";
    std::optional<double> source_protocol;
    std::optional<double> source_frequency;
    bool loaded_from_disk = false;
    bool promoted_from_static = false;

    bool loaded_guess = false;
    if (is_saved) {
      loaded_guess = load_response_vector(world, num_orbitals, pt, x_0);
      if (loaded_guess) {
        restart_source_kind = "same_protocol_archive";
        source_protocol = pt.threshold();
        source_frequency = pt.frequency();
        loaded_from_disk = true;
      }
    }
    if (!loaded_guess && thresh_index > 0) {
      LinearResponsePoint coarser_pt{state_desc, thresh_index - 1, freq_index};
      loaded_guess = load_response_vector(world, num_orbitals, coarser_pt, x_0);
      if (loaded_guess) {
        restart_source_kind = "coarser_protocol_archive";
        source_protocol = coarser_pt.threshold();
        source_frequency = coarser_pt.frequency();
        loaded_from_disk = true;
      }
    }

    if (!loaded_guess) {
      if (!pt.is_static()) {
        if (have_previous_freq_response) {
          x_0 = previous_response;
          restart_source_kind = "previous_frequency_memory";
          source_protocol = pt.threshold();
          if (freq_index > 0) {
            source_frequency = state_desc.frequency(freq_index - 1);
          }
          if (freq_index > 0) {
            LinearResponsePoint prev_freq_pt{state_desc, thresh_index,
                                             freq_index - 1};
            if (prev_freq_pt.is_static()) {
              promote_response_vector(world, previous_response, x_0);
              promoted_from_static = true;
            }
          }
        } else if (freq_index > 0) {
          LinearResponsePoint prev_freq_pt{state_desc, thresh_index,
                                           freq_index - 1};
          if (load_response_vector(world, num_orbitals, prev_freq_pt, x_0)) {
            restart_source_kind = "previous_frequency_archive";
            source_protocol = prev_freq_pt.threshold();
            source_frequency = prev_freq_pt.frequency();
            loaded_from_disk = true;
            world.gop.fence();
            if (prev_freq_pt.is_static()) {
              promote_response_vector(world, x_0, x_0);
              promoted_from_static = true;
            }
          } else {
            x_0 = initialize_guess_vector(world, ground_state, pt);
          }
        } else {
          x_0 = initialize_guess_vector(world, ground_state, pt);
        }
      } else {
        x_0 = initialize_guess_vector(world, ground_state, pt);
      }
    }

    // Run the solver with logging
    logger.start_state(pt);
    const double state_wall_start = madness::wall_time();
    const double state_cpu_start = madness::cpu_time();
    auto max_iter = response_manager.params().maxiter();
    auto conv_thresh = response_manager.params().dconv();
    bool converged = solve_response_vector(world, response_manager,
                                           ground_state, state_desc, pt, x_0,
                                           logger, max_iter, conv_thresh);
    const double state_wall_seconds = madness::wall_time() - state_wall_start;
    const double state_cpu_seconds = madness::cpu_time() - state_cpu_start;
    persistence.record_timing(pt, state_wall_seconds, state_cpu_seconds);
    persistence.record_restart_provenance(
        pt, restart_source_kind, loaded_from_disk, promoted_from_static,
        source_protocol, source_frequency);
    if (world.rank() == 0) {
      logger.print_timing_table(pt);
      logger.print_values_table(pt);
      madness::print("⏱️ State timing for ", pt.perturbationDescription(),
                     " at thresh ", pt.threshold(),
                     " freq ", pt.frequency(),
                     " wall=", state_wall_seconds,
                     "s cpu=", state_cpu_seconds, "s");
      madness::print("↩️ Restart source for ", pt.perturbationDescription(),
                     " at thresh ", pt.threshold(), " freq ", pt.frequency(),
                     ": kind=", restart_source_kind,
                     " loaded_from_disk=", loaded_from_disk,
                     " promoted_from_static=", promoted_from_static);
    }
    world.gop.fence();
    // save and record the response vector
    save_response_vector(world, pt, x_0);
    world.gop.fence();
    persistence.record_status(pt, converged);

    previous_response = x_0;
    have_previous_freq_response = true;
  }
}
