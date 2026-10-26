#pragma once

#include "../molresponse/solvers/function_kinetics.hpp"

#include "NamedStages.hpp"
#include <apps/molresponse/solvers/fd_save_load.hpp>
#include <madness/chem/StageCheckpoint.h>

#include <apps/molresponse/calc/calc_manager.hpp>
#include <madness/chem/LadderPlanner.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <map>
#include <set>

#include <iomanip>
#include <sstream>

namespace named_stages {

// B9: one response stage's convergence settings, built once. The engine rows
// come from the engine's own targets (ConvergencePolicy::effective_for_thresh);
// its gates are per leg and relative, so they are shown as formulas, with the
// per-channel verdicts in the channel records. The written targets come from
// fd_requested_accuracy_targets, the same source the solver stops on.
[[nodiscard]] inline madness::convergence::ConvergenceSettings
response_convergence_settings(const Stage& stage, double final_thresh) {
  namespace cv = madness::convergence;
  using cv::Group;
  using cv::Relation;
  using cv::SetBy;
  const auto& rp = stage.params.get<ResponseParameters>();
  const auto text = [](double value) { return cv::detail::scientific(value); };
  cv::ConvergenceSettings settings;
  settings.stage = stage.name;
  settings.application = "molresponse";
  settings.decider = rp.get<bool>(cv::use_default_criteria_key)
                         ? cv::Decider::engine_check
                         : cv::Decider::requested_criteria;
  settings.decider_set_by_user
      = stage.written.contains(cv::use_default_criteria_key);
  molresponse_v3::ConvergencePolicy policy;
  policy.dconv_user = rp.dconv();
  const auto limits = policy.effective_for_thresh(final_thresh);
  const auto dconv_by
      = stage.written.contains("dconv") ? SetBy::user : SetBy::engine_default;
  const std::string effective = "max(dconv " + text(rp.dconv()) + ", thresh "
                                + text(final_thresh) + ")";
  settings.metrics.push_back(
      {Group::engine_check, "BSH residual per leg", "dconv",
       "engine_bsh_residual", "", "response_amplitude_L2_atomic_units",
       rp.dconv(), dconv_by, Relation::less_than_or_equal, std::nullopt,
       text(limits.bsh_residual) + " + " + text(limits.rtol)
           + " x ||x|| (relative), from 5 x " + effective,
       "per channel; verdicts in the channel records", ""});
  settings.metrics.push_back(
      {Group::engine_check, "density change per leg", "dconv",
       "engine_density_change", "", "kernel_induced_density_L2_atomic_units",
       rp.dconv(), dconv_by, Relation::less_than_or_equal, std::nullopt,
       text(limits.density_residual) + " + " + text(limits.rtol)
           + " x ||rho|| (relative), from 5 x " + effective,
       "per channel; verdicts in the channel records", ""});

  struct Requested {
    const char *name, *keyword, *measurement, *quantity, *units;
    std::string engine_instead;
  };
  const std::vector<Requested> requested = {
      {"induced density change", "dconv", "charge_density_change",
       "kernel_induced_charge_density_change_L2",
       "kernel_induced_density_L2_atomic_units",
       "compares each leg's density change against "
           + text(limits.density_residual) + " + " + text(limits.rtol)
           + " x ||rho||"},
      {"BSH residual (max over functions)", "accuracy.bsh_residual",
       "bsh_residual", "maximum_raw_fixed_point_function_L2",
       "response_amplitude_L2_atomic_units",
       "compares each leg's BSH residual against " + text(limits.bsh_residual)
           + " + " + text(limits.rtol) + " x ||x||"},
      {"kinetic trace change", "accuracy.kinetic_trace_change",
       "delta_kinetic_energy", "absolute_kinetic_trace_change",
       "hartree_times_response_amplitude_scale_squared",
       "has no kinetic term"}};
  const auto targets = molresponse_v3::fd_requested_accuracy_targets(rp);
  for (const auto& row : requested) {
    const bool written = targets.contains(row.measurement);
    settings.metrics.push_back(
        {Group::requested_criteria, row.name, row.keyword, row.measurement,
         row.quantity, row.units,
         written ? std::optional<double>(
                       targets.at(row.measurement).at("value").get<double>())
                 : std::nullopt,
         written ? SetBy::user : SetBy::not_set, Relation::less_than_or_equal,
         std::nullopt, "", "", row.engine_instead});
  }
  settings.metrics.push_back({Group::requested_criteria, "energy change", "",
                              "energy_change", "", "hartree", std::nullopt,
                              SetBy::not_set, Relation::less_than_or_equal,
                              std::nullopt, "",
                              "not applicable: a driven response has no local "
                              "energy; see the reference stage",
                              ""});
  settings.iteration_caps = {
      {static_cast<int>(rp.maxiter())},
      stage.written.contains("maxiter") ? SetBy::user : SetBy::engine_default};
  // A3: no save boundary inside a solve yet (DESIGN §1.4).
  settings.checkpoints.fixed
      = "at every solve end (molresponse has no save boundary inside a "
        "solve); every generation kept";
  settings.warnings = cv::unenforced_target_warnings(settings);
  for (auto& warning : cv::incomplete_request_warnings(
           settings, {"charge_density_change", "bsh_residual"}))
    settings.warnings.push_back(std::move(warning));
  return settings;
}

// The end-of-stage outcome of a response stage (B9): for each written-target
// measurement, the largest saved-state value over the channels and
// frequencies of the final rung, so the limiting channel decides the row.
[[nodiscard]] inline std::vector<madness::convergence::MetricOutcome>
response_saved_outcomes(const json& metadata,
                        const std::string& final_protocol) {
  std::map<std::string, madness::convergence::MetricOutcome> largest;
  // Held in a local: iterating items() of the temporary value() returns
  // would dangle.
  const auto states = metadata.value("fd_states", json::object());
  for (const auto& channel : states.items()) {
    if (!channel.value().contains(final_protocol))
      continue;
    for (const auto& point : channel.value().at(final_protocol).items()) {
      const auto measurements = molresponse_v3::fd_accuracy_measurements(
          point.value()
              .value("saved_state_validation", json::object())
              .value("measurements", json::object()));
      for (const auto& item : measurements.items()) {
        const auto& record = item.value();
        auto& outcome = largest[item.key()];
        outcome.measurement = item.key();
        if (record.value("availability", "") != "measured"
            || !record.at("value").is_number())
          continue;
        const double value = record.at("value").get<double>();
        if (!outcome.value || value > *outcome.value)
          outcome.value = value;
        outcome.availability = "measured";
      }
    }
  }
  std::vector<madness::convergence::MetricOutcome> outcomes;
  outcomes.reserve(largest.size());
  for (auto& item : largest)
    outcomes.push_back(std::move(item.second));
  return outcomes;
}

// B5, B6: how an FD state's rung ended. Entries written since 4f74ffe0c carry
// rung_exit's reason as `stop`; older ones are read through a documented
// mapping: best-effort acceptance and an exhausted rung are an iteration
// limit (or a stall), never convergence.
[[nodiscard]] inline std::string response_stop(const json& entry) {
  if (entry.contains("stop"))
    return entry.at("stop").get<std::string>();
  if (entry.value("converged", false) && !entry.value("accepted", false))
    return "converged";
  if (entry.value("stalled", false))
    return "stalled";
  if (entry.value("accepted", false) || entry.value("rung_exhausted", false))
    return "iteration_limit";
  return "in_progress";
}

// B6: today's engine check as planner rows, judged at the final rung: per
// leg, BSH residual and density change each at most absolute + rtol x norm,
// the solver's own gate (fd_leg_within), from the same policy the solve uses.
[[nodiscard]] inline std::vector<madness::planning::EngineLimit>
response_engine_limits(const Stage& stage, double final_thresh) {
  molresponse_v3::ConvergencePolicy policy;
  policy.dconv_user = stage.params.get<ResponseParameters>().dconv();
  const auto limits = policy.effective_for_thresh(final_thresh);
  return {{"bsh_residual", limits.bsh_residual, limits.rtol, true},
          {"density_change", limits.density_residual, limits.rtol, true}};
}

// B6: a state record (state_records) as the planner reads it. The written
// targets are judged on the saved state's measurements, under the names the
// targets use (fd_accuracy_measurements).
[[nodiscard]] inline madness::planning::SavedState
response_saved_state(const json& record) {
  madness::planning::SavedState state;
  state.identity = record.at("identity").get<std::string>();
  state.reference_generation
      = record.at("reference_generation").get<std::string>();
  state.perturbation = record.at("perturbation").get<std::string>();
  state.frequency = record.at("frequency").get<double>();
  state.rung = {record.at("rung").at("thresh").get<double>(),
                record.at("rung").at("k").get<int>()};
  state.iterations_at_rung = record.at("iterations_at_rung").get<int>();
  state.stop = madness::convergence::stop_reason_from_string(
      record.at("stop").get<std::string>());
  state.diverged = record.at("diverged").get<bool>();
  state.engine_measurements = record.at("engine_measurements");
  state.saved_measurements = molresponse_v3::fd_accuracy_measurements(
      record.at("saved_measurements"));
  return state;
}

// Scientific policy for the named HF/polarizability frontend. The engine still
// owns state serialization, restart guesses, protocol scheduling and analysis.
class ResponseCheckpoint {
  using Manifest = madness::checkpoint::Manifest;
  madness::World& world_;
  const Stage& stage_;
  std::filesystem::path store_, attempt_;
  std::string run_id_;
  json reference_, request_;
  Manifest latest_;
  bool restored_ = false, loaded_ = false;
  std::size_t sequence_ = 0;
  // B6: iterations recorded per state and rung before this attempt (keyed by
  // the engine's channel|protocol|frequency keys), the rungs the planner
  // resumes, and the rungs this attempt solved, so the total across runs is
  // kept.
  std::map<std::string, int> previous_totals_;
  std::set<std::string> resumed_, solved_;
  // B6, D-016: the planner's decision per requested state, for the engine,
  // and its record, for the manifest and results.json.
  molresponse_v3::StatePlans plans_;
  json plan_records_ = json::array();

  static std::string state_key(const std::string& channel,
                               const std::string& protocol,
                               const std::string& frequency) {
    return channel + "|" + protocol + "|" + frequency;
  }

  // B6, the planner's record of every response state in `metadata`, one per
  // (perturbation, frequency, rung): identity, reference generation, rung,
  // iterations at the rung across runs (recorded, never compared, D-007), how
  // it ended, whether it diverged, the deciding mode, and the raw
  // measurements (the deciding sample and the saved state, each with the norm
  // its relative limit uses).
  [[nodiscard]] json state_records(const json& metadata) const {
    json records = json::array();
    // Held in a local: iterating items() of the temporary value() returns
    // would dangle.
    const auto states = metadata.value("fd_states", json::object());
    for (const auto& channel : states.items())
      for (const auto& rung : channel.value().items()) {
        const auto& protocol = metadata.at("protocols").at(rung.key());
        for (const auto& point : rung.value().items()) {
          const auto& e = point.value();
          const auto key = state_key(channel.key(), rung.key(), point.key());
          // The count of the attempt that last solved this rung.
          const int iter = e.value("iter", 0);
          // A rung solved in this attempt adds to what it had when the
          // planner resumed it, else counts this attempt alone; an untouched
          // one keeps its recorded total.
          const auto previous = previous_totals_.find(key);
          const int before
              = previous != previous_totals_.end() ? previous->second : 0;
          int total = iter;
          if (solved_.contains(key))
            total = (resumed_.contains(key) ? before : 0) + iter;
          else if (previous != previous_totals_.end())
            total = before;
          std::ostringstream frequency;
          frequency << std::setprecision(17) << e.at("freq").get<double>();
          records.push_back(
              {{"schema", "madqc.response.state.v1"},
               {"key", key},
               {"identity", channel.key() + "|" + frequency.str()},
               {"perturbation", channel.key()},
               {"frequency", e.at("freq")},
               {"reference_generation", reference_.value("generation", "")},
               {"rung",
                {{"thresh", protocol.at("thresh")}, {"k", protocol.at("k")}}},
               {"iterations_at_rung", total},
               {"stop", response_stop(e)},
               {"diverged", e.value("diverged", false)},
               {"decided_by", e.value("checkpoint_decided_by", "engine_check")},
               {"engine_measurements",
                e.value("convergence_measurements", json::object())},
               {"saved_measurements",
                e.value("saved_state_validation", json::object())
                    .value("measurements", json::object())}});
        }
      }
    return records;
  }

  static void require(bool ok, const std::string& why) {
    if (!ok) throw std::runtime_error("response checkpoint: " + why);
  }

  // Derive the required archive inventory from state definitions, rather than
  // trusting whatever subset happens to be listed in the storage manifest.
  static std::vector<std::string> members(const json& metadata) {
    std::vector<std::string> names;
    require(metadata.at("excited_states").empty() && metadata.at("vbc_states").empty(),
            "unsupported response bundle kind");
    for (const auto& channel : metadata.at("fd_states").items()) {
      require(channel.key() == "dipole_x" || channel.key() == "dipole_y" ||
              channel.key() == "dipole_z", "unsupported channel");
      for (const auto& rung : channel.value().items()) {
        const auto& protocol = metadata.at("protocols").at(rung.key());
        require(molresponse_v3::protocol_key(protocol.at("thresh"), protocol.at("k")) == rung.key(),
                "inconsistent protocol definition");
        for (const auto& point : rung.value().items()) {
          const auto& e = point.value();
          // Named continuation must not silently turn a failed state into a
          // fresh solve. This check runs both before loading and publication,
          // so a diverged solve cannot replace the previous usable generation.
          require(!e.value("diverged", false), "invalid response state: diverged");
          // Older executor verdicts could certify one update without a density
          // comparison. Do not trust that flag on reuse, or certify a stall.
          // Explicit best-effort acceptance remains a separate engine outcome.
          if (e.value("converged", false) && !e.value("accepted", false))
            require(e.value("iter", 0) >= 2 && !e.value("stalled", false),
                    "invalid response state: ineligible strict convergence");
          auto validate_bsh = [&](const json& metric) {
            if (metric.contains("per_function"))
              require(madness::valid_function_residual_record(metric.at("per_function")),
                      "invalid response state: per-function BSH measurement");
          };
          auto validate_kinetic = [&](const json& measurements) {
            if (measurements.contains("bsh_residual")) validate_bsh(measurements.at("bsh_residual"));
            if (measurements.contains("kinetic_energy"))
              require(molresponse_v3::valid_function_kinetic_record(measurements.at("kinetic_energy")),
                      "invalid response state: per-function kinetic measurement");
          };
          if (e.contains("convergence_measurements")) validate_kinetic(e.at("convergence_measurements"));
          if (e.contains("saved_state_validation")) validate_kinetic(e.at("saved_state_validation").at("measurements"));
          if (e.contains("function_bsh_history")) {
            const auto& history = e.at("function_bsh_history");
            require(history.is_array(), "invalid response state: per-function BSH history");
            for (const auto& sample : history) {
              require(sample.contains("per_function") && sample.contains("value") &&
                      sample.at("value").is_number() && sample.value("availability", "") == "measured",
                      "invalid response state: per-function BSH history");
              const double value = sample.at("value").get<double>();
              require(std::isfinite(value) && value >= 0.0,
                      "invalid response state: per-function BSH history aggregate");
              validate_bsh(sample);
            }
          }
          if (e.contains("function_kinetic_history"))
            for (const auto& sample : e.at("function_kinetic_history"))
              require(molresponse_v3::valid_function_kinetic_record(sample),
                      "invalid response state: per-function kinetic history");
          if (e.contains("saved_state_validation"))
            for (const auto& metric : e.at("saved_state_validation").at("measurements").items())
              require(metric.value().value("availability", "") != "invalid",
                      "invalid response state: nonfinite saved-state measurement");
          if (e.contains("convergence_measurements"))
            for (const auto& metric : e.at("convergence_measurements").items())
              require(metric.value().value("availability", "") != "invalid",
                      "invalid response state: nonfinite measurement");
          // Bind progress evidence to the request that produced this attempt,
          // including historical/pending entries retained by later requests.
          const auto origins = metadata.value("checkpoint_requests", json::object());
          const auto history = e.value("function_kinetic_history", json::array());
          const std::string origin = history.empty() ? "" :
              history.at(0).at("sample").value("attempt_directory", "");
          require(!e.contains("checkpoint_progress") || origins.contains(origin),
                  "invalid response state: missing originating progress request");
          if (origins.contains(origin)) {
            const auto& request = origins.at(origin);
            json progress = request.value("progress", json::object());
            if (!progress.empty()) progress["protocol"] = request.at("protocol");
            require(progress.empty() ? !e.contains("checkpoint_progress") :
                        e.value("checkpoint_progress", json::object()) == progress,
                    "invalid response state: originating progress request mismatch");
          }
          require(molresponse_v3::valid_fd_explicit_accuracy(e, protocol.at("thresh"), protocol.at("k")),
                  "invalid response state: explicit accuracy evidence");
          require(molresponse_v3::valid_metric_diagnostics(e),
                  "invalid response state: metric comparison");
          const double frequency = e.at("freq");
          require(std::isfinite(frequency) &&
                  molresponse_v3::ResponseMetadata::freq_key(frequency) == point.key(),
                  "inconsistent frequency definition");
          const auto name = molresponse_v3::response_filename(channel.key(), rung.key(), frequency);
          require(e.at("archive") == name && e.at("backend") == "native" &&
                  e.at("shell") == "closed_shell" &&
                  e.at("type") == (std::abs(frequency) < 1e-12 ? "static" : "full"),
                  "inconsistent archive/channel/model definition");
          require(e.at("checkpoint_dconv").get<double>() > 0.0,
                  "missing response convergence target");
          names.push_back(name);
        }
      }
    }
    std::sort(names.begin(), names.end());
    require(std::adjacent_find(names.begin(), names.end()) == names.end(), "duplicate state member");
    return names;
  }

  // `loaded`: the state keys the planner loaded (D-016); their verdict is the
  // planner's, so only their identity is checked here.
  static bool request_complete(const json& metadata, const json& request,
                               const std::set<std::string>& loaded = {}) {
    const double threshold = request.at("protocol").back();
    const auto key = molresponse_v3::protocol_key(threshold, protocol_k(threshold));
    if (!metadata.at("protocols").contains(key) ||
        metadata.at("protocols").at(key).at("thresh").get<double>() != threshold)
      return false;
    std::set<char> axes;
    for (const auto& c : request.at("components"))
      for (char axis : c.get<std::string>()) axes.insert(axis);
    for (double frequency : request.at("frequencies")) {
      const auto fkey = molresponse_v3::ResponseMetadata::freq_key(frequency);
      for (char axis : axes) {
        const std::string channel = std::string("dipole_") + axis;
        const auto& states = metadata.at("fd_states");
        if (!states.contains(channel) || !states.at(channel).contains(key) ||
            !states.at(channel).at(key).contains(fkey)) return false;
        const auto& e = states.at(channel).at(key).at(fkey);
        if (loaded.contains(state_key(channel, key, fkey))) {
          if (e.at("freq").get<double>() != frequency
              || e.value("diverged", false))
            return false;
          continue;
        }
        // Written targets block completion only when they decide (B8).
        if (request.value("decided_by", "") == "requested_criteria"
            && request.contains("accuracy_targets")) {
          const auto measurements = e.value("saved_state_validation", json::object())
              .value("measurements", json::object());
          if (madness::evaluate_explicit_criteria(molresponse_v3::fd_accuracy_measurements(measurements),
              request.at("accuracy_targets")).at("status") != "satisfied") return false;
        }
        if (e.at("freq").get<double>() != frequency ||
            !e.at("converged").get<bool>() || e.value("accepted", false) ||
            e.value("diverged", false) ||
            e.at("checkpoint_dconv").get<double>() > request.at("dconv").get<double>()) return false;
      }
    }
    return true;
  }

  void validate(const Manifest& m) const {
    require(m.stage == stage_.name && m.method == "TDHF" &&
            m.physical == json{{"schema", "madqc.response.v1"}, {"shell", "closed_shell"}},
            "stage or response model mismatch");
    if (m.numerical.at("request").contains("progress")) {
      const auto origins = m.numerical.at("engine_metadata").value("checkpoint_requests", json::object());
      const auto attempt = m.numerical.value("attempt_directory", "");
      require(origins.contains(attempt) && origins.at(attempt) == m.numerical.at("request"),
              "invalid response state: missing or inconsistent originating request registry");
    }
    auto expected = members(m.numerical.at("engine_metadata"));
    std::vector<std::string> actual;
    for (const auto& a : m.archives) {
      require(a.role == "response" && a.indices.empty(), "invalid response archive role");
      actual.push_back(a.basename);
    }
    std::sort(actual.begin(), actual.end());
    require(!expected.empty() && expected == actual, "required response inventory is missing or inconsistent");
    for (const auto& p : m.properties)
      require(p.at("generation") == m.generation && p.at("name") == "polarizability" &&
              p.at("value") == m.numerical.at("engine_metadata").at("properties"),
              "property/state record is inconsistent");
    if (m.stop_reason == "converged")
      require(m.properties.size() == 1
                  && request_complete(
                      m.numerical.at("engine_metadata"),
                      m.numerical.at("request"),
                      loaded_states(m.numerical.value("plans", json::array()))),
              "completed property lacks converged prerequisites");
  }

  // The state keys a generation's plan records say were loaded and met.
  static std::set<std::string> loaded_states(const json& plan_records) {
    std::set<std::string> loaded;
    for (const auto& record : plan_records)
      if (record.contains("loaded_state"))
        loaded.insert(record.at("loaded_state").get<std::string>());
    return loaded;
  }

  // B6: one plan per requested state, from the generation's state records.
  // `own` is the state's finest saved rung the engine can use. Left out: a
  // rung whose engine key equals a requested rung's but whose threshold
  // differs (keys round it; the engine holds one definition per key), and a
  // rung tighter than the final one (R3): the analysis reads only states at
  // or coarser than the final rung, so such a state is solved again.
  void plan_states(const json& records) {
    std::map<std::string, double> requested_rungs;
    auto request = ladder_request(requested_rungs);
    std::set<int> axes;
    for (const auto& component : stage_.components)
      for (char axis : component)
        axes.insert(axis - 'x');
    for (double frequency : stage_.frequencies)
      for (int axis : axes) {
        const auto perturbation = molresponse_v3::Perturbation::dipole(axis);
        std::ostringstream identity;
        identity << perturbation.description() << '|' << std::setprecision(17)
                 << frequency;
        request.identity = identity.str();
        const auto own = own_state(records, request.identity, requested_rungs,
                                   request.ladder.back());
        adopt(madness::planning::plan_ladder(own, std::nullopt, request), own,
              request, perturbation, frequency);
      }
  }

  // The request every state of this stage is planned against, less its
  // identity; `requested_rungs` maps each requested rung's engine key to its
  // threshold.
  [[nodiscard]] madness::planning::LadderRequest
  ladder_request(std::map<std::string, double>& requested_rungs) const {
    madness::planning::LadderRequest request;
    for (double threshold :
         stage_.params.get<CalculationParameters>().protocol()) {
      request.ladder.push_back({threshold, protocol_k(threshold)});
      requested_rungs.emplace(
          molresponse_v3::protocol_key(threshold, protocol_k(threshold)),
          threshold);
    }
    const double final_thresh = request.ladder.back().thresh;
    request.reference_generation = reference_.value("generation", "");
    request.settings = response_convergence_settings(stage_, final_thresh);
    request.engine_limits = response_engine_limits(stage_, final_thresh);
    request.forced_iterations
        = stage_.params.get<ResponseParameters>().get<int>(
              "diagnostic.iterations_per_rung")
          > 0;
    return request;
  }

  // The record a state is planned from: its finest saved rung the engine can
  // use (see plan_states).
  [[nodiscard]] static std::optional<madness::planning::SavedState>
  own_state(const json& records, const std::string& identity,
            const std::map<std::string, double>& requested_rungs,
            const madness::planning::Rung& final_rung) {
    std::optional<madness::planning::SavedState> own;
    for (const auto& record : records) {
      if (record.at("identity") != identity)
        continue;
      auto saved = response_saved_state(record);
      const auto rung = requested_rungs.find(
          molresponse_v3::protocol_key(saved.rung.thresh, saved.rung.k));
      const bool aliased
          = rung != requested_rungs.end() && rung->second != saved.rung.thresh;
      const bool tighter_than_final
          = saved.rung.thresh < final_rung.thresh
            && !madness::scf_continuation::same_rung(saved.rung, final_rung);
      if (aliased || tighter_than_final)
        continue;
      if (!own || saved.rung.thresh < own->rung.thresh)
        own = std::move(saved);
    }
    return own;
  }

  // One state's plan: the engine's input (plans_), the iteration bookkeeping
  // (resumed_) and the record for results.json and the manifest.
  void adopt(const madness::planning::LadderPlan& plan,
             const std::optional<madness::planning::SavedState>& own,
             const madness::planning::LadderRequest& request,
             const molresponse_v3::Perturbation& perturbation,
             double frequency) {
    namespace planning = madness::planning;
    using Kind = molresponse_v3::StatePlan::Kind;
    require(plan.action != planning::Action::seed,
            "no response seed was offered");
    const std::map<planning::Action, const char*> words
        = {{planning::Action::load, "loaded"},
           {planning::Action::resume, "resumed"},
           {planning::Action::fresh, "fresh"}};
    const auto channel = perturbation.description();
    const auto fkey = molresponse_v3::ResponseMetadata::freq_key(frequency);
    const auto& start = request.ladder.at(plan.start_rung);
    molresponse_v3::StatePlan engine_plan;
    engine_plan.start_thresh = start.thresh;
    json record = {{"state", request.identity},
                   {"perturbation", channel},
                   {"frequency", frequency},
                   {"action", words.at(plan.action)},
                   {"rung",
                    {{"index", plan.start_rung + 1},
                     {"thresh", start.thresh},
                     {"k", start.k}}},
                   {"iterations_done_at_rung", plan.iterations_done},
                   {"why", plan.why},
                   {"unmet", plan.unmet},
                   {"unmet_not_enforceable", plan.unmet_not_enforceable},
                   {"summary", planning::describe(plan, request.ladder)}};
    if (plan.action == planning::Action::load) {
      engine_plan.kind = Kind::load;
      engine_plan.loaded_key
          = molresponse_v3::protocol_key(own->rung.thresh, own->rung.k);
      // R2 loads a stalled final rung with the request unmet: not met.
      engine_plan.met = plan.unmet.empty();
      if (engine_plan.met)
        record["loaded_state"]
            = state_key(channel, engine_plan.loaded_key, fkey);
    } else if (plan.action == planning::Action::resume) {
      engine_plan.kind = Kind::resume;
      if (madness::scf_continuation::same_rung(own->rung, start))
        resumed_.insert(state_key(
            channel, molresponse_v3::protocol_key(start.thresh, start.k),
            fkey));
    } else {
      engine_plan.kind = Kind::fresh;
    }
    plans_[molresponse_v3::fd_node_id(perturbation, frequency)] = engine_plan;
    plan_records_.push_back(std::move(record));
  }

public:
  std::size_t iterations = 0;
  json work_items = json::array();
  std::string work = "solved";
  bool reference_changed = false;

  ResponseCheckpoint(madness::World& world, const Stage& stage,
                     const std::filesystem::path& stage_directory,
                     const std::filesystem::path& attempt, std::string run_id,
                     json reference)
      : world_(world), stage_(stage), store_(stage_directory / "checkpoints"),
        attempt_(attempt), run_id_(std::move(run_id)), reference_(std::move(reference)) {
    const auto& cp = stage.params.get<CalculationParameters>();
    request_ = {{"protocol", cp.protocol()}, {"frequencies", stage.frequencies},
                {"components", stage.components},
                {"dconv", stage.params.get<ResponseParameters>().dconv()}};
    const auto targets = molresponse_v3::fd_requested_accuracy_targets(stage.params.get<ResponseParameters>());
    if (!targets.empty()) request_["accuracy_targets"] = targets;
    // Who decides convergence (B8); reuse and property eligibility follow it.
    request_["decided_by"] = stage.params.get<ResponseParameters>().get<bool>(
                                 madness::convergence::use_default_criteria_key)
                                 ? "engine_check"
                                 : "requested_criteria";
    const auto progress = molresponse_v3::fd_requested_progress(stage.params.get<ResponseParameters>());
    if (!progress.empty()) request_["progress"] = progress;
    const int diagnostic = stage.params.get<ResponseParameters>().get<int>("diagnostic.iterations_per_rung");
    if (diagnostic > 0) request_["diagnostic_iterations_per_rung"] = diagnostic;
  }

  [[nodiscard]] bool targets_decide() const {
    return request_.value("decided_by", "") == "requested_criteria";
  }

  // Return true only when the requested analysis itself is already complete.
  // Otherwise restore prerequisites privately and let the engine schedule work.
  bool prepare(bool continuing) {
    namespace checkpoint = madness::checkpoint;
    if (!continuing) return false;
    bool exists = false;
    checkpoint::detail::on_root(world_, [&] {
      exists = std::filesystem::exists(store_ / "latest");
      return std::string{};
    });
    world_.gop.broadcast(exists, 0);
    if (!exists) return false; // a newly added stage, or no completed first save
    latest_ = checkpoint::read_published(world_, store_);
    validate(latest_);
    reference_changed = latest_.reference != reference_;
    if (reference_changed) return false; // no old descendants enter this attempt
    restored_ = true;
    loaded_ = stage_.params.get<ResponseParameters>().get<int>("diagnostic.iterations_per_rung") == 0 &&
              latest_.stop_reason == "converged" &&
              latest_.numerical.at("request") == request_ && !latest_.properties.empty();
    if (loaded_) { work = "loaded"; return true; }

    auto metadata = latest_.numerical.at("engine_metadata");
    // B6: the generation's state records (or, for one written before them,
    // records derived from its entries): totals so far, and the plans.
    auto records = latest_.numerical.value("states", json::array());
    if (records.empty())
      records = state_records(metadata);
    for (const auto& record : records)
      previous_totals_[record.at("key").get<std::string>()]
          = record.at("iterations_at_rung").get<int>();
    plan_states(records);
    // The engine's state entries are copied as saved (D7, D-016); the plans
    // above steer it. Analysis is derived: old property rows are left out so
    // they cannot certify this request.
    metadata["properties"] = json::object();
    checkpoint::detail::on_root(world_, [&] {
      for (const auto& a : latest_.archives)
        for (const auto& shard : a.shards)
          std::filesystem::copy_file(store_ / latest_.generation / shard.path,
                                     attempt_ / shard.path);
      checkpoint::detail::write_exclusive(attempt_ / "response_metadata.json", metadata.dump(2));
      return std::string{};
    });
    return false;
  }

  const Manifest& latest() const { return latest_; }
  [[nodiscard]] const molresponse_v3::StatePlans& state_plans() const {
    return plans_;
  }
  /// What --continue did per state and why (B6), for results.json.
  [[nodiscard]] json continuation() const {
    return {{"planner", "plan_ladder"}, {"states", plan_records_}};
  }
  json metadata() const { return latest_.numerical.at("engine_metadata"); }
  json properties() const { return metadata().at("properties"); }

  // Called collectively after native FD save and metadata update have finished.
  void observe(madness::World&, const json& event) {
    iterations += event.at("iterations").get<std::size_t>();
    solved_.insert(state_key(event.at("perturbation").get<std::string>(),
                             event.at("protocol").get<std::string>(),
                             molresponse_v3::ResponseMetadata::freq_key(
                                 event.at("frequency").get<double>())));
    auto item = event;
    item["work"] = event.at("seed") == "fd_restart" ? "resumed" : "solved";
    work_items.push_back(item);
    const auto metadata = read_metadata();
    publish(metadata, false, false, event);
  }

  void finish(const json& root_metadata, bool converged) {
    if (loaded_) return;
    std::string text = world_.rank() == 0 ? root_metadata.dump() : std::string{};
    world_.gop.broadcast_serializable(text, 0);
    work = iterations == 0 ? "analyzed" : restored_ ? "resumed" : "solved";
    publish(json::parse(text), true, converged, {{"phase", "property_analysis"}});
  }

private:
  json read_metadata() const {
    const auto text = madness::checkpoint::detail::on_root(world_, [&] {
      return madness::checkpoint::detail::read_json(attempt_ / "response_metadata.json").dump();
    });
    return json::parse(text);
  }

  void publish(json metadata, bool analyzed, bool converged, const json& position) {
    namespace checkpoint = madness::checkpoint;
    if (!analyzed) metadata["properties"] = json::object();
    // Native saves replace entries but preserve their sample attempt identity.
    // Keep the corresponding request registry across intermediate publications.
    auto origins = metadata.value("checkpoint_requests", json::object());
    if (!latest_.generation.empty() && latest_.reference == reference_) {
      const auto previous = latest_.numerical.at("engine_metadata").value("checkpoint_requests", json::object());
      for (const auto& origin : previous.items()) origins[origin.key()] = origin.value();
    }
    origins[attempt_.string()] = request_;
    metadata["checkpoint_requests"] = std::move(origins);
    Manifest m;
    m.stage = stage_.name;
    m.method = "TDHF";
    m.generation = "g-" + stage_.name + "-" + run_id_ + "-" + std::to_string(++sequence_);
    m.physical = {{"schema", "madqc.response.v1"}, {"shell", "closed_shell"}};
    m.reference = reference_;
    m.numerical = {{"request", request_},
                   {"attempt_directory", attempt_.string()},
                   {"engine_metadata", metadata},
                   {"states", state_records(metadata)},
                   {"plans", plan_records_}};
    m.targets = {{"dconv", request_.at("dconv")}};
    m.position = position;
    m.position["iterations_this_run"] = iterations;
    if (!latest_.generation.empty()) m.position["previous_generation"] = latest_.generation;
    m.stop_reason = converged ? "converged" : analyzed ? "unfinished" : "in_progress";
    // Keep analyzed values as diagnostics in engine_metadata, but an unmet
    // explicit request cannot publish them as eligible stage properties.
    if (analyzed && (converged || !targets_decide()))
      m.properties.push_back({{"name", "polarizability"}, {"generation", m.generation},
                              {"converged", converged}, {"value", metadata.at("properties")}});
    const auto directory = store_ / m.generation;
    const auto names = members(metadata);
    if (converged)
      require(
          request_complete(metadata, request_, loaded_states(plan_records_)),
          "analysis lacks converged prerequisites");
    // Writers have already closed collectively. Copy the complete prerequisite
    // set; immutable generations never share writable files with an attempt.
    checkpoint::detail::on_root(world_, [&] {
      std::filesystem::create_directories(store_);
      require(std::filesystem::create_directory(directory), "generation already exists");
      for (const auto& name : names) {
        const int writers = checkpoint::detail::header(attempt_ / (name + ".00000"), true);
        const auto inventory = checkpoint::inspect_archive(attempt_, "response", name, {}, writers);
        for (const auto& shard : inventory.shards)
          std::filesystem::copy_file(attempt_ / shard.path, directory / shard.path);
        m.archives.push_back(checkpoint::inspect_archive(directory, "response", name, {}, writers));
      }
      return std::string{};
    });
    checkpoint::publish(world_, store_, m);
    latest_ = checkpoint::read_published(world_, store_);
  }
};
} // namespace named_stages
