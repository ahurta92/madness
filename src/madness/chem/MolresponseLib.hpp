/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680


  $Id$
 */
#pragma once

#include "ResponseParameters.hpp"
#include "madness_exception.h"
#include <apps/molresponse_v2/DerivedStatePlanner.hpp>
#include <apps/molresponse_v2/ExcitedResponse.hpp>
#include <apps/molresponse_v2/FrequencyLoop.hpp>
#include <apps/molresponse_v2/GroundStateData.hpp>
#include <apps/molresponse_v2/PropertyManager.hpp>
#include <apps/molresponse_v2/ResponseDebugLogger.hpp>
#include <apps/molresponse_v2/ResponseManager.hpp>
#include <apps/molresponse_v2/ResponseRecord.hpp>
#include <apps/molresponse_v2/StateGenerator.hpp>
#include <apps/molresponse_v2/StateParallelPlanner.hpp>
#include <apps/molresponse_v2/VBCMacrotask.hpp>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <chrono>
#include <cctype>
#include <numeric>
#include <optional>
#include <memory>
#include <sstream>
#include <thread>
#include <utility>
#include <vector>
#include <madness/chem/InputWriter.hpp>
#include <madness/chem/ParameterManager.hpp>
#include <madness/chem/Results.h>
#include <madness/mra/macrotaskq.h>

/*
Developer Overview
- Orchestration file for molresponse, structured as three stages:
  1) bootstrap + planning, 2) state solves, 3) property assembly.
- Terminology used in this file:
  - perturbation channel: one perturbation descriptor (e.g. Dipole_x)
  - frequency series: ordered frequencies for one perturbation channel
  - channel point: one (channel, frequency, protocol) solve target
- Stage 1: read ground/checkpoint context, generate linear states, build derived
  request plan, and build state-parallel ownership/execution plan.
- Stage 2: solve linear states over protocol thresholds; fresh runs start
  state-oriented, then later thresholds can fan out by frequency-point.
  Restart runs with complete protocol-0 saved data may promote fanout to ti==0.
  When subgroup mode is enabled, work runs in macrotask subworlds with shard
  metadata/log merge and synchronization barriers between protocol steps.
- Stage 2c: evaluate derived-state dependency gate and execute ready requests
  before property assembly (serial or subgroup lanes).
- Stage 3: compute requested properties (alpha/beta/raman) with optional
  property-group execution and broadcast of assembled outputs.
*/

/// Top-level response workflow orchestrator.
///
/// Responsibilities:
/// - Build runtime context from checkpoint + response inputs.
/// - Plan linear, derived, and excited-state work.
/// - Execute state solve stages (serial or subgroup modes).
/// - Assemble requested properties and return structured JSON payloads.
///
/// Logical sections (in order of appearance):
///   1. Context types          — Results, GroundContext, ExcitedStateBundlePlan,
///                               PlannedStates, SolvedStates, PropertyStageOutput
///   2. State persistence      — JsonStateSolvePersistence
///   3. JSON file I/O          — write_json_file, read_json_file_or_object,
///                               broadcast_json*, with_subworld,
///                               merge_state_metadata_json, merge_debug_log_json,
///                               point_ready_in_metadata,
///                               point_marked_for_frequency_removal_in_metadata
///   4. Manifest/claim helpers — sanitize_manifest_token,
///                               derived_request_{manifest_scope,done_file,
///                               claim_file,done_record_exists,
///                               write_done_record}
///   5. Log capture            — FilteredLineStreambuf, ScopedRankLogRedirect,
///                               group_{shard,console,derived_timing}_file
///   6. Schedule context types — ProtocolExecutionPolicy,
///                               RuntimePointOwnershipPolicy,
///                               StateSolveScheduleContext,
///                               PendingPointWorkItem, PendingProtocolManifest
///                               and associated helper functions
///   7. Schedule context builder — build_state_solve_schedule_context
///                                 and its sub-steps
///   8. Linear state solve     — prepare_protocol_context, log_pending_manifest,
///                               execute_manifest_work, cached_or_built_manifest,
///                               execute_serial_state_solve,
///                               execute_subgroup_state_solve
///   9. Excited-state bundle   — ExcitedExecutionResult,
///                               ensure_excited_protocol_placeholder_node,
///                               build_excited_root_manifest,
///                               execute_excited_state_bundle_stage
///  10. Derived state          — DerivedRequestTiming, run_derived_request,
///                               DerivedExecutionResult,
///                               execute_derived_state_requests
///  11. Property stage         — PropertyContext, PropertyType,
///                               parse_property_name, compute_polarizability,
///                               compute_hyperpolarizability, print_raman_table,
///                               compute_raman, compute_requested_properties,
///                               compute_requested_properties_with_property_group
///  12. Top-level orchestration — solve_all_states, run_response
struct molresponse_lib {
  // ============================================================================
  // SECTION 1: Context types
  //   Public result struct + private planning and solve-stage context bundles.
  //   These are value types passed by move across the Stage 1 → 2 → 3 pipeline.
  // ============================================================================

  /// Structured output returned to workflow drivers (`madqc --wf=response`).
  struct Results {
    nlohmann::json metadata;
    nlohmann::json properties;
    nlohmann::json vibrational_analysis;
    nlohmann::json raman_spectra;
    nlohmann::json debug_log;
  };

  static constexpr const char *label() { return "molresponse"; }

private:
  struct GroundContext {
    // Molecule read from checkpoint; reused across all response stages.
    Molecule molecule;
    // Ground-state object used to build and solve response states.
    GroundStateData ground;
    // Runtime manager configured from response input.
    ResponseManager response_manager;
    // Archive and Fock filenames resolved for this run directory.
    std::string archive_file;
    std::string fock_json_file;
  };

  struct ExcitedStateBundlePlan {
    // Enables the excited-state bundle protocol stage in Stage 2.
    bool enabled = false;
    // Target number of excited states in one coupled bundle solve.
    size_t num_states = 1;
    // Optional TDA mode for the future excited-state solver stage.
    bool tda = false;
    // Guess stage iteration budget.
    size_t guess_max_iter = 5;
    // Main excited-state iteration budget.
    size_t maxiter = 20;
    // Subspace size for iterative diagonalization.
    size_t maxsub = 8;
    // Reserved owner lane for the future excited-state bundle solve.
    size_t owner_group = 0;
    // Protocol thresholds that will index excited-bundle restart metadata.
    std::vector<double> protocols;

    [[nodiscard]] nlohmann::json to_json() const {
      nlohmann::json protocol_keys = nlohmann::json::array();
      for (const auto threshold : protocols) {
        protocol_keys.push_back(ResponseRecord2::protocol_key(threshold));
      }
      return {{"enabled", enabled},
              {"num_states", num_states},
              {"tda", tda},
              {"guess_max_iter", guess_max_iter},
              {"maxiter", maxiter},
              {"maxsub", maxsub},
              {"owner_group", owner_group},
              {"protocols", std::move(protocol_keys)}};
    }
  };

  struct PlannedStates {
    // Linear states generated directly from requested perturbations.
    GeneratedStateData generated_states;
    // Derived-state requests (currently VBC-driven scaffolding/execution).
    DerivedStatePlan derived_state_plan;
    // Excited-state bundle execution plan (protocol-indexed).
    ExcitedStateBundlePlan excited_state_bundle_plan;
    // Ownership and subgroup execution plan for Stage 2.
    StateParallelPlan state_parallel_plan;
  };

  struct SolvedStates {
    // Carries full planning context forward into property stage.
    PlannedStates planned_states;
    // Merged response metadata after linear + derived execution.
    nlohmann::json metadata;
    // Combined debug log (including subgroup shards when used).
    nlohmann::json debug_log;
  };

  struct PropertyStageOutput {
    // Aggregated property JSON (alpha/beta/raman blocks).
    nlohmann::json properties;
    // Vibrational artifacts needed for Raman post-processing/output.
    VibrationalResults vibrational_analysis;
    RamanResults raman_spectra;
  };

  // ============================================================================
  // SECTION 2: State persistence
  //   JsonStateSolvePersistence implements StateSolvePersistence by delegating
  //   to ResponseRecord2 (metadata) and ResponseDebugLogger (debug log).  It is
  //   instantiated once per execution path (serial or per-subgroup) and drives
  //   the FrequencyLoop convergence queries.
  // ============================================================================

  /// JSON-backed implementation of `StateSolvePersistence`.
  ///
  /// Wraps `ResponseRecord2` and `ResponseDebugLogger` so the frequency solver
  /// can remain backend-agnostic while orchestration keeps restart metadata,
  /// diagnostics, and debug logs in stable files.
  class JsonStateSolvePersistence final : public StateSolvePersistence {
  public:
    using ProgressPollFn = std::function<void(bool)>;

    JsonStateSolvePersistence(World &world, const std::string &meta_file,
                              const std::string &debug_file,
                              nlohmann::json baseline_metadata =
                                  nlohmann::json::object(),
                              bool force_retry_removed_frequencies = false,
                              ProgressPollFn progress_poll_fn =
                                  ProgressPollFn{})
        : response_record_(world, meta_file), debug_logger_(debug_file),
          baseline_metadata_(std::move(baseline_metadata)),
          force_retry_removed_frequencies_(force_retry_removed_frequencies),
          progress_poll_fn_(std::move(progress_poll_fn)) {}

    void
    initialize_states(const std::vector<LinearResponseDescriptor> &states) {
      response_record_.initialize_states(states);
    }

    void initialize_excited_bundle(const ExcitedStateBundlePlan &plan) {
      response_record_.initialize_excited_bundle(
          plan.enabled, plan.num_states, plan.tda, plan.guess_max_iter,
          plan.maxiter, plan.maxsub, plan.owner_group, plan.protocols);
    }

    void print_summary() const { response_record_.print_summary(); }

    [[nodiscard]] bool is_saved(const LinearResponsePoint &pt) const override {
      return response_record_.is_saved(pt.perturbationDescription(),
                                       pt.threshold(), pt.frequency()) ||
             point_ready_in_baseline(pt, /*require_converged=*/false);
    }

    [[nodiscard]] bool
    is_converged(const LinearResponsePoint &pt) const override {
      return response_record_.is_converged(pt.perturbationDescription(),
                                           pt.threshold(), pt.frequency()) ||
             point_ready_in_baseline(pt, /*require_converged=*/true);
    }

    [[nodiscard]] bool
    is_removed_from_frequency_set(const LinearResponsePoint &pt) const override {
      return response_record_.is_marked_for_frequency_removal(
                 pt.perturbationDescription(), pt.threshold(), pt.frequency()) ||
             point_marked_for_frequency_removal_in_metadata(baseline_metadata_,
                                                            pt);
    }

    [[nodiscard]] bool force_retry_removed_frequencies() const override {
      return force_retry_removed_frequencies_;
    }

    void record_status(const LinearResponsePoint &pt, bool c) override {
      response_record_.record_status(pt, c);
    }

    void record_solver_diagnostics(const LinearResponsePoint &pt,
                                   const ResponseSolveDiagnostics &d,
                                   bool used_fallback_retry) override {
      response_record_.record_solver_diagnostics(
          pt, d.converged, d.iterations_performed, d.final_residual_norm,
          d.final_density_change, used_fallback_retry, d.final_alpha,
          d.max_consecutive_negative_alpha, d.reached_iteration_limit,
          d.remove_from_frequency_set, d.residual_remove_cutoff,
          d.failure_reason);
    }

    void record_timing(const LinearResponsePoint &pt, double wall_seconds,
                       double cpu_seconds) override {
      response_record_.record_timing(pt, wall_seconds, cpu_seconds);
    }

    void record_restart_provenance(
        const LinearResponsePoint &pt, const std::string &source_kind,
        bool loaded_from_disk, bool promoted_from_static,
        const std::optional<double> &source_protocol,
        const std::optional<double> &source_frequency) override {
      response_record_.record_restart_provenance(
          pt, source_kind, loaded_from_disk, promoted_from_static,
          source_protocol, source_frequency);
    }

    ResponseDebugLogger &logger() override { return debug_logger_; }

    void flush_debug_log(World &world) override {
      if (debug_logger_.enabled() && world.rank() == 0) {
        debug_logger_.write_to_disk();
      }
    }

    void maybe_poll_progress(World &world, bool force = false) override {
      (void)world;
      if (progress_poll_fn_) {
        progress_poll_fn_(force);
      }
    }

    [[nodiscard]] nlohmann::json metadata_json() const {
      return response_record_.to_json();
    }

    [[nodiscard]] nlohmann::json debug_log_json() const {
      return debug_logger_.to_json();
    }

  private:
    [[nodiscard]] bool
    point_ready_in_baseline(const LinearResponsePoint &pt,
                            bool require_converged) const {
      return point_ready_in_metadata(baseline_metadata_, pt,
                                     /*require_saved=*/true,
                                     require_converged);
    }

    ResponseRecord2 response_record_;
    ResponseDebugLogger debug_logger_;
    nlohmann::json baseline_metadata_;
    bool force_retry_removed_frequencies_ = false;
    ProgressPollFn progress_poll_fn_;
  };

  // ============================================================================
  // SECTION 3: Subgroup shard filename generators
  //   Produce per-group file paths for metadata shards, console logs, and
  //   derived-request timing shards.  Used by both the linear-solve and
  //   derived-state subgroup execution paths.
  // ============================================================================

  static std::string group_shard_file(const std::string &filename, size_t gid) {
    std::filesystem::path in(filename);
    const auto parent = in.parent_path();
    const auto stem = in.stem().string();
    const auto ext = in.extension().string();
    const auto grouped = ext.empty()
                             ? stem + ".group" + std::to_string(gid)
                             : stem + ".group" + std::to_string(gid) + ext;
    if (parent.empty()) {
      return grouped;
    }
    return (parent / grouped).string();
  }

  static std::string group_console_file(size_t gid) {
    return group_shard_file("response_console.log", gid);
  }

  static std::string group_derived_timing_file(size_t gid) {
    return group_shard_file("derived_request_timings.json", gid);
  }

  // ============================================================================
  // SECTION 4: Derived-state manifest and claim-file helpers
  //   File-based coordination primitives for derived-request idempotency.
  //   sanitize_manifest_token / derived_request_manifest_scope produce stable
  //   path tokens.  The done/claim file pair implements a lightweight
  //   distributed task-claim protocol used in both subgroup and serial paths.
  // ============================================================================

  static std::string sanitize_manifest_token(const std::string &raw) {
    std::string token;
    token.reserve(raw.size());
    for (const char ch : raw) {
      const bool is_alnum =
          std::isalnum(static_cast<unsigned char>(ch)) != 0;
      token.push_back(is_alnum ? ch : '_');
    }
    if (token.empty()) {
      token = "unnamed";
    }
    return token;
  }

  static std::string derived_request_manifest_scope(size_t final_ti,
                                                    double final_thresh) {
    std::ostringstream oss;
    oss << "ti" << final_ti << "_"
        << sanitize_manifest_token(ResponseRecord2::protocol_key(final_thresh));
    return oss.str();
  }

  static std::string derived_request_done_file(size_t final_ti,
                                               double final_thresh,
                                               const DerivedStateRequest &req,
                                               size_t request_index) {
    std::ostringstream filename;
    filename << "req" << request_index << "_"
             << sanitize_manifest_token(req.derived_state_id) << ".done.json";
    const auto scope = derived_request_manifest_scope(final_ti, final_thresh);
    return (std::filesystem::path("derived_request_manifest") / scope /
            filename.str())
        .string();
  }

  static std::string derived_request_claim_file(const std::string &claim_prefix,
                                                size_t request_index) {
    return claim_prefix + ".req" + std::to_string(request_index) + ".lock";
  }

  static bool derived_request_done_record_exists(size_t final_ti,
                                                 double final_thresh,
                                                 const DerivedStateRequest &req,
                                                 size_t request_index) {
    const auto done_file =
        derived_request_done_file(final_ti, final_thresh, req, request_index);
    if (!std::filesystem::exists(done_file)) {
      return false;
    }
    const auto done_json = read_json_file_or_object(done_file);
    return done_json.value("success", false);
  }

  static void write_derived_request_done_record(size_t final_ti,
                                                double final_thresh,
                                                const DerivedStateRequest &req,
                                                size_t request_index,
                                                size_t subgroup_id,
                                                double wall_seconds,
                                                double cpu_seconds) {
    const auto done_file =
        derived_request_done_file(final_ti, final_thresh, req, request_index);
    std::error_code ec;
    const auto done_dir = std::filesystem::path(done_file).parent_path();
    if (!done_dir.empty()) {
      std::filesystem::create_directories(done_dir, ec);
      if (ec) {
        throw std::runtime_error("Failed to create derived manifest directory '" +
                                 done_dir.string() + "': " + ec.message());
      }
    }
    nlohmann::json done_payload = {
        {"success", true},
        {"request_index", request_index},
        {"derived_state_id", req.derived_state_id},
        {"protocol_threshold", final_thresh},
        {"subgroup_id", subgroup_id},
        {"wall_seconds", wall_seconds},
        {"cpu_seconds", cpu_seconds}};
    write_json_file(done_file, done_payload);
  }

  // ============================================================================
  // SECTION 5: Console log capture (stream-redirect utilities)
  //   FilteredLineStreambuf suppresses noisy MADNESS hung-queue warnings when
  //   redirecting subgroup stdout/stderr to per-group shard log files.
  //   ScopedRankLogRedirect is an RAII guard that activates the redirect for
  //   the lifetime of a subgroup execution block.
  // ============================================================================

  class FilteredLineStreambuf final : public std::streambuf {
  public:
    explicit FilteredLineStreambuf(std::streambuf *destination,
                                   std::string dropped_substring)
        : destination_(destination),
          dropped_substring_(std::move(dropped_substring)) {}

  protected:
    int overflow(int ch) override {
      if (ch == traits_type::eof()) {
        flush_buffered_line();
        return traits_type::not_eof(ch);
      }
      buffered_line_.push_back(static_cast<char>(ch));
      if (ch == '\n') {
        flush_buffered_line();
      }
      return ch;
    }

    int sync() override {
      flush_buffered_line();
      if (destination_ != nullptr) {
        return destination_->pubsync();
      }
      return 0;
    }

  private:
    void flush_buffered_line() {
      if (buffered_line_.empty() || destination_ == nullptr) {
        buffered_line_.clear();
        return;
      }
      const bool should_drop =
          !dropped_substring_.empty() &&
          buffered_line_.find(dropped_substring_) != std::string::npos;
      if (!should_drop) {
        destination_->sputn(buffered_line_.data(),
                            static_cast<std::streamsize>(buffered_line_.size()));
      }
      buffered_line_.clear();
    }

    std::streambuf *destination_ = nullptr;
    std::string dropped_substring_;
    std::string buffered_line_;
  };

  class ScopedRankLogRedirect final {
  public:
    ScopedRankLogRedirect(bool enabled, const std::string &filename,
                          const std::string &header = "")
        : enabled_(enabled) {
      if (!enabled_) {
        return;
      }
      sink_.open(filename, std::ios::out | std::ios::app);
      if (!sink_) {
        return;
      }

      old_cout_ = std::cout.rdbuf(sink_.rdbuf());
      filtered_cerr_ = std::make_unique<FilteredLineStreambuf>(
          sink_.rdbuf(), "!!MADNESS: Hung queue?");
      old_cerr_ = std::cerr.rdbuf(filtered_cerr_.get());
      active_ = true;
      if (!header.empty()) {
        std::cout << "\n=== " << header << " ===\n";
      }
    }

    ~ScopedRankLogRedirect() { restore(); }
    ScopedRankLogRedirect(const ScopedRankLogRedirect &) = delete;
    ScopedRankLogRedirect &operator=(const ScopedRankLogRedirect &) = delete;

  private:
    void restore() {
      if (!active_) {
        return;
      }
      std::cout.flush();
      std::cerr.flush();
      std::cout.rdbuf(old_cout_);
      std::cerr.rdbuf(old_cerr_);
      active_ = false;
    }

    bool enabled_ = false;
    bool active_ = false;
    std::ofstream sink_;
    std::streambuf *old_cout_ = nullptr;
    std::streambuf *old_cerr_ = nullptr;
    std::unique_ptr<FilteredLineStreambuf> filtered_cerr_;
  };

  // ============================================================================
  // SECTION 6: JSON file I/O and broadcast utilities
  //   write_json_file / read_json_file_or_object — safe disk read/write with
  //   in-flight-parse guard (returns empty object on malformed content).
  //   broadcast_json_object / broadcast_json — serialize → gop broadcast →
  //   deserialize pattern for distributing metadata across all MPI ranks.
  //   with_subworld — RAII subworld create/pmap-swap/restore template.
  //   merge_state_metadata_json / merge_debug_log_json — or-merge shard JSONs
  //   into the canonical metadata file at global barrier points.
  //   point_ready_in_metadata / point_marked_for_frequency_removal_in_metadata
  //   — stateless key-path lookups into the metadata JSON used by stages 2-3.
  // ============================================================================

  static void write_json_file(const std::string &filename,
                              const nlohmann::json &json_data) {
    std::ofstream out(filename);
    if (!out) {
      throw std::runtime_error("Cannot open " + filename + " for writing");
    }
    out << std::setw(2) << json_data << "\n";
  }

  static nlohmann::json read_json_file_or_object(const std::string &filename) {
    if (!std::filesystem::exists(filename)) {
      return nlohmann::json::object();
    }
    std::ifstream in(filename);
    if (!in) {
      return nlohmann::json::object();
    }
    std::stringstream buffer;
    buffer << in.rdbuf();
    const std::string content = buffer.str();
    if (content.empty()) {
      return nlohmann::json::object();
    }
    const auto parsed =
        nlohmann::json::parse(content, nullptr, /*allow_exceptions=*/false);
    if (parsed.is_discarded() || !parsed.is_object()) {
      // Polling may read a shard while it is being rewritten.
      return nlohmann::json::object();
    }
    return parsed;
  }

  static nlohmann::json
  broadcast_json_object(World &world, nlohmann::json payload,
                        int root_rank = 0) {
    std::string payload_dump;
    if (world.rank() == root_rank) {
      payload_dump = payload.dump();
    }
    world.gop.broadcast_serializable(payload_dump, root_rank);
    if (payload_dump.empty()) {
      return nlohmann::json::object();
    }
    return nlohmann::json::parse(payload_dump, nullptr,
                                 /*allow_exceptions=*/true);
  }

  /// In-place broadcast of a JSON object from `root` to all ranks.
  /// On rank `root` the object is serialized and sent; on other ranks
  /// the object is replaced with the deserialized result.
  static void broadcast_json(World &world, nlohmann::json &obj, int root = 0) {
    std::string dump;
    if (world.rank() == root) dump = obj.dump();
    world.gop.broadcast_serializable(dump, root);
    if (world.rank() != root)
      obj = dump.empty() ? nlohmann::json::object() : nlohmann::json::parse(dump);
  }

  /// RAII-style helper: create `ngroups` MacroTaskQ subworlds, set the default
  /// pmap to the subworld for the duration of `fn`, then restore the original
  /// pmap even if `fn` throws.
  template <typename Fn>
  static void with_subworld(World &world, size_t ngroups, Fn &&fn) {
    auto subworld_ptr = MacroTaskQ::create_worlds(world, ngroups);
    if (!subworld_ptr) throw std::runtime_error("subworld creation returned null");
    World &subworld = *subworld_ptr;
    auto old_pmap = FunctionDefaults<3>::get_pmap();
    FunctionDefaults<3>::set_default_pmap(subworld);
    try {
      fn(subworld);
      FunctionDefaults<3>::set_pmap(old_pmap);
    } catch (...) {
      FunctionDefaults<3>::set_pmap(old_pmap);
      throw;
    }
  }

  static void merge_state_metadata_json(nlohmann::json &merged,
                                        const nlohmann::json &shard) {
    if (!merged.is_object()) {
      merged = nlohmann::json::object();
    }
    if (!merged.contains("states") || !merged["states"].is_object()) {
      merged["states"] = nlohmann::json::object();
    }
    if (!shard.is_object()) {
      return;
    }

    if (shard.contains("states") && shard["states"].is_object()) {
      for (const auto &[state_id, shard_state] : shard["states"].items()) {
        auto &dst_state = merged["states"][state_id];
        if (!dst_state.is_object()) {
          dst_state = nlohmann::json::object();
        }
        if (!dst_state.contains("protocols") ||
            !dst_state["protocols"].is_object()) {
          dst_state["protocols"] = nlohmann::json::object();
        }

        if (shard_state.contains("final_saved") &&
            shard_state["final_saved"].is_boolean()) {
          const bool lhs = dst_state.contains("final_saved") &&
                           dst_state["final_saved"].is_boolean() &&
                           dst_state["final_saved"].get<bool>();
          const bool rhs = shard_state["final_saved"].get<bool>();
          dst_state["final_saved"] = (lhs || rhs);
        }

        if (!shard_state.contains("protocols") ||
            !shard_state["protocols"].is_object()) {
          continue;
        }

        for (const auto &[protocol_key, shard_proto] :
             shard_state["protocols"].items()) {
          auto &dst_proto = dst_state["protocols"][protocol_key];
          if (!dst_proto.is_object()) {
            dst_proto = nlohmann::json::object();
          }
          if (!dst_proto.contains("saved") || !dst_proto["saved"].is_object()) {
            dst_proto["saved"] = nlohmann::json::object();
          }
          if (!dst_proto.contains("converged") ||
              !dst_proto["converged"].is_object()) {
            dst_proto["converged"] = nlohmann::json::object();
          }
          if (!dst_proto.contains("timings") ||
              !dst_proto["timings"].is_object()) {
            dst_proto["timings"] = nlohmann::json::object();
          }
          if (!dst_proto.contains("solver_diagnostics") ||
              !dst_proto["solver_diagnostics"].is_object()) {
            dst_proto["solver_diagnostics"] = nlohmann::json::object();
          }
          if (!dst_proto.contains("restart_provenance") ||
              !dst_proto["restart_provenance"].is_object()) {
            dst_proto["restart_provenance"] = nlohmann::json::object();
          }

          auto merge_flag_map = [&](const char *name) {
            if (!shard_proto.contains(name) || !shard_proto[name].is_object()) {
              return;
            }
            for (const auto &[freq_key, shard_value] :
                 shard_proto[name].items()) {
              const bool rhs =
                  shard_value.is_boolean() && shard_value.get<bool>();
              const bool lhs = dst_proto[name].contains(freq_key) &&
                               dst_proto[name][freq_key].is_boolean() &&
                               dst_proto[name][freq_key].get<bool>();
              dst_proto[name][freq_key] = (lhs || rhs);
            }
          };

          merge_flag_map("saved");
          merge_flag_map("converged");
          if (shard_proto.contains("timings") &&
              shard_proto["timings"].is_object()) {
            for (const auto &[freq_key, timing_value] :
                 shard_proto["timings"].items()) {
              if (!dst_proto["timings"].contains(freq_key)) {
                dst_proto["timings"][freq_key] = timing_value;
              }
            }
          }
          if (shard_proto.contains("solver_diagnostics") &&
              shard_proto["solver_diagnostics"].is_object()) {
            for (const auto &[freq_key, diagnostics_value] :
                 shard_proto["solver_diagnostics"].items()) {
              if (!dst_proto["solver_diagnostics"].contains(freq_key)) {
                dst_proto["solver_diagnostics"][freq_key] = diagnostics_value;
              }
            }
          }
          if (shard_proto.contains("restart_provenance") &&
              shard_proto["restart_provenance"].is_object()) {
            for (const auto &[freq_key, provenance_value] :
                 shard_proto["restart_provenance"].items()) {
              if (!dst_proto["restart_provenance"].contains(freq_key)) {
                dst_proto["restart_provenance"][freq_key] = provenance_value;
              }
            }
          }
        }
      }
    }

    if (shard.contains("excited_states") && shard["excited_states"].is_object()) {
      if (!merged.contains("excited_states") ||
          !merged["excited_states"].is_object()) {
        merged["excited_states"] = nlohmann::json::object();
      }
      auto &dst_excited = merged["excited_states"];
      const auto &src_excited = shard["excited_states"];

      if (src_excited.contains("plan") && src_excited["plan"].is_object()) {
        if (!dst_excited.contains("plan") || !dst_excited["plan"].is_object()) {
          dst_excited["plan"] = src_excited["plan"];
        } else {
          for (const auto &[key, value] : src_excited["plan"].items()) {
            dst_excited["plan"][key] = value;
          }
        }
      }

      if (src_excited.contains("protocols") &&
          src_excited["protocols"].is_object()) {
        if (!dst_excited.contains("protocols") ||
            !dst_excited["protocols"].is_object()) {
          dst_excited["protocols"] = nlohmann::json::object();
        }
        auto &dst_protocols = dst_excited["protocols"];
        for (const auto &[protocol_key, src_protocol] :
             src_excited["protocols"].items()) {
          auto &dst_protocol = dst_protocols[protocol_key];
          if (!dst_protocol.is_object()) {
            dst_protocol = nlohmann::json::object();
          }
          const bool src_saved =
              src_protocol.contains("saved") && src_protocol["saved"].is_boolean() &&
              src_protocol["saved"].get<bool>();
          const bool src_converged = src_protocol.contains("converged") &&
                                     src_protocol["converged"].is_boolean() &&
                                     src_protocol["converged"].get<bool>();
          const bool dst_saved =
              dst_protocol.contains("saved") &&
              dst_protocol["saved"].is_boolean() &&
              dst_protocol["saved"].get<bool>();
          const bool dst_converged =
              dst_protocol.contains("converged") &&
              dst_protocol["converged"].is_boolean() &&
              dst_protocol["converged"].get<bool>();
          dst_protocol["saved"] = (dst_saved || src_saved);
          dst_protocol["converged"] = (dst_converged || src_converged);

          if (src_protocol.contains("timings") &&
              src_protocol["timings"].is_object()) {
            if (!dst_protocol.contains("timings") ||
                !dst_protocol["timings"].is_object()) {
              dst_protocol["timings"] = src_protocol["timings"];
            } else {
              for (const auto &[k, v] : src_protocol["timings"].items()) {
                const bool take_value =
                    !dst_protocol["timings"].contains(k) ||
                    (dst_protocol["timings"][k].is_number() &&
                     dst_protocol["timings"][k].get<double>() == 0.0);
                if (take_value) {
                  dst_protocol["timings"][k] = v;
                }
              }
            }
          }

          if (src_protocol.contains("energies") &&
              src_protocol["energies"].is_array()) {
            if (!dst_protocol.contains("energies") ||
                !dst_protocol["energies"].is_array() ||
                dst_protocol["energies"].empty()) {
              dst_protocol["energies"] = src_protocol["energies"];
            }
          }
          if (src_protocol.contains("state_names") &&
              src_protocol["state_names"].is_array()) {
            if (!dst_protocol.contains("state_names") ||
                !dst_protocol["state_names"].is_array() ||
                dst_protocol["state_names"].empty()) {
              dst_protocol["state_names"] = src_protocol["state_names"];
            }
          }
          if (src_protocol.contains("roots") && src_protocol["roots"].is_array()) {
            if (!dst_protocol.contains("roots") ||
                !dst_protocol["roots"].is_array() ||
                dst_protocol["roots"].empty()) {
              dst_protocol["roots"] = src_protocol["roots"];
            }
          }
          if (src_protocol.contains("iterations") &&
              src_protocol["iterations"].is_number_unsigned()) {
            const auto src_iterations =
                src_protocol["iterations"].get<size_t>();
            const auto dst_iterations =
                (dst_protocol.contains("iterations") &&
                 dst_protocol["iterations"].is_number_unsigned())
                    ? dst_protocol["iterations"].get<size_t>()
                    : static_cast<size_t>(0);
            dst_protocol["iterations"] =
                std::max(dst_iterations, src_iterations);
          }
          if (src_protocol.contains("residual_norms") &&
              src_protocol["residual_norms"].is_array()) {
            if (!dst_protocol.contains("residual_norms") ||
                !dst_protocol["residual_norms"].is_array() ||
                dst_protocol["residual_norms"].empty()) {
              dst_protocol["residual_norms"] = src_protocol["residual_norms"];
            }
          }
          if (src_protocol.contains("density_change_norms") &&
              src_protocol["density_change_norms"].is_array()) {
            if (!dst_protocol.contains("density_change_norms") ||
                !dst_protocol["density_change_norms"].is_array() ||
                dst_protocol["density_change_norms"].empty()) {
              dst_protocol["density_change_norms"] =
                  src_protocol["density_change_norms"];
            }
          }
          if (src_protocol.contains("relative_residual_norms") &&
              src_protocol["relative_residual_norms"].is_array()) {
            if (!dst_protocol.contains("relative_residual_norms") ||
                !dst_protocol["relative_residual_norms"].is_array() ||
                dst_protocol["relative_residual_norms"].empty()) {
              dst_protocol["relative_residual_norms"] =
                  src_protocol["relative_residual_norms"];
            }
          }
          if (src_protocol.contains("iteration_max_residuals") &&
              src_protocol["iteration_max_residuals"].is_array()) {
            if (!dst_protocol.contains("iteration_max_residuals") ||
                !dst_protocol["iteration_max_residuals"].is_array() ||
                dst_protocol["iteration_max_residuals"].size() <
                    src_protocol["iteration_max_residuals"].size()) {
              dst_protocol["iteration_max_residuals"] =
                  src_protocol["iteration_max_residuals"];
            }
          }
          if (src_protocol.contains("iteration_max_density_changes") &&
              src_protocol["iteration_max_density_changes"].is_array()) {
            if (!dst_protocol.contains("iteration_max_density_changes") ||
                !dst_protocol["iteration_max_density_changes"].is_array() ||
                dst_protocol["iteration_max_density_changes"].size() <
                    src_protocol["iteration_max_density_changes"].size()) {
              dst_protocol["iteration_max_density_changes"] =
                  src_protocol["iteration_max_density_changes"];
            }
          }
          if (src_protocol.contains("iteration_max_relative_residuals") &&
              src_protocol["iteration_max_relative_residuals"].is_array()) {
            if (!dst_protocol.contains("iteration_max_relative_residuals") ||
                !dst_protocol["iteration_max_relative_residuals"].is_array() ||
                dst_protocol["iteration_max_relative_residuals"].size() <
                    src_protocol["iteration_max_relative_residuals"].size()) {
              dst_protocol["iteration_max_relative_residuals"] =
                  src_protocol["iteration_max_relative_residuals"];
            }
          }
          if (src_protocol.contains("convergence_mode") &&
              src_protocol["convergence_mode"].is_string()) {
            if (!dst_protocol.contains("convergence_mode") ||
                !dst_protocol["convergence_mode"].is_string() ||
                dst_protocol["convergence_mode"].get<std::string>().empty()) {
              dst_protocol["convergence_mode"] =
                  src_protocol["convergence_mode"];
            }
          }
          if (src_protocol.contains("accelerator_mode") &&
              src_protocol["accelerator_mode"].is_string()) {
            if (!dst_protocol.contains("accelerator_mode") ||
                !dst_protocol["accelerator_mode"].is_string() ||
                dst_protocol["accelerator_mode"].get<std::string>().empty()) {
              dst_protocol["accelerator_mode"] =
                  src_protocol["accelerator_mode"];
            }
          }
          if (src_protocol.contains("accelerator_subspace") &&
              src_protocol["accelerator_subspace"].is_number_unsigned()) {
            const auto src_subspace =
                src_protocol["accelerator_subspace"].get<size_t>();
            const auto dst_subspace =
                (dst_protocol.contains("accelerator_subspace") &&
                 dst_protocol["accelerator_subspace"].is_number_unsigned())
                    ? dst_protocol["accelerator_subspace"].get<size_t>()
                    : static_cast<size_t>(0);
            dst_protocol["accelerator_subspace"] =
                std::max(dst_subspace, src_subspace);
          }
          if (src_protocol.contains("density_convergence_target") &&
              src_protocol["density_convergence_target"].is_number()) {
            if (!dst_protocol.contains("density_convergence_target") ||
                !dst_protocol["density_convergence_target"].is_number()) {
              dst_protocol["density_convergence_target"] =
                  src_protocol["density_convergence_target"];
            }
          }
          if (src_protocol.contains("relative_convergence_target") &&
              src_protocol["relative_convergence_target"].is_number()) {
            if (!dst_protocol.contains("relative_convergence_target") ||
                !dst_protocol["relative_convergence_target"].is_number()) {
              dst_protocol["relative_convergence_target"] =
                  src_protocol["relative_convergence_target"];
            }
          }
          if (src_protocol.contains("max_rotation") &&
              src_protocol["max_rotation"].is_number()) {
            if (!dst_protocol.contains("max_rotation") ||
                !dst_protocol["max_rotation"].is_number()) {
              dst_protocol["max_rotation"] = src_protocol["max_rotation"];
            }
          }
        }
      }
    }
  }

  static void merge_debug_log_recursive(nlohmann::json &dst,
                                        const nlohmann::json &src) {
    if (!src.is_object()) {
      dst = src;
      return;
    }

    if (!dst.is_object()) {
      dst = nlohmann::json::object();
    }

    for (const auto &[key, src_value] : src.items()) {
      if (!dst.contains(key)) {
        dst[key] = src_value;
        continue;
      }

      auto &dst_value = dst[key];
      if (dst_value.is_object() && src_value.is_object()) {
        merge_debug_log_recursive(dst_value, src_value);
      } else if (dst_value.is_array() && src_value.is_array()) {
        for (const auto &entry : src_value) {
          dst_value.push_back(entry);
        }
      } else {
        dst_value = src_value;
      }
    }
  }

  static void merge_debug_log_json(nlohmann::json &merged,
                                   const nlohmann::json &shard) {
    if (!shard.is_object()) {
      return;
    }
    if (!merged.is_object()) {
      merged = nlohmann::json::object();
    }
    merge_debug_log_recursive(merged, shard);
  }

  static bool point_ready_in_metadata(const nlohmann::json &metadata,
                                      const LinearResponsePoint &pt,
                                      bool require_saved,
                                      bool require_converged) {
    if (!metadata.is_object() || !metadata.contains("states") ||
        !metadata["states"].is_object()) {
      return false;
    }

    const auto state_key = pt.perturbationDescription();
    const auto protocol_key = ResponseRecord2::protocol_key(pt.threshold());
    const auto freq_key = ResponseRecord2::freq_key(pt.frequency());

    const auto states_it = metadata["states"].find(state_key);
    if (states_it == metadata["states"].end()) {
      return false;
    }
    if (!states_it->contains("protocols") ||
        !(*states_it)["protocols"].is_object()) {
      return false;
    }
    const auto protos_it = (*states_it)["protocols"].find(protocol_key);
    if (protos_it == (*states_it)["protocols"].end()) {
      return false;
    }

    auto check_flag = [&](const char *flag_name, bool required) {
      if (!required) {
        return true;
      }
      if (!protos_it->contains(flag_name) ||
          !(*protos_it)[flag_name].is_object()) {
        return false;
      }
      const auto values_it = (*protos_it)[flag_name].find(freq_key);
      if (values_it == (*protos_it)[flag_name].end()) {
        return false;
      }
      return values_it->is_boolean() && values_it->get<bool>();
    };

    return check_flag("saved", require_saved) &&
           check_flag("converged", require_converged);
  }

  static bool
  point_marked_for_frequency_removal_in_metadata(const nlohmann::json &metadata,
                                                 const LinearResponsePoint &pt) {
    if (!metadata.is_object() || !metadata.contains("states") ||
        !metadata["states"].is_object()) {
      return false;
    }

    const auto state_key = pt.perturbationDescription();
    const auto protocol_key = ResponseRecord2::protocol_key(pt.threshold());
    const auto freq_key = ResponseRecord2::freq_key(pt.frequency());

    const auto states_it = metadata["states"].find(state_key);
    if (states_it == metadata["states"].end()) {
      return false;
    }
    if (!states_it->contains("protocols") ||
        !(*states_it)["protocols"].is_object()) {
      return false;
    }
    const auto protos_it = (*states_it)["protocols"].find(protocol_key);
    if (protos_it == (*states_it)["protocols"].end()) {
      return false;
    }
    if (!protos_it->contains("solver_diagnostics") ||
        !(*protos_it)["solver_diagnostics"].is_object()) {
      return false;
    }
    const auto solver_it = (*protos_it)["solver_diagnostics"].find(freq_key);
    if (solver_it == (*protos_it)["solver_diagnostics"].end() ||
        !solver_it->is_object()) {
      return false;
    }
    if (!solver_it->contains("remove_from_frequency_set")) {
      return false;
    }
    return (*solver_it)["remove_from_frequency_set"].is_boolean() &&
           (*solver_it)["remove_from_frequency_set"].get<bool>();
  }

  // ============================================================================
  // SECTION 7: Stage 1 — Planning (ground context + state generation)
  //   make_ground_context  — reads SCF checkpoint, builds GroundStateData and
  //                          ResponseManager for the run directory.
  //   build_excited_state_bundle_plan — translates response knobs into
  //                          ExcitedStateBundlePlan.
  //   plan_required_states — combines StateGenerator, DerivedStatePlanner, and
  //                          StateParallelPlanner outputs into PlannedStates.
  // ============================================================================

  /// Build stage-1 ground/runtime context from SCF checkpoint artifacts.
  ///
  /// Loads molecule data, resolves archive/fock paths, and constructs the
  /// reusable ground-state and response-manager objects used throughout stages 2-3.
  static GroundContext
  make_ground_context(World &world, const CalculationParameters &calc_params,
                      const std::shared_ptr<SCF> &scf_calc,
                      const std::filesystem::path &outdir) {
    auto indir = scf_calc->work_dir;
    auto rel = std::filesystem::relative(indir, outdir);
    auto prox = std::filesystem::proximate(indir, outdir);

    if (world.rank() == 0) {
      print("RUN_CONTEXT outdir=", outdir);
      print("GROUND_ARCHIVE path=", indir);
      print("GROUND_ARCHIVE_RELATIVE path=", rel);
      print("GROUND_ARCHIVE_PROXIMATE path=", prox);
    }

    const auto &prefix = calc_params.prefix();
    std::string archive_name = prefix + ".restartdata";
    std::string fock_json_file = (prox / (prefix + ".fock.json")).string();
    std::string moldft_checkpt = prox / "moldft.calc_info.json";
    auto relative_archive = prox / archive_name;

    if (!std::filesystem::exists(moldft_checkpt)) {
      if (world.rank() == 0) {
        print("ERROR MISSING_GROUND_CHECKPOINT file=", moldft_checkpt);
      }
      throw std::runtime_error("Missing ground-state checkpoint file");
    }

    auto read_molecule = [](const std::string &scf_ckpt) -> Molecule {
      std::ifstream ifs(scf_ckpt);
      nlohmann::json j;
      ifs >> j;
      ifs.close();
      Molecule mol;
      if (j.contains("molecule")) {
        mol.from_json(j["molecule"]);
      } else {
        throw std::runtime_error(
            "Molecule information missing from checkpoint JSON.");
      }
      return mol;
    };

    Molecule molecule = read_molecule(moldft_checkpt);
    //print("Read molecule with ", molecule.natom(), " atoms.");

    GroundStateData ground(world, relative_archive.string(), molecule);
    ResponseManager response_manager(world, calc_params);

    return GroundContext{std::move(molecule), std::move(ground),
                         std::move(response_manager), relative_archive.string(),
                         std::move(fock_json_file)};
  }

  /// Translate response knobs into a protocol-indexed excited-stage plan.
  static ExcitedStateBundlePlan
  build_excited_state_bundle_plan(const CalculationParameters &calc_params,
                                  const ResponseParameters &response_params,
                                  const StateParallelPlan &state_parallel_plan) {
    ExcitedStateBundlePlan plan;
    plan.enabled = response_params.excited_enable();
    plan.num_states = response_params.excited_num_states();
    plan.tda = response_params.excited_tda();
    plan.guess_max_iter = response_params.excited_guess_max_iter();
    plan.maxiter = response_params.excited_maxiter();
    plan.maxsub = response_params.excited_maxsub();
    const size_t max_owner_group =
        state_parallel_plan.mapping_groups > 0
            ? (state_parallel_plan.mapping_groups - 1)
            : 0;
    plan.owner_group =
        std::min(response_params.excited_owner_group(), max_owner_group);
    plan.protocols = calc_params.protocol();
    return plan;
  }

  /// Stage-1 planner entrypoint for linear, derived, and excited requests.
  static PlannedStates
  plan_required_states(World &world, const CalculationParameters &calc_params,
                       const GroundContext &ctx,
                       const ResponseParameters &response_params) {
    StateGenerator state_generator(ctx.molecule, calc_params.protocol(),
                                   ctx.ground.isSpinRestricted(),
                                   response_params);
    auto generated_states = state_generator.generateStates();

    if (world.rank() == 0) {
      GeneratedStateData::print_generated_state_map(generated_states.state_map);
    }
    DerivedStatePlan derived_state_plan =
        DerivedStatePlanner::build_vbc_driven_quadratic_plan(
            response_params, ctx.molecule, ctx.ground.isSpinRestricted(),
            calc_params.protocol());
    StateParallelPlan state_parallel_plan = StateParallelPlanner::build(
        response_params, world.size(), generated_states.states);
    ExcitedStateBundlePlan excited_state_bundle_plan =
        build_excited_state_bundle_plan(calc_params, response_params,
                                        state_parallel_plan);
    if (world.rank() == 0 && !derived_state_plan.requests.empty()) {
      print("PLAN_DERIVED requests=", derived_state_plan.requests.size(),
            " kind=vbc_quadratic");
    }
    if (world.rank() == 0 && excited_state_bundle_plan.enabled) {
      print("PLAN_EXCITED_BUNDLE states=",
            excited_state_bundle_plan.num_states,
            " owner_group=", excited_state_bundle_plan.owner_group,
            " tda=", excited_state_bundle_plan.tda);
    }
    if (world.rank() == 0 && response_params.state_parallel() != "off") {
      print(
          "STATE_PARALLEL_PLAN mode=",
          state_parallel_plan.effective_mode,
          " requested_groups=", state_parallel_plan.requested_groups,
          " mapping_groups=", state_parallel_plan.mapping_groups,
          " channel_owner_groups=", state_parallel_plan.channel_owner_groups,
          " effective_point_groups=", state_parallel_plan.effective_point_groups,
          " point_parallel_start_protocol_index=",
          state_parallel_plan.point_parallel_start_protocol_index,
          " frequency_policy=", state_parallel_plan.frequency_partition_policy,
          " reason=", state_parallel_plan.reason);
    }
    world.gop.fence();
    return PlannedStates{std::move(generated_states),
                         std::move(derived_state_plan),
                         std::move(excited_state_bundle_plan),
                         std::move(state_parallel_plan)};
  }

  // ============================================================================
  // SECTION 8: Stage 2a — Schedule context types and builders
  //   ProtocolExecutionPolicy — per-protocol active-groups / mode descriptor.
  //   RuntimePointOwnershipPolicy — restart-aware runtime policy aggregated
  //                                 from planner output + existing metadata.
  //   StateSolveScheduleContext — full scheduling context passed to the serial
  //                               and subgroup solve paths.
  //   PendingPointWorkItem / PendingProtocolManifest — work lists for one
  //                               protocol step in one owner lane.
  //   build_owner_by_channel_index — maps channel_index → owner lane.
  //   build_protocol_execution_policy — reads existing metadata and computes
  //                               per-protocol active groups and mode.
  //   compute_runtime_point_ownership_policy — wraps restart + point-parallel
  //                               start computation into RuntimePointOwnershipPolicy.
  //   build_local_channel_workset — filters global channel list to subgroup.
  //   use_channel_series_ownership_for_protocol_runtime — inline policy gate.
  //   active_owner_groups_for_protocol_runtime — inline policy lookup.
  //   point_needs_solving / point_needs_solving_from_metadata — convergence
  //                               predicates (persistence-backed and metadata-backed).
  //   any_state_point_needs_solving — short-circuit scan across state indices.
  //   run_protocol_threshold_loop — generic protocol iteration driver.
  //   run_frequency_loop_with_flush — thin wrapper around computeFrequencyLoop.
  //   build_pending_work_manifest / build_pending_manifest_from_metadata —
  //                               construct per-protocol, per-lane work lists.
  //   print_state_solve_execution_mode — logs selected execution strategy.
  //   build_state_solve_schedule_context — top-level schedule context factory.
  // ============================================================================

  struct ProtocolExecutionPolicy {
    bool use_channel_series_mode = true;
    size_t active_groups = 1;
    size_t pending_channels = 0;
    size_t pending_points = 0;

    [[nodiscard]] nlohmann::json to_json() const {
      return {{"use_channel_series_mode", use_channel_series_mode},
              {"active_groups", active_groups},
              {"pending_channels", pending_channels},
              {"pending_points", pending_points}};
    }
  };

  struct RuntimePointOwnershipPolicy {
    size_t point_parallel_start_protocol_index = 0;
    size_t restart_start_protocol_index = 0;
    bool restart_final_protocol_only = false;
    bool force_retry_removed_frequencies = false;
    size_t remaining_final_pending_points = 0;
    size_t new_points_without_history = 0;
    size_t runtime_execution_groups = 1;
    size_t runtime_effective_point_groups = 1;
    std::vector<ProtocolExecutionPolicy> protocol_policies;
  };

  static std::vector<size_t>
  build_owner_by_channel_index(const std::vector<LinearResponseDescriptor> &states,
                             const StateParallelPlan &state_parallel_plan) {
    std::vector<size_t> owner_by_channel_index(states.size(), 0);
    for (const auto &assignment : state_parallel_plan.channel_assignments) {
      if (assignment.channel_index < owner_by_channel_index.size()) {
        owner_by_channel_index[assignment.channel_index] =
            assignment.owner_group;
      }
    }
    return owner_by_channel_index;
  }

  static std::vector<ProtocolExecutionPolicy>
  build_protocol_execution_policy(World &world,
                                  const CalculationParameters &calc_params,
                                  const StateParallelPlan &state_parallel_plan,
                                  const std::vector<LinearResponseDescriptor>
                                      &linear_states,
                                  bool owner_group_schedule,
                                  size_t point_parallel_start_protocol_index,
                                  size_t restart_start_protocol_index,
                                  bool force_retry_removed_frequencies,
                                  size_t runtime_execution_groups,
                                  size_t runtime_effective_point_groups) {
    const auto &protocol = calc_params.protocol();
    std::vector<ProtocolExecutionPolicy> protocol_policies(protocol.size());
    if (protocol.empty()) {
      return protocol_policies;
    }

    if (world.rank() == 0) {
      const nlohmann::json existing_metadata =
          read_json_file_or_object("response_metadata.json");
      const size_t point_owner_groups =
          std::max<size_t>(
              1, std::min(runtime_effective_point_groups, runtime_execution_groups));
      const size_t channel_owner_groups =
          std::max<size_t>(
              1, std::min(state_parallel_plan.channel_owner_groups,
                          runtime_execution_groups));
      for (size_t ti = 0; ti < protocol.size(); ++ti) {
        ProtocolExecutionPolicy policy;
        if (ti < restart_start_protocol_index) {
          policy.use_channel_series_mode = true;
          policy.active_groups = 1;
          protocol_policies[ti] = policy;
          continue;
        }
        const bool at_final_protocol = (ti + 1 == protocol.size());

        policy.use_channel_series_mode =
            !owner_group_schedule || ti < point_parallel_start_protocol_index;

        for (const auto &state : linear_states) {
          bool state_has_pending = false;
          for (size_t freq_idx = 0; freq_idx < state.num_frequencies();
               ++freq_idx) {
            LinearResponsePoint pt{state, ti, freq_idx};
            const bool point_pending = point_needs_solving_from_metadata(
                existing_metadata, pt, at_final_protocol,
                force_retry_removed_frequencies);
            if (point_pending) {
              ++policy.pending_points;
              state_has_pending = true;
            }
          }
          if (state_has_pending) {
            ++policy.pending_channels;
          }
        }

        // Respect requested protocol ownership split:
        // - point_start_protocol == 0 enables channel-point ownership at
        //   protocol-0, exposing frequency-level parallelism immediately.
        // - restart_start_protocol_index (handled above) skips completed
        //   protocol prefixes and starts from earliest pending protocol.

        if (!owner_group_schedule || policy.pending_points == 0) {
          policy.active_groups = 1;
        } else if (policy.use_channel_series_mode) {
          // Keep full state-owner lanes active so fixed owner mapping remains
          // stable across restarts/partial completion.
          policy.active_groups =
              std::max<size_t>(1, std::min(channel_owner_groups,
                                           runtime_execution_groups));
        } else {
          policy.active_groups = std::max<size_t>(
              1, std::min({point_owner_groups, policy.pending_points,
                           runtime_execution_groups}));
        }
        protocol_policies[ti] = policy;
      }
    }

    nlohmann::json payload = nlohmann::json::array();
    if (world.rank() == 0) {
      for (const auto &policy : protocol_policies) {
        payload.push_back(policy.to_json());
      }
    }
    broadcast_json(world, payload);
    if (!payload.is_array()) {
      return protocol_policies;
    }
    protocol_policies.clear();
    protocol_policies.reserve(payload.size());
    for (const auto &entry : payload) {
      ProtocolExecutionPolicy policy;
      if (entry.is_object()) {
        policy.use_channel_series_mode =
            entry.value("use_channel_series_mode", true);
        policy.active_groups = std::max<size_t>(1, entry.value("active_groups", 1));
        policy.pending_channels = entry.value("pending_channels", size_t(0));
        policy.pending_points = entry.value("pending_points", size_t(0));
      }
      protocol_policies.push_back(policy);
    }
    return protocol_policies;
  }

  static RuntimePointOwnershipPolicy compute_runtime_point_ownership_policy(
      World &world, const CalculationParameters &calc_params,
      const StateParallelPlan &state_parallel_plan,
      const std::vector<LinearResponseDescriptor> &linear_states,
      bool owner_group_schedule, const ResponseParameters &response_params) {
    RuntimePointOwnershipPolicy runtime_policy;
    runtime_policy.point_parallel_start_protocol_index =
        state_parallel_plan.point_parallel_start_protocol_index;
    runtime_policy.force_retry_removed_frequencies =
        response_params.force_retry_removed_frequencies();
    runtime_policy.runtime_execution_groups =
        std::max<size_t>(1, state_parallel_plan.execution_groups);
    runtime_policy.runtime_effective_point_groups =
        std::max<size_t>(1, state_parallel_plan.effective_point_groups);

    const bool has_protocol_thresholds = !calc_params.protocol().empty();
    if (has_protocol_thresholds) {
      if (world.rank() == 0) {
        const nlohmann::json existing_metadata =
            read_json_file_or_object("response_metadata.json");
        const size_t protocol_count = calc_params.protocol().size();
        const size_t final_protocol_index = calc_params.protocol().size() - 1;
        size_t final_pending_points = 0;
        size_t new_points_without_history = 0;
        size_t earliest_pending_protocol = protocol_count;
        for (const auto &state : linear_states) {
          for (size_t freq_idx = 0; freq_idx < state.num_frequencies();
               ++freq_idx) {
            if (!runtime_policy.force_retry_removed_frequencies) {
              const LinearResponsePoint protocol0_pt{state, /*thresh_index=*/0,
                                                     freq_idx};
              if (point_marked_for_frequency_removal_in_metadata(
                      existing_metadata, protocol0_pt)) {
                continue;
              }
            }
            bool has_saved_history = false;
            for (size_t ti = 0; ti < calc_params.protocol().size(); ++ti) {
              LinearResponsePoint historical_pt{state, ti, freq_idx};
              if (point_ready_in_metadata(existing_metadata, historical_pt,
                                          /*require_saved=*/true,
                                          /*require_converged=*/false)) {
                has_saved_history = true;
                break;
              }
            }
            if (!has_saved_history) {
              ++new_points_without_history;
            }

            for (size_t ti = 0; ti < protocol_count; ++ti) {
              const bool at_final_protocol = (ti + 1 == protocol_count);
              LinearResponsePoint pt{state, ti, freq_idx};
              const bool needs_solving = point_needs_solving_from_metadata(
                  existing_metadata, pt, at_final_protocol,
                  runtime_policy.force_retry_removed_frequencies);
              if (!needs_solving) {
                continue;
              }
              earliest_pending_protocol = std::min(earliest_pending_protocol, ti);
              if (ti == final_protocol_index) {
                ++final_pending_points;
              }
            }
          }
        }

        runtime_policy.restart_start_protocol_index = earliest_pending_protocol;
        runtime_policy.restart_final_protocol_only =
            (runtime_policy.restart_start_protocol_index == final_protocol_index);
        runtime_policy.remaining_final_pending_points = final_pending_points;
        runtime_policy.new_points_without_history = new_points_without_history;
      }
      world.gop.broadcast_serializable(runtime_policy.restart_start_protocol_index,
                                       0);
      world.gop.broadcast_serializable(runtime_policy.restart_final_protocol_only,
                                       0);
      world.gop.broadcast_serializable(runtime_policy.remaining_final_pending_points,
                                       0);
      world.gop.broadcast_serializable(runtime_policy.new_points_without_history,
                                       0);
    }

    if (has_protocol_thresholds) {
      runtime_policy.point_parallel_start_protocol_index = std::max(
          state_parallel_plan.point_parallel_start_protocol_index,
          runtime_policy.restart_start_protocol_index);
    }
    if (runtime_policy.restart_final_protocol_only && has_protocol_thresholds &&
        runtime_policy.restart_start_protocol_index ==
            calc_params.protocol().size() - 1) {
      const size_t capped_groups = std::max<size_t>(
          1, std::min(state_parallel_plan.execution_groups,
                      runtime_policy.remaining_final_pending_points));
      runtime_policy.runtime_execution_groups = capped_groups;
      runtime_policy.runtime_effective_point_groups = std::max<size_t>(
          1, std::min(state_parallel_plan.effective_point_groups, capped_groups));
    }
    runtime_policy.protocol_policies = build_protocol_execution_policy(
        world, calc_params, state_parallel_plan, linear_states,
        owner_group_schedule, runtime_policy.point_parallel_start_protocol_index,
        runtime_policy.restart_start_protocol_index,
        runtime_policy.force_retry_removed_frequencies,
        runtime_policy.runtime_execution_groups,
        runtime_policy.runtime_effective_point_groups);
    return runtime_policy;
  }

  static void build_local_channel_workset(
      const std::vector<LinearResponseDescriptor> &linear_states,
      const std::vector<size_t> &owner_by_channel_index,
      const PointOwnershipScheduler &point_scheduler, size_t subgroup_id,
      size_t mapping_groups, std::vector<size_t> &local_channel_indices,
      std::vector<LinearResponseDescriptor> &local_channels) {
    local_channel_indices.clear();
    local_channels.clear();
    local_channel_indices.reserve(
        linear_states.size() / std::max<size_t>(1, mapping_groups));
    local_channels.reserve(linear_states.size());

    const size_t point_owner_groups = point_scheduler.owner_groups();
    for (size_t channel_index = 0; channel_index < linear_states.size();
         ++channel_index) {
      const bool owns_channel_at_first_protocol =
          owner_by_channel_index[channel_index] == subgroup_id;
      // Point-mode ownership is computed from protocol-local pending manifests.
      // Keep all channels available in each subgroup once point mode is enabled
      // to avoid dropping channels whose owned points are determined only at
      // runtime.
      const bool owns_any_point_after_first_protocol = (point_owner_groups > 1);
      if (owns_channel_at_first_protocol) {
        local_channel_indices.push_back(channel_index);
      }
      if (owns_channel_at_first_protocol || owns_any_point_after_first_protocol) {
        local_channels.push_back(linear_states[channel_index]);
      }
    }
  }

  struct StateSolveScheduleContext {
    // True when subgroup solve path is enabled and requested.
    bool subgroup_parallel_requested = false;
    // True when deterministic owner-lane scheduling is active.
    bool owner_group_schedule = false;
    // Planner output (group counts, mode, protocol switch threshold).
    const StateParallelPlan &state_parallel_plan;
    // Full set of linear states generated in Stage 1.
    const std::vector<LinearResponseDescriptor> &linear_states;
    // channel_index -> owner lane for protocol ranges using channel-series
    // ownership.
    std::vector<size_t> owner_by_channel_index;
    // Deterministic (channel,freq) -> owner lane mapping for channel-point
    // ownership.
    PointOwnershipScheduler point_scheduler;
    // Runtime protocol index where channel-series ownership switches to
    // channel-point ownership.
    size_t runtime_point_parallel_start_protocol_index = 0;
    // Runtime cap for subgroup lanes used during stage-2 execution.
    size_t runtime_execution_groups = 1;
    // Restart diagnostics propagated to stage metadata.
    size_t restart_start_protocol_index = 0;
    bool restart_final_protocol_only = false;
    bool force_retry_removed_frequencies = false;
    size_t remaining_final_pending_points = 0;
    size_t new_points_without_history = 0;
    // Runtime per-protocol mode + active owner groups after restart analysis.
    std::vector<ProtocolExecutionPolicy> runtime_protocol_policies;

    StateSolveScheduleContext(
        bool subgroup_parallel_requested_, bool owner_group_schedule_,
        const StateParallelPlan &state_parallel_plan_,
        const std::vector<LinearResponseDescriptor> &linear_states_,
        std::vector<size_t> owner_by_channel_index_, size_t effective_point_groups,
        size_t runtime_point_parallel_start_protocol_index_,
        size_t runtime_execution_groups_,
        size_t restart_start_protocol_index_,
        bool restart_final_protocol_only_,
        bool force_retry_removed_frequencies_,
        size_t remaining_final_pending_points_,
        size_t new_points_without_history_,
        std::vector<ProtocolExecutionPolicy> runtime_protocol_policies_)
        : subgroup_parallel_requested(subgroup_parallel_requested_),
          owner_group_schedule(owner_group_schedule_),
          state_parallel_plan(state_parallel_plan_),
          linear_states(linear_states_),
          owner_by_channel_index(std::move(owner_by_channel_index_)),
          point_scheduler(linear_states_, effective_point_groups),
          runtime_point_parallel_start_protocol_index(
              runtime_point_parallel_start_protocol_index_),
          runtime_execution_groups(runtime_execution_groups_),
          restart_start_protocol_index(restart_start_protocol_index_),
          restart_final_protocol_only(restart_final_protocol_only_),
          force_retry_removed_frequencies(force_retry_removed_frequencies_),
          remaining_final_pending_points(remaining_final_pending_points_),
          new_points_without_history(new_points_without_history_),
          runtime_protocol_policies(std::move(runtime_protocol_policies_)) {}

    [[nodiscard]] size_t point_owner_groups() const {
      return point_scheduler.owner_groups();
    }

    [[nodiscard]] ProtocolExecutionPolicy
    protocol_policy(size_t protocol_index) const {
      if (protocol_index < runtime_protocol_policies.size()) {
        return runtime_protocol_policies[protocol_index];
      }
      ProtocolExecutionPolicy fallback;
      fallback.use_channel_series_mode =
          protocol_index < runtime_point_parallel_start_protocol_index;
      fallback.active_groups = fallback.use_channel_series_mode
                                   ? std::max<size_t>(1, std::min(
                                                            state_parallel_plan
                                                                .channel_owner_groups,
                                                            runtime_execution_groups))
                                   : std::max<size_t>(
                                         1, std::min(point_owner_groups(),
                                                     runtime_execution_groups));
      return fallback;
    }
  };

  // Runtime policy gate shared by serial and subgroup paths.
  // For protocol indices below the runtime switch threshold we keep all
  // frequencies of a channel together; after the switch threshold we fan out by
  // independent channel-frequency points.
  static bool
  use_channel_series_ownership_for_protocol_runtime(const StateSolveScheduleContext &ctx,
                                           size_t protocol_index) {
    if (ctx.state_parallel_plan.mapping_groups <= 1) {
      return true;
    }
    return ctx.protocol_policy(protocol_index).use_channel_series_mode;
  }

  static size_t active_owner_groups_for_protocol_runtime(
      const StateSolveScheduleContext &ctx, size_t protocol_index) {
    if (ctx.state_parallel_plan.mapping_groups <= 1) {
      return 1;
    }
    return std::max<size_t>(1, ctx.protocol_policy(protocol_index).active_groups);
  }

  template <typename PersistenceT>
  static bool
  point_needs_solving(PersistenceT &persistence, const LinearResponsePoint &pt,
                      bool at_final_protocol,
                      bool force_retry_removed_frequencies) {
    if (!force_retry_removed_frequencies) {
      const LinearResponsePoint protocol0_pt{pt.desc, /*thresh_index=*/0,
                                             pt.freq_index};
      if (persistence.is_removed_from_frequency_set(protocol0_pt)) {
        return false;
      }
    }
    const bool is_saved = persistence.is_saved(pt);
    return !is_saved || (at_final_protocol && !persistence.is_converged(pt));
  }

  static bool point_needs_solving_from_metadata(
      const nlohmann::json &metadata, const LinearResponsePoint &pt,
      bool at_final_protocol, bool force_retry_removed_frequencies) {
    if (!force_retry_removed_frequencies) {
      const LinearResponsePoint protocol0_pt{pt.desc, /*thresh_index=*/0,
                                             pt.freq_index};
      if (point_marked_for_frequency_removal_in_metadata(metadata,
                                                         protocol0_pt)) {
        return false;
      }
    }
    const bool is_saved = point_ready_in_metadata(metadata, pt,
                                                  /*require_saved=*/true,
                                                  /*require_converged=*/false);
    if (!is_saved) {
      return true;
    }
    if (!at_final_protocol) {
      return false;
    }
    return !point_ready_in_metadata(metadata, pt,
                                    /*require_saved=*/true,
                                    /*require_converged=*/true);
  }

  template <typename PointNeedsFn>
  static bool any_state_point_needs_solving(const std::vector<size_t> &state_indices,
                                            const StateSolveScheduleContext &schedule_ctx,
                                            size_t thresh_index,
                                            PointNeedsFn &&point_needs_solving_fn) {
    for (const auto state_index : state_indices) {
      const auto &state = schedule_ctx.linear_states[state_index];
      for (size_t freq_idx = 0; freq_idx < state.num_frequencies(); ++freq_idx) {
        if (point_needs_solving_fn(state, freq_idx)) {
          return true;
        }
      }
    }
    return false;
  }

  template <typename NeedsFn, typename SkipFn, typename PrepareFn,
            typename ExecuteFn, typename FinalizeFn>
  static void run_protocol_threshold_loop(const std::vector<double> &protocol,
                                          NeedsFn &&needs_solving_at_protocol,
                                          SkipFn &&on_skip_protocol,
                                          PrepareFn &&prepare_protocol,
                                          ExecuteFn &&execute_protocol_work,
                                          FinalizeFn &&finalize_protocol) {
    for (size_t ti = 0; ti < protocol.size(); ++ti) {
      const double thresh = protocol[ti];
      if (!needs_solving_at_protocol(thresh, ti)) {
        on_skip_protocol(thresh, ti);
        continue;
      }

      prepare_protocol(thresh, ti);
      const bool at_final_protocol = (ti + 1 == protocol.size());
      execute_protocol_work(thresh, ti, at_final_protocol);
      finalize_protocol(thresh, ti);
    }
  }

  template <typename PersistenceT, typename StateT>
  static void run_frequency_loop_with_flush(
      World &exec_world, ResponseManager &response_manager,
      const StateT &state_for_compute, size_t thresh_index,
      GroundStateData &ground_state, PersistenceT &persistence,
      bool at_final_protocol) {
    computeFrequencyLoop(exec_world, response_manager, state_for_compute,
                         thresh_index, ground_state, persistence,
                         at_final_protocol);
    persistence.flush_debug_log(exec_world);
  }

  struct PendingPointWorkItem {
    size_t state_index = 0;
    size_t freq_index = 0;
  };

  struct PendingProtocolManifest {
    // Ownership mode used for this protocol index.
    bool use_channel_series_mode = true;
    // State-level work list for state ownership mode.
    std::vector<size_t> pending_channel_indices;
    // Point-level work list for point ownership mode.
    std::vector<PendingPointWorkItem> pending_points;

    [[nodiscard]] bool has_work() const {
      return !pending_channel_indices.empty() || !pending_points.empty();
    }
  };

  template <typename PointNeedsFn>
  static PendingProtocolManifest build_pending_work_manifest(
      const StateSolveScheduleContext &schedule_ctx, size_t thresh_index,
      size_t lane_begin, size_t lane_end,
      PointNeedsFn &&point_needs_solving_fn) {
    PendingProtocolManifest manifest;
    manifest.use_channel_series_mode =
        use_channel_series_ownership_for_protocol_runtime(schedule_ctx, thresh_index);
    const size_t active_groups =
        active_owner_groups_for_protocol_runtime(schedule_ctx, thresh_index);
    const size_t lane_end_clamped = std::min(lane_end, active_groups);
    if (lane_begin >= lane_end_clamped) {
      return manifest;
    }

    std::vector<PendingPointWorkItem> pending_points_all;
    std::vector<size_t> pending_channel_indices_all;
    std::vector<std::vector<char>> pending_mask_by_channel(
        schedule_ctx.linear_states.size());
    pending_points_all.reserve(schedule_ctx.linear_states.size());
    pending_channel_indices_all.reserve(schedule_ctx.linear_states.size());
    for (size_t state_index = 0; state_index < schedule_ctx.linear_states.size();
         ++state_index) {
      const auto &state = schedule_ctx.linear_states[state_index];
      pending_mask_by_channel[state_index].assign(state.num_frequencies(), 0);
      bool state_has_pending = false;
      for (size_t freq_idx = 0; freq_idx < state.num_frequencies(); ++freq_idx) {
        if (point_needs_solving_fn(state, freq_idx)) {
          state_has_pending = true;
          PendingPointWorkItem pending_point{state_index, freq_idx};
          pending_points_all.push_back(pending_point);
          pending_mask_by_channel[state_index][freq_idx] = 1;
        }
      }
      if (state_has_pending) {
        pending_channel_indices_all.push_back(state_index);
      }
    }

    if (manifest.use_channel_series_mode) {
      std::vector<char> owned_states(schedule_ctx.linear_states.size(), 0);
      for (const auto state_index : pending_channel_indices_all) {
        const size_t owner_lane =
            schedule_ctx.owner_by_channel_index[state_index] %
            std::max<size_t>(1, active_groups);
        if (owner_lane >= lane_begin && owner_lane < lane_end_clamped) {
          manifest.pending_channel_indices.push_back(state_index);
          owned_states[state_index] = 1;
        }
      }
      for (const auto &pending_point : pending_points_all) {
        if (owned_states[pending_point.state_index]) {
          manifest.pending_points.push_back(pending_point);
        }
      }
    } else {
      // Point-mode strategy is protocol-dependent:
      // - Protocol 0: channel-partitioned contiguous frequency chains to
      //   maximize in-channel restart reuse.
      // - Later protocols: globally balanced distribution across all pending
      //   points for better wall-time balance near resonances.
      if (thresh_index == 0) {
        struct ChannelFrequencyBlock {
          size_t state_index = 0;
          size_t freq_begin = 0; // inclusive
          size_t freq_end = 0;   // inclusive
          size_t owner_lane = 0;
        };

        const size_t channel_count = pending_channel_indices_all.size();
        if (channel_count == 0) {
          return manifest;
        }

        // If owner lanes are fewer than pending channels, fall back to global
        // balancing to guarantee coverage of all channels.
        if (active_groups < channel_count) {
          std::vector<PendingPointWorkItem> globally_ordered_points;
          size_t max_frequency_count = 0;
          for (const auto state_index : pending_channel_indices_all) {
            max_frequency_count = std::max(
                max_frequency_count,
                pending_mask_by_channel[state_index].size());
          }
          for (size_t freq_idx = 0; freq_idx < max_frequency_count; ++freq_idx) {
            for (const auto state_index : pending_channel_indices_all) {
              if (freq_idx < pending_mask_by_channel[state_index].size() &&
                  pending_mask_by_channel[state_index][freq_idx]) {
                globally_ordered_points.push_back(
                    PendingPointWorkItem{state_index, freq_idx});
              }
            }
          }
          for (size_t point_index = 0; point_index < globally_ordered_points.size();
               ++point_index) {
            const size_t owner_lane = point_index % active_groups;
            if (owner_lane >= lane_begin && owner_lane < lane_end_clamped) {
              manifest.pending_points.push_back(globally_ordered_points[point_index]);
            }
          }
        } else {
          // Allocate point-mode lanes across pending channels (as evenly as
          // possible, weighted by pending-point counts).
          std::vector<size_t> pending_count_by_channel(channel_count, 0);
          for (size_t i = 0; i < channel_count; ++i) {
            const auto state_index = pending_channel_indices_all[i];
            for (const auto flag : pending_mask_by_channel[state_index]) {
              pending_count_by_channel[i] += flag ? 1 : 0;
            }
          }

          std::vector<size_t> groups_per_channel(channel_count, 1);
          size_t assigned_groups = channel_count;
          while (assigned_groups < active_groups) {
            size_t best_channel = 0;
            double best_score = -1.0;
            for (size_t i = 0; i < channel_count; ++i) {
              const double score = static_cast<double>(pending_count_by_channel[i]) /
                                   static_cast<double>(groups_per_channel[i]);
              if (score > best_score) {
                best_score = score;
                best_channel = i;
              }
            }
            ++groups_per_channel[best_channel];
            ++assigned_groups;
          }

          std::vector<size_t> channel_lane_begin(channel_count, 0);
          size_t lane_cursor = 0;
          for (size_t i = 0; i < channel_count; ++i) {
            channel_lane_begin[i] = lane_cursor;
            lane_cursor += groups_per_channel[i];
          }

          std::vector<ChannelFrequencyBlock> blocks;
          blocks.reserve(pending_points_all.size());

          for (size_t i = 0; i < channel_count; ++i) {
            const size_t state_index = pending_channel_indices_all[i];
            const auto &mask = pending_mask_by_channel[state_index];
            std::vector<size_t> pending_freq_indices;
            pending_freq_indices.reserve(mask.size());
            for (size_t freq_idx = 0; freq_idx < mask.size(); ++freq_idx) {
              if (mask[freq_idx]) {
                pending_freq_indices.push_back(freq_idx);
              }
            }
            if (pending_freq_indices.empty()) {
              continue;
            }

            const size_t chain_count =
                std::min(groups_per_channel[i], pending_freq_indices.size());
            const size_t target_chain_size =
                std::max<size_t>(1, (pending_freq_indices.size() + chain_count - 1) /
                                        chain_count);

            size_t run_start_pos = 0;
            size_t local_block_index = 0;
            while (run_start_pos < pending_freq_indices.size()) {
              size_t run_end_pos = run_start_pos;
              while (run_end_pos + 1 < pending_freq_indices.size() &&
                     pending_freq_indices[run_end_pos + 1] ==
                         pending_freq_indices[run_end_pos] + 1) {
                ++run_end_pos;
              }

              const size_t run_first_freq = pending_freq_indices[run_start_pos];
              const size_t run_last_freq = pending_freq_indices[run_end_pos];
              const size_t run_length = run_last_freq - run_first_freq + 1;

              for (size_t offset = 0; offset < run_length;
                   offset += target_chain_size) {
                const size_t block_first_freq = run_first_freq + offset;
                const size_t block_last_freq =
                    std::min(run_last_freq,
                             block_first_freq + target_chain_size - 1);
                const size_t owner_lane =
                    channel_lane_begin[i] +
                    (local_block_index % groups_per_channel[i]);
                blocks.push_back(ChannelFrequencyBlock{state_index,
                                                       block_first_freq,
                                                       block_last_freq,
                                                       owner_lane});
                ++local_block_index;
              }

              run_start_pos = run_end_pos + 1;
            }
          }

          for (const auto &block : blocks) {
            if (block.owner_lane < lane_begin || block.owner_lane >= lane_end_clamped) {
              continue;
            }
            for (size_t freq_idx = block.freq_begin; freq_idx <= block.freq_end;
                 ++freq_idx) {
              manifest.pending_points.push_back(
                  PendingPointWorkItem{block.state_index, freq_idx});
            }
          }
        }
      } else {
        std::vector<PendingPointWorkItem> globally_ordered_points;
        size_t max_frequency_count = 0;
        for (const auto state_index : pending_channel_indices_all) {
          max_frequency_count = std::max(max_frequency_count,
                                         pending_mask_by_channel[state_index].size());
        }
        for (size_t freq_idx = 0; freq_idx < max_frequency_count; ++freq_idx) {
          for (const auto state_index : pending_channel_indices_all) {
            if (freq_idx < pending_mask_by_channel[state_index].size() &&
                pending_mask_by_channel[state_index][freq_idx]) {
              globally_ordered_points.push_back(
                  PendingPointWorkItem{state_index, freq_idx});
            }
          }
        }
        for (size_t point_index = 0; point_index < globally_ordered_points.size();
             ++point_index) {
          const size_t owner_lane = point_index % active_groups;
          if (owner_lane >= lane_begin && owner_lane < lane_end_clamped) {
            manifest.pending_points.push_back(globally_ordered_points[point_index]);
          }
        }
      }

      for (const auto &pending_point : manifest.pending_points) {
        if (manifest.pending_channel_indices.empty() ||
            manifest.pending_channel_indices.back() != pending_point.state_index) {
          manifest.pending_channel_indices.push_back(pending_point.state_index);
        }
      }
    }
    return manifest;
  }

  static PendingProtocolManifest build_pending_manifest_from_metadata(
      const StateSolveScheduleContext &schedule_ctx, size_t thresh_index,
      size_t lane_begin, size_t lane_end, const nlohmann::json &metadata,
      bool at_final_protocol, bool force_retry_removed_frequencies) {
    auto needs_solving = [&](const LinearResponseDescriptor &state,
                             size_t freq_idx) {
      LinearResponsePoint pt{state, thresh_index, freq_idx};
      return point_needs_solving_from_metadata(
          metadata, pt, at_final_protocol, force_retry_removed_frequencies);
    };
    return build_pending_work_manifest(schedule_ctx, thresh_index, lane_begin,
                                       lane_end, needs_solving);
  }

  static void
  print_state_solve_execution_mode(World &world,
                                   const StateParallelPlan &state_parallel_plan,
                                   bool subgroup_parallel_requested,
                                   bool owner_group_schedule) {
    if (world.rank() == 0 && subgroup_parallel_requested) {
      print("STAGE2_EXECUTION_MODE mode=subgroup owner_mapping_groups=",
            state_parallel_plan.mapping_groups,
            " protocol0_owner_groups=",
            state_parallel_plan.channel_owner_groups,
            " strategy=owner_group_parallel_subworld");
    } else if (world.rank() == 0 && owner_group_schedule) {
      print("STAGE2_EXECUTION_MODE mode=serial_lanes owner_mapping_groups=",
            state_parallel_plan.mapping_groups,
            " protocol0_owner_groups=",
            state_parallel_plan.channel_owner_groups,
            " strategy=owner_group_serial_lanes");
    } else if (world.rank() == 0 && state_parallel_plan.mapping_groups > 1) {
      print("STAGE2_EXECUTION_MODE mode=serial_fallback owner_mapping_groups=",
            state_parallel_plan.mapping_groups,
            " strategy=plain_serial_channel_loop");
    }
  }

  /// Build runtime scheduling context for stage-2 linear solves.
  ///
  /// This combines static planner output with restart metadata so the solver can
  /// decide per protocol whether to run channel-series or channel-point mode.
  static StateSolveScheduleContext
  build_state_solve_schedule_context(World &world,
                                     const CalculationParameters &calc_params,
                                     const PlannedStates &planned_states,
                                     const ResponseParameters &response_params) {
    const auto &state_parallel_plan = planned_states.state_parallel_plan;
    const bool subgroup_parallel_requested =
        state_parallel_plan.subgroup_parallel_enabled &&
        state_parallel_plan.execution_groups > 1;
    const bool owner_group_schedule =
        state_parallel_plan.mapping_groups > 1 &&
        state_parallel_plan.effective_mode != "serial";
    print_state_solve_execution_mode(world, state_parallel_plan,
                                     subgroup_parallel_requested,
                                     owner_group_schedule);

    const auto &linear_states = planned_states.generated_states.states;
    auto owner_by_channel_index =
        build_owner_by_channel_index(linear_states, state_parallel_plan);
    const RuntimePointOwnershipPolicy runtime_policy =
        compute_runtime_point_ownership_policy(
            world, calc_params, state_parallel_plan, linear_states,
            owner_group_schedule, response_params);
    if (world.rank() == 0) {
      if (runtime_policy.force_retry_removed_frequencies) {
        print("RESTART_REMOVED_POLICY mode=force_retry");
      } else {
        print("RESTART_REMOVED_POLICY mode=skip_removed_unless_forced");
      }
      const auto protocol_count = calc_params.protocol().size();
      if (runtime_policy.restart_start_protocol_index < protocol_count) {
        print("RESTART_PROTOCOL_START protocol_index=",
              runtime_policy.restart_start_protocol_index,
              " threshold=",
              calc_params.protocol()[runtime_policy.restart_start_protocol_index],
              " strategy=earliest_missing_protocol");
      } else {
        print("RESTART_PROTOCOL_START protocol_index=none "
              "reason=no_pending_points");
      }
      if (runtime_policy.restart_final_protocol_only) {
        print("RESTART_FINAL_ONLY_ACTIVE remaining_final_pending_points=",
              runtime_policy.remaining_final_pending_points,
              " runtime_execution_groups=",
              runtime_policy.runtime_execution_groups,
              " runtime_effective_point_groups=",
              runtime_policy.runtime_effective_point_groups);
      } else if (runtime_policy.new_points_without_history > 0) {
        print("RESTART_NEW_POINTS count=",
              runtime_policy.new_points_without_history,
              " action=full_protocol_progression_for_new_points");
      }
    }
    if (world.rank() == 0 && !runtime_policy.protocol_policies.empty()) {
      for (size_t ti = 0; ti < runtime_policy.protocol_policies.size(); ++ti) {
        const auto &policy = runtime_policy.protocol_policies[ti];
        print("PROTOCOL_POLICY protocol_index=", ti, " mode=",
              policy.use_channel_series_mode ? "channel_series"
                                             : "channel_point",
              " active_groups=", policy.active_groups,
              " pending_channels=", policy.pending_channels,
              " pending_points=", policy.pending_points);
      }
      if (runtime_policy.protocol_policies[0].use_channel_series_mode &&
          runtime_policy.protocol_policies[0].pending_points > 0) {
        print("PROTOCOL0_CHANNEL_OWNERSHIP rows=",
              state_parallel_plan.channel_assignments.size());
        const size_t max_rows_to_print = 64;
        const size_t rows_to_print =
            std::min(max_rows_to_print, state_parallel_plan.channel_assignments.size());
        for (size_t i = 0; i < rows_to_print; ++i) {
          const auto &assignment = state_parallel_plan.channel_assignments[i];
          print("PROTOCOL0_CHANNEL_OWNER channel=", assignment.channel_label,
                " owner_group=", assignment.owner_group);
        }
        if (state_parallel_plan.channel_assignments.size() > rows_to_print) {
          print("PROTOCOL0_CHANNEL_OWNER_OMITTED count=",
                state_parallel_plan.channel_assignments.size() - rows_to_print,
                " reason=row_limit");
        }
      }
    }

    return StateSolveScheduleContext(
        subgroup_parallel_requested, owner_group_schedule, state_parallel_plan,
        linear_states, std::move(owner_by_channel_index),
        runtime_policy.runtime_effective_point_groups,
        runtime_policy.point_parallel_start_protocol_index,
        runtime_policy.runtime_execution_groups,
        runtime_policy.restart_start_protocol_index,
        runtime_policy.restart_final_protocol_only,
        runtime_policy.force_retry_removed_frequencies,
        runtime_policy.remaining_final_pending_points,
        runtime_policy.new_points_without_history,
        std::move(runtime_policy.protocol_policies));
  }

  // ============================================================================
  // SECTION 9: Stage 2b — Linear state solve (serial and subgroup paths)
  //   prepare_protocol_context — sets wavelet order + Coulomb/Fock operators
  //                              for a protocol threshold on a given world.
  //   log_pending_manifest — prints pending-work summary (with subgroup label).
  //   execute_manifest_work — dispatches channel-series or channel-point work
  //                           via caller-supplied solve lambdas.
  //   cached_or_built_manifest — lazy manifest cache to avoid redundant builds.
  //   execute_serial_state_solve — universe-communicator protocol loop;
  //                                initialises JsonStateSolvePersistence once
  //                                and iterates over pending manifests.
  //   execute_subgroup_state_solve — MacroTaskQ subworld protocol loop;
  //                                  each subgroup owns a metadata shard and
  //                                  merges at barriers; includes tail derived
  //                                  polling for subgroup-derived hand-off.
  // ============================================================================

  /// Execute stage-2 linear solves on the universe communicator.
  ///
  /// Used as the deterministic baseline path and as fallback when subgroup mode
  /// is unavailable or fails.
  static void prepare_protocol_context(World &exec_world,
                                       ResponseManager &response_manager,
                                       GroundStateData &ground_state,
                                       double threshold,
                                       const std::string &fock_json_file) {
    response_manager.setProtocol(exec_world, ground_state.getL(), threshold);
    ground_state.prepareOrbitals(exec_world, FunctionDefaults<3>::get_k(),
                                 threshold);
    ground_state.computePreliminaries(exec_world,
                                      *response_manager.getCoulombOp(),
                                      response_manager.getVtol(),
                                      fock_json_file);
  }

  static void log_pending_manifest(World &exec_world, size_t protocol_index,
                                   const PendingProtocolManifest &manifest,
                                   const std::optional<size_t> subgroup_id =
                                       std::nullopt) {
    if (exec_world.rank() != 0) {
      return;
    }
    if (subgroup_id.has_value()) {
      print("PENDING_MANIFEST");
      print("  scope = subgroup");
      print("  subgroup = ", *subgroup_id);
      print("  protocol_index = ", protocol_index);
      print("  mode = ",
            manifest.use_channel_series_mode ? "channel_series"
                                             : "channel_point");
      print("  pending_channels = ",
            manifest.pending_channel_indices.size());
      print("  pending_points = ", manifest.pending_points.size());
      print("------------------------------------------------------------");
      return;
    }

    print("PENDING_MANIFEST");
    print("  scope = world");
    print("  protocol_index = ", protocol_index);
    print("  mode = ",
          manifest.use_channel_series_mode ? "channel_series" : "channel_point");
    print("  pending_channels = ",
          manifest.pending_channel_indices.size());
    print("  pending_points = ", manifest.pending_points.size());
    print("------------------------------------------------------------");
  }

  template <typename SolveStateFn, typename SolvePointFn>
  static void execute_manifest_work(const PendingProtocolManifest &manifest,
                                    SolveStateFn &&solve_state,
                                    SolvePointFn &&solve_state_frequency) {
    if (manifest.use_channel_series_mode) {
      for (const auto state_index : manifest.pending_channel_indices) {
        solve_state(state_index);
      }
      return;
    }
    for (const auto &work_item : manifest.pending_points) {
      solve_state_frequency(work_item.state_index, work_item.freq_index);
    }
  }

  template <typename BuildManifestFn>
  static const PendingProtocolManifest &cached_or_built_manifest(
      std::vector<std::optional<PendingProtocolManifest>> &manifest_cache,
      size_t protocol_index, BuildManifestFn &&build_manifest) {
    if (!manifest_cache[protocol_index].has_value()) {
      manifest_cache[protocol_index] = build_manifest();
    }
    return *manifest_cache[protocol_index];
  }

  static void execute_serial_state_solve(
      World &world, const CalculationParameters &calc_params, GroundContext &ctx,
      const StateSolveScheduleContext &schedule_ctx,
      const ExcitedStateBundlePlan &excited_state_bundle_plan,
      nlohmann::json &state_metadata_json, nlohmann::json &debug_log_json) {
    // Serial execution path:
    // 1) initialize metadata/log persistence,
    // 2) iterate protocol thresholds and solve pending work,
    // 3) emit in-memory metadata/debug JSON.
    JsonStateSolvePersistence persistence(world, "response_metadata.json",
                                          "response_log.json",
                                          nlohmann::json::object(),
                                          schedule_ctx.force_retry_removed_frequencies);
    persistence.initialize_states(schedule_ctx.linear_states);
    persistence.initialize_excited_bundle(excited_state_bundle_plan);

    if (world.rank() == 0) {
      persistence.print_summary();
    }
    world.gop.fence();

    const auto &protocol = calc_params.protocol();
    std::vector<std::optional<PendingProtocolManifest>> pending_manifest_by_ti(
        protocol.size());

    auto needs_solving_at_protocol = [&](double protocol_thresh,
                                         size_t thresh_index) {
      // A protocol threshold is considered "active" when at least one point is
      // missing on disk, or not converged at the final threshold.
      const bool at_final_protocol = protocol_thresh == calc_params.protocol().back();
      if (thresh_index < schedule_ctx.restart_start_protocol_index) {
        return false;
      }
      auto serial_point_needs_solving = [&](const LinearResponseDescriptor &state,
                                            size_t freq_idx) {
        LinearResponsePoint pt{state, thresh_index, freq_idx};
        const bool is_saved = persistence.is_saved(pt);
        const bool should_solve =
            point_needs_solving(persistence, pt, at_final_protocol,
                                schedule_ctx.force_retry_removed_frequencies);
        if (world.rank() == 0) {
          print("POINT_SOLVE_CHECK state=", pt.perturbationDescription(),
                " threshold=", protocol_thresh, " frequency=", pt.frequency(),
                " is_saved=", is_saved,
                " at_final_protocol=", at_final_protocol,
                " should_solve=", should_solve);
        }
        return should_solve;
      };

      if (!schedule_ctx.owner_group_schedule) {
        std::vector<size_t> all_state_indices(schedule_ctx.linear_states.size());
        std::iota(all_state_indices.begin(), all_state_indices.end(), 0);
        return any_state_point_needs_solving(all_state_indices, schedule_ctx,
                                             thresh_index,
                                             serial_point_needs_solving);
      }

      const size_t lane_count =
          active_owner_groups_for_protocol_runtime(schedule_ctx, thresh_index);
      const auto &manifest = cached_or_built_manifest(
          pending_manifest_by_ti, thresh_index, [&]() {
            return build_pending_work_manifest(
                schedule_ctx, thresh_index, 0, lane_count,
                serial_point_needs_solving);
          });
      return manifest.has_work();
    };

    auto solve_state = [&](size_t state_index, size_t thresh_index,
                           bool at_final_protocol) {
      auto &state = schedule_ctx.linear_states[state_index];
      run_frequency_loop_with_flush(world, ctx.response_manager, state,
                                    thresh_index, ctx.ground, persistence,
                                    at_final_protocol);
    };
    auto solve_state_frequency = [&](size_t state_index, size_t freq_index,
                                     size_t thresh_index,
                                     bool at_final_protocol) {
      const auto &state = schedule_ctx.linear_states[state_index];
      computeFrequencyPoint(world, ctx.response_manager, state, thresh_index,
                            freq_index, ctx.ground, persistence,
                            at_final_protocol);
      persistence.flush_debug_log(world);
    };
    run_protocol_threshold_loop(
        protocol, needs_solving_at_protocol,
        [&](double thresh, size_t /*thresh_index*/) {
          if (world.rank() == 0) {
            madness::print("PROTOCOL_SKIP protocol=", thresh,
                           " reason=all_states_converged");
          }
        },
        [&](double thresh, size_t /*thresh_index*/) {
          prepare_protocol_context(world, ctx.response_manager, ctx.ground, thresh,
                                   ctx.fock_json_file);
        },
        [&](double thresh, size_t ti, bool at_final_protocol) {
          (void)thresh;
          if (!schedule_ctx.owner_group_schedule) {
            for (size_t state_index = 0;
                 state_index < schedule_ctx.linear_states.size(); ++state_index) {
              solve_state(state_index, ti, at_final_protocol);
            }
            return;
          }

          const size_t lane_count =
              active_owner_groups_for_protocol_runtime(schedule_ctx, ti);
          const auto &manifest = cached_or_built_manifest(
              pending_manifest_by_ti, ti, [&]() {
                auto serial_point_needs_solving =
                    [&](const LinearResponseDescriptor &state, size_t freq_idx) {
                      LinearResponsePoint pt{state, ti, freq_idx};
                      return point_needs_solving(persistence, pt,
                                                 at_final_protocol,
                                                 schedule_ctx
                                                     .force_retry_removed_frequencies);
                    };
                return build_pending_work_manifest(
                    schedule_ctx, ti, 0, lane_count,
                    serial_point_needs_solving);
              });
          log_pending_manifest(world, ti, manifest);
          execute_manifest_work(
              manifest,
              [&](size_t state_index) {
                solve_state(state_index, ti, at_final_protocol);
              },
              [&](size_t state_index, size_t freq_index) {
                solve_state_frequency(state_index, freq_index, ti,
                                      at_final_protocol);
              });
        },
        [&](double /*thresh*/, size_t /*thresh_index*/) {});

    state_metadata_json = persistence.metadata_json();
    debug_log_json = persistence.debug_log_json();
  }

  /// Execute stage-2 linear solves in macrotask subgroups.
  ///
  /// Returns `true` when subgroup execution completed successfully and produced
  /// merged metadata/debug output. Returns `false` when caller should fall back
  /// to serial execution.
  static bool execute_subgroup_state_solve(
      World &world, const CalculationParameters &calc_params, GroundContext &ctx,
      const StateSolveScheduleContext &schedule_ctx,
      const PlannedStates &planned_states,
      const ExcitedStateBundlePlan &excited_state_bundle_plan,
      nlohmann::json &state_metadata_json, nlohmann::json &debug_log_json) {
    // Subgroup execution path:
    // 1) create macrotask subworlds,
    // 2) solve local ownership shards in each subgroup,
    // 3) merge subgroup metadata/debug shards on world rank 0 and broadcast.
    const auto &derived_state_plan = planned_states.derived_state_plan;
    const auto &generated_states = planned_states.generated_states;
    const size_t runtime_execution_groups =
        std::max<size_t>(1, schedule_ctx.runtime_execution_groups);
    const auto &protocol = calc_params.protocol();
    const bool have_protocol = !protocol.empty();
    const size_t final_protocol_index = have_protocol ? (protocol.size() - 1) : 0;
    const double final_threshold = have_protocol ? protocol.back()
                                                 : FunctionDefaults<3>::get_thresh();
    nlohmann::json baseline_state_metadata = nlohmann::json::object();
    if (world.rank() == 0) {
      baseline_state_metadata = read_json_file_or_object("response_metadata.json");
    }
    broadcast_json(world, baseline_state_metadata);
    std::string tail_derived_claim_prefix;
    if (world.rank() == 0 && !derived_state_plan.requests.empty()) {
      tail_derived_claim_prefix =
          (std::filesystem::path("derived_request_claims") /
           derived_request_manifest_scope(final_protocol_index, final_threshold) /
           ("run_" + iso_timestamp()))
              .string();
    }
    world.gop.broadcast_serializable(tail_derived_claim_prefix, 0);
    const double subgroup_stage_wall_start = madness::wall_time();
    double last_global_manifest_sync_wall = subgroup_stage_wall_start;
    constexpr double k_manifest_sync_poll_interval_s = 60.0;

    auto sync_global_manifests_and_progress = [&](bool force_sync) {
      if (world.rank() != 0) {
        return;
      }
      const double now = madness::wall_time();
      if (!force_sync &&
          (now - last_global_manifest_sync_wall) <
              k_manifest_sync_poll_interval_s) {
        return;
      }

      // Start from the pre-run metadata snapshot so subgroup shard defaults
      // ("saved=false/converged=false" placeholders) do not erase completed
      // points during restart-only runs with no new linear solves.
      nlohmann::json merged_metadata = baseline_state_metadata;
      if (!merged_metadata.is_object()) {
        merged_metadata = nlohmann::json::object();
      }
      for (size_t gid = 0; gid < runtime_execution_groups; ++gid) {
        merge_state_metadata_json(
            merged_metadata,
            read_json_file_or_object(group_shard_file("response_metadata.json",
                                                      gid)));
      }
      write_json_file("response_metadata.json", merged_metadata);

      const size_t protocol_count = protocol.size();
      size_t total_points = 0;
      size_t removed_points = 0;
      size_t active_points = 0;
      size_t final_converged_points = 0;
      if (protocol_count > 0) {
        const size_t final_protocol_index = protocol_count - 1;
        for (const auto &state : schedule_ctx.linear_states) {
          for (size_t freq_idx = 0; freq_idx < state.num_frequencies();
               ++freq_idx) {
            ++total_points;
            LinearResponsePoint protocol0_pt{state, /*thresh_index=*/0,
                                             freq_idx};
            const bool removed = point_marked_for_frequency_removal_in_metadata(
                merged_metadata, protocol0_pt);
            if (removed) {
              ++removed_points;
            }
            if (!schedule_ctx.force_retry_removed_frequencies && removed) {
              continue;
            }
            ++active_points;
            LinearResponsePoint final_pt{state, final_protocol_index, freq_idx};
            if (point_ready_in_metadata(merged_metadata, final_pt,
                                        /*require_saved=*/true,
                                        /*require_converged=*/true)) {
              ++final_converged_points;
            }
          }
        }
      }

      const size_t final_pending_points =
          (active_points > final_converged_points)
              ? (active_points - final_converged_points)
              : 0;

      if (world.rank() == 0) {
        print("STAGE2_PROGRESS_POLL");
        print("  elapsed_s = ", now - subgroup_stage_wall_start);
        print("  total_points = ", total_points);
        print("  active_points = ", active_points);
        print("  removed_points = ", removed_points);
        print("  final_converged_points = ", final_converged_points);
        print("  final_pending_points = ", final_pending_points);
        for (size_t ti = 0; ti < protocol_count; ++ti) {
          const bool at_final_protocol = (ti + 1 == protocol_count);
          const size_t lane_count =
              active_owner_groups_for_protocol_runtime(schedule_ctx, ti);
          const auto manifest = build_pending_manifest_from_metadata(
              schedule_ctx, ti, 0, lane_count, merged_metadata,
              at_final_protocol, schedule_ctx.force_retry_removed_frequencies);
          print("  protocol[", ti, "] threshold=", protocol[ti], " mode=",
                manifest.use_channel_series_mode ? "channel_series"
                                                 : "channel_point",
                " pending_channels=",
                manifest.pending_channel_indices.size(),
                " pending_points=", manifest.pending_points.size());
        }
        print("------------------------------------------------------------");
      }
      last_global_manifest_sync_wall = now;
    };
    const bool previous_console_enabled = TimedValueLogger::console_enabled();
    TimedValueLogger::set_console_enabled(false);
    auto restore_console_setting = [&]() {
      TimedValueLogger::set_console_enabled(previous_console_enabled);
    };
    try {
      auto subworld_ptr = MacroTaskQ::create_worlds(
          world, runtime_execution_groups);
      if (!subworld_ptr) {
        throw std::runtime_error("subworld creation returned null");
      }
      World &subworld = *subworld_ptr;
      const size_t subgroup_id = static_cast<size_t>(
          world.rank() % static_cast<int>(runtime_execution_groups));
      // Redirect every subgroup root (including subgroup 0) into its shard log.
      // This keeps main stdout limited to world-level orchestration messages.
      ScopedRankLogRedirect subgroup_console_redirect(
          subworld.rank() == 0, group_console_file(subgroup_id),
          "stage2-linear subgroup=" + std::to_string(subgroup_id));

      auto old_pmap3 = FunctionDefaults<3>::get_pmap();
      auto restore_pmap = [&]() { FunctionDefaults<3>::set_pmap(old_pmap3); };

      // Molresponse state solves are strictly 3D; only swap the 3D default
      // pmap.
      FunctionDefaults<3>::set_default_pmap(subworld);
      try {
        std::vector<size_t> local_channel_indices;
        std::vector<LinearResponseDescriptor> local_channels;
        build_local_channel_workset(
            schedule_ctx.linear_states, schedule_ctx.owner_by_channel_index,
            schedule_ctx.point_scheduler, subgroup_id, runtime_execution_groups,
            local_channel_indices, local_channels);

        const std::string metadata_shard_file =
            group_shard_file("response_metadata.json", subgroup_id);
        const std::string debug_shard_file =
            group_shard_file("response_log.json", subgroup_id);
        const std::string fock_shard_file =
            group_shard_file(ctx.fock_json_file, subgroup_id);

        JsonStateSolvePersistence local_persistence(
            subworld, metadata_shard_file, debug_shard_file,
            baseline_state_metadata,
            schedule_ctx.force_retry_removed_frequencies,
            [&](bool force_sync) {
              sync_global_manifests_and_progress(force_sync);
            });
        local_persistence.initialize_states(local_channels);
        local_persistence.initialize_excited_bundle(excited_state_bundle_plan);
        if (subworld.rank() == 0) {
          print("SUBGROUP_OWNERSHIP subgroup=", subgroup_id,
                " protocol0_channels=", local_channel_indices.size(),
                " active_channels_all_protocols=", local_channels.size());
          local_persistence.print_summary();
        }
        subworld.gop.fence();

        size_t protocols_with_work = 0;
        size_t dispatched_channel_batches = 0;
        size_t dispatched_point_solves = 0;
        const bool has_local_linear_work = !local_channels.empty();

        if (has_local_linear_work) {
          GroundStateData local_ground(subworld, ctx.archive_file, ctx.molecule);
          ResponseManager local_response_manager(subworld, calc_params);
          std::vector<std::optional<PendingProtocolManifest>>
              local_pending_manifest_by_ti(protocol.size());
          auto build_local_protocol_manifest =
              [&](size_t thresh_index, bool at_final_protocol) {
                auto local_point_needs_solving =
                    [&](const LinearResponseDescriptor &state,
                        size_t freq_idx) {
                      LinearResponsePoint pt{state, thresh_index, freq_idx};
                      return point_needs_solving(
                          local_persistence, pt, at_final_protocol,
                          schedule_ctx.force_retry_removed_frequencies);
                    };
                return build_pending_work_manifest(
                    schedule_ctx, thresh_index, subgroup_id, subgroup_id + 1,
                    local_point_needs_solving);
              };
          auto local_needs_solving_at_protocol = [&](double protocol_thresh,
                                                     size_t thresh_index) {
            const bool at_final_protocol = protocol_thresh == protocol.back();
            if (thresh_index < schedule_ctx.restart_start_protocol_index) {
              return false;
            }
            const auto &manifest = cached_or_built_manifest(
                local_pending_manifest_by_ti, thresh_index, [&]() {
                  return build_local_protocol_manifest(thresh_index,
                                                       at_final_protocol);
                });
            return manifest.has_work();
          };

          run_protocol_threshold_loop(
              protocol, local_needs_solving_at_protocol,
              [&](double thresh, size_t /*thresh_index*/) {
                if (subworld.rank() == 0) {
                  print("SUBGROUP_PROTOCOL_SKIP subgroup=", subgroup_id,
                        " threshold=", thresh,
                        " reason=no_pending_channels");
                }
                // Intentionally no global per-protocol fence: each subgroup can
                // advance independently so stalled channels do not block all
                // other available work.
              },
              [&](double thresh, size_t /*thresh_index*/) {
                prepare_protocol_context(subworld, local_response_manager,
                                         local_ground, thresh, fock_shard_file);
              },
              [&](double thresh, size_t ti, bool at_final_protocol) {
                const auto &manifest = cached_or_built_manifest(
                    local_pending_manifest_by_ti, ti, [&]() {
                      return build_local_protocol_manifest(ti,
                                                           at_final_protocol);
                    });
                log_pending_manifest(subworld, ti, manifest, subgroup_id);
                if (manifest.has_work()) {
                  ++protocols_with_work;
                  dispatched_channel_batches +=
                      manifest.pending_channel_indices.size();
                  dispatched_point_solves += manifest.pending_points.size();
                }
                if (!manifest.use_channel_series_mode && subworld.rank() == 0) {
                  print("SUBGROUP_POINT_MODE");
                  print("  subgroup = ", subgroup_id);
                  print("  threshold = ", thresh);
                  print("  mode = independent_state_frequency_points");
                  print("------------------------------------------------------------");
                }
                execute_manifest_work(
                    manifest,
                    [&](size_t state_index) {
                      auto &state = schedule_ctx.linear_states[state_index];
                      run_frequency_loop_with_flush(
                          subworld, local_response_manager, state, ti,
                          local_ground, local_persistence, at_final_protocol);
                    },
                    [&](size_t state_index, size_t freq_index) {
                      const auto &state = schedule_ctx.linear_states[state_index];
                      computeFrequencyPoint(
                          subworld, local_response_manager, state, ti, freq_index,
                          local_ground, local_persistence, at_final_protocol);
                      local_persistence.flush_debug_log(subworld);
                    });
              },
              [&](double /*thresh*/, size_t /*thresh_index*/) {
                local_persistence.maybe_poll_progress(subworld, /*force=*/true);
              });
        } else if (subworld.rank() == 0) {
          print("SUBGROUP_STAGE2_NO_OWNED_WORK subgroup=", subgroup_id);
        }

        size_t tail_claimed_requests = 0;
        size_t tail_completed_requests = 0;
        size_t tail_failed_requests = 0;
        size_t tail_skipped_done_requests = 0;
        size_t tail_idle_polls = 0;
        if (have_protocol && !derived_state_plan.requests.empty()) {
          std::vector<char> done_seen(derived_state_plan.requests.size(), 0);
          std::unique_ptr<GroundStateData> tail_ground;
          std::unique_ptr<ResponseManager> tail_response_manager;
          std::unique_ptr<SimpleVBCComputer> tail_vbc_computer;
          bool tail_context_ready = false;
          auto ensure_tail_context = [&]() {
            if (tail_context_ready) {
              return;
            }
            tail_ground = std::make_unique<GroundStateData>(
                subworld, ctx.archive_file, ctx.molecule);
            tail_response_manager =
                std::make_unique<ResponseManager>(subworld, calc_params);
            tail_response_manager->setProtocol(subworld, tail_ground->getL(),
                                               final_threshold);
            tail_ground->prepareOrbitals(subworld, FunctionDefaults<3>::get_k(),
                                         final_threshold);
            tail_ground->computePreliminaries(
                subworld, *tail_response_manager->getCoulombOp(),
                tail_response_manager->getVtol(), fock_shard_file);
            tail_vbc_computer =
                std::make_unique<SimpleVBCComputer>(subworld, *tail_ground);
            tail_context_ready = true;
            if (subworld.rank() == 0) {
              print("SUBGROUP_TAIL_DERIVED_CONTEXT_READY subgroup=", subgroup_id,
                    " final_threshold=", final_threshold);
            }
          };

          constexpr size_t k_tail_max_idle_polls = 720;
          constexpr std::chrono::seconds k_tail_poll_sleep{5};

          // TODO(step4): extract execute_subgroup_tail_derived_poll
          while (true) {
            size_t final_pending_linear_points = 0;
            std::vector<size_t> ready_request_indices;
            if (subworld.rank() == 0) {
              // Preserve baseline solved-point status while layering subgroup
              // shard updates for this poll iteration.
              nlohmann::json merged_metadata = baseline_state_metadata;
              if (!merged_metadata.is_object()) {
                merged_metadata = nlohmann::json::object();
              }
              for (size_t gid = 0; gid < runtime_execution_groups; ++gid) {
                merge_state_metadata_json(
                    merged_metadata,
                    read_json_file_or_object(
                        group_shard_file("response_metadata.json", gid)));
              }
              const size_t lane_count = active_owner_groups_for_protocol_runtime(
                  schedule_ctx, final_protocol_index);
              const auto final_manifest = build_pending_manifest_from_metadata(
                  schedule_ctx, final_protocol_index, 0, lane_count,
                  merged_metadata,
                  /*at_final_protocol=*/true,
                  schedule_ctx.force_retry_removed_frequencies);
              final_pending_linear_points = final_manifest.pending_points.size();

              const auto derived_gate = DerivedStatePlanner::evaluate_dependency_gate(
                  derived_state_plan, generated_states, final_protocol_index,
                  [&](const LinearResponsePoint &pt) {
                    return point_ready_in_metadata(merged_metadata, pt,
                                                   /*require_saved=*/true,
                                                   /*require_converged=*/true);
                  });
              ready_request_indices.reserve(derived_gate.ready_requests);
              for (size_t req_index = 0;
                   req_index < derived_state_plan.requests.size(); ++req_index) {
                if (req_index >= derived_gate.entries.size() ||
                    !derived_gate.entries[req_index].ready) {
                  continue;
                }
                const auto &req = derived_state_plan.requests[req_index];
                if (derived_request_done_record_exists(
                        final_protocol_index, final_threshold, req,
                        req_index)) {
                  if (!done_seen[req_index]) {
                    done_seen[req_index] = 1;
                    ++tail_skipped_done_requests;
                  }
                  continue;
                }
                ready_request_indices.push_back(req_index);
              }
            }
            subworld.gop.broadcast_serializable(final_pending_linear_points, 0);
            subworld.gop.broadcast_serializable(ready_request_indices, 0);

            bool claimed_any_request = false;
            for (const auto req_index : ready_request_indices) {
              bool request_claimed = false;
              if (subworld.rank() == 0) {
                const auto &req = derived_state_plan.requests[req_index];
                if (!derived_request_done_record_exists(
                        final_protocol_index, final_threshold, req, req_index)) {
                  std::ostringstream claim_payload;
                  claim_payload << "subgroup=" << subgroup_id
                                << ";request_index=" << req_index
                                << ";derived_state_id=" << req.derived_state_id;
                  request_claimed = try_claim_property_component_task(
                      derived_request_claim_file(tail_derived_claim_prefix,
                                                 req_index),
                      claim_payload.str());
                }
              }
              subworld.gop.broadcast_serializable(request_claimed, 0);
              if (!request_claimed) {
                continue;
              }

              claimed_any_request = true;
              ensure_tail_context();
              long local_completed = 0;
              long local_failed = 0;
              const auto &req = derived_state_plan.requests[req_index];
              const auto timing =
                  run_derived_request(subworld, *tail_ground, req,
                                      *tail_vbc_computer, local_completed,
                                      local_failed);
              if (subworld.rank() == 0) {
                ++tail_claimed_requests;
                if (timing.success) {
                  ++tail_completed_requests;
                  write_derived_request_done_record(
                      final_protocol_index, final_threshold, req, req_index,
                      subgroup_id, timing.wall_seconds, timing.cpu_seconds);
                } else {
                  ++tail_failed_requests;
                }
              }
            }

            if (claimed_any_request) {
              tail_idle_polls = 0;
              continue;
            }

            ++tail_idle_polls;
            const bool final_linear_done = final_pending_linear_points == 0;
            if (final_linear_done) {
              break;
            }
            if (tail_idle_polls >= k_tail_max_idle_polls) {
              if (subworld.rank() == 0) {
                print("WARN SUBGROUP_TAIL_DERIVED_IDLE_TIMEOUT subgroup=",
                      subgroup_id, " idle_polls=", tail_idle_polls,
                      " pending_linear_points=", final_pending_linear_points);
              }
              break;
            }
            std::this_thread::sleep_for(k_tail_poll_sleep);
          }
        }

        if (subworld.rank() == 0) {
          print("SUBGROUP_STAGE2_COMPLETE subgroup=", subgroup_id,
                " protocols_with_work=", protocols_with_work,
                " dispatched_channels=", dispatched_channel_batches,
                " dispatched_points=", dispatched_point_solves,
                " tail_derived_claimed=", tail_claimed_requests,
                " tail_derived_completed=", tail_completed_requests,
                " tail_derived_failed=", tail_failed_requests,
                " tail_derived_done_skipped=", tail_skipped_done_requests);
        }

        local_persistence.maybe_poll_progress(subworld, /*force=*/true);
        local_persistence.flush_debug_log(subworld);
        subworld.gop.fence();
        restore_pmap();
      } catch (...) {
        restore_pmap();
        throw;
      }

      world.gop.fence();
      if (world.rank() == 0) {
        nlohmann::json merged_metadata = baseline_state_metadata;
        if (!merged_metadata.is_object()) {
          merged_metadata = nlohmann::json::object();
        }
        nlohmann::json merged_debug_log = nlohmann::json::object();
        for (size_t gid = 0; gid < runtime_execution_groups; ++gid) {
          const std::string metadata_file =
              group_shard_file("response_metadata.json", gid);
          const std::string debug_file =
              group_shard_file("response_log.json", gid);
          merge_state_metadata_json(merged_metadata,
                                    read_json_file_or_object(metadata_file));
          merge_debug_log_json(merged_debug_log,
                               read_json_file_or_object(debug_file));
        }

        write_json_file("response_metadata.json", merged_metadata);
        write_json_file("response_log.json", merged_debug_log);
        state_metadata_json = std::move(merged_metadata);
        debug_log_json = std::move(merged_debug_log);
      }

      std::string metadata_dump;
      std::string debug_dump;
      if (world.rank() == 0) {
        metadata_dump = state_metadata_json.dump();
        debug_dump = debug_log_json.dump();
      }
      world.gop.broadcast_serializable(metadata_dump, 0);
      world.gop.broadcast_serializable(debug_dump, 0);
      if (world.rank() != 0) {
        state_metadata_json = metadata_dump.empty()
                                  ? nlohmann::json::object()
                                  : nlohmann::json::parse(metadata_dump);
        debug_log_json = debug_dump.empty()
                             ? nlohmann::json::object()
                             : nlohmann::json::parse(debug_dump);
      }
      world.gop.fence();

      restore_console_setting();
      return true;
    } catch (const std::exception &ex) {
      restore_console_setting();
      if (world.rank() == 0) {
        print("ERROR STAGE2_SUBGROUP_FAILED message=", ex.what(),
              " action=fallback_serial_state_loop");
      }
      return false;
    }
  }

  // ============================================================================
  // SECTION 10: Stage 2d — Derived state execution
  //   DerivedRequestTiming — wall/cpu result for one VBC request.
  //   run_derived_request — single-request runner shared by subgroup and
  //                         serial fallback paths; delegates to
  //                         SimpleVBCComputer::compute_and_save.
  //   DerivedExecutionResult — dependency gate + execution summary bundle.
  //   execute_derived_state_requests — evaluates dependency gate, partitions
  //                         ready requests across owner lanes, and runs them in
  //                         subgroup subworlds (with serial fallback).
  // ============================================================================

  struct DerivedRequestTiming {
    bool success = false;
    double wall_seconds = 0.0;
    double cpu_seconds = 0.0;
  };

  static DerivedRequestTiming
  run_derived_request(World &exec_world, const GroundStateData &ground,
                      const DerivedStateRequest &req,
                      SimpleVBCComputer &vbc_computer, long &completed,
                      long &failed) {
    // Single derived-request runner used by both subgroup and serial fallback
    // execution paths.
    const double start_wall = madness::wall_time();
    const double start_cpu = madness::cpu_time();
    DerivedRequestTiming timing;
    try {
      VBCResponseState vbc_state =
          DerivedStatePlanner::make_vbc_state(req, ground.isSpinRestricted());
      vbc_computer.compute_and_save(vbc_state);
      timing.success = true;
      if (exec_world.rank() == 0) {
        ++completed;
      }
    } catch (const std::exception &ex) {
      if (exec_world.rank() == 0) {
        ++failed;
        print("ERROR DERIVED_REQUEST_FAILED derived_state_id=",
              req.derived_state_id, " message=", ex.what());
      }
    }
    timing.wall_seconds = madness::wall_time() - start_wall;
    timing.cpu_seconds = madness::cpu_time() - start_cpu;
    return timing;
  }

  struct DerivedExecutionResult {
    // Dependency readiness report for all derived requests.
    DerivedStateGateReport dependency_gate;
    // Execution summary JSON written into response metadata.
    nlohmann::json execution;
  };

  // ============================================================================
  // SECTION 11: Stage 2c — Excited-state bundle execution
  //   ExcitedExecutionResult — execution summary JSON produced by stage 2c.
  //   FinalProtocolState — final threshold + convergence status forwarded to
  //                        stage 3.
  //   ensure_excited_protocol_placeholder_node — fills missing keys in the
  //                        metadata excited-protocol node with defaults.
  //   build_excited_root_manifest — serialises ExcitedRootDescriptor list.
  //   execute_excited_state_bundle_stage — protocol loop that dispatches
  //                        ExcitedResponse::solve_protocol for pending
  //                        protocols, records results via ResponseRecord2,
  //                        and merges metadata.
  // ============================================================================

  struct ExcitedExecutionResult {
    // Execution summary JSON written into response metadata.
    nlohmann::json execution;
  };

  struct FinalProtocolState {
    double threshold = 0.0;
    size_t threshold_index = 0;
    bool all_linear_points_converged = true;
  };

  static void ensure_excited_protocol_placeholder_node(nlohmann::json &node) {
    static const nlohmann::json kDefaultExcitedProtocolNode = []() {
      ExcitedProtocolResult defaults;
      nlohmann::json j;
      to_json(j, defaults);
      j["stage_status"] = "placeholder_pending_solver";
      j["timings"] = {{"wall_seconds", 0.0}, {"cpu_seconds", 0.0}};
      return j;
    }();
    if (!node.is_object()) node = nlohmann::json::object();
    for (auto &[key, val] : kDefaultExcitedProtocolNode.items()) {
      if (!node.contains(key)) node[key] = val;
    }
  }

  static nlohmann::json build_excited_root_manifest(
      const std::vector<ExcitedRootDescriptor> &roots) {
    return nlohmann::json(roots);
  }

  static ExcitedExecutionResult execute_excited_state_bundle_stage(
      World &world, const PlannedStates &planned_states,
      const GroundContext &ground_ctx, const ResponseParameters &response_params,
      nlohmann::json &state_metadata_json) {
    const auto &plan = planned_states.excited_state_bundle_plan;
    ResponseRecord2 excited_metadata_record(world, "response_metadata.json");
    excited_metadata_record.initialize_excited_bundle(
        plan.enabled, plan.num_states, plan.tda, plan.guess_max_iter,
        plan.maxiter, plan.maxsub, plan.owner_group, plan.protocols);
    state_metadata_json = excited_metadata_record.to_json();

    ExcitedSolverConfig solver_config;
    solver_config.archive_file = ground_ctx.archive_file;
    solver_config.output_prefix = response_params.prefix();
    solver_config.protocols = plan.protocols;
    solver_config.print_level = response_params.print_level();
    ExcitedResponse excited_solver(solver_config);
    const std::string solver_adapter_name = "ExcitedResponse";

    nlohmann::json execution = {
        {"attempted", false},
        {"mode", plan.enabled ? "ExcitedResponse" : "disabled"},
        {"solver_adapter", solver_adapter_name},
        {"enabled", plan.enabled},
        {"owner_group", plan.owner_group},
        {"protocol_count", plan.protocols.size()},
        {"ready_protocol_placeholders", 0},
        {"restart_ready_protocols", 0},
        {"pending_protocols", 0},
        {"skipped_protocols", 0},
        {"completed_protocols", 0},
        {"failed_protocols", 0},
        {"total_wall_seconds", 0.0},
        {"total_cpu_seconds", 0.0},
        {"protocol_events", nlohmann::json::array()}};

    auto protocol_result_to_json = [](const ExcitedProtocolResult &result) {
      return nlohmann::json{
          {"attempted", result.attempted},
          {"saved", result.saved},
          {"converged", result.converged},
          {"failed", result.failed},
          {"skipped", result.skipped},
          {"restart_reused", result.restart_reused},
          {"stage_status", result.stage_status},
          {"response_variant", result.response_variant},
          {"restart_support_mode", result.restart_support_mode},
          {"restart_source", result.restart_source},
          {"snapshot_kind", result.snapshot_kind},
          {"bundle_state_present", result.bundle_state_present},
          {"restart_capable", result.restart_capable},
          {"restart_source_threshold", result.restart_source_threshold},
          {"iterations", result.iterations},
          {"energies", result.energies},
          {"state_names", result.state_names},
          {"roots", result.roots},
          {"slot_permutation", result.slot_permutation},
          {"residual_norms", result.residual_norms},
          {"density_change_norms", result.density_change_norms},
          {"relative_residual_norms", result.relative_residual_norms},
          {"iteration_max_residuals", result.iteration_max_residuals},
          {"iteration_max_density_changes",
           result.iteration_max_density_changes},
          {"iteration_max_relative_residuals",
           result.iteration_max_relative_residuals},
          {"convergence_mode", result.convergence_mode},
          {"accelerator_mode", result.accelerator_mode},
          {"accelerator_subspace", result.accelerator_subspace},
          {"density_convergence_target",
           result.density_convergence_target},
          {"relative_convergence_target",
           result.relative_convergence_target},
          {"max_rotation", result.max_rotation}};
    };
    auto default_roots_from_legacy_fields =
        [](const std::vector<double> &energies,
           const std::vector<std::string> &state_names) {
          std::vector<ExcitedRootDescriptor> roots;
          roots.reserve(energies.size());
          for (size_t i = 0; i < energies.size(); ++i) {
            const std::string display_name =
                (i < state_names.size() && !state_names[i].empty())
                    ? state_names[i]
                    : ("es" + std::to_string(i + 1));
            roots.push_back(
                {std::string("es_root_") +
                     [&]() {
                       std::ostringstream os;
                       os << std::setw(4) << std::setfill('0') << i;
                       return os.str();
                     }(),
                 i,
                 i,
                 energies[i],
                 display_name});
          }
          return roots;
        };
    auto protocol_result_from_json = [&](const nlohmann::json &node) {
      ExcitedProtocolResult result;
      if (!node.is_object()) {
        return result;
      }
      result.attempted = node.value("attempted", false);
      result.saved = node.value("saved", false);
      result.converged = node.value("converged", false);
      result.failed = node.value("failed", false);
      result.skipped = node.value("skipped", false);
      result.restart_reused = node.value("restart_reused", false);
      result.stage_status = node.value("stage_status",
                                       std::string("placeholder_pending_solver"));
      result.response_variant = node.value("response_variant",
                                           std::string("unknown"));
      result.restart_support_mode = node.value(
          "restart_support_mode", std::string("guess_only"));
      result.restart_source =
          node.value("restart_source", std::string("none"));
      result.snapshot_kind =
          node.value("snapshot_kind", std::string("none"));
      result.bundle_state_present = node.value("bundle_state_present", false);
      result.restart_capable = node.value("restart_capable", false);
      result.restart_source_threshold =
          node.value("restart_source_threshold", 0.0);
      result.iterations = node.value("iterations", static_cast<size_t>(0));
      if (node.contains("energies") && node["energies"].is_array()) {
        result.energies = node["energies"].get<std::vector<double>>();
      }
      if (node.contains("state_names") && node["state_names"].is_array()) {
        result.state_names = node["state_names"].get<std::vector<std::string>>();
      }
      if (node.contains("roots") && node["roots"].is_array()) {
        result.roots = node["roots"].get<std::vector<ExcitedRootDescriptor>>();
      }
      if (node.contains("slot_permutation") &&
          node["slot_permutation"].is_array()) {
        result.slot_permutation =
            node["slot_permutation"].get<std::vector<size_t>>();
      }
      if (node.contains("residual_norms") && node["residual_norms"].is_array()) {
        result.residual_norms = node["residual_norms"].get<std::vector<double>>();
      }
      if (node.contains("density_change_norms") &&
          node["density_change_norms"].is_array()) {
        result.density_change_norms =
            node["density_change_norms"].get<std::vector<double>>();
      }
      if (node.contains("relative_residual_norms") &&
          node["relative_residual_norms"].is_array()) {
        result.relative_residual_norms =
            node["relative_residual_norms"].get<std::vector<double>>();
      }
      if (node.contains("iteration_max_residuals") &&
          node["iteration_max_residuals"].is_array()) {
        result.iteration_max_residuals =
            node["iteration_max_residuals"].get<std::vector<double>>();
      }
      if (node.contains("iteration_max_density_changes") &&
          node["iteration_max_density_changes"].is_array()) {
        result.iteration_max_density_changes =
            node["iteration_max_density_changes"].get<std::vector<double>>();
      }
      if (node.contains("iteration_max_relative_residuals") &&
          node["iteration_max_relative_residuals"].is_array()) {
        result.iteration_max_relative_residuals =
            node["iteration_max_relative_residuals"].get<std::vector<double>>();
      }
      result.convergence_mode = node.value("convergence_mode",
                                           std::string("max_residual"));
      result.accelerator_mode = node.value("accelerator_mode",
                                           std::string("none"));
      result.accelerator_subspace =
          node.value("accelerator_subspace", static_cast<size_t>(0));
      result.density_convergence_target =
          node.value("density_convergence_target", 0.0);
      result.relative_convergence_target =
          node.value("relative_convergence_target", 0.0);
      result.max_rotation = node.value("max_rotation", 0.0);
      if (result.roots.empty() && !result.energies.empty()) {
        result.roots = default_roots_from_legacy_fields(result.energies,
                                                        result.state_names);
      }
      if (result.slot_permutation.empty() && !result.roots.empty()) {
        result.slot_permutation.reserve(result.roots.size());
        for (const auto &root : result.roots) {
          result.slot_permutation.push_back(root.stable_index);
        }
      }
      return result;
    };

    size_t ready_protocol_placeholders = 0;
    size_t restart_ready_protocols = 0;
    size_t pending_protocols = 0;
    size_t skipped_protocols = 0;
    size_t completed_protocols = 0;
    size_t failed_protocols = 0;
    size_t attempted_protocols = 0;
    double total_wall_seconds = 0.0;
    double total_cpu_seconds = 0.0;
    nlohmann::json protocol_events = nlohmann::json::array();

    for (size_t protocol_index = 0; protocol_index < plan.protocols.size();
         ++protocol_index) {
      const double threshold = plan.protocols[protocol_index];
      const double protocol_wall_start = madness::wall_time();
      const double protocol_cpu_start = madness::cpu_time();

      nlohmann::json dispatch = nlohmann::json::object();
      if (world.rank() == 0) {
        const std::string protocol_key = ResponseRecord2::protocol_key(threshold);
        auto &excited = state_metadata_json["excited_states"];
        auto &node = excited["protocols"][protocol_key];
        ensure_excited_protocol_placeholder_node(node);
        ++ready_protocol_placeholders;
        dispatch["protocol_key"] = protocol_key;
        dispatch["saved"] = node["saved"].get<bool>();
        dispatch["converged"] = node["converged"].get<bool>();
        dispatch["response_variant"] = node["response_variant"];
        dispatch["restart_support_mode"] = node["restart_support_mode"];
        dispatch["restart_source"] = node["restart_source"];
        dispatch["snapshot_kind"] = node["snapshot_kind"];
        dispatch["bundle_state_present"] = node["bundle_state_present"];
        dispatch["restart_capable"] = node["restart_capable"];
        dispatch["restart_source_threshold"] = node["restart_source_threshold"];
        dispatch["energies"] = node["energies"];
        dispatch["state_names"] = node["state_names"];
        dispatch["roots"] = node["roots"];
        dispatch["slot_permutation"] = node["slot_permutation"];
        dispatch["iterations"] = node["iterations"];
        dispatch["residual_norms"] = node["residual_norms"];
        dispatch["density_change_norms"] = node["density_change_norms"];
        dispatch["relative_residual_norms"] =
            node["relative_residual_norms"];
        dispatch["iteration_max_residuals"] = node["iteration_max_residuals"];
        dispatch["iteration_max_density_changes"] =
            node["iteration_max_density_changes"];
        dispatch["iteration_max_relative_residuals"] =
            node["iteration_max_relative_residuals"];
        dispatch["convergence_mode"] = node["convergence_mode"];
        dispatch["accelerator_mode"] = node["accelerator_mode"];
        dispatch["accelerator_subspace"] = node["accelerator_subspace"];
        dispatch["density_convergence_target"] =
            node["density_convergence_target"];
        dispatch["relative_convergence_target"] =
            node["relative_convergence_target"];
        dispatch["max_rotation"] = node["max_rotation"];
        const bool node_restart_ready =
            dispatch["saved"].get<bool>() && dispatch["converged"].get<bool>() &&
            dispatch["bundle_state_present"].get<bool>() &&
            dispatch["restart_capable"].get<bool>();
        if (!plan.enabled) {
          dispatch["solver_needed"] = false;
          dispatch["stage_status"] = "disabled_skip";
        } else if (node_restart_ready) {
          dispatch["solver_needed"] = false;
          dispatch["stage_status"] = "restart_ready_skip";
          dispatch["restart_source"] = "current_protocol_snapshot";
          dispatch["restart_source_threshold"] = threshold;
        } else {
          dispatch["solver_needed"] = true;
          dispatch["stage_status"] = "placeholder_pending_solver";
        }
      }
      broadcast_json(world, dispatch);

      ExcitedProtocolResult protocol_result;
      const bool solver_needed = dispatch.value("solver_needed", false);
      const bool node_saved = dispatch.value("saved", false);
      const bool node_converged = dispatch.value("converged", false);
      if (!solver_needed) {
        protocol_result.attempted = false;
        protocol_result.saved = node_saved;
        protocol_result.converged = node_converged;
        protocol_result.failed = false;
        protocol_result.skipped = true;
        protocol_result.restart_reused = node_saved && node_converged;
        protocol_result.stage_status = dispatch.value(
            "stage_status", std::string("placeholder_pending_solver"));
        protocol_result.response_variant =
            dispatch.value("response_variant", std::string("unknown"));
        protocol_result.restart_support_mode = dispatch.value(
            "restart_support_mode", std::string("guess_only"));
        protocol_result.restart_source =
            dispatch.value("restart_source", std::string("none"));
        protocol_result.snapshot_kind =
            dispatch.value("snapshot_kind", std::string("none"));
        protocol_result.bundle_state_present =
            dispatch.value("bundle_state_present", false);
        protocol_result.restart_capable =
            dispatch.value("restart_capable", false);
        protocol_result.restart_source_threshold =
            dispatch.value("restart_source_threshold", 0.0);
        protocol_result.iterations =
            dispatch.value("iterations", static_cast<size_t>(0));
        if (dispatch.contains("energies") && dispatch["energies"].is_array()) {
          protocol_result.energies =
              dispatch["energies"].get<std::vector<double>>();
        }
        if (dispatch.contains("state_names") &&
            dispatch["state_names"].is_array()) {
          protocol_result.state_names =
              dispatch["state_names"].get<std::vector<std::string>>();
        }
        if (dispatch.contains("roots") && dispatch["roots"].is_array()) {
          protocol_result.roots =
              dispatch["roots"].get<std::vector<ExcitedRootDescriptor>>();
        }
        if (dispatch.contains("slot_permutation") &&
            dispatch["slot_permutation"].is_array()) {
          protocol_result.slot_permutation =
              dispatch["slot_permutation"].get<std::vector<size_t>>();
        }
        if (dispatch.contains("residual_norms") &&
            dispatch["residual_norms"].is_array()) {
          protocol_result.residual_norms =
              dispatch["residual_norms"].get<std::vector<double>>();
        }
        if (dispatch.contains("density_change_norms") &&
            dispatch["density_change_norms"].is_array()) {
          protocol_result.density_change_norms =
              dispatch["density_change_norms"].get<std::vector<double>>();
        }
        if (dispatch.contains("relative_residual_norms") &&
            dispatch["relative_residual_norms"].is_array()) {
          protocol_result.relative_residual_norms =
              dispatch["relative_residual_norms"].get<std::vector<double>>();
        }
        if (dispatch.contains("iteration_max_residuals") &&
            dispatch["iteration_max_residuals"].is_array()) {
          protocol_result.iteration_max_residuals =
              dispatch["iteration_max_residuals"].get<std::vector<double>>();
        }
        if (dispatch.contains("iteration_max_density_changes") &&
            dispatch["iteration_max_density_changes"].is_array()) {
          protocol_result.iteration_max_density_changes =
              dispatch["iteration_max_density_changes"]
                  .get<std::vector<double>>();
        }
        if (dispatch.contains("iteration_max_relative_residuals") &&
            dispatch["iteration_max_relative_residuals"].is_array()) {
          protocol_result.iteration_max_relative_residuals =
              dispatch["iteration_max_relative_residuals"]
                  .get<std::vector<double>>();
        }
        protocol_result.convergence_mode = dispatch.value(
            "convergence_mode", std::string("max_residual"));
        protocol_result.accelerator_mode = dispatch.value(
            "accelerator_mode", std::string("none"));
        protocol_result.accelerator_subspace = dispatch.value(
            "accelerator_subspace", static_cast<size_t>(0));
        protocol_result.density_convergence_target = dispatch.value(
            "density_convergence_target", 0.0);
        protocol_result.relative_convergence_target = dispatch.value(
            "relative_convergence_target", 0.0);
        protocol_result.max_rotation = dispatch.value("max_rotation", 0.0);
      } else {
        ExcitedProtocolInput protocol_input;
        protocol_input.threshold = threshold;
        protocol_input.dconv = response_params.dconv();
        protocol_input.protocol_index = protocol_index;
        protocol_input.restart_saved = node_saved;
        protocol_input.restart_converged = node_converged;
        protocol_input.owner_group = plan.owner_group;
        protocol_input.tda = plan.tda;
        protocol_input.num_states = plan.num_states;
        protocol_input.guess_max_iter = plan.guess_max_iter;
        protocol_input.maxiter = plan.maxiter;
        protocol_input.maxsub = plan.maxsub;
        protocol_result = excited_solver.solve_protocol(world, protocol_input);
        if (protocol_result.stage_status.empty()) {
          protocol_result.stage_status = "placeholder_pending_solver";
        }
      }

      nlohmann::json protocol_result_json = nlohmann::json::object();
      if (world.rank() == 0) {
        protocol_result_json = protocol_result_to_json(protocol_result);
      }
      broadcast_json(world, protocol_result_json);
      protocol_result = protocol_result_from_json(protocol_result_json);
      protocol_result =
          ResponseRecord2::normalize_excited_protocol_result(
              threshold, std::move(protocol_result));

      const double protocol_wall_seconds =
          madness::wall_time() - protocol_wall_start;
      const double protocol_cpu_seconds =
          madness::cpu_time() - protocol_cpu_start;

      excited_metadata_record.record_excited_protocol_result(
          threshold, plan.owner_group, protocol_result, protocol_wall_seconds,
          protocol_cpu_seconds);
      state_metadata_json = excited_metadata_record.to_json();

      if (world.rank() == 0) {
        const std::string protocol_key = dispatch.value(
            "protocol_key", ResponseRecord2::protocol_key(threshold));
        total_wall_seconds += protocol_wall_seconds;
        total_cpu_seconds += protocol_cpu_seconds;
        if (protocol_result.attempted) {
          ++attempted_protocols;
        }
        if (protocol_result.restart_reused) {
          ++restart_ready_protocols;
        }
        if (protocol_result.skipped) {
          ++skipped_protocols;
        }
        if (!protocol_result.skipped && !protocol_result.failed &&
            !protocol_result.converged) {
          ++pending_protocols;
        }
        if (protocol_result.failed) {
          ++failed_protocols;
        } else {
          ++completed_protocols;
        }

        protocol_events.push_back(
            {{"protocol_key", protocol_key},
             {"threshold", threshold},
             {"saved", protocol_result.saved},
             {"converged", protocol_result.converged},
             {"attempted", protocol_result.attempted},
             {"failed", protocol_result.failed},
             {"skipped", protocol_result.skipped},
             {"restart_reused", protocol_result.restart_reused},
             {"stage_status", protocol_result.stage_status},
             {"response_variant", protocol_result.response_variant},
             {"restart_support_mode", protocol_result.restart_support_mode},
             {"restart_source", protocol_result.restart_source},
             {"snapshot_kind", protocol_result.snapshot_kind},
             {"bundle_state_present", protocol_result.bundle_state_present},
             {"restart_capable", protocol_result.restart_capable},
             {"restart_source_threshold",
              protocol_result.restart_source_threshold},
             {"iterations", protocol_result.iterations},
             {"energies", protocol_result.energies},
             {"state_names", protocol_result.state_names},
             {"roots", build_excited_root_manifest(protocol_result.roots)},
             {"slot_permutation", protocol_result.slot_permutation},
             {"residual_norms", protocol_result.residual_norms},
             {"density_change_norms", protocol_result.density_change_norms},
             {"relative_residual_norms",
              protocol_result.relative_residual_norms},
             {"iteration_max_residuals",
              protocol_result.iteration_max_residuals},
             {"iteration_max_density_changes",
              protocol_result.iteration_max_density_changes},
             {"iteration_max_relative_residuals",
              protocol_result.iteration_max_relative_residuals},
             {"convergence_mode", protocol_result.convergence_mode},
             {"accelerator_mode", protocol_result.accelerator_mode},
             {"accelerator_subspace",
              protocol_result.accelerator_subspace},
             {"density_convergence_target",
              protocol_result.density_convergence_target},
             {"relative_convergence_target",
              protocol_result.relative_convergence_target},
             {"max_rotation", protocol_result.max_rotation},
             {"wall_seconds", protocol_wall_seconds},
             {"cpu_seconds", protocol_cpu_seconds}});
      }
    }

    if (world.rank() == 0) {
      execution["attempted"] = attempted_protocols > 0;
      execution["ready_protocol_placeholders"] = ready_protocol_placeholders;
      execution["restart_ready_protocols"] = restart_ready_protocols;
      execution["pending_protocols"] = pending_protocols;
      execution["skipped_protocols"] = skipped_protocols;
      execution["completed_protocols"] = completed_protocols;
      execution["failed_protocols"] = failed_protocols;
      execution["total_wall_seconds"] = total_wall_seconds;
      execution["total_cpu_seconds"] = total_cpu_seconds;
      execution["protocol_events"] = std::move(protocol_events);

      if (plan.enabled || ready_protocol_placeholders > 0) {
        print("EXCITED_BUNDLE_SUMMARY protocols=",
              ready_protocol_placeholders,
              " pending=", pending_protocols,
              " restart_ready=", restart_ready_protocols,
              " skipped=", skipped_protocols, ".");
      }
    }

    if (world.rank() == 0) {
      state_metadata_json["excited_state_planner"] = {
          {"note",
           "Stage 2c executes the protocol-aware excited-state bundle adapter. "
           "Per-protocol status, timings, restart provenance, and solver "
           "diagnostics are recorded in metadata."},
          {"plan", plan.to_json()},
          {"execution", execution}};
      write_json_file("response_metadata.json", state_metadata_json);
    }

    broadcast_json(world, state_metadata_json);
    broadcast_json(world, execution);
    return ExcitedExecutionResult{std::move(execution)};
  }

  static DerivedExecutionResult execute_derived_state_requests(
      World &world, const CalculationParameters &calc_params, GroundContext &ctx,
      const PlannedStates &planned_states, const nlohmann::json &state_metadata,
      bool subgroup_parallel_requested, bool owner_group_schedule,
      size_t final_ti, double final_thresh) {
    // Derived stage:
    // 1) evaluate dependency gate against final linear-state readiness,
    // 2) execute ready requests in subgroup mode when available,
    // 3) fall back to deterministic serial lanes when needed,
    // 4) return both gate diagnostics and execution summary.
    const auto &state_parallel_plan = planned_states.state_parallel_plan;

    DerivedStateGateReport derived_gate =
        DerivedStatePlanner::evaluate_dependency_gate(
            planned_states.derived_state_plan, planned_states.generated_states,
            final_ti, [&](const LinearResponsePoint &pt) {
              return point_ready_in_metadata(state_metadata, pt,
                                             /*require_saved=*/true,
                                             /*require_converged=*/true);
            });
    if (world.rank() == 0 && derived_gate.total_requests > 0) {
      print("DERIVED_GATE ready=", derived_gate.ready_requests, "/",
            derived_gate.total_requests,
            " blocked=", derived_gate.blocked_requests);
    }

    // Execution summary persisted under metadata["derived_state_planner"].
    nlohmann::json derived_execution = {
        {"attempted", false},
        {"mode", "none_ready"},
        {"execution_groups", 1},
        {"ready_requests", derived_gate.ready_requests},
        {"blocked_requests", derived_gate.blocked_requests},
        {"precompleted_requests", 0},
        {"completed_requests", 0},
        {"failed_requests", 0},
        {"total_wall_seconds", 0.0},
        {"total_cpu_seconds", 0.0},
        {"request_timings", nlohmann::json::array()}};

    // Plan indices that passed the dependency gate at the final threshold.
    std::vector<size_t> ready_request_indices;
    ready_request_indices.reserve(planned_states.derived_state_plan.requests.size());
    size_t precompleted_requests = 0;
    for (size_t i = 0; i < planned_states.derived_state_plan.requests.size();
         ++i) {
      if (i >= derived_gate.entries.size() || !derived_gate.entries[i].ready) {
        continue;
      }
      const auto &req = planned_states.derived_state_plan.requests[i];
      if (derived_request_done_record_exists(final_ti, final_thresh, req, i)) {
        ++precompleted_requests;
        continue;
      }
      ready_request_indices.push_back(i);
    }
    derived_execution["precompleted_requests"] = precompleted_requests;
    if (world.rank() == 0 && precompleted_requests > 0) {
      print("DERIVED_RESTART_SKIP precompleted=", precompleted_requests,
            " remaining_ready=", ready_request_indices.size());
    }

    if (!ready_request_indices.empty()) {
      const size_t derived_owner_groups =
          subgroup_parallel_requested
              ? state_parallel_plan.execution_groups
              : std::max<size_t>(1, state_parallel_plan.mapping_groups);
      // Deterministic lane assignment for ready derived requests.
      std::vector<size_t> owner_by_ready_index(ready_request_indices.size(), 0);
      for (size_t i = 0; i < ready_request_indices.size(); ++i) {
        owner_by_ready_index[i] = i % derived_owner_groups;
      }

      derived_execution["attempted"] = true;
      derived_execution["execution_groups"] = derived_owner_groups;

      bool ran_subgroup_derived = false;
      if (subgroup_parallel_requested && derived_owner_groups > 1) {
        derived_execution["mode"] = "owner_group_subworld";
        try {
          auto subworld_ptr = MacroTaskQ::create_worlds(
              world, state_parallel_plan.execution_groups);
          if (!subworld_ptr) {
            throw std::runtime_error("derived subworld creation returned null");
          }

          World &subworld = *subworld_ptr;
          const size_t subgroup_id = static_cast<size_t>(
              world.rank() %
              static_cast<int>(state_parallel_plan.execution_groups));
          ScopedRankLogRedirect subgroup_console_redirect(
              subworld.rank() == 0, group_console_file(subgroup_id),
              "stage2-derived subgroup=" + std::to_string(subgroup_id));
          const std::string timing_shard_file =
              group_derived_timing_file(subgroup_id);

          auto old_pmap3 = FunctionDefaults<3>::get_pmap();
          auto restore_pmap = [&]() { FunctionDefaults<3>::set_pmap(old_pmap3); };

          long local_completed = 0;
          long local_failed = 0;
          double local_wall_seconds = 0.0;
          double local_cpu_seconds = 0.0;
          nlohmann::json local_request_timings = nlohmann::json::array();
          FunctionDefaults<3>::set_default_pmap(subworld);
          try {
            std::vector<size_t> local_ready_positions;
            local_ready_positions.reserve(
                ready_request_indices.size() /
                std::max<size_t>(1, derived_owner_groups));
            for (size_t pos = 0; pos < ready_request_indices.size(); ++pos) {
              if (owner_by_ready_index[pos] == subgroup_id) {
                local_ready_positions.push_back(pos);
              }
            }

            if (!local_ready_positions.empty()) {
              GroundStateData local_ground(subworld, ctx.archive_file,
                                           ctx.molecule);
              ResponseManager local_response_manager(subworld, calc_params);
              local_response_manager.setProtocol(subworld, local_ground.getL(),
                                                 final_thresh);
              local_ground.prepareOrbitals(
                  subworld, FunctionDefaults<3>::get_k(), final_thresh);

              SimpleVBCComputer local_vbc_computer(subworld, local_ground);
              for (const auto pos : local_ready_positions) {
                const auto &req =
                    planned_states.derived_state_plan
                        .requests[ready_request_indices[pos]];
                const auto timing = run_derived_request(
                    subworld, local_ground, req, local_vbc_computer,
                    local_completed, local_failed);
                if (subworld.rank() == 0) {
                  local_wall_seconds += timing.wall_seconds;
                  local_cpu_seconds += timing.cpu_seconds;
                  local_request_timings.push_back(
                      {{"derived_state_id", req.derived_state_id},
                       {"owner_group", subgroup_id},
                       {"success", timing.success},
                       {"wall_seconds", timing.wall_seconds},
                       {"cpu_seconds", timing.cpu_seconds}});
                }
              }
            } else if (subworld.rank() == 0) {
              print("DERIVED_SUBGROUP_NO_WORK subgroup=", subgroup_id);
            }

            if (subworld.rank() == 0) {
              write_json_file(timing_shard_file, local_request_timings);
            }
            subworld.gop.fence();
            restore_pmap();
          } catch (...) {
            restore_pmap();
            throw;
          }

          if (subworld.rank() != 0) {
            local_completed = 0;
            local_failed = 0;
            local_wall_seconds = 0.0;
            local_cpu_seconds = 0.0;
          }
          world.gop.sum(local_completed);
          world.gop.sum(local_failed);
          world.gop.sum(local_wall_seconds);
          world.gop.sum(local_cpu_seconds);
          derived_execution["completed_requests"] = local_completed;
          derived_execution["failed_requests"] = local_failed;
          derived_execution["total_wall_seconds"] = local_wall_seconds;
          derived_execution["total_cpu_seconds"] = local_cpu_seconds;
          if (world.rank() == 0) {
            print("DERIVED_TIMING_SUMMARY completed=", local_completed,
                  " failed=", local_failed, " wall_s=", local_wall_seconds,
                  " cpu_s=", local_cpu_seconds);
          }

          world.gop.fence();
          nlohmann::json merged_request_timings = nlohmann::json::array();
          if (world.rank() == 0) {
            for (size_t gid = 0; gid < derived_owner_groups; ++gid) {
              const auto shard =
                  read_json_file_or_object(group_derived_timing_file(gid));
              if (!shard.is_array()) {
                continue;
              }
              for (const auto &entry : shard) {
                merged_request_timings.push_back(entry);
              }
            }
          }
          std::string request_timing_dump;
          if (world.rank() == 0) {
            request_timing_dump = merged_request_timings.dump();
          }
          world.gop.broadcast_serializable(request_timing_dump, 0);
          if (world.rank() != 0) {
            merged_request_timings = request_timing_dump.empty()
                                         ? nlohmann::json::array()
                                         : nlohmann::json::parse(
                                               request_timing_dump);
          }
          derived_execution["request_timings"] =
              std::move(merged_request_timings);
          world.gop.fence();
          ran_subgroup_derived = true;
        } catch (const std::exception &ex) {
          if (world.rank() == 0) {
            print("ERROR DERIVED_SUBGROUP_FAILED message=", ex.what(),
                  " action=fallback_serial_derived_loop");
          }
        }
      }

      if (!ran_subgroup_derived) {
        // Deterministic serial fallback keeps a predictable owner-lane order
        // for easier debugging and reproducible logs.
        derived_execution["mode"] =
            (owner_group_schedule && derived_owner_groups > 1)
                ? "owner_group_serial_lanes"
                : "serial";
        if (world.rank() == 0 &&
            derived_execution["mode"].get<std::string>() ==
                "owner_group_serial_lanes") {
          print("DERIVED_LANE_MODE owner_groups=",
                derived_owner_groups,
                " strategy=deterministic_serial_lanes");
        }

        long local_completed = 0;
        long local_failed = 0;
        double local_wall_seconds = 0.0;
        double local_cpu_seconds = 0.0;
        nlohmann::json request_timings = nlohmann::json::array();
        SimpleVBCComputer serial_vbc_computer(world, ctx.ground);
        if (owner_group_schedule && derived_owner_groups > 1) {
          for (size_t gid = 0; gid < derived_owner_groups; ++gid) {
            if (world.rank() == 0) {
              print("DERIVED_LANE_START lane=", gid, " lane_max=",
                    derived_owner_groups - 1);
            }
            for (size_t pos = 0; pos < ready_request_indices.size(); ++pos) {
              if (owner_by_ready_index[pos] != gid) {
                continue;
              }
              const auto &req = planned_states.derived_state_plan
                                    .requests[ready_request_indices[pos]];
              const auto timing = run_derived_request(
                  world, ctx.ground, req, serial_vbc_computer, local_completed,
                  local_failed);
              if (world.rank() == 0) {
                local_wall_seconds += timing.wall_seconds;
                local_cpu_seconds += timing.cpu_seconds;
                request_timings.push_back({{"derived_state_id", req.derived_state_id},
                                           {"owner_group", gid},
                                           {"success", timing.success},
                                           {"wall_seconds", timing.wall_seconds},
                                           {"cpu_seconds", timing.cpu_seconds}});
              }
            }
          }
        } else {
          for (const auto idx : ready_request_indices) {
            const auto &req = planned_states.derived_state_plan.requests[idx];
            const auto timing =
                run_derived_request(world, ctx.ground, req, serial_vbc_computer,
                                    local_completed, local_failed);
            if (world.rank() == 0) {
              local_wall_seconds += timing.wall_seconds;
              local_cpu_seconds += timing.cpu_seconds;
              request_timings.push_back({{"derived_state_id", req.derived_state_id},
                                         {"owner_group", 0},
                                         {"success", timing.success},
                                         {"wall_seconds", timing.wall_seconds},
                                         {"cpu_seconds", timing.cpu_seconds}});
            }
          }
        }

        if (world.rank() != 0) {
          local_completed = 0;
          local_failed = 0;
          local_wall_seconds = 0.0;
          local_cpu_seconds = 0.0;
        }
        world.gop.sum(local_completed);
        world.gop.sum(local_failed);
        world.gop.sum(local_wall_seconds);
        world.gop.sum(local_cpu_seconds);
        derived_execution["completed_requests"] = local_completed;
        derived_execution["failed_requests"] = local_failed;
        derived_execution["total_wall_seconds"] = local_wall_seconds;
        derived_execution["total_cpu_seconds"] = local_cpu_seconds;
        if (world.rank() == 0) {
          print("DERIVED_TIMING_SUMMARY completed=", local_completed,
                " failed=", local_failed, " wall_s=", local_wall_seconds,
                " cpu_s=", local_cpu_seconds);
        }
        derived_execution["request_timings"] = std::move(request_timings);
      }
    }

    std::string derived_execution_dump;
    if (world.rank() == 0) {
      derived_execution_dump = derived_execution.dump();
    }
    world.gop.broadcast_serializable(derived_execution_dump, 0);
    if (world.rank() != 0) {
      derived_execution =
          derived_execution_dump.empty()
              ? nlohmann::json::object()
              : nlohmann::json::parse(derived_execution_dump);
    }

    return DerivedExecutionResult{std::move(derived_gate),
                                  std::move(derived_execution)};
  }

  /// Rebuild final-protocol context and summarize linear convergence status.
  ///
  /// This prepares stage-3 property assembly with protocol-consistent operators
  /// and reports whether any final-threshold points remain unconverged.
  static FinalProtocolState prepare_and_validate_final_protocol_state(
      World &world, const CalculationParameters &calc_params, GroundContext &ctx,
      const PlannedStates &planned_states,
      const nlohmann::json &state_metadata_json) {
    const auto &protocol = calc_params.protocol();
    FinalProtocolState final_state{protocol.back(), protocol.size() - 1};

    ctx.response_manager.setProtocol(world, ctx.ground.getL(),
                                     final_state.threshold);
    ctx.ground.prepareOrbitals(world, FunctionDefaults<3>::get_k(),
                               final_state.threshold);
    ctx.ground.computePreliminaries(world, *ctx.response_manager.getCoulombOp(),
                                    ctx.response_manager.getVtol(),
                                    ctx.fock_json_file);

    bool all_are_converged = true;
    size_t unconverged_points = 0;
    for (const auto &state : planned_states.generated_states.states) {
      for (size_t fi = 0; fi < state.num_frequencies(); ++fi) {
        LinearResponsePoint pt{state, final_state.threshold_index, fi};
        if (!point_ready_in_metadata(state_metadata_json, pt,
                                     /*require_saved=*/false,
                                     /*require_converged=*/true)) {
          all_are_converged = false;
          ++unconverged_points;
        }
      }
    }
    final_state.all_linear_points_converged = all_are_converged;
    if (!all_are_converged && world.rank() == 0) {
      print("WARN FINAL_PROTOCOL_UNCONVERGED_POINTS count=",
            unconverged_points,
            " action=continue_partial_property_execution");
    }
    return final_state;
  }

  /// Assemble stage-2 metadata consumed by stage 3 and restart diagnostics.
  static nlohmann::json build_state_stage_metadata(
      const PlannedStates &planned_states, const nlohmann::json &state_metadata,
      const StateSolveScheduleContext &schedule_ctx,
      const FinalProtocolState &final_state,
      const ExcitedExecutionResult &excited_result,
      const DerivedExecutionResult &derived_result) {
    auto metadata = state_metadata;
    metadata["state_parallel_planner"] =
        planned_states.state_parallel_plan.to_json();
    metadata["state_parallel_runtime"] = {
        {"effective_point_groups", schedule_ctx.point_owner_groups()},
        {"runtime_execution_groups", schedule_ctx.runtime_execution_groups},
        {"restart_start_protocol_index",
         schedule_ctx.restart_start_protocol_index},
        {"effective_point_parallel_start_protocol_index",
         schedule_ctx.runtime_point_parallel_start_protocol_index},
        {"restart_final_protocol_only",
         schedule_ctx.restart_final_protocol_only},
        {"force_retry_removed_frequencies",
         schedule_ctx.force_retry_removed_frequencies},
        {"remaining_final_pending_points",
         schedule_ctx.remaining_final_pending_points},
        {"new_points_without_history", schedule_ctx.new_points_without_history},
        {"all_linear_points_converged", final_state.all_linear_points_converged}};
    nlohmann::json protocol_policy = nlohmann::json::array();
    for (const auto &policy : schedule_ctx.runtime_protocol_policies) {
      protocol_policy.push_back(policy.to_json());
    }
    metadata["state_parallel_runtime"]["protocol_execution_policy"] =
        std::move(protocol_policy);
    if (schedule_ctx.runtime_execution_groups > 1) {
      nlohmann::json group_console_logs = nlohmann::json::array();
      for (size_t gid = 0;
           gid < schedule_ctx.runtime_execution_groups; ++gid) {
        group_console_logs.push_back(group_console_file(gid));
      }
      metadata["state_parallel_runtime"]["group_console_logs"] =
          std::move(group_console_logs);
    }
    metadata["derived_state_planner"] = {
        {"note",
         "Stage 2c: ready VBC-derived requests are executed before property "
         "assembly; blocked requests remain gated."},
        {"plan", planned_states.derived_state_plan.to_json()},
        {"dependency_gate", derived_result.dependency_gate.to_json()},
        {"execution", derived_result.execution}};
    metadata["excited_state_planner"] = {
        {"note",
         "Stage 2c executes the protocol-aware excited-state bundle adapter. "
         "Per-protocol status, timings, restart provenance, and solver "
         "diagnostics are recorded in metadata."},
        {"plan", planned_states.excited_state_bundle_plan.to_json()},
        {"execution", excited_result.execution}};
    return metadata;
  }

  /// Stage-2 orchestrator.
  ///
  /// Sequence:
  /// 1. Build runtime ownership policy.
  /// 2. Solve linear states (subgroup path first, serial fallback second).
  /// 3. Run excited-state bundle adapter stage.
  /// 4. Run dependency-gated derived requests.
  /// 5. Publish merged stage metadata/debug payloads.
  static SolvedStates
  solve_all_states(World &world, const CalculationParameters &calc_params,
                   GroundContext &ctx,
                   const ResponseParameters &response_params,
                   PlannedStates planned_states) {
    // Stage 2a setup: build runtime ownership schedule and restart-aware policy.
    const StateSolveScheduleContext schedule_ctx =
        build_state_solve_schedule_context(world, calc_params, planned_states,
                                           response_params);

    nlohmann::json state_metadata_json = nlohmann::json::object();
    nlohmann::json debug_log_json = nlohmann::json::object();

    // Stage 2b linear solve: subgroup path first, serial fallback second.
    bool ran_subgroup_path = false;
    if (schedule_ctx.subgroup_parallel_requested) {
      ran_subgroup_path = execute_subgroup_state_solve(
          world, calc_params, ctx, schedule_ctx, planned_states,
          planned_states.excited_state_bundle_plan, state_metadata_json,
          debug_log_json);
    }
    if (!ran_subgroup_path) {
      execute_serial_state_solve(world, calc_params, ctx, schedule_ctx,
                                 planned_states.excited_state_bundle_plan,
                                 state_metadata_json, debug_log_json);
    }

    const FinalProtocolState final_state = prepare_and_validate_final_protocol_state(
        world, calc_params, ctx, planned_states, state_metadata_json);

    // Stage 2c excited-state bundle: protocol-aware adapter execution.
    const ExcitedExecutionResult excited_result =
        execute_excited_state_bundle_stage(world, planned_states, ctx,
                                           response_params, state_metadata_json);

    // Stage 2d derived solve: dependency-gated execution + summary metadata.
    const DerivedExecutionResult derived_result = execute_derived_state_requests(
        world, calc_params, ctx, planned_states, state_metadata_json,
        schedule_ctx.subgroup_parallel_requested, schedule_ctx.owner_group_schedule,
        final_state.threshold_index, final_state.threshold);

    // Stage 2 metadata assembly for downstream property stage.
    auto metadata = build_state_stage_metadata(
        planned_states, state_metadata_json, schedule_ctx, final_state, excited_result,
        derived_result);

    return SolvedStates{std::move(planned_states), std::move(metadata),
                        std::move(debug_log_json)};
  }

  // ============================================================================
  // SECTION 12: Stage 3 — Property assembly
  //   PropertyContext — aggregates world, parameters, ground context, solved
  //                     states, and PropertyManager for property dispatch.
  //   PropertyType — internal dispatch enum (Alpha, Beta, Raman).
  //   parse_property_name — maps user string → PropertyType.
  //   compute_polarizability — stage-3 alpha assembly with frequency filter.
  //   compute_hyperpolarizability — stage-3 beta assembly.
  //   print_raman_table — formatted per-mode Raman output on rank 0.
  //   compute_raman — Hessian + Raman tensor assembly.
  //   compute_requested_properties — serial stage-3 dispatch loop.
  //   compute_requested_properties_with_property_group — subgroup-aware path
  //                     with component precompute, property-group execution,
  //                     and broadcast-to-world.
  // ============================================================================

  struct PropertyContext {
    // Universe world used for property assembly.
    World &world;
    // Parsed response inputs (frequencies, directions, requested properties).
    const ResponseParameters &response_params;
    // Ground-state objects and managers prepared in Stage 1.
    const GroundContext &ground_ctx;
    // Solved-state metadata/map produced by Stage 2.
    const SolvedStates &solved_states;
    // Property accumulator/writer.
    PropertyManager &properties;
    // Optional SCF handle needed by Raman/hessian routines.
    std::shared_ptr<SCF> scf_calc;
  };

  enum class PropertyType { Alpha, Beta, Raman };

  /// Parse user-requested property names into internal dispatch enum values.
  inline static PropertyType parse_property_name(const std::string &raw) {
    auto key = raw;
    if (key.size() >= 2 && ((key.front() == '"' && key.back() == '"') ||
                            (key.front() == '\'' && key.back() == '\''))) {
      key = key.substr(1, key.size() - 2);
    }
    if (key == "polarizability")
      return PropertyType::Alpha;
    if (key == "hyperpolarizability")
      return PropertyType::Beta;
    if (key == "raman")
      return PropertyType::Raman;
    MADNESS_EXCEPTION(std::string("Unknown property: " + key).c_str(), 0);
  }

  /// Stage-3 wrapper for alpha assembly.
  inline static void compute_polarizability(PropertyContext &ctx) {
    if (ctx.world.rank() == 0)
      madness::print("PROPERTY_STAGE_START property=polarizability");

    const auto &state_map =
        ctx.solved_states.planned_states.generated_states.state_map;
    const auto &all_frequencies = ctx.response_params.dipole_frequencies();
    const auto dipole_directions = ctx.response_params.dipole_directions();

    std::vector<double> ready_frequencies;
    ready_frequencies.reserve(all_frequencies.size());
    std::vector<double> dropped_unconverged_frequencies;
    std::vector<double> dropped_policy_frequencies;

    for (size_t freq_index = 0; freq_index < all_frequencies.size();
         ++freq_index) {
      bool all_direction_points_ready = true;
      bool any_direction_flagged_for_removal = false;
      for (const char dir : dipole_directions) {
        const std::string state_key =
            std::string("Dipole_") + std::string(1, dir);
        const auto state_it = state_map.find(state_key);
        if (state_it == state_map.end()) {
          all_direction_points_ready = false;
          break;
        }
        const auto &state_desc = state_it->second;
        if (freq_index >= state_desc.num_frequencies()) {
          all_direction_points_ready = false;
          break;
        }
        LinearResponsePoint pt{state_desc, state_desc.thresholds.size() - 1,
                               freq_index};
        if (!point_ready_in_metadata(ctx.solved_states.metadata, pt,
                                     /*require_saved=*/false,
                                     /*require_converged=*/true)) {
          all_direction_points_ready = false;
          break;
        }
        if (point_marked_for_frequency_removal_in_metadata(
                ctx.solved_states.metadata, pt)) {
          any_direction_flagged_for_removal = true;
          break;
        }
      }

      const double omega = all_frequencies[freq_index];
      if (all_direction_points_ready && !any_direction_flagged_for_removal) {
        ready_frequencies.push_back(omega);
      } else if (any_direction_flagged_for_removal) {
        dropped_policy_frequencies.push_back(omega);
      } else {
        dropped_unconverged_frequencies.push_back(omega);
      }
    }

    if (ctx.world.rank() == 0) {
      madness::print(
          "ALPHA_FREQUENCY_FILTER requested=", all_frequencies.size(),
          " ready=", ready_frequencies.size(),
          " dropped_unconverged=", dropped_unconverged_frequencies.size(),
          " dropped_failure_policy=", dropped_policy_frequencies.size());
      if (!dropped_policy_frequencies.empty()) {
        madness::print("ALPHA_FREQUENCY_FILTER_REMOVED frequencies=",
                       dropped_policy_frequencies);
      }
      if (!dropped_unconverged_frequencies.empty()) {
        madness::print("ALPHA_FREQUENCY_FILTER_UNCONVERGED frequencies=",
                       dropped_unconverged_frequencies);
      }
    }

    if (ready_frequencies.empty()) {
      if (ctx.world.rank() == 0) {
        madness::print(
            "WARN ALPHA_ASSEMBLY_SKIPPED reason=no_ready_frequencies");
      }
      return;
    }

    compute_alpha(
        ctx.world, state_map, ctx.ground_ctx.ground, ready_frequencies,
        dipole_directions, ctx.properties);

    ctx.properties.save();
  }

  /// Stage-3 wrapper for beta assembly using configured triplet policy.
  inline static void compute_hyperpolarizability(PropertyContext &ctx) {
    if (ctx.world.rank() == 0)
      madness::print("PROPERTY_STAGE_START property=hyperpolarizability");

    auto dip_dirs = ctx.response_params.dipole_directions();
    ::compute_hyperpolarizability(ctx.world, ctx.ground_ctx.ground,
                                  ctx.response_params.dipole_frequencies(),
                                  dip_dirs, ctx.properties,
                                  ctx.response_params.beta_shg(),
                                  ctx.response_params.beta_or(),
                                  ctx.response_params.beta_all_triplets());

    ctx.properties.save();
  }

  /// Print per-frequency Raman mode table to stdout on rank 0.
  static void print_raman_table(World &world, const RamanResults &raman,
                                int print_level) {
    if (world.rank() != 0) return;

    // 1 a.u. of photon energy → wavelength (nm)
    auto wavelength_nm_from_au = [](double omega_au) {
      const double AU_TO_NM_FACTOR = 45.5633525316;
      return (omega_au > 0.0) ? (AU_TO_NM_FACTOR / omega_au) : 0.0;
    };
    struct ColumnSpec {
      int width;
      int precision;
    };
    static const ColumnSpec COL_FREQ{10, 2};
    static const ColumnSpec COL_FLOAT{11, 6};

    auto print_header = [&](std::ostream &os, double omega_au) {
      const double lambda_nm = wavelength_nm_from_au(omega_au);
      os << "     Raman related properties for freq.  " << std::fixed
         << std::setprecision(6) << std::setw(9) << omega_au
         << " au  = " << std::setw(9) << std::setprecision(2) << lambda_nm
         << " nm\n";
      os << "     "
            "----------------------------------------------------------"
            "----"
            "-\n\n";
      os << " Mode    Freq.     Alpha**2   Beta(a)**2   Pol.Int.   "
            "Depol.Int.  Dep. Ratio \n\n";
    };

    auto print_row = [&](std::ostream &os,
                         const RamanResults::RamanModeRow &r) {
      double pol  = r.pol_int   ? *r.pol_int   : (45.0 * r.alpha2 + 4.0 * r.beta2);
      double depi = r.depol_int ? *r.depol_int : (3.0 * r.beta2);
      double rho;
      if (r.dep_ratio) {
        rho = *r.dep_ratio;
      } else {
        rho = (pol != 0.0) ? (depi / pol) : 0.0;
      }
      os << std::setw(5) << r.mode
         << " "
         << std::fixed << std::setw(COL_FREQ.width)
         << std::setprecision(COL_FREQ.precision) << r.freq_cm1
         << "  "
         << std::setw(COL_FLOAT.width) << std::setprecision(COL_FLOAT.precision)
         << r.alpha2
         << std::setw(COL_FLOAT.width) << std::setprecision(COL_FLOAT.precision)
         << r.beta2
         << std::setw(COL_FLOAT.width) << std::setprecision(COL_FLOAT.precision)
         << pol
         << std::setw(COL_FLOAT.width) << std::setprecision(COL_FLOAT.precision)
         << depi
         << std::setw(COL_FLOAT.width) << std::setprecision(COL_FLOAT.precision)
         << rho << "\n";
    };

    auto print_block =
        [&](std::ostream &os, double omega_au,
            const std::vector<RamanResults::RamanModeRow> &rows) {
          print_header(os, omega_au);
          for (auto r_it = rows.rbegin(); r_it != rows.rend(); ++r_it)
            print_row(os, *r_it);
          os << "\n";
        };

    for (double omega_au : raman.polarization_frequencies) {
      auto it = raman.raman_spectra.find(omega_au);
      if (it != raman.raman_spectra.end()) {
        print_block(std::cout, omega_au, it->second);
      }
    }
  }

  /// Stage-3 Raman wrapper (Hessian + Raman tensor/intensity assembly).
  inline static void compute_raman(PropertyContext &ctx,
                                   VibrationalResults &vib,
                                   RamanResults &raman) {
    if (ctx.world.rank() == 0) {
      madness::print("PROPERTY_STAGE_START property=raman");
    }

    vib = compute_hessian(
        ctx.world, ctx.solved_states.planned_states.generated_states.state_map,
        ctx.ground_ctx.ground, ctx.response_params.dipole_directions(),
        ctx.scf_calc);

    const double csg_factor = 142.9435756;
    Tensor<double> normal_modes = *vib.normalmodes_atomic;
    if (ctx.world.rank() == 0) {
      print("RAMAN_NORMAL_MODES_AU modes=", normal_modes);
    }
    auto vib_freq = *vib.frequencies * constants::au2invcm;

    std::vector<int> mode;
    for (int i = 0; i < vib_freq.size(); ++i) {

      if (abs(vib_freq(i)) > 1e-2) {
        mode.push_back(i);
        raman.vibrational_frequencies.push_back(vib_freq(i));
      }
    }
    auto nnmodes = normal_modes(_, Slice(mode[0], -1, 1));

    raman.normal_modes = nnmodes;

    if (ctx.world.rank() == 0) {
      print("RAMAN_ACTIVE_MODE_INDICES indices=", mode);
      print("RAMAN_ACTIVE_MODE_RANGE first=", mode[0],
            " last=", mode[mode.size() - 1]);
    }
    auto alpha_derivatives =
        compute_Raman(ctx.world, ctx.ground_ctx.ground,
                      ctx.response_params.dipole_frequencies(),
                      ctx.response_params.dipole_directions(),
                      ctx.response_params.nuclear_directions(), ctx.properties);
    raman.polarization_frequencies = ctx.response_params.dipole_frequencies();
    ctx.properties.save();
    ctx.world.gop.fence();
    auto compute_alpha2 = [](const Tensor<double> &alpha) {
      auto alpha_mean = 0.0;
      for (int i = 0; i < 3; ++i) {
        alpha_mean += alpha(i, i);
      }
      alpha_mean *= (1.0 / 3.0);
      return alpha_mean * alpha_mean;
    };
    auto compute_beta2 = [](const Tensor<double> &alpha) {
      auto beta2 = 0.0;
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          beta2 +=
              0.5 * (3 * alpha(i, j) * alpha(i, j) - alpha(i, i) * alpha(j, j));
        }
      }
      return beta2;
    };

    // Raman Intensity linearly polarized light
    auto RamanIntensityL = [](const double alpha2, const double beta2) {
      return 45 * alpha2 + 4 * beta2; // in a.u.
    };
    auto DepolarizationRatio = [](const double alpha2, const double beta2) {
      return 3 * beta2 / (45 * alpha2 + 4 * beta2);
    };
    bool debug = ctx.world.rank() == 0 && ctx.response_params.print_level() > 1;
    if (debug) {
      for (size_t freq_idx = 0;
           freq_idx < raman.polarization_frequencies.size(); ++freq_idx) {
        print("RAMAN_ALPHA_DERIVATIVES frequency=",
              raman.polarization_frequencies[freq_idx], " (a.u.): \n",
              alpha_derivatives[freq_idx]);
      };
    }
    // Compute Normal Mode
    for (size_t freq_idx = 0; freq_idx < raman.polarization_frequencies.size();
         ++freq_idx) {

      const auto &alpha_dxyz = alpha_derivatives[freq_idx];
      double pol_freq = raman.polarization_frequencies[freq_idx];
      // Save polarizability derivatives to Results
      raman.polarizability_derivatives.push_back(alpha_dxyz);
      auto alpha_qi = inner(alpha_dxyz, nnmodes);
      raman.polarizability_derivatives_normal_modes.push_back(alpha_qi);
      if (debug) {
        print("RAMAN_ALPHA_DERIVATIVES_NORMAL_MODES frequency=",
              pol_freq, " (a.u.): \n", alpha_qi);
      }
    }
    using RamanModeRow = RamanResults::RamanModeRow;

    for (size_t freq_idx = 0; freq_idx < raman.polarization_frequencies.size();
         ++freq_idx) {
      auto alpha_qi = raman.polarizability_derivatives_normal_modes[freq_idx];
      vector<RamanModeRow> raman_rows;
      for (size_t i = 0; i < mode.size(); ++i) {
        RamanModeRow row;
        double vib_freq_i = vib_freq(mode[i]);
        row.mode = static_cast<int>(i) + 1;
        row.freq_cm1 = vib_freq_i;
        auto alpha_i = copy(alpha_qi(_, i));
        auto alpha = alpha_i.reshape(3, 3);
        auto alpha2_au = compute_alpha2(alpha);
        row.alpha2 = alpha2_au * csg_factor;
        auto beta2_au = compute_beta2(alpha);
        row.beta2 = beta2_au * csg_factor;
        auto depol_int_au = RamanIntensityL(alpha2_au, beta2_au);
        row.pol_int = depol_int_au * csg_factor;
        row.dep_ratio = DepolarizationRatio(alpha2_au, beta2_au);
        row.depol_int = (*row.dep_ratio) * (*row.pol_int);
        ctx.world.gop.fence();
        raman_rows.push_back(row);
      }
      raman.raman_spectra[raman.polarization_frequencies[freq_idx]] =
          raman_rows;
      ctx.world.gop.fence();
    }
    print_raman_table(ctx.world, raman, ctx.response_params.print_level());
  }

  /// Serial stage-3 property dispatch over requested property list.
  ///
  /// This is the baseline assembly path and is also used as fallback when any
  /// subgroup property path fails.
  static PropertyStageOutput compute_requested_properties(
      World &world, const ResponseParameters &response_params,
      const GroundContext &ground_ctx, const SolvedStates &solved_states,
      const std::shared_ptr<SCF> &scf_calc) {
    VibrationalResults vib;
    RamanResults raman;
    PropertyManager properties(world, "properties.json");

    PropertyContext prop_ctx{world,         response_params, ground_ctx,
                             solved_states, properties,      scf_calc};

    for (const std::string &prop : response_params.requested_properties()) {
      auto prop_type = parse_property_name(prop);
      switch (prop_type) {
      case PropertyType::Alpha:
        compute_polarizability(prop_ctx);
        break;
      case PropertyType::Beta:
        compute_hyperpolarizability(prop_ctx);
        break;
      case PropertyType::Raman:
        compute_raman(prop_ctx, vib, raman);
        break;
      }
    }

    return PropertyStageOutput{properties.to_json(), std::move(vib),
                               std::move(raman)};
  }

  /// Subgroup-aware property stage with robust fallbacks.
  ///
  /// When subgroup mode is active this can precompute beta/Raman components in
  /// distributed shards, then execute final property assembly on a selected
  /// property subgroup and broadcast results back to the world.
  static PropertyStageOutput compute_requested_properties_with_property_group(
      World &world, const CalculationParameters &calc_params,
      const ResponseParameters &response_params,
      const GroundContext &ground_ctx, const SolvedStates &solved_states,
      const std::shared_ptr<SCF> &scf_calc) {
    const auto &state_parallel_plan =
        solved_states.planned_states.state_parallel_plan;
    const bool subgroup_property_mode =
        state_parallel_plan.subgroup_parallel_enabled &&
        state_parallel_plan.execution_groups > 1;

    if (!subgroup_property_mode) {
      return compute_requested_properties(world, response_params, ground_ctx,
                                          solved_states, scf_calc);
    }

    const size_t property_group =
        response_params.state_parallel_property_group();
    if (property_group >= state_parallel_plan.execution_groups) {
      if (world.rank() == 0) {
        print("WARN PROPERTY_SUBGROUP_INVALID subgroup=", property_group,
              " execution_groups=",
              state_parallel_plan.execution_groups,
              " action=fallback_world_property_stage");
      }
      return compute_requested_properties(world, response_params, ground_ctx,
                                          solved_states, scf_calc);
    }

    const auto requested_has = [&](PropertyType needle) {
      for (const auto &prop : response_params.requested_properties()) {
        if (parse_property_name(prop) == needle) {
          return true;
        }
      }
      return false;
    };

    const bool request_beta_components = requested_has(PropertyType::Beta);
    const bool request_raman_components = requested_has(PropertyType::Raman);

    if (request_beta_components || request_raman_components) {
      if (world.rank() == 0) {
        print("PROPERTY_COMPONENT_PRECOMPUTE enabled=true beta=",
              request_beta_components ? "on" : "off",
              ", raman=", request_raman_components ? "on" : "off",
              " execution_groups=", state_parallel_plan.execution_groups);
      }

      std::string component_claim_prefix;
      if (world.rank() == 0) {
        component_claim_prefix =
            (std::filesystem::path("property_component_claims") /
             ("run_" + iso_timestamp()))
                .string();
      }
      world.gop.broadcast_serializable(component_claim_prefix, 0);

      bool local_component_stage_failed = false;
      std::string local_component_stage_error;
      try {
        auto component_subworld_ptr = MacroTaskQ::create_worlds(
            world, state_parallel_plan.execution_groups);
        if (!component_subworld_ptr) {
          throw std::runtime_error(
              "property component subworld creation returned null");
        }

        World &component_subworld = *component_subworld_ptr;
        const size_t component_subgroup_id = static_cast<size_t>(
            world.rank() %
            static_cast<int>(state_parallel_plan.execution_groups));

        auto old_component_pmap3 = FunctionDefaults<3>::get_pmap();
        auto restore_component_pmap = [&]() {
          FunctionDefaults<3>::set_pmap(old_component_pmap3);
        };

        bool local_subgroup_component_failed = false;
        std::string local_subgroup_component_error;

        FunctionDefaults<3>::set_default_pmap(component_subworld);
        try {
          const auto &protocol = calc_params.protocol();
          const double final_thresh = protocol.empty()
                                          ? FunctionDefaults<3>::get_thresh()
                                          : protocol.back();
          const std::string component_fock_file =
              group_shard_file(ground_ctx.fock_json_file, component_subgroup_id);
          const std::string component_shard_file =
              group_shard_file("properties_components.json", component_subgroup_id);

          if (component_subworld.rank() == 0) {
            std::error_code ec;
            std::filesystem::remove(component_shard_file, ec);
          }
          component_subworld.gop.fence();

          GroundStateData component_ground(component_subworld,
                                           ground_ctx.archive_file,
                                           ground_ctx.molecule);
          ResponseManager component_response_manager(component_subworld,
                                                     calc_params);
          component_response_manager.setProtocol(component_subworld,
                                                 component_ground.getL(),
                                                 final_thresh);
          component_ground.prepareOrbitals(component_subworld,
                                          FunctionDefaults<3>::get_k(),
                                          final_thresh);
          component_ground.computePreliminaries(
              component_subworld, *component_response_manager.getCoulombOp(),
              component_response_manager.getVtol(), component_fock_file);

          PropertyManager component_properties(component_subworld,
                                              component_shard_file);

          if (request_beta_components) {
            ::compute_hyperpolarizability(
                component_subworld, component_ground,
                response_params.dipole_frequencies(),
                response_params.dipole_directions(), component_properties,
                response_params.beta_shg(), response_params.beta_or(),
                response_params.beta_all_triplets(),
                component_claim_prefix + ".beta", component_subgroup_id);
          }
          if (request_raman_components) {
            ::compute_Raman_components(
                component_subworld, component_ground,
                response_params.dipole_frequencies(),
                response_params.dipole_directions(),
                response_params.nuclear_directions(), component_properties,
                component_claim_prefix + ".raman", component_subgroup_id);
          }
          component_properties.save();
        } catch (const std::exception &ex) {
          local_subgroup_component_failed = true;
          local_subgroup_component_error = ex.what();
        } catch (...) {
          local_subgroup_component_failed = true;
          local_subgroup_component_error =
              "unknown exception during component precompute";
        }

        long subgroup_component_failed_flag =
            local_subgroup_component_failed ? 1 : 0;
        component_subworld.gop.max(subgroup_component_failed_flag);
        restore_component_pmap();

        if (subgroup_component_failed_flag != 0) {
          local_component_stage_failed = true;
          if (local_subgroup_component_failed &&
              !local_subgroup_component_error.empty()) {
            local_component_stage_error =
                std::move(local_subgroup_component_error);
          }
        }
      } catch (const std::exception &ex) {
        local_component_stage_failed = true;
        local_component_stage_error = ex.what();
      } catch (...) {
        local_component_stage_failed = true;
        local_component_stage_error =
            "unknown exception during component stage setup";
      }

      long any_component_failure = local_component_stage_failed ? 1 : 0;
      world.gop.max(any_component_failure);
      if (any_component_failure != 0) {
        if (world.rank() == 0) {
          if (!local_component_stage_error.empty()) {
            print("ERROR PROPERTY_COMPONENT_STAGE_FAILED message=",
                  local_component_stage_error,
                  " action=fallback_world_property_stage");
          } else {
            print("ERROR PROPERTY_COMPONENT_STAGE_FAILED "
                  "message=rank_failure action=fallback_world_property_stage");
          }
        }
        return compute_requested_properties(world, response_params, ground_ctx,
                                            solved_states, scf_calc);
      }

      if (world.rank() == 0) {
        nlohmann::json merged_components = nlohmann::json::array();
        const auto baseline = read_json_file_or_object("properties.json");
        if (baseline.is_array()) {
          for (const auto &row : baseline) {
            merged_components.push_back(row);
          }
        }
        for (size_t gid = 0; gid < state_parallel_plan.execution_groups; ++gid) {
          const auto shard = read_json_file_or_object(
              group_shard_file("properties_components.json", gid));
          if (!shard.is_array()) {
            continue;
          }
          for (const auto &row : shard) {
            merged_components.push_back(row);
          }
        }
        write_json_file("properties.json", merged_components);
      }
      world.gop.fence();
    }

    if (world.rank() == 0) {
      print("PROPERTY_SUBGROUP_EXECUTE subgroup=", property_group,
            " subgroup_max=", state_parallel_plan.execution_groups - 1);
    }

    std::string properties_dump;
    std::string vib_dump;
    std::string raman_dump;
    bool local_property_stage_failed = false;
    std::string local_property_stage_error;

    try {
      auto subworld_ptr = MacroTaskQ::create_worlds(
          world, state_parallel_plan.execution_groups);
      if (!subworld_ptr) {
        throw std::runtime_error("property subworld creation returned null");
      }

      World &subworld = *subworld_ptr;
      const size_t subgroup_id = static_cast<size_t>(
          world.rank() %
          static_cast<int>(state_parallel_plan.execution_groups));
      const bool in_property_group = (subgroup_id == property_group);

      auto old_pmap3 = FunctionDefaults<3>::get_pmap();
      auto restore_pmap = [&]() { FunctionDefaults<3>::set_pmap(old_pmap3); };

      bool local_subgroup_failed = false;
      std::string local_subgroup_error;

      // Property assembly uses 3D response objects; only swap the 3D pmap.
      FunctionDefaults<3>::set_default_pmap(subworld);
      try {
        if (in_property_group) {
          const auto &protocol = calc_params.protocol();
          const double final_thresh = protocol.empty()
                                          ? FunctionDefaults<3>::get_thresh()
                                          : protocol.back();
          const std::string property_fock_file =
              group_shard_file(ground_ctx.fock_json_file, property_group);

          GroundStateData local_ground(subworld, ground_ctx.archive_file,
                                       ground_ctx.molecule);
          ResponseManager local_response_manager(subworld, calc_params);
          local_response_manager.setProtocol(subworld, local_ground.getL(),
                                             final_thresh);
          local_ground.prepareOrbitals(subworld, FunctionDefaults<3>::get_k(),
                                       final_thresh);
          local_ground.computePreliminaries(
              subworld, *local_response_manager.getCoulombOp(),
              local_response_manager.getVtol(), property_fock_file);

          GroundContext local_ground_ctx{
              ground_ctx.molecule, std::move(local_ground),
              std::move(local_response_manager), ground_ctx.archive_file,
              property_fock_file};

          PropertyStageOutput subgroup_output = compute_requested_properties(
              subworld, response_params, local_ground_ctx, solved_states,
              scf_calc);
          if (subworld.rank() == 0) {
            properties_dump = subgroup_output.properties.dump();
            vib_dump = subgroup_output.vibrational_analysis.to_json().dump();
            raman_dump = subgroup_output.raman_spectra.to_json().dump();
          }
        }
      } catch (const std::exception &ex) {
        local_subgroup_failed = true;
        local_subgroup_error = ex.what();
      } catch (...) {
        local_subgroup_failed = true;
        local_subgroup_error =
            "unknown exception during subgroup property execution";
      }

      long subgroup_failed_flag = local_subgroup_failed ? 1 : 0;
      subworld.gop.max(subgroup_failed_flag);
      restore_pmap();

      if (subgroup_failed_flag != 0) {
        local_property_stage_failed = true;
        if (local_subgroup_failed && !local_subgroup_error.empty()) {
          local_property_stage_error = std::move(local_subgroup_error);
        }
      }
    } catch (const std::exception &ex) {
      local_property_stage_failed = true;
      local_property_stage_error = ex.what();
    } catch (...) {
      local_property_stage_failed = true;
      local_property_stage_error =
          "unknown exception during property subgroup setup";
    }

    long any_property_failure = local_property_stage_failed ? 1 : 0;
    world.gop.max(any_property_failure);
    if (any_property_failure != 0) {
      if (world.rank() == 0) {
        if (!local_property_stage_error.empty()) {
          print("ERROR PROPERTY_SUBGROUP_FAILED message=",
                local_property_stage_error,
                " action=fallback_world_property_stage");
        } else {
          print("ERROR PROPERTY_SUBGROUP_FAILED message=rank_failure "
                "action=fallback_world_property_stage");
        }
      }
      return compute_requested_properties(world, response_params, ground_ctx,
                                          solved_states, scf_calc);
    }

    const int property_world_root = static_cast<int>(property_group);
    world.gop.broadcast_serializable(properties_dump, property_world_root);
    world.gop.broadcast_serializable(vib_dump, property_world_root);
    world.gop.broadcast_serializable(raman_dump, property_world_root);

    PropertyStageOutput stage_output;
    stage_output.properties = properties_dump.empty()
                                  ? nlohmann::json::object()
                                  : nlohmann::json::parse(properties_dump);
    if (!vib_dump.empty()) {
      stage_output.vibrational_analysis.from_json(
          nlohmann::json::parse(vib_dump));
    }
    if (!raman_dump.empty()) {
      stage_output.raman_spectra.from_json(nlohmann::json::parse(raman_dump));
    }
    world.gop.fence();
    return stage_output;
  }

  // ============================================================================
  // SECTION 13: Top-level orchestration
  //   solve_all_states — Stage 2 master function: builds schedule context,
  //                      runs linear solve (subgroup or serial), calls excited-
  //                      bundle stage, calls derived-state stage, and assembles
  //                      the SolvedStates bundle consumed by stage 3.
  //   prepare_and_validate_final_protocol_state — rebuilds final-threshold
  //                      context and reports unconverged-point count.
  //   build_state_stage_metadata — assembles the combined stage-2 metadata
  //                      JSON from planner/runtime/excited/derived outputs.
  //   run_response (public) — entry point called by WorkflowBuilders; drives
  //                      stages 1 → 2 → 3, finalises stats, returns Results.
  // ============================================================================

public:
  /**
   * @brief Run the full molecular response & property workflow.
   *
   * 1. Builds contexts and plans states
   * 2. Solves all linear and derived states
   * 3. Computes requested properties
   *
   * @param world      The MADNESS world communicator
   * @param params     Unified parameters containing response and molecule info
   * @param scf_calc   Optional SCF handle used by Raman/Hessian pathways
   * @param outdir     Directory where all outputs will be written
   * @return Results   Structured JSON fragments: metadata + properties
   */
  inline static Results run_response(World &world, const Params &params,
                                     const std::shared_ptr<SCF> &scf_calc,
                                     const std::filesystem::path &outdir) {
    const auto &calc_params = params.get<CalculationParameters>();
    const auto &rp_copy = params.get<ResponseParameters>();

    if (world.rank() == 0) {
      json response_input_json = {};
      response_input_json["response"] =
          rp_copy.to_json_if_precedence("defined");
      print("RESPONSE_INPUT_JSON json=", response_input_json.dump(4));
      std::ofstream ofs("response.in");
      write_json_to_input_file(response_input_json, {"response"}, ofs);
      ofs.close();
    }
    world.gop.fence();
    commandlineparser parser;
    parser.set_keyval("input", "response.in");
    if (world.rank() == 0) {
      ::print("RESPONSE_INPUT_FILE path=", parser.value("input"));
    }

    auto response_params = ResponseParameters(world, parser);
    GroundContext ground_ctx =
        make_ground_context(world, calc_params, scf_calc, outdir);

    PlannedStates planned_states =
        plan_required_states(world, calc_params, ground_ctx, response_params);
    SolvedStates solved_states =
        solve_all_states(world, calc_params, ground_ctx, response_params,
                         std::move(planned_states));
    PropertyStageOutput property_stage =
        compute_requested_properties_with_property_group(
            world, calc_params, response_params, ground_ctx, solved_states,
            scf_calc);

    // finalize & stats
    if (world.rank() == 0) {
      madness::print("WORKFLOW_COMPLETE stage=molresponse_property");
    }
    world.gop.fence();
    world.gop.fence();
    print_stats(world);

    // aggregate JSON results
    Results results;
    results.metadata = std::move(solved_states.metadata);
    results.properties = std::move(property_stage.properties);
    results.debug_log = std::move(solved_states.debug_log);
    results.vibrational_analysis =
        property_stage.vibrational_analysis.to_json();
    results.raman_spectra = property_stage.raman_spectra.to_json();
    return results;
  }
}; // namespace molresponse_lib
