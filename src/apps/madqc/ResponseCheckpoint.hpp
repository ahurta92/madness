#pragma once

#include "../molresponse/solvers/function_kinetics.hpp"

#include "NamedStages.hpp"
#include <apps/molresponse/solvers/fd_save_load.hpp>
#include <madness/chem/StageCheckpoint.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <map>
#include <set>

namespace named_stages {

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

  static bool request_complete(const json& metadata, const json& request) {
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
      require(m.properties.size() == 1 &&
              request_complete(m.numerical.at("engine_metadata"), m.numerical.at("request")),
              "completed property lacks converged prerequisites");
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
    loaded_ = latest_.stop_reason == "converged" &&
              latest_.numerical.at("request") == request_ && !latest_.properties.empty();
    if (loaded_) { work = "loaded"; return true; }

    auto metadata = latest_.numerical.at("engine_metadata");
    // Analysis is derived: never let old property rows certify this request.
    metadata["properties"] = json::object();
    prepare_states(metadata);
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
  json metadata() const { return latest_.numerical.at("engine_metadata"); }
  json properties() const { return metadata().at("properties"); }

  // Called collectively after native FD save and metadata update have finished.
  void observe(madness::World&, const json& event) {
    iterations += event.at("iterations").get<std::size_t>();
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
  // Engine filenames are rounded keys. Compare exact definitions before giving
  // any entry to its restart planner, including sources at coarser rungs.
  void prepare_states(json& metadata) const {
    std::map<std::string, double> frequencies, protocols;
    std::set<std::string> channels;
    for (double f : stage_.frequencies)
      frequencies.emplace(molresponse_v3::ResponseMetadata::freq_key(f), f);
    for (double t : stage_.params.get<CalculationParameters>().protocol())
      protocols.emplace(molresponse_v3::protocol_key(t, protocol_k(t)), t);
    for (const auto& component : stage_.components)
      for (char axis : component) channels.insert(std::string("dipole_") + axis);
    const double target = stage_.params.get<ResponseParameters>().dconv();
    for (auto& channel : metadata["fd_states"].items()) {
      for (auto rung = channel.value().begin(); rung != channel.value().end();) {
        const auto wanted = protocols.find(rung.key());
        if (wanted != protocols.end() &&
            metadata.at("protocols").at(rung.key()).at("thresh").get<double>() != wanted->second) {
          // The engine protocol registry has one definition per key. An alias
          // cannot coexist in this attempt; its immutable old generation remains.
          rung = channel.value().erase(rung);
          continue;
        }
        for (auto point = rung.value().begin(); point != rung.value().end();) {
          auto& e = point.value();
          const auto frequency = frequencies.find(point.key());
          const bool requested = channels.count(channel.key()) && frequency != frequencies.end();
          if (requested && e.at("freq").get<double>() != frequency->second) {
            point = rung.value().erase(point);
            continue;
          }
          // A tighter subset must not erase independent convergence evidence.
          if (requested && e.at("checkpoint_dconv").get<double>() > target)
            e["converged"] = false;
          if (!e.at("converged").get<bool>()) {
            e["iteration_at_checkpoint"] = e.at("iter");
            e["iter"] = 0;
            e["stalled"] = false;
          }
          ++point;
        }
        ++rung;
      }
    }
    // save_fd_state registers only new keys, so replace an aliased registry
    // definition here after every old state at that key has been removed.
    for (const auto& [key, threshold] : protocols)
      if (metadata["protocols"].contains(key)) {
        metadata["protocols"][key]["thresh"] = threshold;
        metadata["protocols"][key]["k"] = protocol_k(threshold);
      }
  }

  json read_metadata() const {
    const auto text = madness::checkpoint::detail::on_root(world_, [&] {
      return madness::checkpoint::detail::read_json(attempt_ / "response_metadata.json").dump();
    });
    return json::parse(text);
  }

  void publish(json metadata, bool analyzed, bool converged, const json& position) {
    namespace checkpoint = madness::checkpoint;
    if (!analyzed) metadata["properties"] = json::object();
    Manifest m;
    m.stage = stage_.name;
    m.method = "TDHF";
    m.generation = "g-" + stage_.name + "-" + run_id_ + "-" + std::to_string(++sequence_);
    m.physical = {{"schema", "madqc.response.v1"}, {"shell", "closed_shell"}};
    m.reference = reference_;
    m.numerical = {{"request", request_}, {"engine_metadata", metadata}};
    m.targets = {{"dconv", request_.at("dconv")}};
    m.position = position;
    m.position["iterations_this_run"] = iterations;
    if (!latest_.generation.empty()) m.position["previous_generation"] = latest_.generation;
    m.stop_reason = converged ? "converged" : analyzed ? "unfinished" : "in_progress";
    if (analyzed)
      m.properties.push_back({{"name", "polarizability"}, {"generation", m.generation},
                              {"converged", converged}, {"value", metadata.at("properties")}});
    const auto directory = store_ / m.generation;
    const auto names = members(metadata);
    if (converged)
      require(request_complete(metadata, request_), "analysis lacks converged prerequisites");
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
