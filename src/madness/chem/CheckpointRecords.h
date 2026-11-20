#pragma once

// Shared task evidence only. Native archive inventory/validation belongs to
// StageCheckpoint; scientific definitions and numerical state belong to apps.
#include <madness/chem/CheckpointRetention.h>
#include <madness/chem/IterationMetrics.h>
#include <nlohmann/json.hpp>
#include <set>

namespace madness::checkpoint::records {
using json = nlohmann::json;
namespace it = madness::iteration;
inline void require(bool valid, const std::string &why) {
  if (!valid)
    throw std::runtime_error("checkpoint record: " + why);
}
inline void identifier(const std::string &id) {
  require(!id.empty() && id.front() != '.' &&
              id.find('/') == std::string::npos &&
              id.find('\0') == std::string::npos && id != "latest",
          "unsafe identity: " + id);
}
inline std::string id(const json &value) {
  auto text = value.get<std::string>();
  identifier(text);
  return text;
}
// Walk the small, already parsed JSON record.
// NOLINTNEXTLINE(misc-no-recursion)
inline void finite_numbers(const json &record) {
  if (record.is_number_float())
    require(std::isfinite(record.get<double>()), "nonfinite JSON number");
  if (record.is_structured())
    for (const auto &item : record)
      finite_numbers(item);
}
inline std::size_t count(const json &value) {
  require(value.is_number_integer(), "count is not an integer");
  if (!value.is_number_unsigned())
    require(value.get<std::int64_t>() >= 0, "negative count");
  return value.get<std::size_t>();
}
inline it::Protocol read_protocol(const json &value) {
  const auto k = count(value.at("k"));
  require(value.at("threshold").is_number() &&
              std::isfinite(value.at("threshold").get<double>()) &&
              value.at("threshold").get<double>() > 0 && k > 0 &&
              k <= static_cast<std::size_t>(std::numeric_limits<int>::max()),
          "invalid protocol");
  return {value.at("threshold").get<double>(), static_cast<int>(k)};
}
inline json protocol_record(const it::Protocol &value) {
  return {{"threshold", value.thresh}, {"k", value.k}};
}
inline it::Comparison read_comparison(const json &value) {
  return {id(value.at("task_id")), id(value.at("attempt_id")),
          read_protocol(value.at("protocol")), count(value.at("iteration")),
          id(value.at("subject"))};
}
inline json comparison_record(const it::Comparison &value) {
  return {{"task_id", value.task_id},
          {"attempt_id", value.attempt_id},
          {"protocol", protocol_record(value.protocol)},
          {"iteration", value.iteration},
          {"subject", value.subject}};
}
inline GenerationRef reference(const json &entry) {
  return {id(entry.at("calculation_id")), id(entry.at("task_id")),
          id(entry.at("generation"))};
}
inline std::vector<GenerationRef> inputs(const json &list) {
  require(list.is_array(), "inputs must be an array");
  std::vector<GenerationRef> result;
  for (const auto &entry : list) {
    require(entry.at("role").is_string() &&
                !entry.at("role").get<std::string>().empty(),
            "missing input role");
    result.push_back(reference(entry));
  }
  return result;
}
inline void outcome(const json &record) {
  require(record.at("complete").is_boolean() &&
              record.at("reason").is_string() &&
              !record.at("reason").get<std::string>().empty(),
          "invalid outcome");
}
inline std::set<std::string> subjects(const json &definition) {
  std::set<std::string> result;
  if (definition.contains("subjects")) {
    require(!definition.contains("subject") &&
                definition.at("subjects").is_array(),
            "ambiguous subjects");
    for (const auto &item : definition.at("subjects"))
      require(result.insert(id(item)).second, "duplicate subject");
  } else
    result.insert(id(definition.at("subject")));
  require(!result.empty(), "empty subjects");
  return result;
}
inline void validate(const json &record) {
  finite_numbers(record);
  require(record.at("version") == 2, "unsupported task record version");
  id(record.at("calculation_id"));
  id(record.at("attempt_id"));
  id(record.at("generation"));
  const auto &stage = record.at("stage");
  require(stage.at("application").is_string() &&
              !stage.at("application").get<std::string>().empty(),
          "missing application");
  require(stage.at("settings").is_object() &&
              stage.at("requested_properties").is_array(),
          "invalid stage snapshot");
  inputs(stage.at("inputs"));
  outcome(stage.at("outcome"));
  const auto &task = record.at("task");
  id(task.at("task_id"));
  inputs(task.at("inputs"));
  outcome(task.at("outcome"));
  require(task.at("definition").is_object(), "invalid task definition");
  const auto &definitions = record.at("metric_definitions");
  require(definitions.is_object(), "invalid metric definitions");
  if (task.at("kind") == "direct") {
    require(!task.contains("iteration"),
            "direct task cannot have iteration evidence");
    return;
  }
  require(task.at("kind") == "iterative", "unknown task kind");
  const auto &definition = task.at("definition");
  require(!definition.at("context_id").get<std::string>().empty(),
          "missing scientific context identity");
  const auto protocol = read_protocol(definition.at("protocol"));
  const auto expected = subjects(definition);
  const auto &iteration = task.at("iteration");
  require(!iteration.contains("sample_id"), "unsupported sample_id field");
  require(count(iteration.at("updates_this_attempt")) <=
              count(iteration.at("updates_at_protocol")),
          "inconsistent update counts");
  require(iteration.at("protocol_complete").is_boolean(),
          "invalid protocol completion");
  const auto &metrics = iteration.at("metrics");
  require(metrics.is_array(), "metrics must be an array");
  std::set<std::pair<std::string, std::string>> pairs;
  std::set<std::string> present;
  for (const auto &metric : metrics) {
    const auto name = metric.at("metric").get<std::string>();
    const auto kind = it::metric_from_name(name);
    const auto subject = metric.contains("subject") ? id(metric.at("subject"))
                         : expected.size() == 1     ? *expected.begin()
                                                    : std::string{};
    require(expected.contains(subject) && pairs.emplace(subject, name).second,
            "unknown subject or duplicate metric");
    present.insert(subject);
    const auto definition_id = metric.at("definition").get<std::string>();
    const auto &meaning = definitions.at(definition_id);
    for (const auto *field :
         {"metric", "meaning", "units", "fields", "weights", "normalization",
          "reduction", "comparison", "scale"})
      require(meaning.contains(field) && !meaning.at(field).is_null(),
              "incomplete metric definition");
    require(meaning.at("metric") == name &&
                !meaning.at("units").get<std::string>().empty(),
            "metric definition mismatch");
    for (const auto *forbidden : {"absolute", "relative", "limit", "tolerance",
                                  "satisfied", "result", "target"})
      require(!metric.contains(forbidden),
              "request-dependent check in checkpoint metrics");
    const auto status = metric.at("status").get<std::string>();
    if (status == "measured") {
      require(metric.at("value").is_number() &&
                  it::finite_nonnegative(metric.at("value").get<double>()),
              "invalid measured value");
      if (kind != it::Metric::bsh_residual)
        require(metric.contains("comparison"), "missing comparison");
    } else {
      require(status == "unavailable" || status == "not_applicable" ||
                  status == "invalid",
              "unknown metric status");
      require(!metric.contains("value") &&
                  !metric.at("reason").get<std::string>().empty(),
              "invalid unavailable metric");
    }
    if (metric.contains("comparison")) {
      const auto comparison = read_comparison(metric.at("comparison"));
      require(comparison.task_id == task.at("task_id").get<std::string>(),
              "comparison belongs to another task");
      require(comparison.attempt_id !=
                      record.at("attempt_id").get<std::string>() ||
                  comparison.protocol != protocol ||
                  comparison.iteration <
                      count(iteration.at("updates_this_attempt")),
              "comparison is not an earlier iteration");
    }
    if (metric.contains("scale")) {
      const auto &scale = metric.at("scale");
      require(status != "invalid" && status != "not_applicable",
              "inapplicable metric has a scale");
      if (scale.at("status") == "measured")
        require(it::finite_nonnegative(scale.at("value").get<double>()),
                "invalid scale");
      else
        require(scale.at("status") == "unavailable" &&
                    !scale.contains("value") &&
                    !scale.at("reason").get<std::string>().empty(),
                "invalid unavailable scale");
    }
  }
  require(present == expected,
          "checkpoint omits a required subject's evidence");
}
inline json measurements(const it::IterationMetrics &metrics,
                         const std::map<it::Metric, std::string> &definitions) {
  json entries = json::array();
  for (const auto &subject : metrics.states)
    for (const auto &[kind, measurement] : subject.metrics) {
      json entry = {{"metric", it::metric_name(kind)},
                    {"subject", subject.id},
                    {"definition", definitions.at(kind)}};
      switch (measurement.status) {
      case it::Availability::measured:
        entry["status"] = "measured";
        require(measurement.value.has_value(), "missing measured value");
        entry["value"] = *measurement.value;
        break;
      case it::Availability::unavailable:
        entry["status"] = "unavailable";
        break;
      case it::Availability::not_applicable:
        entry["status"] = "not_applicable";
        break;
      case it::Availability::invalid:
        entry["status"] = "invalid";
        break;
      }
      if (measurement.status != it::Availability::measured)
        entry["reason"] = measurement.reason;
      const auto &comparison =
          measurement.comparison ? measurement.comparison : subject.comparison;
      if (comparison && kind != it::Metric::bsh_residual)
        entry["comparison"] = comparison_record(*comparison);
      if (measurement.status == it::Availability::measured ||
          measurement.status == it::Availability::unavailable) {
        if (measurement.scale.status == it::Availability::measured &&
            measurement.scale.value &&
            it::finite_nonnegative(*measurement.scale.value))
          entry["scale"] = {{"status", "measured"},
                            {"value", *measurement.scale.value}};
        else
          entry["scale"] = {{"status", "unavailable"},
                            {"reason", measurement.scale.reason.empty()
                                           ? "invalid_scale"
                                           : measurement.scale.reason}};
      }
      const auto &evidence = measurement.evidence;
      json values = json::object();
      if (evidence.current)
        values["current"] = *evidence.current;
      if (evidence.previous)
        values["previous"] = *evidence.previous;
      if (evidence.current_norm)
        values["current_norm"] = *evidence.current_norm;
      if (evidence.previous_norm)
        values["previous_norm"] = *evidence.previous_norm;
      if (!evidence.components.empty()) {
        values["components"] = json::array();
        for (const auto &component : evidence.components) {
          json item = {{"id", component.id},
                       {"residual_norm", component.residual_norm}};
          if (component.function_norm)
            item["function_norm"] = *component.function_norm;
          values["components"].push_back(std::move(item));
        }
      }
      if (evidence.limiting_id)
        values["limiting_id"] = *evidence.limiting_id;
      if (!values.empty())
        entry["evidence"] = std::move(values);
      entries.push_back(std::move(entry));
    }
  return entries;
}
} // namespace madness::checkpoint::records
