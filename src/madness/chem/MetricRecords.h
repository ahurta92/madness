#pragma once

// JSON records of iteration measurements, shared by the node archive and the
// history file. Pure: no World, MPI, archives or files.
#include <madness/chem/IterationMetrics.h>
#include <nlohmann/json.hpp>
#include <cmath>
#include <limits>
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
// The earlier state a change was measured against: its task, level, update
// number and subject.
inline json comparison_record(const it::Comparison &value) {
  return {{"task_id", value.task_id},
          {"protocol", protocol_record(value.protocol)},
          {"iteration", value.iteration},
          {"subject", value.subject}};
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
