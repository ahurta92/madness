#include "../solvers/iterate.hpp"
#include <madness/chem/ConvergenceSettings.h>
#include <stdexcept>
#include <vector>
#include <sstream>

namespace it = madness::iteration;
struct Solver {
  struct State {
    int iter = 0;
    bool diverged = false;
    double x = 0.0;
  };
  int evaluations = 0, updates = 0;
  bool proposal = false;
  void reset_shared() { proposal = false; }
  it::IterationMetrics evaluate(State &state, const it::IterationInfo &info) {
    ++evaluations;
    proposal = true;
    it::StateMetrics metrics;
    metrics.id = "response";
    metrics.metrics[it::Metric::bsh_residual] =
        it::MetricValue::measured(1.0 - state.x, 1.0);
    return {info, {metrics}};
  }
  void advance(State &state) {
    if (!proposal)
      throw std::runtime_error("raw work rebuilt or discarded");
    proposal = false;
    ++updates;
    ++state.iter;
    state.x += 0.75;
  }
};
void require(bool value) {
  if (!value)
    throw std::runtime_error("shared iteration contract failed");
}
struct ManualClock {
  using duration = std::chrono::milliseconds;
  using time_point = std::chrono::time_point<ManualClock>;
  static inline time_point current{};
  static time_point now() { return current; }
};
int main() {
  it::ElapsedTimer<ManualClock> timer;
  ManualClock::current += std::chrono::milliseconds(125);
  require(timer.lap() == 0.125);
  ManualClock::current += std::chrono::milliseconds(250);
  require(timer.lap() == 0.25 && timer.seconds() == 0.375);

  const it::ExpectedTask expected{"fd", "context", {"response"}};
  const it::IterationInfo info{"fd", "attempt", {1e-4, 6}, 0, "context"};
  it::Request request;
  request.criteria = {it::Metric::bsh_residual};
  request.max_updates = 2;
  request.tolerances[it::Metric::bsh_residual] = {0.5, 0.0};
  const it::ModelMetrics model{
      {it::Metric::bsh_residual, {it::Support::supported, {}, {}}}};
  Solver solver;
  Solver::State state;
  std::vector<double> samples;
  auto observe = [&](const Solver::State &s, const it::IterationMetrics &m,
                     const it::Decision &decision, const it::Progress &) {
    require(*m.states[0].metrics.at(it::Metric::bsh_residual).value ==
            1.0 - s.x);
    if (decision.checkpoint)
      samples.push_back(s.x);
  };
  std::vector<it::IterationTiming> timings;
  const auto timed = [&](const it::IterationInfo& current, const it::IterationTiming& timing) {
    require(current.iteration == timings.size());
    require(timing.seconds >= 0.0);
    for (const auto& step : timing.steps) require(step.second >= 0.0);
    timings.push_back(timing);
  };
  auto decision = molresponse_v3::solvers::iterate_shared(
      solver, state, it::resolve(request, model, info.protocol, true), expected,
      info, 0, observe, timed);
  require(decision.action == it::Action::complete && state.iter == 1);
  require(solver.evaluations == 2 && solver.updates == 1 &&
          samples == std::vector<double>{0.75});
  require(timings.size() == 2);
  require(timings[0].steps.back().first == "advance");
  require(timings[1].steps.back().first == "checkpoint/report");
  std::ostringstream text;
  it::write_iteration_timing(text, "response", info, timings.back());
  require(text.str().find("Solver update") == std::string::npos);
  require(text.str().find("seconds, driver rank") != std::string::npos);
  {
    auto request_with_missing = request;
    request_with_missing.criteria = {it::Metric::bsh_residual, it::Metric::density_change,
                                    it::Metric::energy_change};
    auto partial_model = model;
    partial_model[it::Metric::density_change] = {it::Support::supported, 0.1, {}};
    partial_model[it::Metric::energy_change] = {it::Support::not_applicable, {}, {}};
    it::StateMetrics partial;
    partial.id = "response";
    partial.metrics[it::Metric::bsh_residual] = it::MetricValue::measured(0.0);
    partial.metrics[it::Metric::energy_change] = it::MetricValue::measured(5.0);
    const it::IterationMetrics frame{info, {partial}};
    const auto controls = it::resolve(request_with_missing, partial_model, info.protocol, true);
    const auto result = it::decide(controls, expected, {0, 0, true}, frame);
    std::ostringstream printed;
    it::write_iteration_metrics(printed, "response", frame, result);
    require(printed.str().find("not measured") != std::string::npos);
    require(printed.str().find("model not applicable") != std::string::npos);
    require(printed.str().find("5.0000e+00") == std::string::npos);
    require(printed.str().find("1.0000e-01") != std::string::npos);
    require(printed.str().find("unavailable (required)") != std::string::npos);
    auto applicable_request = request;
    applicable_request.criteria.push_back(it::Metric::energy_change);
    const auto applicable_controls = it::resolve(applicable_request, partial_model, info.protocol, true);
    const auto complete = it::decide(applicable_controls, expected, {0, 0, true}, frame);
    require(complete.action == it::Action::complete);
    std::ostringstream completion, reuse;
    it::write_iteration_metrics(completion, "response", frame, complete);
    it::write_restart_summary(reuse, "response", frame, info.protocol,
        it::assess(applicable_controls, expected, frame), true);
    require(completion.str().find("applicable selected criteria") != std::string::npos);
    require(reuse.str().find("applicable selected tolerances") != std::string::npos);
    // A relative-only bound cannot be known without its scale; it must remain
    // visibly required even when no effective limit can be printed.
    request_with_missing.tolerances[it::Metric::density_change] = {0.0, 0.2};
    const auto relative_controls = it::resolve(request_with_missing, partial_model, info.protocol, true);
    const auto pending = it::decide(relative_controls, expected, {0, 0, true}, frame);
    std::ostringstream relative;
    it::write_iteration_metrics(relative, "response", frame, pending);
    require(relative.str().find("unavailable (required)") != std::string::npos);
  }
  {
    // Reporting must not hide a measured diagnostic just because it is not a gate.
    it::StateMetrics measured;
    measured.id = "a_very_long_response_state_name";
    measured.metrics[it::Metric::bsh_residual] = it::MetricValue::measured(1e-110, 2.0);
    measured.metrics[it::Metric::kinetic_change] = it::MetricValue::measured(1.25e110);
    auto relative_request = request;
    relative_request.tolerances[it::Metric::bsh_residual] = {0.1, 0.2};
    const it::IterationMetrics frame{info, {measured}};
    const it::ExpectedTask subject{"fd", "context", {measured.id}};
    const auto checked = it::decide(it::resolve(relative_request, model, info.protocol, true),
                                     subject, {0, 0, true}, frame);
    std::ostringstream printed;
    printed << std::fixed << std::setprecision(1);
    const auto flags = printed.flags();
    it::write_iteration_metrics(printed, "response", frame, checked);
    require(printed.flags() == flags && printed.precision() == 1);
    const auto output = printed.str();
    require(output.find("reported only") != std::string::npos);
    require(output.find("1.2500e+110") != std::string::npos);
    require(output.find("5.0000e-01") != std::string::npos); // absolute + relative * scale
    require(output.find("Measured") != std::string::npos);
    std::istringstream lines(output);
    std::string line;
    std::size_t limit_end = 0;
    while (std::getline(lines, line)) {
      if (line.find("Measured") != std::string::npos) limit_end = line.find("Limit") + 5;
      if (line.find("5.0000e-01") != std::string::npos)
        require(line.find("5.0000e-01") + 10 == limit_end);
    }
    require(limit_end > 0);
  }
  {
    it::StateMetrics first, second;
    first.id = "a long coupled root identifier";
    second.id = "root2";
    first.metrics[it::Metric::bsh_residual] = it::MetricValue::measured(1e-5);
    second.metrics[it::Metric::bsh_residual] = it::MetricValue::measured(1e-2);
    first.metrics[it::Metric::density_change] = it::MetricValue::unavailable("no_comparison");
    second.metrics[it::Metric::kinetic_change] = it::MetricValue::invalid("nonfinite_quantity");
    const it::IterationMetrics frame{info, {first, second}};
    const auto checked = it::decide(it::resolve(request, model, info.protocol, false),
        {"fd", "context", {first.id, second.id}}, {0, 0, true}, frame);
    const it::MetricLabels labels{{it::Metric::bsh_residual, "Maximum regularized per-function BSH residual"}};
    std::ostringstream printed;
    it::write_iteration_metrics(printed, "coupled roots", frame, checked, labels, false);
    const auto output = printed.str();
    require(output.find("intermediate protocol") != std::string::npos);
    require(output.find("All selected criteria met") == std::string::npos);
    require(output.find("no comparison") != std::string::npos);
    require(output.find("invalid") != std::string::npos);
    const auto first_value = output.find("1.0000e-05");
    const auto second_value = output.find("1.0000e-02");
    require(first_value != std::string::npos && second_value != std::string::npos);
    require(first_value - output.rfind('\n', first_value) ==
            second_value - output.rfind('\n', second_value));
  }
  Solver fresh;
  Solver::State initial;
  request.tolerances[it::Metric::bsh_residual] = {2.0, 0.0};
  decision = molresponse_v3::solvers::iterate_shared(
      fresh, initial, it::resolve(request, model, info.protocol, true),
      expected, info, 0, observe);
  require(fresh.updates == 0 && fresh.evaluations == 1 &&
          decision.action == it::Action::complete);
  Solver capped;
  Solver::State capped_state;
  request.tolerances[it::Metric::bsh_residual] = {0.0, 0.0};
  request.max_updates = 1;
  decision = molresponse_v3::solvers::iterate_shared(
      capped, capped_state, it::resolve(request, model, info.protocol, true),
      expected, info, 37, observe);
  require(capped.evaluations == 2 && capped.updates == 1 &&
          decision.action == it::Action::stop);
}
