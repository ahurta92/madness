// Who decides convergence (DESIGN B8) and the settings table (B9): the pure
// functions in ConvergenceSettings.h. No MPI, no files.

#include <madness/chem/ConvergenceSettings.h>

#include <iostream>
#include <string>
#include <vector>

namespace {
using namespace madness::convergence;

int failures = 0;
void check(bool ok, const std::string& what) {
  if (!ok) {
    ++failures;
    std::cerr << "FAILED: " << what << '\n';
  }
}

/// A moldft-like stage: H2O (3 atoms), dconv 1e-4, final thresh 1e-6.
ConvergenceSettings scf_like(Decider decider, bool user_targets) {
  ConvergenceSettings settings;
  settings.stage = "gs";
  settings.application = "moldft";
  settings.decider = decider;
  settings.decider_set_by_user = decider == Decider::requested_criteria;
  const SetBy dconv_by = user_targets ? SetBy::user : SetBy::engine_default;
  settings.metrics = {
      {Group::engine_check, "spin density change (max)", "dconv",
       "density_change", "", "a0^(-3/2)", 1e-4, dconv_by, Relation::less_than,
       5e-4, "max(dconv 1.0e-04, thresh 1.0e-06) x max(5, natom 3)", "", ""},
      {Group::engine_check, "BSH residual (max)", "dconv", "bsh_residual", "",
       "dimensionless", 1e-4, dconv_by, Relation::less_than, 5e-4,
       "5 x max(dconv, thresh)", "", ""},
      {Group::requested_criteria, "charge density change", "dconv",
       "charge_density_change", "", "a0^(-3/2)",
       user_targets ? std::optional<double>(1e-4) : std::nullopt,
       user_targets ? SetBy::user : SetBy::not_set,
       Relation::less_than_or_equal, std::nullopt, "", "",
       "compares spin-density change against 5.0e-04"},
      {Group::requested_criteria, "energy change", "econv", "energy_change", "",
       "hartree", user_targets ? std::optional<double>(1e-6) : std::nullopt,
       user_targets ? SetBy::user : SetBy::not_set,
       Relation::less_than_or_equal, std::nullopt, "",
       user_targets ? "" : "engine default econv 1.0e-05 is not a target (D8)",
       "has no energy term"}};
  settings.iteration_caps = {{30}, SetBy::engine_default};
  return settings;
}

} // namespace

int main() {
  // B8 case table (criteria-modes.md), both modes.
  struct Row {
    bool engine, criteria, engine_mode_converged, criteria_mode_converged;
  };
  for (const auto& row : std::vector<Row>{{true, true, true, true},
                                          {true, false, true, false},
                                          {false, true, false, true},
                                          {false, false, false, false}}) {
    check(converged(Decider::engine_check, row.engine, row.criteria)
              == row.engine_mode_converged,
          "engine mode follows the engine check only");
    check(converged(Decider::requested_criteria, row.engine, row.criteria)
              == row.criteria_mode_converged,
          "criteria mode follows the requested criteria only");
  }

  // A default deck: the engine check is enforced, nothing is requested.
  const auto plain = scf_like(Decider::engine_check, false);
  check(enforced(plain, plain.metrics[0]) && enforced(plain, plain.metrics[1]),
        "engine rows are enforced in engine mode");
  check(!enforced(plain, plain.metrics[2])
            && !enforced(plain, plain.metrics[3]),
        "rows without a target are reported only");
  check(!has_requested_criteria(plain), "a default deck requests nothing");
  check(requested_targets(plain).empty(), "no targets from defaults (B3, D8)");
  check(unenforced_target_warnings(plain).empty(),
        "no warning when the user wrote no target");

  // User targets in engine mode: reported, not enforced, with a warning each.
  auto engine = scf_like(Decider::engine_check, true);
  const auto warnings = unenforced_target_warnings(engine);
  check(warnings.size() == 2, "one warning per written target");
  check(!warnings.empty()
            && warnings.front()
                   == "dconv 1.0e-04 is reported, not enforced: the engine "
                      "check compares spin-density change against 5.0e-04. "
                      "Set convergence.use_default_criteria false to enforce "
                      "your targets.",
        "warning text names the keyword, its value and what the engine does");
  engine.warnings = warnings; // the stage's settings builder stores them
  const auto targets = requested_targets(engine);
  check(targets.size() == 2 && targets.at("energy_change").at("value") == 1e-6
            && targets.at("charge_density_change").at("units") == "a0^(-3/2)",
        "written targets become evaluate_explicit_criteria targets");

  // Criteria mode: only the requested rows are enforced; a missing BSH target
  // is warned about.
  auto criteria = scf_like(Decider::requested_criteria, true);
  check(!enforced(criteria, criteria.metrics[0])
            && enforced(criteria, criteria.metrics[2]),
        "criteria mode enforces requested rows only");
  const auto incomplete = incomplete_request_warnings(
      criteria, {"charge_density_change", "bsh_residual"});
  check(incomplete.size() == 1
            && incomplete.front().find("bsh_residual") != std::string::npos,
        "a request without a BSH target is warned about");
  check(incomplete_request_warnings(engine, {"bsh_residual"}).empty(),
        "no incomplete-request warning in engine mode");

  // Reuse (B6, B8).
  const nlohmann::json unmet_energy
      = {{"criteria",
          {{"charge_density_change", {{"status", "satisfied"}}},
           {"energy_change", {{"status", "unmet"}}}}}};
  const auto engine_reuse = decide_reuse(engine, true, unmet_energy);
  check(engine_reuse.action == ReuseAction::load
            && engine_reuse.unmet_not_enforceable
                   == std::vector<std::string>{"energy_change"},
        "engine mode loads and reports unenforceable unmet targets");
  check(decide_reuse(engine, false, unmet_energy).action == ReuseAction::resume,
        "engine mode resumes when the saved engine check fails");
  check(decide_reuse(criteria, true, unmet_energy).action
            == ReuseAction::resume,
        "criteria mode resumes when a requested target is unmet");
  const nlohmann::json all_met
      = {{"criteria", {{"energy_change", {{"status", "satisfied"}}}}}};
  check(
      decide_reuse(criteria, false, all_met).action == ReuseAction::load,
      "criteria mode loads an engine-mode generation that meets every target");

  // B9 one source: the JSON record and the printed table agree row by row.
  const auto record = to_json(engine);
  const auto table = render_settings_table(engine);
  check(record.at("decided_by") == "engine_check"
            && table.find("Decided by: engine check") != std::string::npos,
        "deciding mode in record and table");
  for (std::size_t i = 0; i < engine.metrics.size(); ++i) {
    const auto& row = record.at("metrics").at(i);
    const auto& metric = engine.metrics[i];
    check(row.at("enforced") == enforced(engine, metric)
              && row.at("set_by") == to_string(metric.set_by)
              && table.find(metric.name) != std::string::npos,
          "row " + metric.name + " agrees between record and table");
  }
  check(table.find("max(dconv 1.0e-04, thresh 1.0e-06) x max(5, natom 3)")
            != std::string::npos,
        "derived limits show their formula");
  check(render_settings_table(plain).find("engine default econv 1.0e-05")
            != std::string::npos,
        "an engine default is shown as one, not as a target");
  check(table.find("Iteration cap (maxiter): 30 at every rung (engine default)")
            != std::string::npos,
        "the iteration cap is shown with who set it");
  check(table.find("WARNING: dconv 1.0e-04 is reported") != std::string::npos,
        "warnings are printed with the table");

  // End of stage: verdicts and the limiting metric.
  const std::vector<MetricOutcome> outcomes
      = {{"density_change", 2e-4, "measured"},
         {"bsh_residual", 4.9e-4, "measured"},
         {"charge_density_change", 3e-4, "measured"},
         {"energy_change", std::nullopt, "unavailable"}};
  const auto outcome = outcome_to_json(engine, outcomes);
  check(outcome.at("metrics").at(0).at("verdict") == "met"
            && outcome.at("metrics").at(2).at("verdict") == "unmet"
            && outcome.at("metrics").at(3).at("verdict") == "unavailable",
        "verdicts: met, unmet (reported only), unavailable");
  check(outcome.at("limiting_metric") == "BSH residual (max)",
        "the enforced row closest to its limit is the limiting metric");
  check(render_outcome_table(engine, outcomes, "converged")
                .find("Limiting metric: BSH residual (max)")
            != std::string::npos,
        "the outcome table names the limiting metric");

  // D-095: a row judged per leg, each leg on its own limit; a row whose note
  // says it does not apply; the record's limit, units, phase and legs.
  {
    ConvergenceSettings legs = engine;
    MetricSetting per_leg;
    per_leg.group = Group::engine_check;
    per_leg.name = "BSH residual per leg";
    per_leg.measurement = "engine_bsh_residual";
    per_leg.units = "a.u.";
    per_leg.limit_formula = "5e-4 + rtol x ||x||";
    MetricSetting energy;
    energy.name = "energy change";
    energy.measurement = "energy_change";
    energy.note = "not applicable: see the reference stage";
    legs.metrics = {per_leg, energy};
    legs.scope = "limits of the one protocol step (thresh 1.0e-04, k 6)";
    MetricOutcome a{
        "engine_bsh_residual", 1e-4, "measured", 1e-3, "a.u.", "solver_step",
        "dipole_z f0.00000"};
    MetricOutcome b{
        "engine_bsh_residual", 2e-3, "measured", 1.5e-3, "a.u.", "solver_step",
        "dipole_z f0.02000"};
    check(row_verdict(per_leg, {a}) == "met"
              && row_verdict(per_leg, {a, b}) == "unmet"
              && verdict(per_leg,
                         MetricOutcome{"engine_bsh_residual", 1e-4, "measured"})
                     == "reported"
              && verdict(energy, MetricOutcome{"energy_change"})
                     == "not applicable",
          "verdicts (D-095): per leg on each leg's limit; not applicable; "
          "reported without a limit");
    const auto record = outcome_record(legs, {a, b});
    const auto& row = record.at("metrics").at(0);
    check(row.at("legs").size() == 2 && row.at("value") == 2e-3
              && row.at("verdict") == "unmet" && row.at("units") == "a.u."
              && row.at("phase") == "solver_step"
              && row.at("legs").at(1).at("verdict") == "unmet"
              && record.at("scope") == legs.scope
              && !to_json(legs).contains("scope"),
          "outcome_record (D-095): the legs, units, phase and scope; the "
          "settings record has no scope");
    const auto table = render_outcome_table(legs, {a, b}, "converged");
    check(table.find("Scope: limits of the one protocol step")
                  != std::string::npos
              && table.find("Engine check (values: solver step)")
                     != std::string::npos
              && table.find("dipole_z f0.02000") != std::string::npos
              && table.find("not applicable") != std::string::npos
              && table.find("Limiting metric: BSH residual per leg")
                     != std::string::npos,
          "outcome table (D-095): groups, phase, legs, scope, limiting leg");
  }

  // B5: rung_exit orders the exits once. Each row: decider, engine check,
  // criteria met, stall, iterations this run, cap, forced, expected exit.
  constexpr auto E = Decider::engine_check;
  constexpr auto C = Decider::requested_criteria;
  constexpr auto none = StopReason::in_progress;
  constexpr auto conv = StopReason::converged;
  constexpr auto stall = StopReason::stalled;
  constexpr auto cap = StopReason::iteration_limit;
  struct ExitCase {
    const char* name;
    RungSample sample;
    StopReason expected;
  };
  const std::vector<ExitCase> exits = {
      {"nothing yet", {E, false, false, false, 3, 30, 0}, none},
      {"engine check passes", {E, true, false, false, 3, 30, 0}, conv},
      {"engine mode ignores met criteria",
       {E, false, true, false, 3, 30, 0},
       none},
      {"criteria met", {C, false, true, false, 3, 30, 0}, conv},
      {"criteria mode ignores the engine check",
       {C, true, false, false, 3, 30, 0},
       none},
      {"cap reached", {E, false, false, false, 30, 30, 0}, cap},
      {"beyond the cap", {E, false, false, false, 31, 30, 0}, cap},
      {"stall", {C, false, false, true, 5, 30, 0}, stall},
      {"converged and cap: converged",
       {E, true, false, false, 30, 30, 0},
       conv},
      {"converged and stall: converged",
       {E, true, false, true, 5, 30, 0},
       conv},
      {"cap and stall: iteration limit",
       {C, false, false, true, 30, 30, 0},
       cap},
      {"all three: converged", {C, false, true, true, 30, 30, 0}, conv},
      {"forced: no exit before N even if converged",
       {E, true, false, false, 2, 30, 5},
       none},
      {"forced: no exit before N on a stall",
       {C, false, false, true, 2, 30, 5},
       none},
      {"forced: the cap does not end it early",
       {E, false, false, false, 3, 3, 5},
       none},
      {"forced: converged at N", {E, true, false, false, 5, 30, 5}, conv},
      {"forced: unmet at N is the iteration limit",
       {E, false, false, false, 5, 30, 5},
       cap},
      {"forced: stall ignored at N", {C, false, true, true, 5, 30, 5}, conv},
      {"forced: unmet and stalled at N",
       {C, false, false, true, 5, 30, 5},
       cap},
  };
  for (const auto& row : exits)
    check(rung_exit(row.sample) == row.expected,
          std::string("rung_exit: ") + row.name);
  // Over every verdict combination: never converged unless the decider's
  // verdict holds, and a stall never wins over the cap or convergence.
  for (const auto decider : {E, C})
    for (int bits = 0; bits < 16; ++bits) {
      const bool engine = (bits & 1) != 0;
      const bool criteria = (bits & 2) != 0;
      const bool stalled = (bits & 4) != 0;
      const bool at_cap = (bits & 8) != 0;
      const RungSample sample{decider,         engine, criteria, stalled,
                              at_cap ? 30 : 3, 30,     0};
      const auto exit = rung_exit(sample);
      const bool decided = converged(decider, engine, criteria);
      check((exit == conv) == decided,
            "rung_exit: converged exactly when the decider's verdict holds");
      check(exit != stall || (!decided && !at_cap && stalled),
            "rung_exit: stalled only below the cap and unconverged");
    }
  check(std::string(to_string(stall)) == "stalled"
            && stop_reason_from_string("iteration_limit") == cap
            && stop_reason_from_string("in_progress") == none
            && stop_reason_from_string("") == none,
        "stop reasons keep their saved spelling");

  std::cout << "convergence settings failures: " << failures << '\n';
  return failures ? 1 : 0;
}
