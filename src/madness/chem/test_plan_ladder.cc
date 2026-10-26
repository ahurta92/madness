// Table of cases for plan_ladder (DESIGN B6, the planner): what --continue does
// with one saved state. Row numbers follow planner-design.md section 3. Pure:
// no MPI, no files.

#include <madness/chem/LadderPlanner.h>

#include <iostream>
#include <optional>
#include <string>
#include <vector>

namespace {
using namespace madness::planning;
using madness::convergence::ConvergenceSettings;
using madness::convergence::Decider;
using madness::convergence::Group;
using madness::convergence::MetricSetting;
using madness::convergence::Relation;
using madness::convergence::SetBy;
using json = nlohmann::json;

int failures = 0;
void check(bool ok, const std::string& what) {
  if (!ok) {
    ++failures;
    std::cerr << "FAILED: " << what << '\n';
  }
}

const Rung coarse{1e-4, 6};
const Rung middle{1e-5, 7};
const Rung fine{1e-6, 8};
const std::vector<Rung> ladder = {coarse, fine};

json measured(double value, const char* units,
              std::optional<double> norm = std::nullopt) {
  json record
      = {{"value", value}, {"availability", "measured"}, {"units", units}};
  if (norm)
    record["normalization_scale"] = *norm;
  return record;
}

/// The deciding sample of an SCF stage: density and BSH changes.
json scf_sample(double density, double bsh) {
  return {{"density_change", measured(density, "a0^(-3/2)")},
          {"bsh_residual", measured(bsh, "dimensionless")}};
}

/// Measurements on the saved orbitals, with an energy change.
json saved_values(double energy_change) {
  return {{"energy_change", measured(energy_change, "hartree")}};
}

/// SCF's engine check at the final rung, H2 with dconv 1e-5: limits
/// 5e-6 (density) and 5e-6 (BSH), compared with <.
std::vector<EngineLimit> scf_limits(double dconv = 1e-5) {
  const double effective = std::max(dconv, fine.thresh);
  return {{"density_change", effective * 5, 0.0, false},
          {"bsh_residual", effective * 5, 0.0, false}};
}

ConvergenceSettings settings(Decider decider,
                             std::optional<double> energy_target) {
  ConvergenceSettings s;
  s.stage = "gs";
  s.application = "moldft";
  s.decider = decider;
  if (energy_target)
    s.metrics.push_back({Group::requested_criteria, "energy change", "econv",
                         "energy_change", "", "hartree", energy_target,
                         SetBy::user, Relation::less_than_or_equal,
                         std::nullopt, "", "", "has no energy term"});
  return s;
}

LadderRequest request(std::vector<Rung> rungs = ladder,
                      Decider decider = Decider::engine_check,
                      std::optional<double> energy_target = std::nullopt,
                      double dconv = 1e-5) {
  LadderRequest r;
  r.identity = "H2";
  r.ladder = std::move(rungs);
  r.settings = settings(decider, energy_target);
  r.engine_limits = scf_limits(dconv);
  return r;
}

/// A saved state of H2 at `rung` after `iterations`, stopped by `stop`, whose
/// deciding sample passes (true) or fails today's default engine check.
SavedState state(Rung rung, int iterations, StopReason stop, bool passes,
                 double energy_change = 1e-8) {
  SavedState s;
  s.identity = "H2";
  s.rung = rung;
  s.iterations_at_rung = iterations;
  s.stop = stop;
  s.engine_measurements
      = passes ? scf_sample(1e-6, 1e-6) : scf_sample(1e-3, 1e-3);
  s.saved_measurements = saved_values(energy_change);
  return s;
}

struct Case {
  std::string name;
  std::optional<SavedState> own;
  std::optional<SavedState> seed;
  LadderRequest request;
  Action action;
  std::size_t start_rung;
  int iterations_done;
  std::size_t unmet = 0;
  std::size_t unmet_not_enforceable = 0;
};
} // namespace

int main() {
  using S = StopReason;
  auto forced = request();
  forced.forced_iterations = true;
  auto other_k = state(Rung{1e-6, 9}, 4, S::converged, true);
  auto other_reference = state(fine, 5, S::converged, true);
  other_reference.reference_generation = "g-old";
  auto other_problem = state(fine, 5, S::converged, true);
  other_problem.identity = "H2 stretched";
  auto diverged = state(fine, 3, S::in_progress, false);
  diverged.diverged = true;
  auto coarse_done = state(coarse, 5, S::converged, true);

  const std::vector<Case> cases = {
      {"1 final rung converged, meets: load",
       state(fine, 5, S::converged, true),
       {},
       request(),
       Action::load,
       1,
       0},
      {"2 inside rung 2: resume it, 3 recorded",
       state(fine, 3, S::in_progress, false),
       {},
       request(),
       Action::resume,
       1,
       3},
      {"3 inside rung 1: resume it",
       state(coarse, 4, S::in_progress, false),
       {},
       request(),
       Action::resume,
       0,
       4},
      {"4 rung 1 converged: never redone, resume rung 2",
       coarse_done,
       {},
       request(),
       Action::resume,
       1,
       0},
      {"5 rung 1 ended at its cap: resume rung 2",
       state(coarse, 30, S::iteration_limit, false),
       {},
       request(),
       Action::resume,
       1,
       0},
      {"6 rung 1 stalled: resume rung 2",
       state(coarse, 9, S::stalled, false),
       {},
       request(),
       Action::resume,
       1,
       0},
      {"7 final rung ended at its cap, unmet: resume, fresh maxiter (D-007)",
       state(fine, 30, S::iteration_limit, false),
       {},
       request(),
       Action::resume,
       1,
       30},
      {"8 final rung ended at its cap, meets today: load",
       state(fine, 30, S::iteration_limit, true),
       {},
       request(),
       Action::load,
       1,
       0},
      {"9 killed inside the final rung: resume it",
       state(fine, 30, S::in_progress, true),
       {},
       request(),
       Action::resume,
       1,
       30},
      {"10 final rung converged, tighter written target: resume (criteria)",
       state(fine, 5, S::converged, true, 1e-6),
       {},
       request(ladder, Decider::requested_criteria, 1e-9),
       Action::resume,
       1,
       5},
      {"11 a new middle rung after a finished coarse rung",
       coarse_done,
       {},
       request({coarse, middle, fine}),
       Action::resume,
       1,
       0},
      {"12 coarser requested rungs are skipped",
       state(middle, 2, S::in_progress, false),
       {},
       request({coarse, middle, fine}),
       Action::resume,
       1,
       2},
      {"13 saved rung absent from the request: first tighter one",
       state(middle, 2, S::in_progress, false),
       {},
       request(),
       Action::resume,
       1,
       0},
      {"14 tighter than every requested rung, meets: load (R3)",
       state(Rung{1e-7, 9}, 4, S::converged, true),
       {},
       request(),
       Action::load,
       1,
       0},
      {"15 engine mode: engine check meets, a written target does not: load, "
       "reported not enforceable (B8)",
       state(fine, 5, S::converged, true, 1e-3),
       {},
       request(ladder, Decider::engine_check, 1e-6),
       Action::load,
       1,
       0,
       0,
       1},
      {"16 engine check fails today: resume",
       state(fine, 5, S::converged, false),
       {},
       request(),
       Action::resume,
       1,
       5},
      {"17 criteria mode, a written target unmet: resume",
       state(fine, 5, S::converged, true, 1e-3),
       {},
       request(ladder, Decider::requested_criteria, 1e-6),
       Action::resume,
       1,
       5},
      {"18 an engine-mode state meeting every target, criteria mode: load",
       state(fine, 5, S::converged, true, 1e-8),
       {},
       request(ladder, Decider::requested_criteria, 1e-6),
       Action::load,
       1,
       0},
      {"19 no saved state, no seed: fresh",
       std::nullopt,
       {},
       request(),
       Action::fresh,
       0,
       0},
      {"20 diverged at rung 2, rung 1 converged: seed from it, restart rung 2 "
       "(R6)",
       diverged, coarse_done, request(), Action::seed, 1, 0},
      {"21 diverged, no coarser rung: fresh (R6)",
       diverged,
       {},
       request(),
       Action::fresh,
       1,
       0},
      {"22 tighter dconv the saved values still meet: load (re-judged, B6)",
       state(fine, 5, S::converged, true),
       {},
       request(ladder, Decider::engine_check, std::nullopt, 1e-8),
       Action::load,
       1,
       0},
      {"23 forced iterations always iterate",
       state(fine, 5, S::converged, true),
       {},
       forced,
       Action::resume,
       1,
       5},
      {"24 same thresh, different k: not the same rung",
       other_k,
       {},
       request(),
       Action::resume,
       1,
       0},
      {"25 one-rung ladder, converged, meets: load",
       state(fine, 5, S::converged, true),
       {},
       request({fine}),
       Action::load,
       0,
       0},
      {"26 tighter than every requested rung, does not meet: resume final (R3)",
       state(Rung{1e-7, 9}, 4, S::converged, false),
       {},
       request(),
       Action::resume,
       1,
       0},
      {"27 final rung stalled, unmet: load and report (R2)",
       state(fine, 9, S::stalled, false),
       {},
       request(),
       Action::load,
       1,
       0,
       1},
      {"28 stalled final rung, a finer rung added: resume there (R2)",
       state(fine, 9, S::stalled, false),
       {},
       request({coarse, fine, Rung{1e-7, 9}}),
       Action::resume,
       2,
       0},
      {"29 stalled final rung, loosened target now met: load (R2, re-judged)",
       state(fine, 9, S::stalled, true, 1e-8),
       {},
       request(ladder, Decider::requested_criteria, 1e-6),
       Action::load,
       1,
       0},
      {"30 a different reference generation: fresh (R1)",
       other_reference,
       {},
       request(),
       Action::fresh,
       0,
       0},
      {"31 a different problem: fresh (R1)",
       other_problem,
       {},
       request(),
       Action::fresh,
       0,
       0},
      {"32 a different problem, the deck says seed: seed, rung 1 (R1)",
       other_problem, other_problem, request(), Action::seed, 0, 0},
  };
  for (const auto& c : cases) {
    const auto plan = plan_ladder(c.own, c.seed, c.request);
    const bool ok
        = plan.action == c.action && plan.start_rung == c.start_rung
          && plan.iterations_done == c.iterations_done
          && plan.unmet.size() == c.unmet
          && plan.unmet_not_enforceable.size() == c.unmet_not_enforceable
          && !plan.why.empty();
    check(ok, c.name + " -> " + to_string(plan.action) + " rung "
                  + std::to_string(plan.start_rung) + " done "
                  + std::to_string(plan.iterations_done) + " why: " + plan.why);
  }

  // Row 22's companion: a tighter dconv lowers today's limit below the saved
  // values, so the state resumes. (Above, the saved values already meet it:
  // re-judged, a load is right, B6.)
  auto weak = state(fine, 5, S::converged, true);
  weak.engine_measurements = scf_sample(2e-5, 1e-6);
  check(plan_ladder(weak, {}, request()).action == Action::load,
        "22a values meet the dconv 1e-5 limit (5e-5): load");
  {
    auto tight = request(ladder, Decider::engine_check, std::nullopt, 1e-6);
    check(
        plan_ladder(weak, {}, tight).action == Action::resume,
        "22b dconv 1e-6 lowers the limit to 5e-6 below the saved 2e-5: resume");
  }

  // Row 37 (D-005.4): a relative limit with no recorded norm cannot be
  // re-judged: resume, never load.
  {
    auto response = request();
    response.engine_limits = {{"bsh_residual", 5e-4, 5e-4, true}};
    auto with_norm = state(fine, 5, S::converged, true);
    with_norm.engine_measurements
        = {{"bsh_residual", measured(1e-4, "x", 2.0)}};
    check(plan_ladder(with_norm, {}, response).action == Action::load,
          "37a relative limit with its norm: judged (1e-4 <= 5e-4 + 5e-4 x 2): "
          "load");
    auto without_norm = with_norm;
    without_norm.engine_measurements = {{"bsh_residual", measured(1e-4, "x")}};
    const auto plan = plan_ladder(without_norm, {}, response);
    check(plan.action == Action::resume
              && plan.why.find("cannot be re-judged") != std::string::npos,
          "37b no norm recorded: resume, never load (D-005.4)");
  }

  // R4, R7: the nearest saved frequency of the same perturbation and
  // reference; the lower on a tie; diverged states never seed.
  const auto response_state
      = [](const char* perturbation, double frequency, const char* reference) {
          SavedState s = state(fine, 5, S::converged, true);
          s.perturbation = perturbation;
          s.frequency = frequency;
          s.reference_generation = reference;
          return s;
        };
  {
    const std::vector<SavedState> saved
        = {response_state("dipole_z", 0.02, "g-1"),
           response_state("dipole_z", 0.1, "g-1"),
           response_state("dipole_x", 0.05, "g-1"),
           response_state("dipole_z", 0.05, "g-2")};
    const auto seed = nearest_seed(saved, "dipole_z", "g-1", 0.05);
    check(seed && seed->frequency == 0.02,
          "33 seed from the nearest frequency of the same perturbation and "
          "reference (R4)");
    const std::vector<SavedState> tie
        = {response_state("dipole_z", 0.07, "g-1"),
           response_state("dipole_z", 0.03, "g-1")};
    const auto lower = nearest_seed(tie, "dipole_z", "g-1", 0.05);
    check(lower && lower->frequency == 0.03,
          "34 a tie: the lower frequency (R7)");
    check(!nearest_seed(saved, "dipole_y", "g-1", 0.05),
          "35 no state of the same perturbation and reference: none (R4)");
    auto bad = response_state("dipole_z", 0.05, "g-1");
    bad.diverged = true;
    check(nearest_seed({bad}, "dipole_z", "g-1", 0.05) == std::nullopt,
          "a diverged state never seeds");
  }

  // R6: the finest converged rung coarser than the diverged one.
  {
    const std::vector<SavedState> history
        = {state(coarse, 5, S::converged, true),
           state(middle, 4, S::converged, true),
           state(fine, 3, S::in_progress, false)};
    const auto seed = coarser_seed(history, fine);
    check(seed && seed->rung.thresh == middle.thresh,
          "R6: the finest converged coarser rung");
    check(!coarser_seed(history, coarse), "R6: nothing coarser than rung 1");
  }

  // maxiter in one place (B5).
  check(maxiter_at_rung({}, 30, 0, 1) == 30, "maxiter: one number");
  check(maxiter_at_rung({30, 10}, 25, 0, 1) == 10, "maxiter: one per rung");
  check(maxiter_at_rung({12}, 25, 0, 1) == 12,
        "maxiter: one entry for every rung");
  check(maxiter_at_rung({30, 10}, 25, 5, 1) == 5,
        "maxiter: forced iterations win");

  // The record line.
  const auto resumed
      = plan_ladder(state(fine, 3, S::in_progress, false), {}, request());
  check(
      describe(resumed, ladder)
          == "resumed at rung 2 (thresh 1e-06, k 8), iteration 4 of that rung",
      "describe: " + describe(resumed, ladder));

  std::cout << "plan ladder failures: " << failures << '\n';
  return failures ? 1 : 0;
}
