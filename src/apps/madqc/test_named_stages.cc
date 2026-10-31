#include "MoldftCheckpoint.hpp"
#include "NamedStages.hpp"
#include "ResponseCheckpoint.hpp"
#include "RunNamedStages.hpp"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <unistd.h>

int main(int argc, char** argv) {
    auto& world = madness::initialize(argc, argv);
    int failures = 0;
    auto check = [&](bool ok, const std::string& why) {
        if (!ok) { ++failures; if (world.rank() == 0) std::cerr << why << '\n'; }
    };
    const std::string geometry = "molecule H2\n units atomic\n H 0 0 -0.7\n H 0 0 0.7\nend\n";
    const std::string gs = "moldft GS\n molecule H2\n xc hf\n protocol [1e-4,1e-6]\n dconv 1e-5\n requested_properties [energy,dipole]\nend\n";
    const std::string response = "molresponse optical\n reference GS\n protocol [1e-4]\n dconv 2e-4\n requested_properties [polarizability]\n polarizability.frequencies [0.0,0.02]\n polarizability.frequency_units au\n polarizability.components [zz]\nend\n";
    try {
        auto stages = named_stages::parse(world, geometry + gs + response +
            "moldft other\n molecule H2\n xc hf\n dconv 3e-5\n requested_properties [energy]\nend\n");
        check(stages.size() == 3, "three independently named stages");
        check(stages[1].reference == "GS", "case-sensitive reference");
        check(stages[0].params.get<CalculationParameters>().dconv() == 1e-5, "ground owns dconv");
        check(stages[1].params.get<ResponseParameters>().dconv() == 2e-4, "response owns dconv");
        check(stages[1].params.get<CalculationParameters>().protocol() == std::vector<double>{1e-4}, "response owns protocol");
        stages[2].params.get<CalculationParameters>().set_user_defined_value("dconv", 9e-5);
        check(stages[0].params.get<CalculationParameters>().dconv() == 1e-5, "instances do not alias");
        check(stages[0].params.get<Molecule>().natom() == 2, "named molecule bound");
        CalculationParameters paths;
        paths.set_user_defined_value("prefix", std::string("/tmp/CaseRun/GS/state"));
        check(paths.prefix() == "/tmp/CaseRun/GS/state", "programmatic prefix preserves case");
        named_stages::read_value(world, paths, {1, "prefix /tmp/CaseRun/GS/state"});
        check(paths.prefix() == "/tmp/CaseRun/GS/state", "text prefix preserves case");
    } catch (const std::exception& e) { check(false, e.what()); }
    auto rejects_with
        = [&](const std::function<void()>& call, const std::string& message,
              const std::string& why) {
            try {
              call();
              check(false, "accepted: " + why);
            } catch (const std::exception& e) {
              check(std::string(e.what()).find(message) != std::string::npos,
                    why + ": " + e.what());
            }
          };
    auto rejects = [&](const std::string& text, const std::string& message) {
        try { named_stages::parse(world, text); check(false, "accepted: " + message); }
        catch (const std::exception& e) { check(std::string(e.what()).find(message) != std::string::npos, e.what()); }
    };
    rejects(geometry + gs + gs, "duplicate");
    rejects(geometry + "nemo foo\nend\n", "unsupported stage");
    rejects(geometry + response + gs, "reference");
    rejects(geometry + gs + "molresponse bad\n reference GS\n maxsub 4\nend\n", "maxsub");
    rejects(geometry + gs + "molresponse bad\nend\n", "reference");
    rejects(geometry + gs + "molresponse bad\n reference H2\nend\n", "reference");
    // DESIGN A1: a stage outside the deck is read from its published generation.
    try {
        auto alone = response;
        alone.replace(alone.find("GS"), 2, "outside");
        const auto stages = named_stages::parse(world, alone);
        check(stages.size() == 1 && stages[0].reference == "outside", "reference outside the deck");
        check(stages[0].params.get<Molecule>().natom() == 0, "nothing copied from a reference");
    } catch (const std::exception& e) { check(false, e.what()); }
    rejects(geometry + "moldft bad\n molecule absent\nend\n", "molecule");
    rejects(geometry + "moldft ../bad\nend\n", "name");
    rejects(geometry + "moldft bad\n molecule H2\n dconv 2e-4junk\nend\n", "line");
    rejects(geometry + "moldft bad\n molecule H2\n mystery 1\nend\n", "mystery");
    rejects(geometry + "moldft bad\n molecule H2\n dconv 1e-4\n dconv 2e-4\nend\n", "duplicate");
    rejects(geometry + gs + "molresponse bad\n reference GS\n protocol [1e-8]\n requested_properties [polarizability]\n polarizability.frequencies [0]\n polarizability.frequency_units au\n polarizability.components [zz]\nend\n", "tighter");
    rejects("molecule broken\n H 0 0\nend\n" + gs, "coordinates");
    rejects(geometry + "moldft unfinished\n molecule H2\n", "unterminated");
    auto collisions = response;
    collisions.replace(collisions.find("[0.0,0.02]"), std::string("[0.0,0.02]").size(), "[0.02,0.020001]");
    rejects(geometry + gs + collisions, "frequency-key precision");
    collisions = response;
    collisions.replace(collisions.find("[1e-4]"), std::string("[1e-4]").size(), "[1.1e-4,1e-4]");
    rejects(geometry + gs + collisions, "protocol-key precision");
    // D8, B8, B9: written values are targets; the settings table shows the
    // engine defaults and whether each target is enforced.
    {
      using nlohmann::json;
      namespace cv = madness::convergence;
      const auto saved = [](json value, const char* availability) {
        return json{{"energy_change",
                     {{"value", value},
                      {"availability", availability},
                      {"units", "hartree"}}}};
      };
      CalculationParameters defaults;
      check(!named_stages::requested_econv(defaults),
            "engine default econv is not a request");
      const auto plain
          = named_stages::moldft_convergence_settings("gs", defaults, 3);
      check(plain.decider == cv::Decider::engine_check
                && cv::requested_targets(plain).empty(),
            "a default deck: engine check decides, nothing requested");
      check(plain.metrics.at(0).set_by == cv::SetBy::engine_default
                && plain.metrics.at(0).effective_limit
                       == 5.0 * std::max(1e-4, defaults.dconv()),
            "engine density row: engine default dconv and its derived limit");
      check(plain.metrics.at(4).note.find("engine default econv")
                != std::string::npos,
            "the engine default econv is shown as one, not as a target");
      CalculationParameters requested;
      requested.set_user_defined_value("econv", 1e-6);
      check(named_stages::requested_econv(requested) == 1e-6,
            "user econv is the request");
      const auto engine_mode
          = named_stages::moldft_convergence_settings("gs", requested, 3);
      check(engine_mode.warnings.size() == 1
                && engine_mode.warnings.front().find("econv 1.0e-06")
                       != std::string::npos,
            "engine mode warns about the unenforced econv");
      // The report annotates the shared verdict's evaluation (S4d2).
      const auto judged = [](const json& measurements,
                             const cv::ConvergenceSettings& settings) {
        return named_stages::requested_criteria_report(
            madness::evaluate_explicit_criteria(
                measurements, cv::requested_targets(settings)),
            settings);
      };
      const auto met = judged(saved(5e-7, "measured"), engine_mode);
      check(met.at("status") == "satisfied" && met.at("enforced") == false,
            "met by evidence, reported in engine mode");
      const auto unmet = judged(saved(2e-6, "measured"), engine_mode);
      check(unmet.at("status") == "unmet"
                && unmet.at("criteria").at("energy_change").at("outcome")
                       == "unmet, not enforceable in engine mode",
            "unmet evidence is reported, not enforced, in engine mode");
      check(judged(json::object(), engine_mode).at("status") == "unavailable",
            "missing evidence does not meet a request");
      requested.set_user_defined_value(
          std::string(cv::use_default_criteria_key), false);
      const auto criteria_mode
          = named_stages::moldft_convergence_settings("gs", requested, 3);
      check(criteria_mode.decider == cv::Decider::requested_criteria
                && judged(saved(2e-6, "measured"), criteria_mode).at("enforced")
                       == true,
            "criteria mode enforces the written targets");
      check(criteria_mode.warnings.size() == 2,
            "criteria mode warns about missing density and BSH targets");
    }
    auto criteria_without_targets = gs;
    criteria_without_targets.replace(
        criteria_without_targets.find(" dconv 1e-5\n"),
        std::string(" dconv 1e-5\n").size(),
        " convergence.use_default_criteria false\n");
    rejects(geometry + criteria_without_targets, "needs at least one written");
    auto progress_engine = gs;
    progress_engine.replace(progress_engine.find(" dconv 1e-5\n"), 0,
                            " accuracy.kinetic_trace_change 1e-6\n"
                            " progress.window 2\n"
                            " progress.minimum_relative_improvement 0.1\n");
    rejects(geometry + progress_engine, "progress.* stall rules");
    auto old_energy = gs;
    old_energy.replace(old_energy.find(" dconv 1e-5\n"), 0,
                       " accuracy.energy_change 1e-6\n");
    rejects(geometry + old_energy, "is now econv");
    // B5: maxiter takes one number or one per rung (moldft only).
    try {
      auto capped = gs;
      capped.replace(capped.find(" dconv 1e-5\n"), 0, " maxiter [30, 10]\n");
      const auto stages = named_stages::parse(world, geometry + capped);
      check(stages[0].params.get<CalculationParameters>().get<std::vector<int>>(
                "maxiter_per_rung")
                == std::vector<int>{30, 10},
            "per-rung maxiter list is kept per rung");
    } catch (const std::exception& e) {
      check(false, e.what());
    }
    auto wrong_length = gs;
    wrong_length.replace(wrong_length.find(" dconv 1e-5\n"), 0,
                         " maxiter [30, 10, 5]\n");
    rejects(geometry + wrong_length, "one per protocol rung");
    auto response_list = response;
    response_list.replace(response_list.find(" dconv 2e-4\n"), 0,
                          " maxiter [30]\n");
    rejects(geometry + gs + response_list, "moldft only");
    rejects(geometry + "moldft bad\n molecule H2\n restart.start_rung 1\nend\n",
            "restart.start_rung");
    // A3: checkpoint cadence and retention are stage keys (moldft only).
    try {
        auto cadence = gs;
        cadence.replace(
            cadence.find(" dconv 1e-5\n"), 0,
            " checkpoint.every_iterations 3\n checkpoint.keep_all true\n");
        const auto stages = named_stages::parse(world, geometry + cadence);
        check(stages[0].checkpoint_every_iterations == 3
                  && stages[0].checkpoint_keep_all,
              "checkpoint keys are kept on the moldft stage");
        const auto defaults = named_stages::parse(world, geometry + gs);
        check(defaults[0].checkpoint_every_iterations == 0
                  && !defaults[0].checkpoint_keep_all,
              "default: decision points only, retention on");
    } catch (const std::exception& e) {
        check(false, e.what());
    }
    auto zero_cadence = gs;
    zero_cadence.replace(zero_cadence.find(" dconv 1e-5\n"), 0,
                         " checkpoint.every_iterations 0\n");
    rejects(geometry + zero_cadence, "must be at least 1");
    {
      auto response_cadence = response;
      response_cadence.replace(response_cadence.find(" dconv 2e-4\n"), 0,
                               " checkpoint.every_iterations 2\n");
      rejects(geometry + gs + response_cadence, "publishes at solve ends only");
    }
    // D-070: checkpoint.keep_all is every application's.
    try {
      auto keep = response;
      keep.replace(keep.find(" dconv 2e-4\n"), 0,
                   " checkpoint.keep_all true\n");
      const auto stages = named_stages::parse(world, geometry + gs + keep);
      check(stages[1].checkpoint_keep_all,
            "checkpoint.keep_all is accepted for molresponse (D-070)");
    } catch (const std::exception& e) {
      check(false, std::string("molresponse keep_all: ") + e.what());
    }
    // The states a response request's values are built from (D-019): each leg
    // and source once, nuclear legs at 0, "every atom" from the reference.
    {
        const json request = {{"protocol", {1e-4}}, {"properties", {
            {"polarizability", {{"frequencies", {0.0}}, {"components", {"zz"}}}},
            {"hyperpolarizability", {{"frequencies", {0.0, 0.05}}, {"components", {"zzz", "xzz"}}}},
            {"raman", {{"frequencies", {0.05}}, {"components", {"xz"}}, {"atoms", json::array()},
                       {"directions", "z"}, {"full_set", false}}}}}};
        const auto needs = named_stages::response_needs(request, 2);
        std::vector<std::string> legs;
        for (const auto& [p, f] : needs.legs) legs.push_back(molresponse_v3::fd_node_id(p, f));
        const std::vector<std::string> expected_legs = {
            "fd:dipole_z@f0.00000", "fd:dipole_x@f0.00000", "fd:dipole_z@f0.10000",
            "fd:dipole_z@f0.05000", "fd:dipole_x@f0.10000", "fd:dipole_x@f0.05000",
            "fd:nuc_0_z@f0.00000", "fd:nuc_1_z@f0.00000"};
        check(legs == expected_legs, "response_needs: the legs, each once, in request order");
        const std::vector<std::string> expected_sources = {
            "vbc:dipole_z__dipole_z@f0.00000_f0.00000", "vbc:dipole_z__dipole_z@f0.05000_f0.05000",
            "vbc:dipole_z__nuc_0_z@f0.05000_f0.00000", "vbc:dipole_z__nuc_1_z@f0.05000_f0.00000"};
        check(needs.sources == expected_sources, "response_needs: the B x C sources");
        const json legacy = {{"protocol", {1e-4}}, {"frequencies", {0.0, 0.02}}, {"components", {"zz"}}};
        check(named_stages::response_needs(legacy, 2).legs.size() == 2 &&
                  named_stages::response_needs(legacy, 2).sources.empty(),
              "response_needs: a request saved before it named properties is polarizability");
    }
    // S3a (D-036): the AppSpec rows, one table per field.
    try {
      namespace cv = madness::convergence;
      using madness::checkpoint::Manifest;
      // The rows without their run binding: this test links no engine.
      const auto moldft = named_stages::moldft_fields();
      const auto molresponse = named_stages::response_fields();
      auto cadence = gs;
      cadence.replace(cadence.find(" dconv 1e-5\n"), 0,
                      " checkpoint.every_iterations 3\n");
      const auto stages = named_stages::parse(
          world, geometry + cadence + response
                     + "moldft other\n molecule H2\n xc hf\n "
                       "requested_properties [energy]\nend\n");
      const auto& ground = stages[0];
      const auto& optical = stages[1];
      const auto& other = stages[2];
      // name, validate (empty until S6, D-047); run is bound by <app>_spec().
      check(moldft.name == "moldft" && molresponse.name == "molresponse",
            "spec: names are the deck's block kinds");
      check(!moldft.validate && !molresponse.validate,
            "spec: validate is left empty until S6");
      check(!moldft.run && !molresponse.run,
            "spec: the fields leave the run binding to <app>_spec()");
      // settings: the stage's table plus moldft's cadence, who set each.
      const auto settings = moldft.settings(ground, 2);
      const auto plain = named_stages::moldft_convergence_settings(
          "GS", ground.params.get<CalculationParameters>(), 2);
      check(
          cv::to_json(settings).at("metrics")
                  == cv::to_json(plain).at("metrics")
              && settings.decider == plain.decider,
          "spec settings (moldft): the metrics of moldft_convergence_settings");
      check(settings.checkpoints.every_iterations == 3
                && settings.checkpoints.every_iterations_set_by
                       == cv::SetBy::user
                && !settings.checkpoints.keep_all
                && settings.checkpoints.keep_all_set_by
                       == cv::SetBy::engine_default,
            "spec settings (moldft): written cadence by the user, retention by "
            "default");
      check(cv::to_json(molresponse.settings(optical, 2))
                == cv::to_json(
                    named_stages::response_convergence_settings(optical, 1e-4)),
            "spec settings (response): response_convergence_settings at the "
            "final rung");
      // members: the files a generation must hold.
      check(moldft.members(json::object())
                == std::vector<std::string>{"state.restartdata"},
            "spec members (moldft): one orbital archive");
      const auto key
          = molresponse_v3::protocol_key(1e-4, named_stages::protocol_k(1e-4));
      const json bundle = {{"type", "tda"},
                           {"shell", "closed_shell"},
                           {"backend", "native"},
                           {"bundle_dir", "es__" + key},
                           {"n_roots", 1},
                           {"converged", false},
                           {"roots",
                            {{{"slot", 0},
                              {"omega", 0.5},
                              {"density_residual", nullptr},
                              {"omega_residual", nullptr}}}}};
      const json metadata = {
          {"fd_states", json::object()},
          {"vbc_states", json::object()},
          {"protocols",
           {{key, {{"thresh", 1e-4}, {"k", named_stages::protocol_k(1e-4)}}}}},
          {"excited_states", {{key, bundle}}}};
      const std::string root
          = "es__" + key + "/" + molresponse_v3::detail_save_load::root_file(0);
      check(molresponse.members(metadata) == std::vector<std::string>{root},
            "spec members (response): one archive per root of a bundle");
      auto diverged = metadata;
      diverged["excited_states"][key]["diverged"] = true;
      rejects_with([&] { (void)molresponse.members(diverged); }, "diverged",
                   "spec members (response): a diverged bundle is refused");
      // validate_generation: the checks on the generation --continue reads.
      Manifest scf;
      scf.stage = "GS";
      scf.method = "HF";
      scf.numerical = {{"representation", "mo"}};
      scf.archives.push_back({"orbitals", "state.restartdata", {}, 1, {}});
      moldft.validate_generation(ground, 2, scf);
      check(true, "spec validate_generation (moldft): no report, no controls");
      auto refused = scf;
      refused.stage = "other";
      rejects_with([&] { moldft.validate_generation(ground, 2, refused); },
                   "stage/method mismatch",
                   "spec validate_generation (moldft): another stage's "
                   "generation");
      refused = scf;
      refused.method = "TDHF";
      rejects_with([&] { moldft.validate_generation(ground, 2, refused); },
                   "stage/method mismatch",
                   "spec validate_generation (moldft): another method");
      refused = scf;
      refused.archives.clear();
      rejects_with([&] { moldft.validate_generation(ground, 2, refused); },
                   "inventory",
                   "spec validate_generation (moldft): a missing orbital "
                   "archive");
      refused = scf;
      refused.numerical["representation"] = "ao";
      rejects_with([&] { moldft.validate_generation(ground, 2, refused); },
                   "unsupported representation",
                   "spec validate_generation (moldft): another representation");
      refused = scf;
      refused.metrics
          = {{"convergence_report",
              {{"measurements",
                {{"energy_change", {{"availability", "invalid"}}}}}}}};
      rejects_with([&] { moldft.validate_generation(ground, 2, refused); },
                   "energy_change",
                   "spec validate_generation (moldft): an invalid measurement "
                   "is refused");
      Manifest es;
      es.stage = "optical";
      es.method = "TDHF";
      es.physical
          = {{"schema", "madqc.response.v1"}, {"shell", "closed_shell"}};
      es.numerical = {
          {"request", {{"protocol", {1e-4}}, {"properties", json::object()}}},
          {"engine_metadata", metadata},
          {"bundle_indices",
           {{key,
             {{"roots",
               {{{"file", molresponse_v3::detail_save_load::root_file(0)}}}},
              {"type", "tda"},
              {"shell", "closed_shell"},
              {"n_roots", 1},
              {"backend", "native"},
              {"writer_nproc", 1}}}}}};
      es.archives.push_back({"response", root, {}, 1, {}});
      es.stop_reason = "in_progress";
      molresponse.validate_generation(optical, 2, es);
      check(
          true,
          "spec validate_generation (response): a complete bundle generation");
      auto renamed = es;
      renamed.stage = "elsewhere";
      rejects_with(
          [&] { molresponse.validate_generation(optical, 2, renamed); },
          "stage or response model mismatch",
          "spec validate_generation (response): another stage's generation");
      auto missing = es;
      missing.archives.clear();
      rejects_with(
          [&] { molresponse.validate_generation(optical, 2, missing); },
          "inventory",
          "spec validate_generation (response): a missing archive");
      // properties: the requested records from the engine's output.
      const json evaluated
          = {{"properties",
              {{{"name", "energy"}, {"value", -1.1}, {"units", "hartree"}},
               {{"name", "dipole"},
                {"value", {0.0, 0.0, 0.1}},
                {"units", "e*a0"}}}}};
      check(moldft.properties(other, 2, evaluated)
                == json{{"energy", {{"value", -1.1}, {"units", "hartree"}}}},
            "spec properties (moldft): only the requested ones, without their "
            "names");
      const json alpha_rows = {{{"omega", 0.0},
                                {"directions", "z"},
                                {"alpha", {{6.45}}},
                                {"converged", true}},
                               {{"omega", 0.02},
                                {"directions", "z"},
                                {"alpha", {{6.46}}},
                                {"converged", true}}};
      const json output = {{"properties", {{"alpha", {{key, alpha_rows}}}}},
                           {"metadata", metadata},
                           {"vibrational_analysis", nullptr},
                           {"raman_spectra", nullptr},
                           {"transition", nullptr}};
      const auto direct = [&]() -> json {
        try {
          return named_stages::select_response_properties(
              optical, output.at("properties"), output.at("metadata"), nullptr,
              nullptr, 2, nullptr);
        } catch (const std::exception& e) {
          return e.what();
        }
      }();
      const auto field = [&]() -> json {
        try {
          return molresponse.properties(optical, 2, output);
        } catch (const std::exception& e) {
          return e.what();
        }
      }();
      check(field == direct
                && field.at("polarizability")
                           .at("values")
                           .at(1)
                           .at("components")
                           .at("zz")
                       == 6.46,
            "spec properties (response): select_response_properties on the "
            "output");
      // S3b (D-048): the response stage's read side, as functions of their
      // arguments. The request: what the stage's states are planned against.
      const auto request = named_stages::response_request(optical);
      check(request.at("protocol") == json{1e-4} && request.at("dconv") == 2e-4
                && request.at("decided_by") == "engine_check"
                && request.at("properties").at("polarizability")
                       == json{{"frequencies", {0.0, 0.02}},
                               {"components", {"zz"}}},
            "response_request: the stage's protocol, dconv, decider and "
            "properties");
      const auto legs = named_stages::response_needs(request, 2).legs;
      check(legs.size() == 2 && legs[0].second == 0.0 && legs[1].second == 0.02,
            "response_request: zz at two frequencies needs two z legs");
      // response_state_records: one record per state, with its iterations
      // at the rung across runs (D-007).
      const auto channel = legs[0].first.description();
      const auto fkey = molresponse_v3::ResponseMetadata::freq_key(0.0);
      const auto state = channel + "|" + key + "|" + fkey;
      json fd = metadata;
      fd.erase("excited_states");
      fd["fd_states"]
          = {{channel, {{key, {{fkey, {{"iter", 3}, {"freq", 0.0}}}}}}}};
      const auto iterations = [&](const named_stages::IterationTotals& totals) {
        return named_stages::response_state_records(fd, "g-ref", totals)
            .at(0)
            .at("iterations_at_rung")
            .get<int>();
      };
      const auto records
          = named_stages::response_state_records(fd, "g-ref", {});
      check(records.size() == 1 && records[0].at("key") == state
                && records[0].at("identity") == channel + "|0"
                && records[0].at("reference_generation") == "g-ref"
                && records[0].at("rung")
                       == json{{"thresh", 1e-4},
                               {"k", named_stages::protocol_k(1e-4)}}
                && records[0].at("stop") == "in_progress"
                && records[0].at("diverged") == false,
            "response_state_records: one FD state, unfinished");
      struct TotalsRow {
        named_stages::IterationTotals totals;
        int expected;
        const char* why;
      };
      for (const auto& row : std::vector<TotalsRow>{
               {{}, 3, "first run: this attempt's count"},
               {{{{state, 7}}, {}, {}}, 7, "untouched: the recorded total"},
               {{{{state, 7}}, {state}, {state}},
                10,
                "resumed and solved: the recorded total plus this attempt"},
               {{{{state, 7}}, {}, {state}},
                3,
                "solved again from the start: this attempt alone"}})
        check(iterations(row.totals) == row.expected,
              std::string("response_state_records: ") + row.why);
      const auto bundles
          = named_stages::response_state_records(metadata, "g-ref", {});
      check(bundles.size() == 1 && bundles[0].at("kind") == "excited_states"
                && bundles[0].at("identity")
                       == molresponse_v3::es_node_id(true, 1)
                && bundles[0].at("key")
                       == molresponse_v3::es_node_id(true, 1) + "|" + key
                && bundles[0].at("roots").size() == 1,
            "response_state_records: a bundle is one state per rung");
      // plan_response_states: one plan per leg the request needs.
      using Kind = molresponse_v3::StatePlan::Kind;
      const auto node0 = molresponse_v3::fd_node_id(legs[0].first, 0.0);
      const auto node1 = molresponse_v3::fd_node_id(legs[1].first, 0.02);
      const auto fresh = named_stages::response_engine_plans(
          named_stages::plan_response_states(optical, request, 2, "g-ref",
                                             json::array()));
      check(fresh.engine.size() == 2
                && fresh.engine.at(node0).kind == Kind::fresh
                && fresh.engine.at(node1).kind == Kind::fresh
                && fresh.records.size() == 2
                && fresh.records[0].at("action") == "fresh"
                && fresh.records[1].at("action") == "fresh"
                && fresh.resumed.empty(),
            "plan_response_states: no records, every leg fresh");
      const auto resumed = named_stages::response_engine_plans(
          named_stages::plan_response_states(optical, request, 2, "g-ref",
                                             records));
      check(resumed.engine.at(node0).kind == Kind::resume
                && resumed.records[0].at("action") == "resumed"
                && resumed.records[0].at("iterations_done_at_rung") == 3
                && resumed.resumed == std::set<std::string>{state},
            "plan_response_states: an unfinished rung resumes, and its key "
            "is resumed");
      check(resumed.engine.at(node1).kind == Kind::seed
                && resumed.engine.at(node1).seed_frequency == 0.0
                && resumed.records[1].at("action") == "seeded"
                && resumed.records[1].at("seed_frequency") == 0.0,
            "plan_response_states: a leg with no record starts from the "
            "nearest saved frequency (R4)");
      const auto changed = named_stages::response_engine_plans(
          named_stages::plan_response_states(optical, request, 2, "g-new",
                                             records));
      check(changed.engine.at(node0).kind == Kind::fresh
                && changed.engine.at(node1).kind == Kind::fresh
                && changed.resumed.empty(),
            "plan_response_states: records of another reference are not used "
            "(R1)");
      // S3c1 (D-045): observe, one engine sample as its record.
      const auto residuals = [](double first, double second) {
        return madness::function_residual_record(
            madness::FunctionResidualBlocks{{"alpha", {first, second}}});
      };
      json first_sample = {{"run_id", "r1"},
                           {"iteration", 0},
                           {"threshold", 1e-4},
                           {"k", 6},
                           {"has_previous", false},
                           {"input_energy", -1.1},
                           {"density_alpha", 0.0},
                           {"density_beta", 0.0},
                           {"density_charge", 0.0},
                           {"density_spin", 0.0},
                           {"previous_bsh_residual", 0.0},
                           {"previous_function_residuals", nullptr},
                           {"input_kinetic", 1.2},
                           {"input_kinetic_diagonal", residuals(0.5, 0.25)},
                           {"energy_change", 0.0},
                           {"kinetic_change", 0.0},
                           {"explicit_accuracy", json::object()},
                           {"configured_progress", json::object()},
                           {"engine_check_passed", false},
                           {"dconv", 1e-5},
                           {"conv_only_dens", false},
                           {"econv", 1e-6},
                           {"requested_econv", nullptr},
                           {"diagnostic_iterations_per_rung", 0},
                           {"decided_by", "engine_check"},
                           {"earlier_runs_iteration_history", json::array()}};
      json carried = json::object();
      const auto first = moldft.observe(ground, 2, first_sample, carried);
      const auto& first_evidence
          = first.at("convergence_report").at("measurements");
      check(first.at("energy").at("value") == -1.1
                && !first.contains("density_change")
                && first_evidence.at("density_change").at("availability")
                       == "unavailable"
                && first_evidence.at("density_change").at("value").is_null()
                && first.at("kinetic_trace").at("value") == 0.75
                && first.at("kinetic_trace").at("sample").at("stage") == "GS"
                && first.at("convergence_report").at("iteration_history").size()
                       == 1
                && !first.at("convergence_report")
                        .contains("earlier_runs_iteration_history"),
            "spec observe (moldft): a first sample has no changes, one row of "
            "history");
      check(first_evidence.at("energy_change").at("set_by")
                    == "engine_default_not_a_target"
                && first_evidence.at("energy_change").at("engine_default")
                       == 1e-6,
            "spec observe (moldft): econv not written is the engine's, not a "
            "target");
      auto second_sample = first_sample;
      second_sample["iteration"] = 1;
      second_sample["has_previous"] = true;
      second_sample["density_alpha"] = 2e-5;
      second_sample["density_beta"] = 3e-5;
      second_sample["previous_bsh_residual"] = 6e-4;
      second_sample["previous_function_residuals"] = residuals(1e-4, 5e-5);
      second_sample["input_kinetic_diagonal"] = residuals(0.5, 0.5);
      second_sample["energy_change"] = 1e-7;
      second_sample["requested_econv"] = 1e-6;
      second_sample["earlier_runs_iteration_history"]
          = json::array({{{"x", 1}}});
      const auto second = moldft.observe(ground, 2, second_sample, carried);
      const auto& evidence = second.at("convergence_report").at("measurements");
      // The density limit is max(dconv, thresh) x max(5, natom); the BSH
      // residual's is 5 x max(dconv, thresh), sampled one iteration back.
      check(
          second.at("density_change").at("beta") == 3e-5
              && evidence.at("density_change").at("value") == 3e-5
              && evidence.at("density_change").at("effective_limit") == 1e-4 * 5
              && evidence.at("density_change").at("satisfied") == true
              && evidence.at("bsh_residual").at("effective_limit") == 5.0 * 1e-4
              && evidence.at("bsh_residual").at("satisfied") == false
              && evidence.at("bsh_residual").at("sample").at("iteration") == 0
              && evidence.at("bsh_residual").at("per_function")
                     == residuals(1e-4, 5e-5),
          "spec observe (moldft): changes measured against their limits");
      check(evidence.at("energy_change").at("set_by") == "user"
                && evidence.at("energy_change").at("relation")
                       == "less_than_or_equal"
                && evidence.at("energy_change").at("satisfied") == true,
            "spec observe (moldft): a written econv is reported against, <=");
      check(
          evidence.at("delta_kinetic_energy").at("value") == 0.25
              && second.at("convergence_report").at("iteration_history").size()
                     == 2
              && second.at("convergence_report")
                         .at("earlier_runs_iteration_history")
                     == json::array({{{"x", 1}}})
              && carried.at("previous_kinetic_blocks")
                     == residuals(0.5, 0.5).at("blocks"),
          "spec observe (moldft): the kinetic trace change from the carried "
          "blocks, the history across samples and earlier runs");
      json none;
      for (const auto& [seed, work] :
           std::vector<std::pair<std::string, std::string>>{
               {"fd_restart", "resumed"},
               {"es_restart", "resumed"},
               {"nearest_frequency", "seeded"},
               {"fresh", "solved"}}) {
        const json event = {{"seed", seed}, {"iterations", 4}};
        auto expected = event;
        expected["work"] = work;
        std::string why = "spec observe (response): seed ";
        why += seed;
        why += " is ";
        why += work;
        check(molresponse.observe(optical, 2, event, none) == expected, why);
      }
      // consumer_context (A1): what a consuming stage reads, from the
      // engine's description; response has none.
      const json engine
          = {{"molecule_parameters",
              {{"units", "angstrom"},
               {"no_orient", false},
               {"source_name", "h2.in"},
               {"source_type", "file"},
               {"eprec", 1e-4}}},
             {"atoms",
              {{{"atomic_number", 1}, {"x", 0.0}, {"y", 0.0}, {"z", -0.7}},
               {{"atomic_number", 1}, {"x", 0.0}, {"y", 0.0}, {"z", 0.7}}}},
             {"charge", 0.0},
             {"nalpha", 1},
             {"nbeta", 1},
             {"L", 50.0},
             {"functional", "hf"},
             {"spin_restricted", true},
             {"thresh", 1e-6},
             {"k", 8}};
      const auto context = moldft.consumer_context(engine);
      check(context.at("molecule").at("symbols") == json{"H", "H"}
                && context.at("molecule").at("geometry")
                       == json{{0.0, 0.0, -0.7}, {0.0, 0.0, 0.7}}
                && context.at("molecule").at("parameters")
                       == json{{"units", "atomic"},
                               {"no_orient", true},
                               {"eprec", 1e-4}}
                && context.at("multiplicity") == 1 && context.at("L") == 50.0
                && context.at("thresh") == 1e-6 && context.at("k") == 8,
            "spec consumer_context (moldft): bohr, not reoriented, no file "
            "provenance, the multiplicity from the electrons");
      check(!molresponse.consumer_context,
            "spec consumer_context (response): none");
      // D-087: another stage's manifest that exists but cannot be read is
      // reported; a generation directory with no manifest references nothing.
      {
        namespace fs = std::filesystem;
        const auto root
            = fs::temp_directory_path()
              / ("test_named_stages_refs_" + std::to_string(::getpid()));
        fs::remove_all(root);
        const auto own = root / "gs" / "checkpoints";
        const auto other = root / "optical" / "checkpoints";
        fs::create_directories(own);
        fs::create_directories(other / "g-optical-1-1");
        fs::create_directories(other / "g-optical-1-2");
        std::ofstream(other / "g-optical-1-1" / "manifest.json") << "{\"half";
        const auto refs
            = named_stages::referenced_generations(world, own, "gs");
        check(refs.generations.empty() && refs.unreadable.size() == 1
                  && refs.unreadable[0].find("g-optical-1-1")
                         != std::string::npos,
              "referenced_generations (D-087): an unreadable manifest is "
              "reported, a directory without one is not");
        fs::remove_all(root);
      }
      // S3c2 (D-052, D-055): the skeleton's publish refuses a generation
      // before it is published; the refused one stays unpublished.
      {
        namespace fs = std::filesystem;
        const auto store
            = fs::temp_directory_path()
              / ("test_named_stages_publish_" + std::to_string(::getpid()))
              / "checkpoints";
        fs::remove_all(store.parent_path());
        // S4c1 (D-058): the skeleton names each attempt, refused or not.
        named_stages::GenerationSequence names{ground.name, "1"};
        std::vector<std::string> warnings;
        const auto name = [&](int n) {
          return "g-" + ground.name + "-1-" + std::to_string(n);
        };
        // S4c3: the skeleton assembles the manifest with the AppSpec's
        // manifest field. A stand-in field here returns the manifest each
        // row shapes (`assembled`) and keeps the publish sample it was given.
        Manifest assembled;
        json given;
        auto assembling = moldft;
        assembling.manifest
            = [&](const named_stages::Stage&, int, const json& sample) {
                given = sample;
                return assembled;
              };
        const auto publish = [&](named_stages::GenerationWrite write,
                                 const Manifest& manifest) {
          assembled = manifest;
          return named_stages::publish_generation(world, store, assembling,
                                                  ground, 2, names, warnings,
                                                  std::move(write));
        };
        const auto generation = [&] {
          auto m = scf;
          m.archives.clear();
          m.numerical["threshold"] = 1e-4;
          m.numerical["k"] = 6;
          m.numerical["states"] = json::array(
              {named_stages::moldft_state_record(m, "engine_check")});
          return m;
        };
        const auto writing = [&](const char* row) {
          named_stages::GenerationWrite write;
          write.publish_sample = {{"row", row}};
          write.write_archives = [](const fs::path&, Manifest& m) {
            // The inventory, without files: the checks before publish read
            // the manifest only.
            m.archives = {{"orbitals", "state.restartdata", {}, 1, {}}};
            return std::string{};
          };
          return write;
        };
        auto refused_early = writing("refused early");
        auto early = generation();
        early.generation = "g-set-by-the-application";
        early.properties = json::array({{{"name", "energy"}, {"value", -1.0}}});
        early.working_orbitals = {{"ordering", "archive_order"}};
        refused_early.iterations_this_run = 7;
        Manifest seen;
        refused_early.check_before_write = [&](const Manifest& m) {
          seen = m;
          throw std::runtime_error("not before writing");
        };
        rejects_with([&] { (void)publish(refused_early, early); },
                     "not before writing",
                     "publish_generation: the application's check before "
                     "writing refuses");
        check(!fs::exists(store / name(1)),
              "publish_generation: nothing written after that refusal");
        check(given == json{{"row", "refused early"}},
              "publish_generation (S4c3): the manifest field reads the "
              "application's publish sample");
        check(seen.generation == name(1) && seen.stage == ground.name
                  && seen.position.at("iterations_this_run") == 7
                  && seen.properties.at(0).at("generation") == name(1)
                  && seen.working_orbitals.at("state_id") == name(1),
              "publish_generation (S4c1): the shared fields are filled before "
              "the check before writing, the skeleton's name replacing the "
              "application's");
        auto unlisted = writing("unlisted");
        unlisted.write_archives
            = [](const fs::path&, Manifest&) { return std::string{}; };
        rejects_with([&] { (void)publish(unlisted, generation()); },
                     "does not hold its members' archives",
                     "publish_generation: archives that are not the members");
        auto invalid = generation();
        invalid.numerical["representation"] = "ao";
        rejects_with([&] { (void)publish(writing("invalid"), invalid); },
                     "unsupported representation",
                     "publish_generation (D-052): a generation that fails "
                     "validate_generation");
        check(fs::is_directory(store / name(3))
                  && !fs::exists(store / name(3) / "manifest.json")
                  && !fs::exists(store / "latest"),
              "publish_generation: the refused generation stays unpublished");
        auto unrecorded = generation();
        unrecorded.numerical.erase("states");
        rejects_with([&] { (void)publish(writing("unrecorded"), unrecorded); },
                     "has no state records",
                     "publish_generation (S4a): a generation without state "
                     "records");
        check(!fs::exists(store / name(4) / "manifest.json")
                  && !fs::exists(store / "latest"),
              "publish_generation: that generation stays unpublished");
        // S4d2 (D-077): a final write's sample takes the shared verdict of
        // the state records it assembles to, before the manifest is
        // assembled for publishing. Each row stops at the check before
        // writing.
        for (const bool at_final_rung : {true, false}) {
          auto last = generation();
          json state = named_stages::moldft_state_record(last, "engine_check");
          state["stop"] = "converged";
          if (at_final_rung)
            state["rung"] = moldft.requested_states(ground, 2)[0].at("rung");
          last.numerical["states"] = json::array({state});
          auto write = writing("final");
          write.final = true;
          write.check_before_write
              = [](const Manifest&) { throw std::runtime_error("stop here"); };
          rejects_with([&] { (void)publish(write, last); }, "stop here",
                       "publish_generation (S4d2): a final write");
          check(given.at("converged") == at_final_rung,
                at_final_rung ? "publish_generation (D-077): the final sample "
                                "takes the verdict converged"
                              : "publish_generation (D-077): a state off the "
                                "final rung makes the final sample "
                                "unconverged");
        }
        // D-072: the attempt's copies of the latest generation's shards go;
        // its text output and records stay.
        {
          const auto attempt = store.parent_path() / "attempt-test";
          fs::create_directories(attempt / "es_bundle");
          for (const auto* name :
               {"state.restartdata.00000", "state.restartdata.00001",
                "resolved.in", "es_bundle/root_0.00000",
                "es_bundle/roots.json"})
            std::ofstream(attempt / name) << "x";
          Manifest latest;
          latest.archives = {{"orbitals",
                              "state.restartdata",
                              {},
                              2,
                              {{"state.restartdata.00000", 1, ""},
                               {"state.restartdata.00001", 1, ""}}},
                             {"response",
                              "es_bundle/root_0",
                              {},
                              1,
                              {{"es_bundle/root_0.00000", 1, ""}}}};
          const auto warning
              = named_stages::delete_attempt_copies(world, attempt, latest);
          check(warning.empty()
                    && !fs::exists(attempt / "state.restartdata.00000")
                    && !fs::exists(attempt / "state.restartdata.00001")
                    && !fs::exists(attempt / "es_bundle/root_0.00000")
                    && fs::exists(attempt / "resolved.in")
                    && fs::exists(attempt / "es_bundle/roots.json"),
                "delete_attempt_copies (D-072): the shards go, the records "
                "stay");
          check(named_stages::delete_attempt_copies(world, attempt, latest)
                    .empty(),
                "delete_attempt_copies: nothing left to delete is no warning");
        }
        fs::remove_all(store.parent_path());
      }
      // S4a (D-056, D-059, D-063): one per-state record for both
      // applications, one reader, one check on both sides (D-052).
      {
        namespace planning = madness::planning;
        // A moldft generation as written before S4a: its rung record and its
        // report, no numerical.states.
        Manifest before;
        before.physical = {{"schema", "madqc.moldft.v1"}, {"molecule", "H2"}};
        before.numerical = {{"threshold", 1e-6}, {"k", 8}};
        const auto saved_at = [&](double thresh, int k, const char* stop) {
          auto m = before;
          m.position = {{"rung",
                         {{"thresh", thresh},
                          {"k", k},
                          {"iterations_at_rung", 5},
                          {"stop", stop}}}};
          m.metrics
              = {{"convergence_report",
                  {{"measurements",
                    {{"density_change",
                      {{"value", 2e-5},
                       {"phase", "successive_iteration_input_densities"}}}}},
                   {"saved_state_validation",
                    {{"measurements",
                      {{"density_change",
                        {{"value", 3e-5},
                         {"phase", "saved_orbitals_before_archive"}}}}}}}}}};
          return m;
        };
        planning::LadderRequest request;
        request.identity = before.physical.dump();
        request.ladder = {{1e-4, 6}, {1e-6, 8}};
        request.engine_limits = {{"density_change", 1e-4, 0.0, false}};
        auto elsewhere = request;
        elsewhere.identity = "another problem";
        struct PlanRow {
          Manifest manifest;
          planning::LadderRequest request;
          planning::Action action;
          const char* why;
        };
        for (const auto& row : std::vector<PlanRow>{
                 {saved_at(1e-6, 8, "converged"), request,
                  planning::Action::load, "loaded"},
                 {saved_at(1e-4, 6, "iteration_limit"), request,
                  planning::Action::resume, "resumed"},
                 {saved_at(1e-6, 8, "converged"), elsewhere,
                  planning::Action::fresh, "fresh"}}) {
          auto recorded = row.manifest;
          recorded.numerical["states"] = json::array(
              {named_stages::moldft_state_record(recorded, "engine_check")});
          // Before S4a the planner read the manifest (saved_state); now it
          // reads the record, written or (S4b) rebuilt by state_records.
          const auto old_plan
              = planning::plan_ladder(named_stages::saved_state(row.manifest),
                                      std::nullopt, row.request);
          const auto new_plan = planning::plan_ladder(
              named_stages::saved_state_from_record(
                  moldft.state_records(ground, 2, recorded).at(0)),
              std::nullopt, row.request);
          const auto rebuilt = planning::plan_ladder(
              named_stages::saved_state_from_record(
                  moldft.state_records(ground, 2, row.manifest).at(0)),
              std::nullopt, row.request);
          check(rebuilt.action == old_plan.action
                    && rebuilt.start_rung == old_plan.start_rung
                    && rebuilt.why == old_plan.why,
                std::string("spec state_records (moldft, S4b): a generation "
                            "before S4a plans as before, ")
                    + row.why);
          check(old_plan.action == row.action
                    && new_plan.action == old_plan.action
                    && new_plan.start_rung == old_plan.start_rung
                    && new_plan.iterations_done == old_plan.iterations_done
                    && new_plan.why == old_plan.why,
                std::string("moldft state record (S4a): plans as the "
                            "generation before it, ")
                    + row.why);
        }
        // The record on both sides (D-052): absent, accepted on the read
        // side (publish_generation refuses it, above); present, it must be
        // the one the generation's own rung record and report give.
        auto gs = scf;
        gs.numerical["threshold"] = 1e-4;
        gs.numerical["k"] = 6;
        moldft.validate_generation(ground, 2, gs);
        check(true, "spec validate_generation (moldft): a generation without "
                    "state records, written before S4a, is read");
        gs.numerical["states"]
            = json::array({named_stages::moldft_state_record(gs, "planner")});
        moldft.validate_generation(ground, 2, gs);
        check(gs.numerical.at("states").at(0).at("kind") == "ground_state"
                  && gs.numerical.at("states").at(0).at("schema")
                         == "madqc.state.v1",
              "spec validate_generation (moldft): its state record");
        auto wrong = gs;
        wrong.numerical["states"][0]["rung"]["k"] = 7;
        rejects_with([&] { moldft.validate_generation(ground, 2, wrong); },
                     "state record disagrees with its generation",
                     "spec validate_generation (moldft): a record with "
                     "another rung");
        wrong = gs;
        wrong.numerical["states"][0]["schema"] = "madqc.other.v1";
        rejects_with([&] { moldft.validate_generation(ground, 2, wrong); },
                     "unknown schema",
                     "spec validate_generation (moldft): an unknown schema");
        wrong = gs;
        wrong.numerical["states"].push_back(wrong.numerical["states"][0]);
        rejects_with([&] { moldft.validate_generation(ground, 2, wrong); },
                     "needs one",
                     "spec validate_generation (moldft): two records");
        // Response: FD records get kind fd; every measurement row a phase,
        // null for a row with no sample (D-063); bundle rows the engine's.
        auto sampled = fd;
        sampled["fd_states"][channel][key][fkey]["convergence_measurements"]
            = {{"bsh_residual",
                {{"value", 1e-5},
                 {"phase", "solver_step_input_before_acceleration"},
                 {"sample", {{"iteration", 2}}}}},
               {"energy_change",
                {{"value", nullptr}, {"availability", "not_applicable"}}}};
        const auto fd_records
            = named_stages::response_state_records(sampled, "g-ref", {});
        const auto& rows = fd_records.at(0).at("engine_measurements");
        check(fd_records.at(0).at("schema") == "madqc.state.v1"
                  && fd_records.at(0).at("kind") == "fd"
                  && rows.at("bsh_residual").at("phase")
                         == "solver_step_input_before_acceleration"
                  && rows.at("energy_change").at("phase").is_null(),
              "response_state_records (S4a): kind fd; a sampled row keeps "
              "its phase, a row with no sample has phase null");
        auto measured = metadata;
        measured["excited_states"][key]["roots"][0]["density_residual"] = 1e-5;
        measured["excited_states"][key]["roots"][0]["omega_residual"] = 1e-6;
        const auto bundle_rows
            = named_stages::response_state_records(measured, "g-ref", {})
                  .at(0)
                  .at("engine_measurements");
        check(bundle_rows.at("density_residual").at("phase")
                      == "solver_step_input_before_acceleration"
                  && bundle_rows.at("omega_residual").at("phase")
                         == "solver_step_output"
                  && bundles.at(0)
                         .at("engine_measurements")
                         .at("density_residual")
                         .at("phase")
                         .is_null(),
              "response_state_records (S4a): a bundle's rows carry the "
              "engine's phase, null with no value");
        named_stages::validate_state_records(fd_records);
        named_stages::validate_state_records(bundles);
        check(true, "validate_state_records: response's FD and bundle "
                    "records");
        // The old response name is still read, without kind or phases.
        auto old = fd_records;
        old[0]["schema"] = "madqc.response.state.v1";
        old[0].erase("kind");
        old[0]["engine_measurements"]["energy_change"].erase("phase");
        named_stages::validate_state_records(old);
        const auto fd_reader = molresponse_v3::fd_accuracy_measurements;
        const auto read_old
            = named_stages::saved_state_from_record(old[0], fd_reader);
        const auto read_new
            = named_stages::saved_state_from_record(fd_records[0], fd_reader);
        check(read_old.identity == read_new.identity
                  && read_old.rung.thresh == read_new.rung.thresh
                  && read_old.rung.k == read_new.rung.k
                  && read_old.reference_generation == "g-ref"
                  && read_old.iterations_at_rung == read_new.iterations_at_rung
                  && read_old.stop == read_new.stop,
              "saved_state_from_record: the pre-S4a response name is read");
        auto unknown = fd_records;
        unknown[0]["schema"] = "madqc.other.v1";
        rejects_with(
            [&] { (void)named_stages::saved_state_from_record(unknown[0]); },
            "unknown schema", "saved_state_from_record: an unknown schema");
        rejects_with([&] { named_stages::validate_state_records(unknown); },
                     "unknown schema",
                     "validate_state_records: an unknown schema");
        auto phaseless = fd_records;
        phaseless[0]["engine_measurements"]["energy_change"].erase("phase");
        rejects_with([&] { named_stages::validate_state_records(phaseless); },
                     "has no phase",
                     "validate_state_records: a current record's row without "
                     "a phase");
        auto es_wrong = es;
        es_wrong.numerical["states"] = unknown;
        rejects_with(
            [&] { molresponse.validate_generation(optical, 2, es_wrong); },
            "unknown schema",
            "spec validate_generation (response): a present, wrong record");
        // S4b (D-048, D-064): the state_records and plan fields.
        auto reused = saved_at(1e-6, 8, "converged");
        reused.position["reused_from"] = "g-0-1";
        check(moldft.state_records(ground, 2, reused).at(0).at("decided_by")
                      == "planner"
                  && moldft.state_records(ground, 2,
                                          saved_at(1e-6, 8, "converged"))
                             .at(0)
                             .at("decided_by")
                         == "engine_check",
              "spec state_records (moldft): a rebuilt record says planner "
              "when the generation was published again with its verdict");
        auto with_record = saved_at(1e-4, 6, "iteration_limit");
        with_record.numerical["states"] = json::array(
            {named_stages::moldft_state_record(with_record, "engine_check")});
        check(moldft.state_records(ground, 2, with_record)
                  == with_record.numerical.at("states"),
              "spec state_records (moldft): the generation's own records");
        const named_stages::PlanContext problem{before.physical.dump(), {}};
        const auto ground_request
            = named_stages::moldft_ladder_request(ground, 2, problem);
        const auto& gcp = ground.params.get<CalculationParameters>();
        check(ground_request.identity == problem.identity
                  && ground_request.ladder.size() == gcp.protocol().size()
                  && ground_request.ladder.back().thresh
                         == gcp.protocol().back()
                  && ground_request.engine_limits.at(0).measurement
                         == "density_change",
              "moldft_ladder_request: the stage's ladder and the engine "
              "check's limits, the problem from the context");
        const auto ground_plans = moldft.plan(
            ground, 2, problem, moldft.state_records(ground, 2, with_record));
        const auto direct_plan
            = planning::plan_ladder(named_stages::saved_state(with_record),
                                    std::nullopt, ground_request);
        check(ground_plans.states.size() == 1
                  && ground_plans.states[0].plan.action == direct_plan.action
                  && ground_plans.states[0].plan.start_rung
                         == direct_plan.start_rung
                  && ground_plans.states[0].plan.why == direct_plan.why
                  && ground_plans.states[0].identity == "ground_state"
                  && ground_plans.states[0].own_rung->thresh == 1e-4
                  && ground_plans.continuation.at("planner") == "plan_ladder"
                  && ground_plans.continuation.at("why") == direct_plan.why
                  && ground_plans.continuation.at("rung").at("index")
                         == direct_plan.start_rung + 1,
              "spec plan (moldft): the planner's plan for the ground state, "
              "and its continuation record");
        const auto field_plans
            = molresponse.plan(optical, 2, {{}, "g-ref"}, records);
        const auto direct = named_stages::plan_response_states(
            optical, named_stages::response_request(optical), 2, "g-ref",
            records);
        check(field_plans.continuation == direct.continuation
                  && field_plans.states.size() == 2
                  && field_plans.states[0].fields.at("node") == node0
                  && field_plans.states[0].fields.at("key_at_start") == state,
              "spec plan (response): plan_response_states on the stage's "
              "request, with the engine names of each state");
        Manifest unrecorded_response;
        unrecorded_response.reference
            = {{"stage", "GS"}, {"generation", "g-ref"}};
        unrecorded_response.numerical = {{"engine_metadata", fd}};
        auto recorded_response = unrecorded_response;
        recorded_response.numerical["states"] = records;
        check(molresponse.state_records(optical, 2, unrecorded_response)
                      == named_stages::response_state_records(fd, "g-ref", {})
                  && molresponse.state_records(optical, 2, recorded_response)
                         == records,
              "spec state_records (response): records derived from the "
              "engine metadata against the generation's own reference");
      }
      // S4c1 (D-058): one name form, and the fields every manifest shares.
      {
        named_stages::GenerationSequence names{"excited-tda", "20261028-1234"};
        check(names.next() == "g-excited-tda-20261028-1234-1"
                  && names.next() == "g-excited-tda-20261028-1234-2",
              "GenerationSequence: g-<stage>-<run>-<n>, n from 1 in publish "
              "order");
        Manifest m;
        m.position = {{"phase", "property_analysis"},
                      {"previous_generation", "g-20261027-99-4"}};
        m.properties = json::array({{{"name", "alpha"}, {"value", 1.0}},
                                    {{"name", "beta"}, {"value", 2.0}}});
        named_stages::set_shared_fields(m, "optical", "g-optical-r-3", 12);
        check(m.stage == "optical" && m.generation == "g-optical-r-3"
                  && m.position
                         == json{{"phase", "property_analysis"},
                                 {"previous_generation", "g-20261027-99-4"},
                                 {"iterations_this_run", 12}}
                  && m.properties.at(0).at("generation") == "g-optical-r-3"
                  && m.properties.at(1).at("generation") == "g-optical-r-3"
                  && m.working_orbitals.is_null(),
              "set_shared_fields: stage, name, this run's iterations and the "
              "property records' generation; the previous generation is the "
              "application's; no working orbitals stay none (response)");
        m.working_orbitals
            = {{"state_id", "g-old"}, {"ordering", "archive_order"}};
        named_stages::set_shared_fields(m, "gs", "g-gs-r-1", 0);
        check(m.working_orbitals
                  == json{{"state_id", "g-gs-r-1"},
                          {"ordering", "archive_order"}},
              "set_shared_fields (A1): the working orbitals' state id is the "
              "generation (moldft, a republished state included)");
      }
      // S4c2 (D-054): publish_now, and its second output, the rung record.
      // `ground` has protocol [1e-4, 1e-6], maxiter 25 and
      // checkpoint.every_iterations 3.
      {
        const auto ladder = named_stages::requested_ladder(
            ground.params.get<CalculationParameters>());
        const auto sample
            = [&](std::size_t rung, int iteration, const char* stop) {
                return json{{"threshold", ladder.at(rung).thresh},
                            {"k", ladder.at(rung).k},
                            {"iteration", iteration},
                            {"rung_stop", stop}};
              };
        const json none = json::object();
        const auto at = [&](const named_stages::Stage& stage, const json& s,
                            const json& carried) {
          return moldft.publish_now(stage, 2, s, carried);
        };
        auto decision = at(ground, sample(1, 3, "converged"), none);
        check(decision.publish
                  && decision.rung
                         == json{{"index", 2},
                                 {"thresh", ladder.at(1).thresh},
                                 {"k", ladder.at(1).k},
                                 {"iterations_at_rung", 4},
                                 {"iterations_this_run", 4},
                                 {"stop", "converged"}},
              "spec publish_now (moldft): a rung stop publishes; the rung "
              "record by value, index from 1");
        check(!at(ground, sample(0, 0, "in_progress"), none).publish
                  && !at(ground, sample(0, 1, "in_progress"), none).publish
                  && at(ground, sample(0, 2, "in_progress"), none).publish
                  && !at(ground, sample(0, 3, "in_progress"), none).publish,
              "spec publish_now (moldft): within a rung, every K-th "
              "iteration at the rung (K = 3)");
        const json resumed
            = {{"resumed_rung", {{"start_rung", 0}, {"iterations_done", 2}}}};
        decision = at(ground, sample(0, 0, "in_progress"), resumed);
        check(decision.publish && decision.rung.at("iterations_at_rung") == 3
                  && decision.rung.at("iterations_this_run") == 1,
              "spec publish_now (moldft): a resumed rung counts the "
              "planner's iterations_done toward K (A3: the cadence holds "
              "over a restart)");
        decision = at(ground, sample(1, 0, "in_progress"), resumed);
        check(!decision.publish && decision.rung.at("iterations_at_rung") == 1,
              "spec publish_now (moldft): only the resumed rung carries the "
              "earlier iterations");
        auto capped = ground;
        capped.checkpoint_every_iterations = 0;
        capped.params.get<CalculationParameters>().set_user_defined_value(
            "maxiter_per_rung", std::vector<int>{30, 2});
        check(!at(capped, sample(1, 0, "in_progress"), none).publish
                  && at(capped, sample(1, 1, "in_progress"), none).publish
                  && !at(capped, sample(0, 1, "in_progress"), none).publish,
              "spec publish_now (moldft): this run's iterations reaching the "
              "rung's cap publish (maxiter [30, 2]), at that rung only");
        auto unknown = sample(0, 0, "in_progress");
        unknown["threshold"] = 1e-5;
        rejects_with([&] { (void)at(ground, unknown, none); },
                     "sample rung is not in the protocol",
                     "spec publish_now (moldft): a sample off the ladder");
        decision = molresponse.publish_now(optical, 2, json::object(), none);
        check(decision.publish && decision.rung.is_null(),
              "spec publish_now (response): every save event, no rung");
      }
      // S4c3 (D-054, D-067): the manifest field against recorded publish
      // samples, in test_named_stages_publish_samples.json beside this file:
      // H2 runs (the record check's deck; the reuse case from
      // test_named_reference_from_generation), each with its deck and stage.
      // The expected manifest is the one the class assembled for that sample
      // (the S4c3 guard), before the shared fields; the histories the fields
      // copy unchanged are cut to their first entry on both sides.
      {
        const auto path = std::filesystem::path(__FILE__).parent_path()
                          / "test_named_stages_publish_samples.json";
        std::ifstream in(path);
        check(in.good(), "publish samples: cannot read " + path.string());
        const auto cases = json::parse(in);
        check(cases.size() == 5, "publish samples: five cases");
        for (const auto& [name, recorded] : cases.items()) {
          const bool is_moldft = name.starts_with("moldft");
          const auto stages = named_stages::parse(
              world, recorded.at("deck").get<std::string>());
          const auto stage
              = std::find_if(stages.begin(), stages.end(), [&](const auto& s) {
                  return s.name == recorded.at("stage").get<std::string>();
                });
          const auto& spec = is_moldft ? moldft : molresponse;
          const auto assembled
              = spec.manifest(*stage, 2, recorded.at("sample"));
          check(madness::checkpoint::to_json_record(assembled)
                    == recorded.at("manifest"),
                "spec manifest (" + name
                    + "): the class's manifest for the recorded sample");
        }
        auto unknown = cases.at("moldft_save").at("sample");
        unknown["kind"] = "other";
        const auto stages = named_stages::parse(
            world, cases.at("moldft_save").at("deck").get<std::string>());
        rejects_with([&] { (void)moldft.manifest(stages[0], 2, unknown); },
                     "unknown publish sample kind",
                     "spec manifest (moldft): a sample of no known kind");
        // R-039: the reuse publish takes the skeleton's verdict, as the
        // other final publishes do; without it the loaded stop reason stands.
        auto reuse = cases.at("moldft_reuse").at("sample");
        const auto reuse_stages = named_stages::parse(
            world, cases.at("moldft_reuse").at("deck").get<std::string>());
        reuse["converged"] = false;
        check(moldft.manifest(reuse_stages[0], 2, reuse).stop_reason
                  == reuse.at("loaded").at("stop_reason").get<std::string>(),
              "spec manifest (moldft reuse): an unconverged verdict keeps the "
              "loaded stop reason");
        reuse["converged"] = true;
        check(moldft.manifest(reuse_stages[0], 2, reuse).stop_reason
                  == "converged",
              "spec manifest (moldft reuse): the verdict converged");
      }
      // S4d1 (D-046, D-057): the shared verdict from the state records.
      {
        const auto needed = moldft.requested_states(ground, 2);
        check(needed.size() == 1 && needed[0].at("key") == "ground_state"
                  && needed[0].at("rung").at("thresh") == 1e-6,
              "spec requested_states (moldft): the ground state at the final "
              "rung");
        const auto final_rung = needed[0].at("rung");
        const json measured = {
            {"availability", "measured"}, {"units", "a.u."}, {"value", 1e-7}};
        const auto record = [&](const char* stop, const char* decided_by) {
          return json{{"key", "ground_state"},
                      {"rung", final_rung},
                      {"stop", stop},
                      {"diverged", false},
                      {"decided_by", decided_by},
                      {"saved_measurements",
                       {{"bsh_residual", measured},
                        {"charge_density_change",
                         {{"availability", "measured"},
                          {"units", "a0^(-3/2)"},
                          {"value", 1e-6}}}}}};
        };
        const auto verdict
            = [&](const json& records,
                  const madness::convergence::ConvergenceSettings& settings) {
                Manifest published;
                published.numerical["states"] = records;
                return named_stages::shared_verdict(moldft, ground, 2, settings,
                                                    published);
              };
        auto settings = moldft.settings(ground, 2);
        settings.decider = madness::convergence::Decider::engine_check;
        check(verdict(json::array({record("converged", "engine_check")}),
                      settings)
                  .at("converged"),
              "shared verdict: the needed state stopped converged at its rung");
        check(!verdict(json::array({record("iteration_limit", "engine_check")}),
                       settings)
                   .at("converged"),
              "shared verdict: a state that hit its iteration limit");
        check(verdict(json::array({record("iteration_limit", "planner")}),
                      settings)
                  .at("converged"),
              "shared verdict (D-016): a state loaded under a met request");
        auto coarse = record("converged", "engine_check");
        coarse["rung"]["thresh"] = 1e-4;
        check(!verdict(json::array({coarse}), settings).at("converged"),
              "shared verdict: a state converged at another rung");
        auto diverged = record("converged", "engine_check");
        diverged["diverged"] = true;
        check(!verdict(json::array({diverged}), settings).at("converged"),
              "shared verdict: a diverged state");
        // D-085, R3: a finished state finer than the final rung counts, its
        // record keeping its true rung; a diverged one does not.
        auto finer = record("converged", "engine_check");
        finer["rung"]["thresh"] = 1e-7;
        const auto beyond = verdict(json::array({finer}), settings);
        check(beyond.at("converged")
                  && !beyond.at("states")[0].at("at_rung").get<bool>()
                  && beyond.at("states")[0].at("beyond_rung").get<bool>(),
              "shared verdict (D-085): a converged state finer than the final "
              "rung, beyond_rung");
        auto finer_diverged = finer;
        finer_diverged["diverged"] = true;
        check(!verdict(json::array({finer_diverged}), settings).at("converged"),
              "shared verdict (D-085): a diverged state finer than the final "
              "rung");
        // D-091: an R2 load leaves the stalled generation's own record.
        const auto stalled = verdict(
            json::array({record("stalled", "engine_check")}), settings);
        check(!stalled.at("converged")
                  && stalled.at("states")[0].contains("requested_accuracy"),
              "shared verdict (D-091): a stalled final rung judged by the "
              "engine check is not converged, its requested_accuracy "
              "reported");
        check(named_stages::republish_on_load("iteration_limit", {})
                  && !named_stages::republish_on_load("converged", {})
                  && !named_stages::republish_on_load(
                      "stalled", {"stalled at rung 1 with the request unmet"}),
              "republish_on_load (D-091): only a met load of a generation "
              "that did not stop converged is published again; an R2 load is "
              "not");
        const auto missing = verdict(json::array(), settings);
        check(!missing.at("converged")
                  && !missing.at("states")[0].at("found").get<bool>(),
              "shared verdict: a needed state without a record");
        auto criteria = settings;
        criteria.decider = madness::convergence::Decider::requested_criteria;
        for (auto& metric : criteria.metrics)
          if (metric.measurement == "bsh_residual") {
            metric.group = madness::convergence::Group::requested_criteria;
            metric.target = 1e-8;
            metric.units = "a.u.";
            metric.quantity.clear();
          }
        const auto unmet = verdict(
            json::array({record("converged", "engine_check")}), criteria);
        check(!unmet.at("converged")
                  && unmet.at("states")[0].at("requested_accuracy").at("status")
                         == "unmet",
              "shared verdict (B8): under the written targets, a converged "
              "stop whose saved state misses them, the unmet target listed");
        check(verdict(json::array({record("iteration_limit", "planner")}),
                      criteria)
                  .at("converged"),
              "shared verdict (D-016): a loaded state's met request stands "
              "under the written targets too");
        const auto response_needed = molresponse.requested_states(optical, 2);
        check(response_needed.size() == 2
                  && response_needed[0].at("key")
                         == "dipole_z|1e-04_k6|f0.00000",
              "spec requested_states (response): one FD state per requested "
              "frequency at the final protocol");
      }
      // D-082, D-061: each state judged on the written targets its kind has
      // rows for; the bundle on econv (its engine omega row).
      {
        namespace cv = madness::convergence;
        const auto es_deck = [&](const std::string& properties,
                                 const std::string& extra) {
          return named_stages::parse(
              world,
              geometry + gs
                  + "molresponse es\n reference GS\n protocol [1e-4]\n"
                    " dconv 1e-4\n requested_properties ["
                  + properties + "]\n"
                  + (properties.find("polarizability") == std::string::npos
                         ? ""
                         : " polarizability.frequencies [0.0]\n"
                           " polarizability.frequency_units au\n"
                           " polarizability.components [zz]\n")
                  + " excited_states.roots 1\n excited_states.method tda\n"
                  + extra + "end\n")[1];
        };
        const auto bundle_record = [&](const named_stages::Stage& stage,
                                       double omega, bool engine_rows = true) {
          const auto needed = molresponse.requested_states(stage, 2).back();
          json record = {{"key", needed.at("key")},
                         {"kind", "excited_states"},
                         {"rung", needed.at("rung")},
                         {"stop", "converged"},
                         {"diverged", false},
                         {"decided_by", "engine_check"},
                         {"saved_measurements", json::object()}};
          if (engine_rows)
            record["engine_measurements"]
                = {{"omega_residual",
                    {{"value", omega},
                     {"availability", "measured"},
                     {"root_id", "r0"},
                     {"phase", "solver_step_output"}}}};
          return record;
        };
        const auto judge = [&](const named_stages::Stage& stage,
                               const json& records) {
          Manifest published;
          published.numerical["states"] = records;
          return named_stages::shared_verdict(
              molresponse, stage, 2, molresponse.settings(stage, 2), published);
        };
        const std::string criteria
            = " convergence.use_default_criteria false\n";
        // D-092: criteria mode on an excited-states-only stage needs econv,
        // so the stage whose dconv no needed state measures is in engine mode
        // and the bundle judged on its stop has dipole legs beside it.
        const auto es_dconv = es_deck("excited_states", "");
        const auto mixed_dconv
            = es_deck("polarizability, excited_states", criteria);
        const auto needed = molresponse.requested_states(es_dconv, 2);
        check(needed.size() == 1
                  && needed[0].at("measurements")
                         == json::array({"omega_residual"}),
              "requested_states (D-082): the bundle's rows are omega_residual");
        const auto alone
            = judge(es_dconv, json::array({bundle_record(es_dconv, 1e-5)}));
        check(alone.at("converged")
                  && alone.at("states")[0].at("requested_accuracy").at("status")
                         == "not_requested"
                  && alone.at("unmeasured_targets")
                         == json::array({"charge_density_change"}),
              "shared verdict (D-082): dconv, which no needed state measures, "
              "listed once, not a failed state");
        const auto mixed_needed_dconv
            = molresponse.requested_states(mixed_dconv, 2);
        const json planner_leg = {{"key", mixed_needed_dconv[0].at("key")},
                                  {"kind", "fd"},
                                  {"rung", mixed_needed_dconv[0].at("rung")},
                                  {"stop", "iteration_limit"},
                                  {"diverged", false},
                                  {"decided_by", "planner"},
                                  {"saved_measurements", json::object()}};
        const auto on_stop = judge(
            mixed_dconv,
            json::array({planner_leg, bundle_record(mixed_dconv, 1e-5)}));
        check(
            on_stop.at("converged")
                && on_stop.at("states")[1].at("requested_accuracy").at("status")
                       == "not_requested"
                && on_stop.at("unmeasured_targets").empty(),
            "shared verdict (D-082): a bundle under written targets none of "
            "its kind's is judged on its stop");
        const auto es_econv
            = es_deck("excited_states", criteria + " econv 1e-6\n");
        const auto met
            = judge(es_econv, json::array({bundle_record(es_econv, 1e-7)}));
        const auto missed
            = judge(es_econv, json::array({bundle_record(es_econv, 1e-5)}));
        check(
            met.at("converged")
                && met.at("states")[0].at("requested_accuracy").at("status")
                       == "satisfied"
                && !missed.at("converged")
                && missed.at("states")[0].at("requested_accuracy").at("status")
                       == "unmet",
            "shared verdict (D-061, D-082): the bundle judged on econv, met "
            "and unmet");
        check(
            !judge(es_econv, json::array({bundle_record(es_econv, 0.0, false)}))
                 .at("converged"),
            "shared verdict (D-082): a bundle record without engine rows (a "
            "replay capture) is unavailable on econv, not an error");
        const auto es_engine = es_deck("excited_states", " econv 1e-6\n");
        const auto reported
            = judge(es_engine, json::array({bundle_record(es_engine, 1e-5)}));
        check(reported.at("converged")
                  && reported.at("states")[0]
                             .at("requested_accuracy")
                             .at("status")
                         == "unmet",
              "shared verdict (D-061): in engine mode an unmet econv is "
              "reported, the bundle converged");
        // A stage with legs and a bundle: each kind on its own rows.
        const auto mixed = es_deck("polarizability, excited_states",
                                   criteria + " econv 1e-6\n");
        const auto mixed_needed = molresponse.requested_states(mixed, 2);
        json leg = {{"key", mixed_needed[0].at("key")},
                    {"kind", "fd"},
                    {"rung", mixed_needed[0].at("rung")},
                    {"stop", "iteration_limit"},
                    {"diverged", false},
                    {"decided_by", "planner"},
                    {"saved_measurements", json::object()}};
        const auto both
            = judge(mixed, json::array({leg, bundle_record(mixed, 1e-7)}));
        check(both.at("converged") && both.at("unmeasured_targets").empty()
                  && !both.at("states")[0]
                          .at("requested_accuracy")
                          .at("criteria")
                          .contains("omega_residual")
                  && both.at("states")[0]
                         .at("requested_accuracy")
                         .at("criteria")
                         .contains("charge_density_change")
                  && both.at("states")[1]
                             .at("requested_accuracy")
                             .at("criteria")
                             .size()
                         == 1,
              "shared verdict (D-082): the leg judged on dconv, the bundle on "
              "econv alone, nothing unmeasured");
        // The settings row and the planner's per-kind sets.
        const auto settings = molresponse.settings(mixed, 2);
        const auto row = std::find_if(
            settings.metrics.begin(), settings.metrics.end(),
            [](const auto& m) { return m.measurement == "omega_residual"; });
        check(row != settings.metrics.end() && row->target == 1e-6
                  && std::none_of(settings.metrics.begin(),
                                  settings.metrics.end(),
                                  [](const auto& m) {
                                    return m.measurement == "energy_change";
                                  }),
              "settings (D-061): on an excited-state stage the econv row "
              "replaces the energy change row");
        const auto targets_of = [&](const char* kind) {
          return cv::requested_targets(named_stages::targets_of_kind(
              settings, named_stages::kind_measurements(kind)));
        };
        check(targets_of("fd").contains("charge_density_change")
                  && !targets_of("fd").contains("omega_residual")
                  && targets_of("excited_states")
                         == json{{"omega_residual",
                                  {{"value", 1e-6}, {"units", "hartree"}}}},
              "planner (D-082): the legs plan on their rows, the bundle on "
              "econv (the reset keeps the econv row)");
        check(named_stages::response_request(mixed).at("excited_state_targets")
                      == json{{"omega_residual",
                               {{"value", 1e-6}, {"units", "hartree"}}}}
                  && !named_stages::response_request(es_dconv).contains(
                      "excited_state_targets"),
              "response_request (D-061): econv beside accuracy_targets, only "
              "when written");
        // On --continue: a saved bundle that passes the engine check but not
        // econv is loaded and reported in engine mode (B8's reuse exception)
        // and resumed under the written targets.
        const auto es_key = molresponse_v3::protocol_key(
            1e-4, named_stages::protocol_k(1e-4));
        const json es_meta
            = {{"fd_states", json::object()},
               {"vbc_states", json::object()},
               {"protocols",
                {{es_key,
                  {{"thresh", 1e-4}, {"k", named_stages::protocol_k(1e-4)}}}}},
               {"excited_states",
                {{es_key,
                  {{"type", "tda"},
                   {"n_roots", 1},
                   {"iter", 6},
                   {"converged", true},
                   {"stop_reason", "converged"},
                   {"roots",
                    {{{"slot", 0},
                      {"root_id", "r0"},
                      {"omega", 0.5},
                      {"density_residual", 1e-6},
                      {"omega_residual", 1e-5}}}}}}}}};
        const auto es_records
            = named_stages::response_state_records(es_meta, "g-ref", {});
        const auto bundle_plan = [&](const named_stages::Stage& stage) {
          return named_stages::plan_response_states(
                     stage, named_stages::response_request(stage), 2, "g-ref",
                     es_records)
              .continuation.at("states")
              .back();
        };
        const auto engine_plan = bundle_plan(es_engine);
        const auto criteria_plan = bundle_plan(es_econv);
        check(engine_plan.at("action") == "loaded"
                  && engine_plan.at("unmet_not_enforceable")
                         == json::array({"omega_residual"})
                  && criteria_plan.at("action") != "loaded",
              "planner (D-061, D-082): a saved bundle missing econv, loaded "
              "and reported in engine mode, not loaded under written targets");
        check(bundle_plan(mixed_dconv).at("action") == "loaded",
              "planner (D-082): a bundle under written targets none of its "
              "kind's (dconv only) is loaded on its stop");
      }
      // S5: the Dalton seed's member kind.
      {
        namespace fs = std::filesystem;
        check(named_stages::sha256_hex("")
                      == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca4959"
                         "91b7852b855"
                  && named_stages::sha256_hex("abc")
                         == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb4"
                            "10ff61f20015ad"
                  && named_stages::sha256_hex(std::string(1000000, 'a'))
                         == "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e04"
                            "6d39ccc7112cd0",
              "sha256: the FIPS 180-4 test vectors");
        const auto dalton
            = fs::path(__FILE__).parent_path().parent_path().parent_path()
              / "examples/qc/response_h2o_tpa_seeded/dalton";
        const auto block = named_stages::dalton_seed_block(dalton.string());
        std::map<std::string, std::string> sums;
        for (const auto& f : block.at("files"))
          sums[f.at("path").get<std::string>()]
              = f.at("sha256").get<std::string>();
        check(sums.size() == 2
                  && sums.at("seed/molden.inp")
                         == "86cbac809c3cbe67e79443980fa118f3475aea493afd692264"
                            "863a78655ae3d3"
                  && sums.at("seed/RSPVEC")
                         == "a521dc7815608bea58b991a47ccdbdd34ee071a28873737556"
                            "59e87e03cd2091",
              "dalton_seed_block: the files the reader reads (molden.inp, "
              "RSPVEC; no .out here), by size and SHA-256");
        const auto generation
            = fs::temp_directory_path()
              / ("test_named_stages_seed_" + std::to_string(::getpid()));
        fs::remove_all(generation);
        fs::create_directories(generation / "seed");
        for (const auto& f : block.at("files"))
          fs::copy_file(f.at("source").get<std::string>(),
                        generation / f.at("path").get<std::string>());
        check(named_stages::check_seed_files(generation, block).empty(),
              "check_seed_files: intact copies");
        std::ofstream(generation / "seed/RSPVEC", std::ios::app) << "x";
        check(named_stages::check_seed_files(generation, block)
                      .find("seed/RSPVEC altered")
                  != std::string::npos,
              "check_seed_files: an altered copy");
        fs::remove(generation / "seed/molden.inp");
        check(named_stages::check_seed_files(generation, block)
                      .find("seed/molden.inp missing")
                  != std::string::npos,
              "check_seed_files: a missing copy");
        fs::remove_all(generation);
        // The record a driven leg is seeded from (dalton_fd_guess's choice,
        // nearest_dalton_record), on the qctest's RSPVEC: its ZDIPLEN
        // records are at 0.1609303480 and 0.1916703618 au; the qctest's
        // 2PA legs sit at 0.16125 and 0.19169.
        {
          const auto rsp = read_rspvec((dalton / "RSPVEC").string());
          double d = -1;
          const auto* near = molresponse_v3::nearest_dalton_record(
              rsp.second, "ZDIPLEN", 0.16125, 0.01, &d);
          check(
              near != nullptr && std::abs(near->freq1 - 0.1609303480) < 1e-9
                  && std::abs(d - (0.16125 - 0.1609303480)) < 1e-9,
              "nearest_dalton_record: the ZDIPLEN record within seed.freq_tol "
              "of the leg's frequency");
          const auto* far = molresponse_v3::nearest_dalton_record(
              rsp.second, "ZDIPLEN", 0.19169, 0.01);
          check(far != nullptr && std::abs(far->freq1 - 0.1916703618) < 1e-9,
                "nearest_dalton_record: the nearer of two records");
          check(molresponse_v3::nearest_dalton_record(rsp.second, "ZDIPLEN",
                                                      0.05, 0.01)
                        == nullptr
                    && molresponse_v3::nearest_dalton_record(
                           rsp.second, "ZDIPLEN", 0.16125, 1e-4)
                           == nullptr,
                "nearest_dalton_record: none outside seed.freq_tol (the leg "
                "starts unseeded)");
          check(molresponse_v3::nearest_dalton_record(rsp.second, "ZDIPLEN",
                                                      0.3218606960, 0.01)
                    == nullptr,
                "nearest_dalton_record: an excitation record (EXCITLAB) at the "
                "frequency is not a dipole record");
        }
        // The deck keys (B3, DESIGN v3.27).
        const std::string seeded_es
            = "molresponse es\n reference GS\n protocol [1e-4]\n dconv 1e-4\n"
              " requested_properties [excited_states]\n excited_states.roots "
              "1\n"
              " excited_states.method tda\n seed.dalton_dir "
              + dalton.string()
              + "\n"
                "end\n";
        try {
          const auto stages
              = named_stages::parse(world, geometry + gs + seeded_es);
          check(stages[1].seed_dalton_dir == dalton.lexically_normal().string()
                    && stages[1].seed_freq_tol == 0.0
                    && stages[1]
                               .params.get<ResponseParameters>()
                               .get<std::string>("dalton.dir")
                           == stages[1].seed_dalton_dir,
                "seed.dalton_dir: a molresponse stage's, handed to the engine");
          const auto settings = molresponse.settings(stages[1], 2);
          check(
              settings.seed.find("Dalton directory") != std::string::npos
                  && madness::convergence::render_settings_table(settings).find(
                         "  Seed: ")
                         != std::string::npos
                  && madness::convergence::to_json(settings).contains("seed"),
              "settings (S5): the seed shown and recorded");
          check(!madness::convergence::to_json(molresponse.settings(optical, 2))
                     .contains("seed"),
                "settings (S5): an unseeded stage's record is unchanged");
        } catch (const std::exception& e) {
          check(false, std::string("seeded stage: ") + e.what());
        }
        auto with = [&](std::string block_text, const std::string& lines) {
          block_text.replace(block_text.find("end\n"), 0, lines);
          return block_text;
        };
        rejects(geometry + gs + with(seeded_es, " seed.freq_tol 0.01\n"),
                "seed.freq_tol does nothing here");
        rejects(
            geometry + gs
                + with(response, " seed.dalton_dir " + dalton.string() + "\n"),
            "seed.freq_tol is required");
        rejects(geometry + gs
                    + with(response, " seed.dalton_dir " + dalton.string()
                                         + "\n seed.freq_tol -1\n"),
                "seed.freq_tol must be positive");
        rejects(geometry
                    + with(gs, " seed.dalton_dir " + dalton.string() + "\n"),
                "is a molresponse stage's");
        rejects(geometry + gs
                    + with(response, " dalton.dir " + dalton.string() + "\n"),
                "now seed.dalton_dir");
        rejects(geometry + gs + with(response, " seed.es_warmup true\n"),
                "unsupported parameter seed.es_warmup");
      }
      // release (B5): only a converged outcome is a reference.
      named_stages::StageOutcome outcome;
      outcome.converged = true;
      check(moldft.release(outcome) && molresponse.release(outcome),
            "spec release: a converged outcome binds");
      outcome.converged = false;
      check(!moldft.release(outcome) && !molresponse.release(outcome),
            "spec release: an unconverged outcome does not");
    } catch (const std::exception& e) {
      check(false, std::string("spec tables: ") + e.what());
    }
    world.gop.sum(failures);
    if (world.rank() == 0) std::cout << "named-stage failures: " << failures << '\n';
    madness::finalize();
    return failures ? 1 : 0;
}
