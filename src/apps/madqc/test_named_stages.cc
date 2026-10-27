#include "MoldftCheckpoint.hpp"
#include "NamedStages.hpp"
#include "ResponseCheckpoint.hpp"
#include <iostream>
#include <stdexcept>

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
      const auto met = named_stages::requested_criteria_report(
          saved(5e-7, "measured"), engine_mode);
      check(met.at("status") == "satisfied" && met.at("enforced") == false,
            "met by evidence, reported in engine mode");
      const auto unmet = named_stages::requested_criteria_report(
          saved(2e-6, "measured"), engine_mode);
      check(unmet.at("status") == "unmet"
                && unmet.at("criteria").at("energy_change").at("outcome")
                       == "unmet, not enforceable in engine mode",
            "unmet evidence is reported, not enforced, in engine mode");
      check(named_stages::requested_criteria_report(json::object(), engine_mode)
                    .at("status")
                == "unavailable",
            "missing evidence does not meet a request");
      requested.set_user_defined_value(
          std::string(cv::use_default_criteria_key), false);
      const auto criteria_mode
          = named_stages::moldft_convergence_settings("gs", requested, 3);
      check(criteria_mode.decider == cv::Decider::requested_criteria
                && named_stages::requested_criteria_report(
                       saved(2e-6, "measured"), criteria_mode)
                           .at("enforced")
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
    for (const auto* key :
         {" checkpoint.every_iterations 2\n", " checkpoint.keep_all true\n"}) {
        auto response_cadence = response;
        response_cadence.replace(response_cadence.find(" dconv 2e-4\n"), 0,
                                 key);
        rejects(geometry + gs + response_cadence,
                "publishes at solve ends only");
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
    world.gop.sum(failures);
    if (world.rank() == 0) std::cout << "named-stage failures: " << failures << '\n';
    madness::finalize();
    return failures ? 1 : 0;
}
