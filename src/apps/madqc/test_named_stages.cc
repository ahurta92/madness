#include "MoldftCheckpoint.hpp"
#include "NamedStages.hpp"
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
    rejects(geometry + gs + "molresponse bad\n reference absent\nend\n", "reference");
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
    // D8: only a user-written econv is a request; it is judged on saved
    // evidence and, in engine mode, reported but not enforced.
    {
      using nlohmann::json;
      const auto saved = [](json value, const char* availability) {
        return json{{"energy_change",
                     {{"value", value},
                      {"availability", availability},
                      {"units", "hartree"}}}};
      };
      const auto status = [](const json& report) {
        return report.at("status").get<std::string>();
      };
      CalculationParameters defaults;
      check(!named_stages::requested_econv(defaults),
            "engine default econv is not a request");
      check(named_stages::energy_change_target(defaults).empty(),
            "no econv, no energy target");
      check(status(named_stages::energy_change_request_report(json::object(),
                                                              defaults))
                == "not_requested",
            "no request needs no evidence");
      CalculationParameters requested;
      requested.set_user_defined_value("econv", 1e-6);
      check(named_stages::requested_econv(requested) == 1e-6,
            "user econv is the request");
      const auto met = named_stages::energy_change_request_report(
          saved(5e-7, "measured"), requested);
      check(status(met) == "satisfied" && !met.contains("outcome"),
            "met by evidence");
      check(status(named_stages::energy_change_request_report(
                saved(1e-6, "measured"), requested))
                == "satisfied",
            "met at equality (<=)");
      const auto unmet = named_stages::energy_change_request_report(
          saved(2e-6, "measured"), requested);
      check(status(unmet) == "unmet" && unmet.at("enforced") == false
                && unmet.at("outcome")
                       == "unmet, not enforceable in engine mode",
            "unmet evidence is reported, not enforced");
      check(status(named_stages::energy_change_request_report(
                saved(nullptr, "unavailable"), requested))
                == "unavailable",
            "unavailable evidence does not meet a request");
      check(status(named_stages::energy_change_request_report(json::object(),
                                                              requested))
                == "unavailable",
            "missing evidence does not meet a request");
    }
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
    world.gop.sum(failures);
    if (world.rank() == 0) std::cout << "named-stage failures: " << failures << '\n';
    madness::finalize();
    return failures ? 1 : 0;
}
