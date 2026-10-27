// The named property language (D-017) and the selectors that read the
// engine's rows: pure functions, table-tested on synthetic rows. The parser
// cases go through named_stages::parse; the mapping cases call
// translate_response_properties directly, because parse still refuses
// hyperpolarizability and raman at its end (D-017.2) until the response
// checkpoint carries quadratic sources.
#include "NamedStages.hpp"
#include "ResponseProperties.hpp"
#include <iostream>
#include <stdexcept>

namespace {
using named_stages::Stage;
using nlohmann::json;

// A molresponse stage with one rung, as parse would build it, without a deck.
Stage bare_stage(madness::World& world) {
  named_stages::Block block{"molresponse", "optical", {}};
  Stage stage(world, block);
  named_stages::read_value(world, stage.params.get<CalculationParameters>(),
                           {1, "protocol [1e-4]"});
  return stage;
}
// Read control lines and translate them, as parse does for a molresponse stage.
Stage translated(madness::World& world, const std::vector<std::string>& lines) {
  Stage stage = bare_stage(world);
  named_stages::Controls controls;
  std::size_t number = 1;
  for (const auto& line : lines)
    named_stages::read_value(world, controls, {number++, line});
  stage.properties
      = controls.get<std::vector<std::string>>("requested_properties");
  named_stages::translate_response_properties(
      "optical", controls, stage, stage.params.get<ResponseParameters>());
  return stage;
}
json quadratic_row(const std::string& a, const std::string& b,
                   const std::string& c, double freq_b, double freq_c,
                   double value, bool converged) {
  return {{"A", a},
          {"B", b},
          {"C", c},
          {"freq_b", freq_b},
          {"freq_c", freq_c},
          {"beta", value},
          {"converged", converged}};
}
} // namespace

int main(int argc, char** argv) {
  auto& world = madness::initialize(argc, argv);
  int failures = 0;
  auto check = [&](bool ok, const std::string& why) {
    if (!ok) {
      ++failures;
      if (world.rank() == 0)
        std::cerr << why << '\n';
    }
  };
  auto throws = [&](const std::function<void()>& action,
                    const std::string& message) {
    try {
      action();
      check(false, "accepted: " + message);
    } catch (const std::exception& e) {
      check(std::string(e.what()).find(message) != std::string::npos, e.what());
    }
  };

  // ----- The component tables ---------------------------------------------
  {
    const auto pairs = named_stages::every_component(2);
    const auto triples = named_stages::every_component(3);
    check(pairs.size() == 9 && pairs.front() == "xx" && pairs[1] == "xy"
              && pairs.back() == "zz",
          "every pair, x before y before z in each place");
    check(triples.size() == 27 && triples.front() == "xxx"
              && triples.back() == "zzz",
          "every triple");
    check(named_stages::checked_components("polarizability", {}).size() == 9,
          "omitted polarizability components mean all 9");
    check(named_stages::checked_components("hyperpolarizability", {}).size()
              == 27,
          "omitted hyperpolarizability components mean all 27");
    check(
        named_stages::checked_components("hyperpolarizability", {"zzz", "xzz"})
            == std::vector<std::string>{"zzz", "xzz"},
        "written triples are kept in order");
    throws(
        [&] {
          named_stages::checked_components("hyperpolarizability", {"zz"});
        },
        "triples");
    throws([&] { named_stages::checked_components("raman", {"zzz"}); },
           "pairs");
    throws([&] { named_stages::checked_components("raman", {"zz", "zz"}); },
           "distinct");
    throws([&] { named_stages::checked_components("polarizability", {"zq"}); },
           "Cartesian");
    throws([&] { named_stages::checked_frequencies("raman", {}, "au"); },
           "explicit frequencies");
    throws([&] { named_stages::checked_frequencies("raman", {0.0}, ""); },
           "frequency_units au");
    throws([&] { named_stages::checked_frequencies("raman", {-0.1}, "au"); },
           "nonnegative");
    throws(
        [&] {
          named_stages::checked_frequencies("raman", {0.02, 0.020001}, "au");
        },
        "frequency-key precision");
  }

  // ----- Mapping onto the engine's knobs ----------------------------------
  try {
    const auto beta = translated(
        world,
        {"requested_properties [polarizability, hyperpolarizability]",
         "polarizability.frequencies [0.0, 0.02]",
         "polarizability.frequency_units au", "polarizability.components [zz]",
         "hyperpolarizability.frequencies [0.0, 0.05]",
         "hyperpolarizability.frequency_units au",
         "hyperpolarizability.components [xzz]"});
    const auto& rp = beta.params.get<ResponseParameters>();
    check(beta.frequencies == std::vector<double>{0.0, 0.02}
              && beta.components == std::vector<std::string>{"zz"},
          "the polarizability request stays on stage.frequencies/components");
    check(beta.hyperpolarizability.frequencies == std::vector<double>{0.0, 0.05}
              && beta.hyperpolarizability.components
                     == std::vector<std::string>{"xzz"},
          "the hyperpolarizability request");
    check(rp.dipole_frequencies() == std::vector<double>{0.0, 0.02, 0.05},
          "dipole.frequencies is the union, polarizability first");
    check(rp.dipole_directions() == "xz",
          "dipole.directions is the union of every axis named");
    check(rp.beta_frequencies() == std::vector<double>{0.0, 0.05},
          "beta.frequencies are the driver frequencies");
    check(rp.requested_properties()
              == std::vector<std::string>{"polarizability",
                                          "hyperpolarizability"},
          "requested_properties as written");
    check(!rp.is_user_defined("raman.nuc_atom")
              && !rp.is_user_defined("raman.atoms")
              && !rp.is_user_defined("raman.directions")
              && !rp.is_user_defined("raman.frequencies"),
          "no Raman knobs without raman");

    const auto all = translated(world, {"requested_properties [polarizability]",
                                        "polarizability.frequencies [0.0]",
                                        "polarizability.frequency_units au"});
    check(all.components.size() == 9
              && all.params.get<ResponseParameters>().dipole_directions()
                     == "xyz",
          "omitted polarizability components: all 9, every axis");

    const auto full = translated(world, {"requested_properties [raman]",
                                         "raman.frequencies [0.0]",
                                         "raman.frequency_units au"});
    check(full.raman_full_set && full.raman.components.size() == 9
              && full.params.get<ResponseParameters>().raman_nuc_atom() == -1
              && !full.params.get<ResponseParameters>().is_user_defined(
                  "raman.atoms")
              && !full.params.get<ResponseParameters>().is_user_defined(
                  "raman.directions")
              && full.params.get<ResponseParameters>().raman_frequencies()
                     == std::vector<double>{0.0},
          "raman with atoms and directions omitted: the engine's full set "
          "(raman.nuc_atom -1, no lists) at raman.frequencies");
    const auto one = translated(
        world, {"requested_properties [raman]", "raman.frequencies [0.0]",
                "raman.frequency_units au", "raman.components [zz]",
                "raman.atoms [1]", "raman.directions [z]"});
    check(!one.raman_full_set && one.raman_atoms == std::vector<int>{1}
              && one.raman_directions == "z"
              && one.params.get<ResponseParameters>().raman_atoms()
                     == std::vector<int>{1}
              && one.params.get<ResponseParameters>().raman_directions() == "z"
              && !one.params.get<ResponseParameters>().is_user_defined(
                  "raman.nuc_atom")
              && !one.params.get<ResponseParameters>().is_user_defined(
                  "nuclear.directions")
              && one.params.get<ResponseParameters>().dipole_directions()
                     == "z",
          "one coordinate: raman.atoms and raman.directions as written, "
          "neither raman.nuc_atom nor nuclear.directions");
    const auto every_atom_z = translated(
        world, {"requested_properties [raman]", "raman.frequencies [0.0]",
                "raman.frequency_units au", "raman.directions [x, z]"});
    check(
        every_atom_z.raman_atoms.empty()
            && every_atom_z.raman_directions == "xz"
            && !every_atom_z.params.get<ResponseParameters>().is_user_defined(
                "raman.atoms")
            && !every_atom_z.params.get<ResponseParameters>().is_user_defined(
                "raman.nuc_atom")
            && every_atom_z.params.get<ResponseParameters>().raman_directions()
                   == "xz",
        "every atom times written directions: raman.directions alone");
    const auto atom_xyz = translated(
        world, {"requested_properties [raman]", "raman.frequencies [0.0]",
                "raman.frequency_units au", "raman.atoms [0]"});
    check(!atom_xyz.raman_full_set
              && atom_xyz.params.get<ResponseParameters>().raman_atoms()
                     == std::vector<int>{0}
              && !atom_xyz.params.get<ResponseParameters>().is_user_defined(
                  "raman.directions"),
          "one atom, directions omitted: raman.atoms alone");
    const auto two_atoms = translated(
        world, {"requested_properties [raman]", "raman.frequencies [0.0, 0.05]",
                "raman.frequency_units au", "raman.atoms [1, 0]",
                "raman.directions [z, x]"});
    check(two_atoms.raman_atoms == std::vector<int>{1, 0}
              && two_atoms.raman_directions == "zx"
              && two_atoms.params.get<ResponseParameters>().raman_atoms()
                     == std::vector<int>{1, 0}
              && two_atoms.params.get<ResponseParameters>().raman_directions()
                     == "zx"
              && two_atoms.params.get<ResponseParameters>().raman_frequencies()
                     == std::vector<double>{0.0, 0.05},
          "a list of atoms and directions passes through in deck order");
  } catch (const std::exception& e) {
    check(false, e.what());
  }
  throws(
      [&] {
        translated(world, {"requested_properties [polarizability]",
                           "polarizability.frequencies [0.0]",
                           "polarizability.frequency_units au",
                           "hyperpolarizability.frequencies [0.0]"});
      },
      "hyperpolarizability.frequencies needs hyperpolarizability in "
      "requested_properties");
  throws(
      [&] {
        translated(world,
                   {"requested_properties [polarizability]",
                    "polarizability.frequencies [0.0]",
                    "polarizability.frequency_units au", "raman.atoms [0]"});
      },
      "raman.atoms needs raman");
  throws(
      [&] {
        translated(world, {"requested_properties [hyperpolarizability]"});
      },
      "hyperpolarizability needs explicit frequencies");
  throws([&] { translated(world, {"requested_properties [magnetizability]"}); },
         "unsupported or repeated");
  throws(
      [&] {
        translated(world,
                   {"requested_properties [raman, raman]",
                    "raman.frequencies [0.0]", "raman.frequency_units au"});
      },
      "unsupported or repeated");
  throws([&] { translated(world, {"requested_properties []"}); },
         "at least one");
  throws(
      [&] {
        translated(world,
                   {"requested_properties [raman]", "raman.frequencies [0.0]",
                    "raman.frequency_units au", "raman.atoms [0, 0]"});
      },
      "distinct atom indices");
  throws(
      [&] {
        translated(world,
                   {"requested_properties [raman]", "raman.frequencies [0.0]",
                    "raman.frequency_units au", "raman.atoms [-1]"});
      },
      "numbered from 0");
  throws(
      [&] {
        translated(world,
                   {"requested_properties [raman]", "raman.frequencies [0.0]",
                    "raman.frequency_units au", "raman.directions [q]"});
      },
      "distinct axes");
  throws(
      [&] {
        translated(world,
                   {"requested_properties [raman]", "raman.frequencies [0.0]",
                    "raman.frequency_units au", "raman.directions [z, z]"});
      },
      "distinct axes");

  // ----- Through parse: accepted language, rejected engine spellings -------
  const std::string geometry
      = "molecule LiH\n units angstrom\n Li 0 0 0\n H 0 0 1.5949\nend\n";
  const std::string gs = "moldft gs\n molecule LiH\n xc hf\n protocol [1e-4]\n "
                         "requested_properties [energy]\nend\n";
  auto rejects = [&](const std::string& text, const std::string& message) {
    throws([&] { named_stages::parse(world, text); }, message);
  };
  const std::string beta_deck
      = "molresponse quadratic\n reference gs\n protocol [1e-4]\n"
        " requested_properties [polarizability, hyperpolarizability]\n"
        " polarizability.frequencies [0.0]\n polarizability.frequency_units "
        "au\n"
        " polarizability.components [zz]\n hyperpolarizability.frequencies "
        "[0.0]\n"
        " hyperpolarizability.frequency_units au\n "
        "hyperpolarizability.components [zzz]\nend\n";
  // D-017.2: a valid beta deck passes every check and stops at the last one.
  rejects(geometry + gs + beta_deck, "wait for the response checkpoint");
  rejects(geometry + gs
              + "molresponse r\n reference gs\n protocol [1e-4]\n "
                "requested_properties [raman]\n"
                " raman.frequencies [0.0]\n raman.frequency_units au\nend\n",
          "wait for the response checkpoint");
  try {
    const auto stages = named_stages::parse(
        world, geometry + gs
                   + "molresponse optical\n reference gs\n protocol [1e-4]\n"
                     " requested_properties [polarizability]\n "
                     "polarizability.frequencies [0.0]\n"
                     " polarizability.frequency_units au\nend\n");
    check(stages[1].components.size() == 9,
          "parse: omitted components are all 9 pairs");
  } catch (const std::exception& e) {
    check(false, e.what());
  }
  for (const auto* engine_key :
       {" quadratic true\n", " beta.frequencies [0.0]\n", " raman.nuc_atom 1\n",
        " raman.nuc_axis 2\n", " nuclear.directions z\n",
        " dipole.directions z\n"}) {
    auto deck = beta_deck;
    deck.replace(deck.find(" protocol"), 0, engine_key);
    rejects(geometry + gs + deck, "unsupported parameter");
  }
  rejects(geometry
              + "moldft bad\n molecule LiH\n xc hf\n requested_properties "
                "[energy]\n raman.atoms [0]\nend\n",
          "raman belongs to molresponse");

  // ----- Selectors over synthetic engine rows -----------------------------
  try {
    Stage stage = bare_stage(world);
    const auto pkey = named_stages::final_protocol_key(stage);
    stage.properties = {"polarizability", "hyperpolarizability"};
    stage.frequencies = {0.0};
    stage.components = {"zz"};
    stage.hyperpolarizability.frequencies = {0.0};
    stage.hyperpolarizability.components = {"zzz"};
    json properties;
    properties["alpha"][pkey] = json::array({{{"omega", 0.0},
                                              {"directions", "z"},
                                              {"alpha", {{22.15}}},
                                              {"converged", true}}});
    properties["beta"][pkey] = json::array(
        {quadratic_row("z", "dipole_z", "dipole_z", 0.0, 0.0, 315.73, true),
         quadratic_row("x", "dipole_z", "dipole_z", 0.0, 0.0, 1.5, false)});
    const auto beta = named_stages::selected_beta(stage, properties, json());
    check(beta.size() == 1 && beta[0].at("components").at("zzz") == 315.73
              && beta[0].at("converged") == true
              && beta[0].at("process") == "static"
              && beta[0].at("engine_evidence").contains("zzz"),
          "static beta_zzz selected with its row as evidence");
    auto selected = named_stages::select_response_properties(
        stage, properties, json(), json(), json(), 2);
    check(selected.contains("polarizability")
              && selected.contains("hyperpolarizability")
              && selected.at("polarizability")
                         .at("values")[0]
                         .at("components")
                         .at("zz")
                     == 22.15
              && named_stages::all_requested_converged(selected),
          "both properties selected, stage converged");
    stage.hyperpolarizability.components = {"zzz", "xzz"};
    selected = named_stages::select_response_properties(
        stage, properties, json(), json(), json(), 2);
    check(!named_stages::all_requested_converged(selected),
          "one unconverged row: the stage is not");
    // D-024: an element the engine wrote no row for is published without a
    // value, with the engine's reason from run_summary.dropped_work, and the
    // record is unconverged.
    json metadata;
    metadata["run_summary"]["dropped_work"] = json::array(
        {{{"id", "vbc:dipole_y__dipole_y@f0.00000_f0.00000"},
          {"top_protocol_key", pkey},
          {"reason",
           "prerequisites never converged (gated out of every wave)"}},
         {{"id", "vbc:dipole_z__dipole_z@f0.00000_f0.00000"},
          {"top_protocol_key", pkey},
          {"reason", "stalled (quarantined by the no-progress guard)"}}});
    stage.hyperpolarizability.components = {"zzz", "yyy"};
    auto partial = named_stages::selected_beta(stage, properties, metadata);
    const auto reason_of = [](const json& record, const char* name) {
      return record.at("components_unavailable").at(name).get<std::string>();
    };
    check(partial.size() == 1 && partial[0].at("components").at("zzz") == 315.73
              && partial[0].at("components").at("yyy").is_null()
              && partial[0].at("converged") == false
              && reason_of(partial[0], "yyy").find("vbc:dipole_y__dipole_y")
                     != std::string::npos
              && reason_of(partial[0], "yyy").find("prerequisites never")
                     != std::string::npos
              && !partial[0].at("engine_evidence").contains("yyy"),
          "an element without a row: null value, the engine's reason, "
          "unconverged");
    stage.hyperpolarizability.frequencies = {0.02};
    partial = named_stages::selected_beta(stage, properties, metadata);
    check(partial[0].at("components").at("zzz").is_null()
              && reason_of(partial[0], "zzz").find("does not name")
                     != std::string::npos,
          "a frequency with no rows and no audit entry: unavailable, the "
          "audit's silence recorded");
    stage.hyperpolarizability.frequencies = {0.0};
    stage.hyperpolarizability.components = {"zzz"};
    stage.properties = {"hyperpolarizability"};
    selected = named_stages::select_response_properties(
        stage, json::object(), metadata, json(), json(), 2);
    check(selected.at("hyperpolarizability")
                  .at("values")[0]
                  .at("components")
                  .at("zzz")
                  .is_null()
              && reason_of(selected.at("hyperpolarizability").at("values")[0],
                           "zzz")
                         .find("stalled")
                     != std::string::npos
              && !named_stages::all_requested_converged(selected),
          "the engine assembled no beta at all: a record without a value, "
          "the stage unconverged");
    stage.properties = {"polarizability", "hyperpolarizability"};
    properties["beta"][pkey].push_back(
        quadratic_row("z", "dipole_z", "dipole_z", 0.0, 0.0, 315.73, true));
    throws([&] { named_stages::selected_beta(stage, properties, json()); },
           "duplicate hyperpolarizability row");
    properties["beta"][pkey] = json::array(
        {quadratic_row("z", "dipole_z", "dipole_z", 0.0, 0.0,
                       std::numeric_limits<double>::quiet_NaN(), true)});
    throws([&] { named_stages::selected_beta(stage, properties, json()); },
           "nonfinite hyperpolarizability");

    // Raman: a written subset, then the full set with its completeness record.
    stage.properties = {"raman"};
    stage.raman.frequencies = {0.0};
    stage.raman.components = {"zz"};
    json raman_rows = json::array();
    for (int atom : {0, 1})
      for (char axis : std::string("xyz"))
        raman_rows.push_back(quadratic_row(
            "z", "dipole_z", "nuc_" + std::to_string(atom) + "_" + axis, 0.0,
            0.0, 0.5 * atom + 0.1 * (axis - 'x'), true));
    properties["raman"][pkey] = raman_rows;
    stage.raman_atoms = {1};
    stage.raman_directions = "z";
    stage.raman_full_set = false;
    auto raman = named_stages::selected_raman(stage, properties, json(), json(),
                                              json(), 2);
    check(raman.at("values").size() == 1
              && raman.at("values")[0].at("coordinate") == "nuc_1_z"
              && raman.at("values")[0].at("atom") == 1
              && raman.at("values")[0].at("axis") == "z"
              && raman.at("values")[0].at("components").at("zz") == 0.7
              && raman.at("coordinates") == "written"
              && raman.at("normal_mode_analysis") == false
              && raman.contains("normal_mode_analysis_unavailable"),
          "one written coordinate: its row, no normal modes, and why");
    stage.raman_atoms = {};
    raman = named_stages::selected_raman(stage, properties, json(), json(),
                                         json(), 2);
    check(raman.at("values").size() == 2
              && raman.at("values")[0].at("coordinate") == "nuc_0_z"
              && raman.at("values")[1].at("coordinate") == "nuc_1_z",
          "atoms omitted: every atom of the reference along the written "
          "direction");
    stage.raman_atoms = {3};
    throws(
        [&] {
          named_stages::selected_raman(stage, properties, json(), json(),
                                       json(), 2);
        },
        "reference has 2 atom(s)");
    stage.raman_atoms = {};
    stage.raman_directions = "";
    stage.raman_full_set = true;
    const json vib = {{"n_atoms", 2}};
    json spectra = {{"completeness",
                     {{"f0.00000",
                       {{"rows_found", 54},
                        {"rows_expected", 54},
                        {"all_rows_converged", true}}}}}};
    raman = named_stages::selected_raman(stage, properties, json(), vib,
                                         spectra, 2);
    check(raman.at("values").size() == 6 && raman.at("coordinates") == "all"
              && raman.at("normal_mode_analysis") == true
              && raman.at("spectrum_complete") == true
              && raman.at("vibrational_analysis") == vib
              && raman.at("raman_spectra") == spectra,
          "the full set: 3N rows, the analysis passed through, complete");
    selected = named_stages::select_response_properties(
        stage, properties, json(), vib, spectra, 2);
    check(named_stages::all_requested_converged(selected),
          "full set complete and converged");
    spectra["completeness"]["f0.00000"]["all_rows_converged"] = false;
    selected = named_stages::select_response_properties(
        stage, properties, json(), vib, spectra, 2);
    check(selected.at("raman").at("spectrum_complete") == false
              && !named_stages::all_requested_converged(selected),
          "an incomplete spectrum makes the stage unconverged");
    throws(
        [&] {
          named_stages::selected_raman(stage, properties, json(), vib, json(),
                                       2);
        },
        "no completeness record");
    // D-024 for a written coordinate: no row at this frequency.
    stage.raman.frequencies = {0.02};
    stage.raman_atoms = {1};
    stage.raman_directions = "z";
    stage.raman_full_set = false;
    json raman_metadata;
    raman_metadata["run_summary"]["dropped_work"] = json::array(
        {{{"id", "vbc:dipole_z__nuc_1_z@f0.02000_f0.00000"},
          {"top_protocol_key", pkey},
          {"reason",
           "prerequisites never converged (gated out of every wave)"}}});
    raman = named_stages::selected_raman(stage, properties, raman_metadata,
                                         json(), json(), 2);
    check(raman.at("values").size() == 1
              && raman.at("values")[0].at("coordinate") == "nuc_1_z"
              && raman.at("values")[0].at("components").at("zz").is_null()
              && raman.at("values")[0].at("converged") == false
              && raman.at("values")[0]
                         .at("components_unavailable")
                         .at("zz")
                         .get<std::string>()
                         .find("vbc:dipole_z__nuc_1_z@f0.02000_f0.00000")
                     != std::string::npos,
          "a Raman element without a row: null value, the undelivered "
          "source named, unconverged");
    selected = named_stages::select_response_properties(
        stage, properties, raman_metadata, json(), json(), 2);
    check(!named_stages::all_requested_converged(selected),
          "the stage is unconverged with it");
  } catch (const std::exception& e) {
    check(false, e.what());
  }

  world.gop.sum(failures);
  if (world.rank() == 0)
    std::cout << "response-property failures: " << failures << '\n';
  madness::finalize();
  return failures ? 1 : 0;
}
