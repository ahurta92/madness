// The named property language (D-017) and the selectors that read the
// engine's rows: pure functions, table-tested on synthetic rows. The parser
// cases go through named_stages::parse; the mapping cases call
// translate_response_properties directly.
#include "NamedStages.hpp"
#include "ResponseProperties.hpp"
#include <iostream>
#include <madness/chem/ExplicitCriteria.h>
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

// The selected records' own flags: every requested row converged and a
// Raman spectrum complete. The stage's verdict read this until S4d2; it now
// comes from the state records (shared_verdict, D-077).
static bool every_row_converged(const nlohmann::json& selected) {

  for (const auto& [name, record] : selected.items()) {
    if (record.contains("values")) {
      for (const auto& row : record.at("values"))
        if (!row.at("converged").get<bool>())
          return false;
    } else if (!record.value("converged", false)) {
      return false;
    }
    if (record.contains("spectrum_complete")
        && !record.at("spectrum_complete").get<bool>())
      return false;
  }
  return true;
}

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

    // D-025: the excited-state request onto the engine's excited.* knobs.
    const auto excited = translated(
        world, {"requested_properties [excited_states, polarizability]",
                "excited_states.roots 2", "excited_states.method rpa",
                "polarizability.frequencies [0.0]",
                "polarizability.frequency_units au",
                "polarizability.components [zz]"});
    const auto& erp = excited.params.get<ResponseParameters>();
    check(excited.excited_states.roots == 2
              && excited.excited_states.method == "rpa"
              && excited.excited_states.guess.empty() && erp.excited_enable()
              && erp.excited_num_states() == 2 && !erp.excited_tda()
              && !erp.is_user_defined("excited.guess")
              && !erp.is_user_defined("excited.guess_basis")
              && !erp.is_user_defined("excited.maxiter")
              && !erp.is_user_defined("excited.guess_max_iter")
              && !erp.is_user_defined("excited.maxsub")
              && erp.dipole_frequencies() == std::vector<double>{0.0}
              && erp.dipole_directions() == "z",
          "excited_states: roots and method onto excited.enable, "
          "excited.num_states and excited.tda; no cap of its own; the "
          "dipole legs beside it");
    const auto tda_only = translated(
        world,
        {"requested_properties [excited_states]", "excited_states.roots 1",
         "excited_states.method tda", "excited_states.guess virtual_ao",
         "excited_states.guess_basis aug-cc-pvtz"});
    const auto& trp = tda_only.params.get<ResponseParameters>();
    check(trp.excited_enable() && trp.excited_tda()
              && trp.excited_num_states() == 1
              && trp.excited_guess() == "virtual_ao"
              && trp.excited_guess_basis() == "aug-cc-pvtz"
              && trp.dipole_frequencies().empty()
              && tda_only.frequencies.empty(),
          "excited_states alone with a written guess: excited.tda, the "
          "guess knobs, no dipole grid");
    check(erp.is_user_defined("excited.enable")
              && !all.params.get<ResponseParameters>().is_user_defined(
                  "excited.enable"),
          "no excited.* knobs without excited_states");
  } catch (const std::exception& e) {
    check(false, e.what());
  }
  throws(
      [&] {
        translated(world, {"requested_properties [excited_states]",
                           "excited_states.roots 1"});
      },
      "excited_states.method (tda or rpa); neither has a default");
  throws(
      [&] {
        translated(world, {"requested_properties [excited_states]",
                           "excited_states.method tda"});
      },
      "neither has a default");
  throws(
      [&] {
        translated(world,
                   {"requested_properties [excited_states]",
                    "excited_states.roots 0", "excited_states.method tda"});
      },
      "at least 1");
  throws(
      [&] {
        translated(world,
                   {"requested_properties [excited_states]",
                    "excited_states.roots 1", "excited_states.method tda",
                    "excited_states.guess_basis aug-cc-pvtz"});
      },
      "guess_basis applies to the virtual_ao guess only");
  throws(
      [&] {
        translated(world,
                   {"requested_properties [excited_states]",
                    "excited_states.roots 1", "excited_states.method cis"});
      },
      "excited_states.method");
  throws(
      [&] {
        translated(world, {"requested_properties [polarizability]",
                           "polarizability.frequencies [0.0]",
                           "polarizability.frequency_units au",
                           "excited_states.roots 1"});
      },
      "excited_states.roots needs excited_states in requested_properties");
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
  // The response checkpoint carries quadratic sources (D-019): valid beta
  // and Raman decks parse.
  try {
    const auto beta = named_stages::parse(world, geometry + gs + beta_deck);
    check(beta[1].properties
              == std::vector<std::string>{"polarizability",
                                          "hyperpolarizability"},
          "parse: a beta deck is accepted");
    const auto raman = named_stages::parse(
        world, geometry + gs
                   + "molresponse r\n reference gs\n protocol [1e-4]\n "
                     "requested_properties [raman]\n"
                     " raman.frequencies [0.0]\n raman.frequency_units "
                     "au\nend\n");
    check(raman[1].raman_full_set, "parse: a Raman deck is accepted");
  } catch (const std::exception& e) {
    check(false, e.what());
  }
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
  const std::string es_deck
      = "molresponse excited\n reference gs\n protocol [1e-4]\n"
        " requested_properties [excited_states]\n excited_states.roots 1\n"
        " excited_states.method tda\nend\n";
  // D-025: the response checkpoint takes the bundle, so a valid
  // excited-state deck parses.
  try {
    const auto es = named_stages::parse(world, geometry + gs + es_deck);
    check(es[1].properties == std::vector<std::string>{"excited_states"}
              && es[1].excited_states.roots == 1
              && es[1].excited_states.method == "tda",
          "parse: an excited-state deck is accepted");
  } catch (const std::exception& e) {
    check(false, std::string("parse: an excited-state deck: ") + e.what());
  }
  {
    auto deck = es_deck;
    deck.replace(deck.find(" protocol"), 0, " excited.enable true\n");
    rejects(geometry + gs + deck, "unsupported parameter");
  }
  rejects(geometry
              + "moldft bad\n molecule LiH\n xc hf\n requested_properties "
                "[energy]\n excited_states.roots 1\nend\n",
          "excited_states belongs to molresponse");
  // D-061: econv on an excited-state stage is the target of the excitation
  // energy change per root (D8: written only); on any other molresponse
  // stage it is refused, a driven response having no local energy.
  try {
    auto last = es_deck;
    last.insert(last.rfind("end\n"), " econv 1e-6\n");
    const auto es = named_stages::parse(world, geometry + gs + last);
    check(es[1].excited_states.omega_target == 1e-6
              && es[1].written.contains("econv"),
          "parse: econv on an excited-state stage is its omega target");
    auto first = es_deck;
    first.insert(first.find(" protocol"), " econv 1e-6\n");
    check(named_stages::parse(world, geometry + gs + first)[1]
                  .excited_states.omega_target
              == 1e-6,
          "parse: econv written before requested_properties");
    const auto plain = named_stages::parse(world, geometry + gs + es_deck);
    check(!plain[1].excited_states.omega_target
              && !plain[1].written.contains("econv"),
          "parse: without econv there is no omega target");
  } catch (const std::exception& e) {
    check(false,
          std::string("parse: econv on an excited-state stage: ") + e.what());
  }
  rejects(geometry + gs
              + "molresponse optical\n reference gs\n protocol [1e-4]\n"
                " requested_properties [polarizability]\n"
                " polarizability.frequencies [0.0]\n"
                " polarizability.frequency_units au\n"
                " polarizability.components [zz]\n econv 1e-6\nend\n",
          "econv is not applicable: a driven response has no local energy");
  // D-092: criteria mode on an excited-state stage. A stage that requests
  // only excited_states needs econv; a driven-response property needs one
  // of the FD targets; a stage with both needs the FD target, econv optional.
  {
    auto only = es_deck;
    only.insert(only.rfind("end\n"),
                " econv 1e-6\n convergence.use_default_criteria false\n");
    try {
      const auto st = named_stages::parse(world, geometry + gs + only);
      check(st[1].excited_states.omega_target == 1e-6,
            "parse: criteria mode on an excited-state-only stage with econv");
    } catch (const std::exception& e) {
      check(false, std::string("parse: criteria mode with econv: ") + e.what());
    }
    auto bare = es_deck;
    bare.insert(bare.rfind("end\n"),
                " convergence.use_default_criteria false\n");
    rejects(geometry + gs + bare, "requests only excited_states needs econv");
    const std::string mixed
        = "molresponse excited\n reference gs\n protocol [1e-4]\n"
          " requested_properties [excited_states, polarizability]\n"
          " excited_states.roots 1\n excited_states.method tda\n"
          " polarizability.frequencies [0.0]\n"
          " polarizability.frequency_units au\n"
          " polarizability.components [zz]\n"
          " convergence.use_default_criteria false\n";
    rejects(geometry + gs + mixed + " econv 1e-6\nend\n",
            "written target for the driven-response legs");
    try {
      const auto st = named_stages::parse(world, geometry + gs + mixed
                                                     + " dconv 1e-4\nend\n");
      check(!st[1].excited_states.omega_target && st[1].properties.size() == 2,
            "parse: criteria mode on a mixed stage with dconv alone");
    } catch (const std::exception& e) {
      check(false, std::string("parse: criteria mode, mixed: ") + e.what());
    }
  }
  // D-099 (B3): stall.window and stall.ratio, the bundle's plateau detector.
  {
    auto stall = es_deck;
    stall.insert(stall.rfind("end\n"), " stall.window 10\n stall.ratio 0.2\n");
    try {
      const auto st = named_stages::parse(world, geometry + gs + stall);
      const auto rows = named_stages::bundle_stall_rows(st[1]);
      check(rows.size() == 2 && rows[0].keyword == "stall.window"
                && rows[0].target == 10.0 && rows[0].effective_limit == 10.0
                && rows[0].units == "iterations"
                && rows[0].set_by == madness::convergence::SetBy::user
                && rows[0].group == madness::convergence::Group::engine_check
                && rows[1].keyword == "stall.ratio" && rows[1].target == 0.2
                && rows[1].set_by == madness::convergence::SetBy::user,
            "stall keys written: two engine rows, set by the user");
      const auto plain = named_stages::parse(world, geometry + gs + es_deck);
      const auto defaults = named_stages::bundle_stall_rows(plain[1]);
      check(defaults[0].target == 6.0
                && defaults[0].set_by
                       == madness::convergence::SetBy::engine_default
                && defaults[1].target == 0.1
                && defaults[1].set_by
                       == madness::convergence::SetBy::engine_default,
            "stall keys unwritten: the engine's defaults, 6 and 0.1");
      auto off = es_deck;
      off.insert(off.rfind("end\n"), " stall.window 0\n");
      const auto disabled = named_stages::bundle_stall_rows(
          named_stages::parse(world, geometry + gs + off)[1]);
      check(disabled[0].target == 0.0
                && disabled[0].note.find("detector off") != std::string::npos,
            "stall.window 0: the detector is off, the row says so");
    } catch (const std::exception& e) {
      check(false, std::string("parse: stall keys: ") + e.what());
    }
    rejects(geometry + gs
                + "molresponse optical\n reference gs\n protocol [1e-4]\n"
                  " requested_properties [polarizability]\n"
                  " polarizability.frequencies [0.0]\n"
                  " polarizability.frequency_units au\n"
                  " polarizability.components [zz]\n stall.window 10\nend\n",
            "drives the excited-state bundle's plateau detector");
  }
  {
    auto zero = es_deck;
    zero.insert(zero.rfind("end\n"), " econv 0\n");
    rejects(geometry + gs + zero, "positive finite econv");
  }

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
    check(beta[0].at("n_b") == 0 && beta[0].at("n_c") == 0,
          "beta row (D-116): n_b and n_c, its position in the written list");
    auto selected = named_stages::select_response_properties(
        stage, properties, json(), json(), json(), 2);
    check(selected.contains("polarizability")
              && selected.contains("hyperpolarizability")
              && selected.at("polarizability")
                         .at("values")[0]
                         .at("components")
                         .at("zz")
                     == 22.15
              && every_row_converged(selected)
              && selected.at("polarizability").at("values")[0].at("omega_index") == 0,
          "both properties selected, stage converged; the alpha row's "
          "omega_index (D-116)");
    stage.hyperpolarizability.components = {"zzz", "xzz"};
    selected = named_stages::select_response_properties(
        stage, properties, json(), json(), json(), 2);
    check(!every_row_converged(selected),
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
              && !every_row_converged(selected),
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
              && raman.at("values")[0].at("units")
                     == "a.u. (e^2 bohr^2 / Eh per bohr)"
              && raman.at("coordinates") == "written"
              && raman.at("normal_mode_analysis") == false
              && raman.contains("normal_mode_analysis_unavailable")
              && raman.at("values")[0].at("omega_index") == 0,
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
    check(every_row_converged(selected), "full set complete and converged");
    spectra["completeness"]["f0.00000"]["all_rows_converged"] = false;
    selected = named_stages::select_response_properties(
        stage, properties, json(), vib, spectra, 2);
    check(selected.at("raman").at("spectrum_complete") == false
              && !every_row_converged(selected),
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
    check(!every_row_converged(selected), "the stage is unconverged with it");
  } catch (const std::exception& e) {
    check(false, e.what());
  }

  // ----- The excited-state selector over a synthetic bundle (D-025) -------
  try {
    Stage stage = bare_stage(world);
    const auto pkey = named_stages::final_protocol_key(stage);
    stage.properties = {"excited_states"};
    stage.excited_states.roots = 2;
    stage.excited_states.method = "tda";
    auto root = [](int slot, int stable, double omega, bool converged,
                   bool stalled, bool locked) {
      return json{{"slot", slot},
                  {"stable_index", stable},
                  {"root_id", "es_root_000" + std::to_string(stable)},
                  {"display_name", "root " + std::to_string(slot)},
                  {"omega", omega},
                  {"bsh_residual", 1e-5},
                  {"density_residual", 2e-6},
                  {"omega_residual", 3e-7},
                  {"density_target", 5e-4},
                  {"omega_target", 5e-4},
                  {"converged", converged},
                  {"stalled", stalled},
                  {"locked", locked},
                  {"file", "root_0.bin"},
                  {"coeffs", 12},
                  {"bytes", 96}};
    };
    json metadata;
    metadata["excited_states"][pkey]
        = {{"type", "tda"},
           {"shell", "closed_shell"},
           {"n_roots", 2},
           {"bundle_dir", "es__" + pkey},
           {"converged", true},
           {"diverged", false},
           {"stalled", false},
           {"stop_reason", "converged"},
           {"iter", 7},
           {"slot_permutation", {1, 0}},
           {"kinetic_history", json::array({1.0, 2.0})},
           {"roots", json::array({root(0, 1, 0.42, true, false, true),
                                  root(1, 0, 0.47, true, false, false)})}};
    json transition
        = {{"analysis", "tda"},
           {"protocol_key", pkey},
           {"units", {{"transition_dipole", "a.u. (e a0)"}}},
           {"roots", json::array({{{"slot", 1},
                                   {"root_id", "es_root_0000"},
                                   {"oscillator_strength", 0.3}},
                                  {{"slot", 0},
                                   {"root_id", "es_root_0001"},
                                   {"oscillator_strength", 0.1}}})}};
    auto es
        = named_stages::selected_excited_states(stage, metadata, transition);
    check(es.at("method") == "tda" && es.at("roots_requested") == 2
              && es.at("guess") == "solid_harmonics"
              && es.at("guess_written") == false
              && es.at("protocol_step_key") == pkey
              && es.at("stop_reason") == "converged"
              && es.at("iterations_at_protocol_step") == 7
              && es.at("decided_by") == "engine_check"
              && es.at("converged") == true && es.at("roots").size() == 2
              && es.at("transition_units").at("transition_dipole")
                     == "a.u. (e a0)",
          "a converged bundle: the request, the exit, converged, the "
          "engine's transition units once (D-104)");
    const auto& first = es.at("roots")[0];
    check(first.at("root_id") == "es_root_0001" && first.at("slot") == 0
              && first.at("omega") == 0.42 && first.at("locked") == true
              && first.at("omega_units") == "au"
              && first.at("decided_by") == "engine_check"
              && !first.contains("file") && !first.contains("coeffs")
              && first.at("transition").at("oscillator_strength") == 0.1
              && es.at("roots")[1].at("transition").at("oscillator_strength")
                     == 0.3,
          "each root keeps the engine's fields, loses the archive "
          "bookkeeping, and takes its transition analysis by root_id");
    check(!es.at("engine_evidence").contains("roots")
              && !es.at("engine_evidence").contains("kinetic_history")
              && es.at("engine_evidence").at("bundle_dir") == "es__" + pkey,
          "the bundle entry travels as evidence without the roots");
    // D-061: the written econv beside the engine's gate, in the record: the
    // request and the evidence, never a status (that is the shared verdict's,
    // convergence_report.requested_accuracy under the bundle's state key).
    const auto& none = es.at("requested_omega_target");
    check(none.at("target").is_null() && none.at("set_by") == "not_set"
              && none.at("units") == "hartree"
              && none.at("relation") == "less_than_or_equal"
              && none.at("worst").at("value") == 3e-7
              && none.at("worst").at("root_id") == "es_root_0001"
              && none.at("verdict_in").at("path")
                     == "convergence_report.requested_accuracy"
              && none.at("verdict_in").at("state_key") == "es:tda_n2|" + pkey
              && none.at("verdict_in").at("criterion")
                     == "excitation_energy_change"
              && !none.contains("status") && !none.contains("enforced")
              && !none.contains("outcome")
              && !first.contains("excitation_energy_change_within_requested"),
          "no econv: target null, not_set, the worst root still named, no "
          "per-root field, no status of its own");
    stage.excited_states.omega_target = 1e-6;
    es = named_stages::selected_excited_states(stage, metadata, transition);
    check(
        es.at("requested_omega_target").at("target") == 1e-6
            && es.at("requested_omega_target").at("set_by") == "user"
            && es.at("roots")[0].at("excitation_energy_change_within_requested")
                   == true
            && es.at("roots")[1].at("excitation_energy_change_within_requested")
                   == true
            && es.at("requested_omega_target").at("worst").at("value") == 3e-7,
        "econv 1e-6 over 3e-7 on both roots: within on both");
    stage.excited_states.omega_target = 1e-7;
    es = named_stages::selected_excited_states(stage, metadata, transition);
    check(
        es.at("requested_omega_target").at("target") == 1e-7
            && es.at("roots")[0].at("excitation_energy_change_within_requested")
                   == false
            && es.at("roots")[1].at("excitation_energy_change_within_requested")
                   == false
            && es.at("requested_omega_target").at("worst").at("root_id")
                   == "es_root_0001"
            && es.at("converged") == true,
        "econv 1e-7 under 3e-7: not within on either root; the bundle's "
        "converged is the engine's and stays true");
    metadata["excited_states"][pkey]["roots"][1]["omega_residual"] = nullptr;
    es = named_stages::selected_excited_states(stage, metadata, transition);
    check(
        es.at("requested_omega_target").at("worst").is_null()
            && es.at("roots")[0].at("excitation_energy_change_within_requested")
                   == false
            && es.at("roots")[1]
                   .at("excitation_energy_change_within_requested")
                   .is_null(),
        "a root without an eigenvalue step: worst null, that root null");
    const auto absent
        = named_stages::selected_excited_states(stage, json::object(), json());
    check(absent.at("requested_omega_target").at("target") == 1e-7
              && absent.at("requested_omega_target").at("worst").is_null()
              && absent.at("requested_omega_target")
                         .at("verdict_in")
                         .at("state_key")
                     == "es:tda_n2|" + pkey,
          "no bundle at the final rung: the request is still described, "
          "worst null");
    stage.excited_states.omega_target.reset();
    metadata["excited_states"][pkey]["roots"][1]
        = root(1, 0, 0.47, true, false, false);
    es = named_stages::selected_excited_states(stage, metadata, transition);
    auto selected = named_stages::select_response_properties(
        stage, json::object(), metadata, json(), json(), 2, transition);
    check(selected.contains("excited_states") && every_row_converged(selected),
          "the stage is converged with the bundle");
    // A stalled root: the bundle ended converged by the engine's flag, but
    // a stalled root is never relabelled converged (B5).
    metadata["excited_states"][pkey]["roots"][1]
        = root(1, 0, 0.47, false, true, false);
    es = named_stages::selected_excited_states(stage, metadata, transition);
    check(es.at("converged") == false
              && es.at("roots")[1].at("stalled") == true,
          "one stalled root: the record is unconverged");
    // An iteration-limit exit: unconverged, the transition analysis absent.
    metadata["excited_states"][pkey]["converged"] = false;
    metadata["excited_states"][pkey]["stop_reason"] = "iteration_limit";
    es = named_stages::selected_excited_states(stage, metadata, json());
    check(es.at("converged") == false && !es.contains("transition_units")
              && es.at("stop_reason") == "iteration_limit"
              && es.at("roots")[0]
                         .at("transition_unavailable")
                         .get<std::string>()
                         .find("bundle not converged")
                     != std::string::npos,
          "an iteration-limit exit: unconverged, transition unavailable with "
          "the reason");
    selected = named_stages::select_response_properties(
        stage, json::object(), metadata, json(), json(), 2);
    check(!every_row_converged(selected), "the stage is unconverged with it");
    // The planner's verdict on a met load travels through (D-016).
    metadata["excited_states"][pkey]["decided_by"] = "planner";
    es = named_stages::selected_excited_states(stage, metadata, json());
    check(es.at("decided_by") == "planner"
              && es.at("roots")[0].at("decided_by") == "planner",
          "decided_by from the bundle entry, inherited by its roots");
    // No bundle at the final rung: a record without roots, with the reason.
    es = named_stages::selected_excited_states(stage, json::object(), json());
    check(es.at("converged") == false && es.at("roots").empty()
              && es.at("unavailable").get<std::string>().find(pkey)
                     != std::string::npos,
          "no bundle at the final rung: unavailable, unconverged");
    // Another bundle at the rung than the one requested is an error.
    stage.excited_states.method = "rpa";
    throws(
        [&] { named_stages::selected_excited_states(stage, metadata, json()); },
        "not the requested rpa with ");
    stage.excited_states.method = "tda";
    stage.excited_states.roots = 3;
    throws(
        [&] { named_stages::selected_excited_states(stage, metadata, json()); },
        "not the requested tda with 3 roots");
    // D-113: two_photon's record from the engine's tpa rows: per root the
    // tensor and strengths, its own convergence in the bundle, its legs; no
    // cross sections in GM.
    stage.excited_states.roots = 2;
    auto tpa_row = [](int root, double omega, bool converged) {
      return json{{"es_root_id", root},
                  {"omega", omega},
                  {"S", {{1.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 3.0}}},
                  {"Df", 1.2}, {"Dg", 0.4667}, {"D_linear", 4.2668},
                  {"D_circular", 0.4002}, {"R", 0.1},
                  {"sigma_linear_gm", 9.0}, {"sigma_circular_gm", 8.0},
                  {"row_accuracy", json::array({{{"converged", converged}},
                                                {{"converged", true}},
                                                {{"converged", true}}})},
                  {"converged", converged},
                  {"max_bsh_residual", 1e-5}};
    };
    json tpa_metadata;
    tpa_metadata["excited_states"][pkey]
        = {{"type", "full"}, {"n_roots", 2}, {"converged", true},
           {"stop_reason", "converged"},
           {"roots", json::array({root(0, 0, 0.42, true, false, false),
                                  root(1, 1, 0.47, true, false, false)})}};
    json tpa_properties;
    tpa_properties["tpa"][pkey]
        = json::array({tpa_row(0, 0.42, true), tpa_row(1, 0.47, true)});
    auto tp = named_stages::selected_two_photon(stage, tpa_properties, tpa_metadata);
    const auto& tp_root = tp.at("roots").at(0);
    check(tp.at("converged") == true && tp.at("units") == "au"
              && tp.at("roots").size() == 2 && tp_root.at("root") == 0
              && tp_root.at("omega") == 0.42 && tp_root.at("photon_frequency") == 0.21
              && tp_root.at("D_linear") == 4.2668 && tp_root.at("tensor")[2][2] == 3.0
              && tp_root.at("legs").size() == 3 && tp_root.at("legs_converged") == true
              && tp_root.at("root_convergence").at("density_residual") == 2e-6
              && tp_root.at("root_convergence").at("excitation_energy_change") == 3e-7
              && !tp_root.contains("sigma_linear_gm"),
          "selected_two_photon (D-113): per root the tensor, D, the root's own "
          "convergence and its legs; no GM cross sections");
    tpa_properties["tpa"][pkey] = json::array({tpa_row(0, 0.42, false), tpa_row(1, 0.47, true)});
    check(named_stages::selected_two_photon(stage, tpa_properties, tpa_metadata)
                  .at("converged") == false,
          "selected_two_photon: a leg not converged, the record not converged");
    tpa_properties["tpa"][pkey] = json::array({tpa_row(0, 0.42, true)});
    check(named_stages::selected_two_photon(stage, tpa_properties, tpa_metadata)
                  .at("converged") == false,
          "selected_two_photon: a root without a row, not converged");
    check(named_stages::selected_two_photon(stage, json::object(), tpa_metadata)
              .contains("unavailable"),
          "selected_two_photon: no rows, unavailable");
  } catch (const std::exception& e) {
    check(false, e.what());
  }

  // ----- D-061: the econv row of the settings table and the bundle's judged
  // row, through the shared judge (B1) ----------------------------------------
  try {
    namespace cv = madness::convergence;
    Stage stage = bare_stage(world);
    stage.properties = {"excited_states"};
    stage.excited_states.roots = 2;
    stage.excited_states.method = "tda";
    const std::string instead
        = "stops the bundle when each root's eigenvalue change is below "
          "5.0e-04 (5 x max(dconv, thresh)), strictly, in both modes; a "
          "written econv is judged by the stage verdict";
    const auto unset
        = named_stages::excited_state_omega_target_row(stage, 5e-4);
    check(unset.group == cv::Group::requested_criteria
              && unset.name == "excitation energy change per root"
              && unset.keyword == "econv"
              && unset.measurement == "excitation_energy_change"
              && unset.quantity.empty() && unset.units == "hartree"
              && !unset.target && unset.set_by == cv::SetBy::not_set
              && unset.relation == cv::Relation::less_than_or_equal
              && !unset.effective_limit && unset.limit_formula.empty()
              && unset.note.find("the engine's default is never a target")
                     != std::string::npos
              && unset.engine_check_instead == instead,
          "econv row: unset, the engine's gate named as what is enforced "
          "instead");
    stage.excited_states.omega_target = 1e-6;
    const auto row = named_stages::excited_state_omega_target_row(stage, 5e-4);
    check(row.target == 1e-6 && row.set_by == cv::SetBy::user
              && row.note.empty() && row.engine_check_instead == instead,
          "econv row: a written econv is the target, set by the user");
    cv::ConvergenceSettings settings;
    settings.stage = "excited";
    settings.decider = cv::Decider::engine_check;
    settings.metrics.push_back(row);
    const auto warnings = cv::unenforced_target_warnings(settings);
    check(warnings.size() == 1
              && warnings.front()
                     == "econv 1.0e-06 is reported, not enforced: the engine "
                        "check "
                            + instead,
          "econv row: the engine-mode warning at stage start");
    const auto targets = cv::requested_targets(settings);
    check(targets.size() == 1
              && targets.at("excitation_energy_change").at("value") == 1e-6
              && targets.at("excitation_energy_change").at("units") == "hartree"
              && !targets.at("excitation_energy_change").contains("quantity"),
          "econv row: the requested target the shared judge takes");
    // The bundle's engine rows as the state record carries them
    // (bundle_measurements): the worst root per measurement.
    const json engine = {{"density_residual",
                          {{"value", 2e-6},
                           {"availability", "measured"},
                           {"root_id", "es_root_0001"},
                           {"phase", "solver_step_input_before_acceleration"}}},
                         {"excitation_energy_change",
                          {{"value", 3e-7},
                           {"availability", "measured"},
                           {"root_id", "es_root_0001"},
                           {"phase", "solver_step_output"}}}};
    const auto saved = named_stages::bundle_saved_measurements(engine);
    check(saved.size() == 1
              && saved.at("excitation_energy_change").at("value") == 3e-7
              && saved.at("excitation_energy_change").at("units") == "hartree"
              && saved.at("excitation_energy_change").at("availability")
                     == "measured"
              && saved.at("excitation_energy_change").at("root_id")
                     == "es_root_0001"
              && saved.at("excitation_energy_change").at("phase")
                     == "solver_step_output",
          "the bundle's judged row: the worst root's eigenvalue change with "
          "its units");
    const auto met = madness::evaluate_explicit_criteria(saved, targets);
    check(met.at("status") == "satisfied"
              && met.at("criteria").at("excitation_energy_change").at("value")
                     == 3e-7,
          "the shared judge: 3e-7 meets econv 1e-6");
    stage.excited_states.omega_target = 1e-7;
    settings.metrics
        = {named_stages::excited_state_omega_target_row(stage, 5e-4)};
    check(madness::evaluate_explicit_criteria(saved,
                                              cv::requested_targets(settings))
                  .at("status")
              == "unmet",
          "the shared judge: 3e-7 misses econv 1e-7");
    json none = engine;
    none["excitation_energy_change"] = {{"value", nullptr},
                                        {"availability", "unavailable"},
                                        {"phase", nullptr}};
    check(madness::evaluate_explicit_criteria(
              named_stages::bundle_saved_measurements(none),
              cv::requested_targets(settings))
                  .at("status")
              == "unavailable",
          "the shared judge: a root without an eigenvalue step yet");
    stage.excited_states.omega_target.reset();
    settings.metrics
        = {named_stages::excited_state_omega_target_row(stage, 5e-4)};
    check(cv::requested_targets(settings).empty()
              && madness::evaluate_explicit_criteria(
                     saved, cv::requested_targets(settings))
                         .at("status")
                     == "not_requested",
          "no econv: nothing requested, nothing judged");
    check(named_stages::bundle_saved_measurements(json::object()).empty(),
          "no engine rows: no judged row");
  } catch (const std::exception& e) {
    check(false, std::string("econv row: ") + e.what());
  }

  // D-095: the bundle's engine checks as per-root outcomes under the two
  // engine rows, so the outcome table shows the evidence it judged on.
  try {
    namespace cv = madness::convergence;
    using cv::Group;
    using cv::Relation;
    const json bundle = {
        {"roots",
         json::array({{{"root_id", "es_root_0000"},
                       {"density_residual", 2e-4},
                       {"omega_residual", 3e-7},
                       {"density_target", 5e-4},
                       {"omega_target", 5e-4}},
                      {{"root_id", "es_root_0001"},
                       {"density_residual", 6e-4},
                       {"omega_residual", nullptr},
                       {"omega_residual_unavailable", "no eigenvalue step yet"},
                       {"density_target", nullptr},
                       {"omega_target", nullptr},
                       {"targets_unavailable", "no convergence check yet"}}})}};
    const auto legs = named_stages::bundle_root_outcomes(bundle);
    check(legs.size() == 4 && legs[0].measurement == "engine_es_density_change"
              && legs[0].subject == "es_root_0000" && legs[0].value == 2e-4
              && legs[0].limit == 5e-4
              && legs[0].units == "kernel_induced_density_L2_atomic_units"
              && legs[0].phase == "solver_step_input_before_acceleration"
              && legs[0].availability == "measured"
              && legs[1].subject == "es_root_0001" && legs[1].value == 6e-4
              && !legs[1].limit
              && legs[2].measurement == "engine_es_omega_change"
              && legs[2].subject == "es_root_0000" && legs[2].value == 3e-7
              && legs[2].limit == 5e-4 && legs[2].units == "hartree"
              && legs[2].phase == "solver_step_output"
              && legs[3].subject == "es_root_0001" && !legs[3].value
              && !legs[3].limit && legs[3].availability == "unavailable",
          "per-root outcomes: two rows times two roots, each root's own "
          "value and recorded target, the engine's phases and units");
    check(named_stages::bundle_root_outcomes(json::object()).empty(),
          "no roots: no per-root outcomes");
    // The shared record and verdict over those legs: the two engine rows as
    // the response settings table writes them (ResponseCheckpoint.hpp).
    cv::ConvergenceSettings settings;
    settings.decider = cv::Decider::engine_check;
    settings.metrics.push_back({Group::engine_check, "density change per root",
                                "dconv", "engine_es_density_change", "",
                                "kernel_induced_density_L2_atomic_units", 1e-4,
                                cv::SetBy::user, Relation::less_than, 5e-4,
                                "5 x max(dconv, thresh)",
                                "per root; verdicts in the bundle record", ""});
    settings.metrics.push_back(
        {Group::engine_check, "excitation energy change per root", "dconv",
         "engine_es_omega_change", "", "hartree", 1e-4, cv::SetBy::user,
         Relation::less_than, 5e-4, "5 x max(dconv, thresh)",
         "per root; verdicts in the bundle record", ""});
    const auto record = cv::outcome_record(settings, legs);
    const auto& density = record.at("metrics")[0];
    const auto& omega = record.at("metrics")[1];
    check(density.at("legs").size() == 2
              && density.at("legs")[0].at("subject") == "es_root_0000"
              && density.at("legs")[0].at("value") == 2e-4
              && density.at("legs")[0].at("limit") == 5e-4
              && density.at("legs")[0].at("verdict") == "met"
              && density.at("legs")[1].at("subject") == "es_root_0001"
              && density.at("legs")[1].at("value") == 6e-4
              && density.at("legs")[1].at("limit").is_null()
              && density.at("legs")[1].at("verdict") == "unmet"
              && density.at("value") == 6e-4 && density.at("verdict") == "unmet"
              && density.at("phase") == "solver_step_input_before_acceleration",
          "the record's density row: one leg per root, the second root "
          "judged on the row's limit and unmet, the row unmet");
    check(omega.at("legs").size() == 2
              && omega.at("legs")[0].at("verdict") == "met"
              && omega.at("legs")[1].at("value").is_null()
              && omega.at("legs")[1].at("availability") == "unavailable"
              && omega.at("legs")[1].at("verdict") == "unavailable"
              && omega.at("value") == 3e-7
              && omega.at("availability") == "unavailable"
              && omega.at("verdict") != "met",
          "the record's omega row: a root without an eigenvalue step is "
          "unavailable, so the row is not met");
  } catch (const std::exception& e) {
    check(false, std::string("per-root outcomes: ") + e.what());
  }

  world.gop.sum(failures);
  if (world.rank() == 0)
    std::cout << "response-property failures: " << failures << '\n';
  madness::finalize();
  return failures ? 1 : 0;
}
