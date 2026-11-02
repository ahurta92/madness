#pragma once

// Load MRA function declarations before ResponseProtocol includes its operators.
#include <madness/chem/ParameterManager.hpp>
#include <madness/chem/SCF.h>

#include <apps/molresponse/ResponseProtocol.hpp>
#include <apps/molresponse/solvers/fd_explicit_accuracy.hpp>
#include <apps/molresponse/solvers/response_metadata.hpp>
#include <array>
#include <charconv>
#include <cmath>
#include <filesystem>
#include <functional>
#include <madness/chem/LadderPlanner.h>
#include <madness/chem/StageCheckpoint.h>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>

/// Deliberately small, flat input surface for fresh moldft -> molresponse runs.
/// Scientific values are parsed by the existing typed parameter groups.
/// Flow: input text -> named blocks -> typed stages with resolved references.
/// This file does no SCF or response calculation; RunNamedStages.hpp executes
/// the returned stages. The HF/alpha restrictions below are frontend policy.
namespace named_stages {

// Keep source line numbers so invalid values can be traced back to the deck.
struct Line {
  std::size_t number;
  std::string text;
};
struct Block {
  std::string kind, name;
  std::vector<Line> lines;
};

inline std::string clean(std::string s) {
  s = s.substr(0, s.find('#')); // remove the comment
  const auto first =
      s.find_first_not_of(" \t\r"); // not space,tab, or cariage return
  return first == std::string::npos
             ? ""
             : s.substr(first, s.find_last_not_of(" \t\r") - first + 1);
}
inline std::string key(const std::string &s) {
  auto word = s.substr(0, s.find_first_of(" \t="));
  std::transform(word.begin(), word.end(), word.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return word;
}

// report why the line number fails
inline void fail(const Line &line, const std::string &why) {
  throw std::runtime_error("named stages line " + std::to_string(line.number) +
                           ": " + why);
}
// Select the named-input path; full syntax validation happens in blocks().
inline bool named_input(const std::string &text) {
  std::istringstream input(text);
  std::string line;
  while (std::getline(input, line)) {
    line = clean(line);
    if (line.empty())
      continue;
    std::istringstream words(line);
    std::string kind, name;
    words >> kind >> name;
    return !name.empty() && line.front() != '{';
  }
  return false;
}
// Read only block structure: "moldft gs", its body, then "end".
// Parameter types and reference bindings are handled later by parse().
inline std::vector<Block> blocks(const std::string &text) {

  std::vector<Block> result;   // the blocks it finds
  std::set<std::string> names; // the block names to reject duplicates

  std::istringstream input(text);
  std::string textline;
  bool inside = false; // am i header? or am I inside the block
  std::size_t n = 0;

  while (std::getline(input, textline)) { // for each line of text
    Line line{++n, clean(textline)};      // store the line number and text
    if (line.text.empty())
      continue;
    if (!inside) {
      std::istringstream words(line.text); // stringstream the cleaned line
      Block b; // create a Block, kind, and clean lines
      std::string extra;
      // The assumption is that the first line will be a header
      words >> b.kind >> b.name;
      // So we are checking if it's one of the supported block types
      if (b.kind != "molecule" && b.kind != "moldft" && b.kind != "molresponse")
        fail(line, "unsupported stage type " + b.kind);
      // Header contains more than kind and name not safe
      if (b.name.empty() || words >> extra ||
          !std::isalpha(static_cast<unsigned char>(b.name.front())) ||
          b.name.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQ"
                                   "RSTUVWXYZ0123456789_-") !=
              std::string::npos)
        fail(line, "expected a safe stage name after " + b.kind);
      // if the name was found again, you failed... no duplicate names!
      if (!names.insert(b.name).second)
        fail(line, "duplicate name " + b.name);
      result.push_back(std::move(b)); // safe move from inner to outerscope
      inside = true;                  // now we are inside
    } else if (line.text == "end") {
      inside = false;
    } else {
      if (key(line.text) == "end")
        fail(line, "unexpected text after end");
      // A stage header inside a block: the block above was not closed.
      if (const auto k = key(line.text); k == "moldft" || k == "molresponse")
        fail(line, "a " + k + " block starts here, but "
                       + result.back().kind + " " + result.back().name
                       + " has no end line before it");
      result.back().lines.push_back(std::move(line));
      // to the last block, push_back to the lines of that block
    }
  }
  if (inside)
    throw std::runtime_error("unterminated named block " + result.back().name);
  return result;
}

// Keywords owned by the stage interface. Solver settings such as dconv and
// maxiter continue to use CalculationParameters and ResponseParameters.
struct Controls : madness::QCCalculationParametersBase {
  Controls() {
    initialize<std::string>("molecule", "", "named molecule");
    initialize<std::string>("reference", "", "earlier moldft stage");
    initialize<std::vector<std::string>>("requested_properties", {},
                                         "requested properties");
    // D-017: one block per property; frequencies are always written (B3),
    // components mean every Cartesian tuple when omitted.
    initialize<std::vector<double>>("polarizability.frequencies", {},
                                    "explicit frequencies");
    initialize<std::string>("polarizability.frequency_units", "",
                            "frequency units", {"", "au"});
    initialize<std::vector<std::string>>("polarizability.components", {},
                                         "Cartesian pairs; omitted: all 9");
    initialize<std::vector<double>>("hyperpolarizability.frequencies", {},
                                    "explicit driver frequencies (SHG: "
                                    "output at 2w; static at 0)");
    initialize<std::string>("hyperpolarizability.frequency_units", "",
                            "frequency units", {"", "au"});
    initialize<std::vector<std::string>>("hyperpolarizability.components", {},
                                         "Cartesian triples; omitted: all 27");
    initialize<std::vector<double>>("raman.frequencies", {},
                                    "explicit optical frequencies");
    initialize<std::string>("raman.frequency_units", "", "frequency units",
                            {"", "au"});
    initialize<std::vector<std::string>>("raman.components", {},
                                         "Cartesian pairs; omitted: all 9");
    initialize<std::vector<int>>("raman.atoms", {},
                                 "displaced atoms, numbered from 0; omitted: "
                                 "every atom");
    initialize<std::vector<std::string>>("raman.directions", {},
                                         "displacement axes; omitted: x, y, z");
    // D-025: the excited-state bundle; roots and method are always written.
    initialize<int>("excited_states.roots", 0,
                    "how many roots, from the lowest (required)");
    initialize<std::string>("excited_states.method", "",
                            "tda or rpa (required)", {"", "tda", "rpa"});
    initialize<std::string>("excited_states.guess", "",
                            "initial guess; omitted: the engine's "
                            "solid_harmonics",
                            {"", "solid_harmonics", "virtual_ao", "random"});
    initialize<std::string>("excited_states.guess_basis", "",
                            "AO basis of the virtual_ao guess; omitted: the "
                            "engine's aug-cc-pvdz");
    // DESIGN A3: when the stage's observer publishes, and what it keeps.
    initialize<int>("checkpoint.every_iterations", 0,
                    "also publish every K iterations within a rung (K >= 1); "
                    "default: only at rung ends and stops");
    initialize<bool>("checkpoint.keep_all", false,
                     "keep every generation (debugging); default: delete "
                     "older mid-rung generations");
    // S5, DESIGN B3: the Dalton seed (a starting guess for states with no
    // saved state), of a molresponse stage and (D-112) a moldft stage; no
    // defaults.
    initialize<std::string>("seed.dalton_dir", "",
                            "a Dalton output directory seeding the stage: a "
                            "moldft stage's orbitals from molden.inp; a "
                            "molresponse stage's excited states and driven "
                            "legs from RSPVEC");
    initialize<double>("seed.freq_tol", 0.0,
                       "the largest frequency distance (au) at which a "
                       "driven leg takes a Dalton RSPVEC record as its seed");
  }
  [[nodiscard]] std::string get_tag() const override { return "stage"; }
  [[nodiscard]] std::set<std::string> case_sensitive_keys() const override {
    return {"molecule", "reference"};
  }
};
// One optical property request of a molresponse stage: explicit frequencies
// and the Cartesian components (every tuple when the deck wrote none).
struct PropertyRequest {
  std::vector<double> frequencies;
  std::vector<std::string> components;
};
// One application instance. Its Params is independent of every other stage,
// including stages of the same application type.
struct Stage {
  std::string kind, name, reference;
  Params params;
  std::vector<std::string> properties;
  // The polarizability request (ResponseCheckpoint.hpp reads these two).
  std::vector<std::string> components;
  std::vector<double> frequencies;
  // D-017: the hyperpolarizability and Raman requests. Raman coordinates are
  // atoms x directions; empty atoms means every atom of the reference
  // (numbered from 0), empty directions means x, y and z. Both omitted is the
  // engine's full set with the normal-mode analysis (raman_full_set); a
  // written subset gives derivative rows only.
  PropertyRequest hyperpolarizability, raman;
  std::vector<int> raman_atoms;
  std::string raman_directions;
  bool raman_full_set = false;
  // D-025: the excited-state request: one bundle of the `roots` lowest roots
  // by `method` (tda or rpa). guess and guess_basis as written, empty for
  // the engine's defaults (the record names them).
  struct ExcitedStatesRequest {
    int roots = 0;
    std::string method, guess, guess_basis;
    // D-061: a written econv, the target of the excitation energy change per
    // root (D8); empty when the deck wrote none (the engine's default econv
    // is not a request).
    std::optional<double> omega_target;
  };
  ExcitedStatesRequest excited_states;
  // A3: checkpoint cadence and retention (moldft); 0 is "decision points only".
  int checkpoint_every_iterations = 0;
  bool checkpoint_keep_all = false;
  // S5, B3, D-112: the Dalton seed directory (absolute; empty: unseeded) and
  // the driven legs' frequency tolerance (0 when the stage has no driven
  // legs, always for moldft).
  std::string seed_dalton_dir;
  double seed_freq_tol = 0.0;
  // The keywords this stage's deck wrote, so "who set it" (B9) does not
  // depend on values the parser fills in for the engine.
  std::set<std::string> written;
  Stage(madness::World &world, const Block &b)
      : kind(b.kind), name(b.name), params(world) {}
};
// Wrap one line in a temporary block to reuse the existing typed value parser.
// No temporary file is created.
inline void read_value(madness::World &world,
                       madness::QCCalculationParametersBase &group,
                       const Line &line) {
  std::string text = "stage\n" + line.text + "\nend\n";
  try {
    group.read_internal(world, text, "stage");
  } catch (const std::exception &e) {
    fail(line, e.what());
  }
}
inline void valid_protocol(const std::vector<double> &protocol) {
  if (protocol.empty())
    throw std::runtime_error("protocol must not be empty");
  double previous = std::numeric_limits<double>::infinity();
  for (double t : protocol) {
    if (!std::isfinite(t) || t <= 0 || t >= previous)
      throw std::runtime_error(
          "protocol must contain positive, strictly decreasing thresholds");
    previous = t;
  }
}
inline int protocol_k(double t) {
  return molresponse_v3::default_k_for_thresh(t);
}

// D-107: the nuclear coordinates a Raman request covers, in its written
// order: the written atoms or every atom of the reference, each with the
// written directions or x, y, z. Callers validate the atoms.
[[nodiscard]] inline std::vector<std::pair<int, char>>
raman_coordinates(std::vector<int> atoms, std::string directions, int natom) {
  if (atoms.empty())
    for (int atom = 0; atom < natom; ++atom)
      atoms.push_back(atom);
  if (directions.empty())
    directions = "xyz";
  std::vector<std::pair<int, char>> coordinates;
  for (int atom : atoms)
    for (char axis : directions)
      coordinates.emplace_back(atom, axis);
  return coordinates;
}

// The named property language of a molresponse stage (D-017). Each optical
// property in requested_properties has its own block: <property>.frequencies
// and <property>.frequency_units are required (B3: a number is never
// defaulted), <property>.components is optional and means every Cartesian
// tuple when omitted. Raman adds raman.atoms and raman.directions (see
// Stage). excited_states (D-025) has excited_states.roots and
// excited_states.method, both required, and the optional guess.
inline const std::vector<std::string>& named_response_properties() {
  static const std::vector<std::string> names
      = {"polarizability", "hyperpolarizability", "raman", "excited_states"};
  return names;
}
// The tensor rank of a property's Cartesian components: beta_ABC is a
// triple, alpha_AB and d(alpha_AB)/dQ are pairs.
inline std::size_t component_rank(const std::string& property) {
  return property == "hyperpolarizability" ? 3 : 2;
}
// Every Cartesian tuple of the given rank, x before y before z in each place
// (xx, xy, xz, yx, ... for pairs).
inline std::vector<std::string> every_component(std::size_t rank) {
  std::vector<std::string> tuples = {""};
  for (std::size_t place = 0; place < rank; ++place) {
    std::vector<std::string> longer;
    for (const auto& prefix : tuples)
      for (char axis : std::string("xyz"))
        longer.push_back(prefix + axis);
    tuples = std::move(longer);
  }
  return tuples;
}
// Frequencies must be finite, nonnegative and distinct at the engine's
// frequency-key precision: the engine rounds frequencies in state filenames,
// and two that round alike would share a stored state.
inline std::vector<double>
checked_frequencies(const std::string& property,
                    const std::vector<double>& frequencies,
                    const std::string& units) {
  if (frequencies.empty() || units != "au")
    throw std::runtime_error(property
                             + " needs explicit frequencies and "
                               "frequency_units au");
  std::set<std::string> keys;
  for (double f : frequencies)
    if (!std::isfinite(f) || f < 0
        || !keys.insert(molresponse_v3::ResponseMetadata::freq_key(f)).second)
      throw std::runtime_error(
          property
          + " frequencies must be finite, nonnegative and distinct "
            "at engine frequency-key precision");
  return frequencies;
}
// Components are distinct Cartesian tuples of the property's rank; none
// written means every one.
inline std::vector<std::string>
checked_components(const std::string& property,
                   const std::vector<std::string>& components) {
  const auto rank = component_rank(property);
  if (components.empty())
    return every_component(rank);
  std::set<std::string> distinct;
  for (const auto& c : components)
    if (c.size() != rank || c.find_first_not_of("xyz") != std::string::npos
        || !distinct.insert(c).second)
      throw std::runtime_error(property
                               + " components must be distinct "
                                 "Cartesian "
                               + (rank == 3 ? "triples" : "pairs"));
  return components;
}
// Read the property blocks of a molresponse stage into its requests and map
// them onto the engine's knobs. The engine takes axes and one frequency grid,
// not components: it gets the union of every axis and frequency named (it may
// compute extra tensor elements), and the runner selects the requested
// components from its rows (ResponseProperties.hpp). The raw engine spellings
// (dipole.*, quadratic, beta.*, raman.nuc_*, nuclear.*) are not stage
// controls, so parse rejects them: one spelling per thing.
inline void translate_response_properties(const std::string& stage_name,
                                          const Controls& controls,
                                          Stage& stage,
                                          ResponseParameters& rp) {
  const auto& known = named_response_properties();
  if (stage.properties.empty())
    throw std::runtime_error(
        stage_name
        + ": molresponse requested_properties must name at least "
          "one of polarizability, hyperpolarizability, raman, "
          "excited_states");
  std::set<std::string> requested;
  for (const auto& p : stage.properties)
    if (std::find(known.begin(), known.end(), p) == known.end()
        || !requested.insert(p).second)
      throw std::runtime_error(
          stage_name + ": unsupported or repeated molresponse property " + p
          + " (polarizability, hyperpolarizability, raman, excited_states)");
  // A property's block without the property is a mistake, not a default.
  for (const auto& p : known) {
    if (requested.contains(p))
      continue;
    for (const auto* suffix :
         {".frequencies", ".frequency_units", ".components", ".atoms",
          ".directions", ".roots", ".method", ".guess", ".guess_basis"})
      if (controls.parameter_exists(p + suffix)
          && controls.is_user_defined(p + suffix))
        throw std::runtime_error(stage_name + ": " + p + suffix + " needs " + p
                                 + " in requested_properties");
  }
  auto request = [&](const std::string& p) {
    PropertyRequest r;
    r.frequencies = checked_frequencies(
        p, controls.get<std::vector<double>>(p + ".frequencies"),
        controls.get<std::string>(p + ".frequency_units"));
    r.components = checked_components(
        p, controls.get<std::vector<std::string>>(p + ".components"));
    return r;
  };
  // The engine's grid and axes: the union over the properties requested, the
  // polarizability's first so an alpha-only deck maps as it always did.
  std::vector<double> grid;
  std::set<std::string> grid_keys;
  std::set<char> axes;
  auto gather = [&](const PropertyRequest& r) {
    for (double f : r.frequencies)
      if (grid_keys.insert(molresponse_v3::ResponseMetadata::freq_key(f))
              .second)
        grid.push_back(f);
    for (const auto& c : r.components)
      axes.insert(c.begin(), c.end());
  };
  if (requested.contains("polarizability")) {
    const auto r = request("polarizability");
    stage.frequencies = r.frequencies;
    stage.components = r.components;
    gather(r);
  }
  if (requested.contains("hyperpolarizability")) {
    stage.hyperpolarizability = request("hyperpolarizability");
    gather(stage.hyperpolarizability);
  }
  if (requested.contains("raman")) {
    stage.raman = request("raman");
    gather(stage.raman);
    stage.raman_atoms = controls.get<std::vector<int>>("raman.atoms");
    std::set<int> distinct_atoms;
    for (int atom : stage.raman_atoms)
      if (atom < 0 || !distinct_atoms.insert(atom).second)
        throw std::runtime_error(
            stage_name
            + ": raman.atoms must be distinct atom indices, numbered from 0");
    // The written atoms are checked against the reference's atom count when
    // the stage runs (the engine knows the molecule; the parser does not).
    std::set<char> distinct_axes;
    for (const auto& d :
         controls.get<std::vector<std::string>>("raman.directions")) {
      if (d.size() != 1 || std::string("xyz").find(d[0]) == std::string::npos
          || !distinct_axes.insert(d[0]).second)
        throw std::runtime_error(
            stage_name
            + ": raman.directions must be distinct axes out of "
              "x, y, z");
      stage.raman_directions += d[0];
    }
    stage.raman_full_set
        = stage.raman_atoms.empty() && stage.raman_directions.empty();
  }
  rp.set_user_defined_value("requested_properties", stage.properties);
  rp.set_user_defined_value("dipole.frequencies", grid);
  rp.set_user_defined_value("dipole.directions",
                            std::string(axes.begin(), axes.end()));
  if (requested.contains("hyperpolarizability"))
    // SHG legs at w and 2w for these driver frequencies only (static at 0).
    rp.set_user_defined_value("beta.frequencies",
                              stage.hyperpolarizability.frequencies);
  if (requested.contains("raman")) {
    // The engine's raman.atoms and raman.directions take the stage's lists as
    // written (D-020): the engine plans their product, every atom or x, y, z
    // for an empty one, derivative rows only. The full set (both omitted)
    // is the engine's raman.nuc_atom -1, its only full-set switch until the
    // unnamed path is ported; the lists stay unwritten then.
    rp.set_user_defined_value("raman.frequencies", stage.raman.frequencies);
    if (stage.raman_full_set)
      rp.set_user_defined_value("raman.nuc_atom", -1);
    if (!stage.raman_atoms.empty())
      rp.set_user_defined_value("raman.atoms", stage.raman_atoms);
    if (!stage.raman_directions.empty())
      rp.set_user_defined_value("raman.directions", stage.raman_directions);
  }
  if (requested.contains("excited_states")) {
    // D-025: one bundle of the `roots` lowest roots by `method`. Both are
    // written: a root count is a number (B3) and TDA against the full (RPA)
    // response is a method, not a default. The guess is a method choice
    // with the engine's default. The stage's maxiter caps the bundle's rungs
    // as it caps the dipole legs' (B5, one cap): excited.maxiter stays
    // unwritten, so the engine's bundle budget is maxiter
    // (ExecutorSettings::es_iter_budget).
    auto& es = stage.excited_states;
    if (!controls.is_user_defined("excited_states.roots")
        || !controls.is_user_defined("excited_states.method"))
      throw std::runtime_error(
          stage_name
          + ": excited_states needs excited_states.roots (how many, from the "
            "lowest) and excited_states.method (tda or rpa); neither has a "
            "default");
    es.roots = controls.get<int>("excited_states.roots");
    if (es.roots < 1)
      throw std::runtime_error(stage_name
                               + ": excited_states.roots must be at least 1");
    es.method = controls.get<std::string>("excited_states.method");
    if (es.method != "tda" && es.method != "rpa")
      throw std::runtime_error(stage_name
                               + ": excited_states.method must be tda or rpa");
    es.guess = controls.get<std::string>("excited_states.guess");
    es.guess_basis = controls.get<std::string>("excited_states.guess_basis");
    if (!es.guess_basis.empty() && es.guess != "virtual_ao")
      throw std::runtime_error(
          stage_name
          + ": excited_states.guess_basis applies to the virtual_ao guess "
            "only");
    rp.set_user_defined_value("excited.enable", true);
    rp.set_user_defined_value("excited.num_states",
                              static_cast<std::size_t>(es.roots));
    rp.set_user_defined_value("excited.tda", es.method == "tda");
    if (!es.guess.empty())
      rp.set_user_defined_value("excited.guess", es.guess);
    if (!es.guess_basis.empty())
      rp.set_user_defined_value("excited.guess_basis", es.guess_basis);
  }
}

// D-102: what a refusal says about a reference's saved generation under
// `root`: " ; <ref>'s saved generation <g> is at thresh <t>, k <k>", or empty
// when there is none (or no root). Collective: rank 0 reads, every rank gets
// the same text.
inline std::string saved_reference_rung(madness::World& world,
                                        const std::filesystem::path& root,
                                        const std::string& reference) {
  if (root.empty())
    return {};
  return madness::checkpoint::detail::on_root(world, [&]() -> std::string {
    try {
      const auto store = root / reference / "checkpoints";
      if (!std::filesystem::exists(store / "latest"))
        return {};
      const auto m = madness::checkpoint::read_manifest(
          store
          / madness::checkpoint::detail::read_json(store / "latest")
                .at("generation")
                .get<std::string>());
      std::ostringstream text;
      text << "; " << reference << "'s saved generation " << m.generation
           << " is at thresh "
           << madness::convergence::detail::scientific(
                  m.numerical.at("threshold").get<double>())
           << ", k " << m.numerical.at("k").get<int>();
      return text.str();
    } catch (const std::exception&) {
      return {};
    }
  });
}

// Build and validate the whole sequence before any calculation starts.
// Names enter these lookup tables in file order, so references must point back.
// `output_root` (D-102): where a refusal looks for a reference's saved
// generation to name its rung; empty: nowhere.
inline std::vector<Stage> parse(madness::World& world, const std::string& text,
                                const std::filesystem::path& output_root = {}) {
  std::map<std::string, madness::Molecule> molecules;
  std::map<std::string, std::size_t> indices;
  std::vector<Stage> stages;
  // This bounded allowlist prevents accepting settings whose outputs or
  // scientific semantics have not yet been wired through the named adapter.
  const std::set<std::string> scf_keys
      = {"xc",
         "protocol",
         "dconv",
         "econv",
         "maxiter",
         "l",
         "k",
         "charge",
         "localize",
         "print_level",
         "maxsub",
         "maxrotn",
         "diagnostic.iterations_per_rung",
         "progress.window",
         "progress.minimum_relative_improvement",
         "accuracy.bsh_residual",
         "accuracy.kinetic_trace_change",
         madness::convergence::use_default_criteria_key};
  // B8: dconv and econv are the density- and energy-change targets of a
  // moldft stage; the older keywords point there.
  const std::map<std::string, std::string> replaced_scf_keys
      = {{"accuracy.density_change", "dconv (the density-change target)"},
         {"accuracy.energy_change", "econv (the energy-change target)"}};
  // stall.window and stall.ratio (D-099, B3): the engine's own keys for the
  // excited-state bundle's plateau detector; refused below on a stage
  // without excited_states.
  const std::set<std::string> response_keys
      = {"dconv",
         "maxiter",
         "kain",
         "diagnostic.iterations_per_rung",
         "progress.window",
         "progress.minimum_relative_improvement",
         "accuracy.bsh_residual",
         "accuracy.kinetic_trace_change",
         madness::convergence::use_default_criteria_key,
         "print_level",
         "stall.window",
         "stall.ratio"};
  // B8: dconv is the induced-density-change target of a molresponse stage; a
  // driven response has no local energy, so no energy target.
  const std::map<std::string, std::string> replaced_response_keys
      = {{"accuracy.density_change", "dconv (the density-change target)"},
         {"dalton.dir", "now seed.dalton_dir (the Dalton seed, S5)"},
         {"accuracy.energy_change",
          "not applicable: a driven response has no local energy (see the "
          "reference stage's energy change)"}};
  // DESIGN A1: a reference names an earlier moldft stage of this deck, or a
  // stage outside it whose published generation is read at run time. A name
  // used later in the deck would be a forward reference.
  std::map<std::string, std::string> deck_kinds;
  for (const auto& b : blocks(text))
    deck_kinds.emplace(b.name, b.kind);
  for (const auto &b : blocks(text)) {
    std::set<std::string> seen;

    // Molecule blocks define reusable geometry, not executable stages.
    if (b.kind == "molecule") {
      madness::Molecule mol;
      std::string body = "molecule\n";
      const std::set<std::string> geometry_keys = {"units", "eprec",
                                                   "no_orient"};
      for (const auto &line : b.lines) {
        const auto k = key(line.text);
        if (geometry_keys.count(k)) {
          if (!seen.insert(k).second)
            fail(line, "duplicate parameter " + k);
          read_value(world, mol.parameters, line);
        } else {
          std::istringstream atom(line.text);
          std::string symbol, extra;
          double x = 0, y = 0, z = 0;
          if (!(atom >> symbol >> x >> y >> z) || atom >> extra ||
              !std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
            fail(line, "expected atom and three finite coordinates");
          if (madness::symbol_to_atomic_number(symbol) <= 0)
            fail(line, "expected an ordinary atom in coordinates");
          body += line.text + "\n";
        }
      }
      body += "end\n";
      std::istringstream stream(body);
      mol.read(stream);
      if (!mol.natom())
        throw std::runtime_error("empty molecule " + b.name);
      molecules.emplace(b.name, std::move(mol));
      continue;
    }

    // Route each keyword to the interface controls or its solver group.
    Stage stage(world, b);
    Controls controls;
    // cp holds SCF settings and the protocol used by the response adapter;
    // rp holds this stage's response iteration settings.
    auto &cp = stage.params.get<CalculationParameters>();
    auto &rp = stage.params.get<ResponseParameters>();
    for (const auto &line : b.lines) {
      const auto k = key(line.text);
      if (!seen.insert(k).second)
        fail(line, "duplicate parameter " + k);
      if (b.kind == "moldft" && replaced_scf_keys.contains(k))
        fail(line, k + " is now " + replaced_scf_keys.at(k));
      if (b.kind == "molresponse" && replaced_response_keys.contains(k))
        fail(line, k + " is " + (k == "accuracy.density_change" ? "now " : "")
                       + replaced_response_keys.at(k));
      if (k == "maxiter" && line.text.find('[') != std::string::npos) {
        // B5: "maxiter [30, 10]" gives each rung its own cap (moldft only).
        if (b.kind != "moldft")
          fail(line, "a per-rung maxiter list is supported for moldft only "
                     "(molresponse follows with the reuse planner); "
                     "give one number");
        read_value(world, cp,
                   {line.number,
                    "maxiter_per_rung"
                        + line.text.substr(std::string("maxiter").size())});
      } else if (controls.parameter_exists(k))
        read_value(world, controls, line);
      else if (k == "protocol" || (b.kind == "moldft" && scf_keys.count(k)))
        read_value(world, cp, line);
      else if (b.kind == "molresponse" && k == "econv")
        // D-061: the excitation energy change target of an excited-state
        // stage (D8), read where the key lives; refused below on a stage
        // without excited_states.
        read_value(world, cp, line);
      else if (b.kind == "molresponse" && response_keys.count(k) &&
               rp.parameter_exists(k))
        read_value(world, rp, line);
      else
        fail(line, "unsupported parameter " + k + " in " + b.kind);
    }

    // A3: a stage publishes inside a solve only where it has a save boundary
    // there; molresponse publishes at solve ends only (DESIGN §1.4).
    // checkpoint.keep_all is every application's (D-070).
    if (seen.contains("checkpoint.every_iterations") && b.kind != "moldft")
      throw std::runtime_error(
          b.name
          + ": checkpoint.every_iterations is not available for molresponse, "
            "which publishes at solve ends only");
    // S5, B3: the Dalton seed of a molresponse stage; seed.freq_tol is
    // required when the stage requests driven-response states (the legs it
    // seeds) and refused otherwise: a key that does nothing is refused.
    // D-112: a moldft stage takes seed.dalton_dir too (its orbitals start
    // from the Dalton molden.inp); it has no driven legs, so no
    // seed.freq_tol.
    if (seen.contains("seed.freq_tol") && b.kind != "molresponse")
      throw std::runtime_error(b.name
                               + ": seed.freq_tol is a molresponse stage's "
                                 "(it seeds driven legs); a moldft stage "
                                 "takes seed.dalton_dir alone");
    if (b.kind == "moldft" && seen.contains("seed.dalton_dir"))
      // Relative to where madqc was launched, as for a response stage.
      stage.seed_dalton_dir
          = std::filesystem::absolute(
                controls.get<std::string>("seed.dalton_dir"))
                .lexically_normal()
                .string();
    if (b.kind == "molresponse") {
      const auto requested
          = controls.get<std::vector<std::string>>("requested_properties");
      const bool driven = std::any_of(
          requested.begin(), requested.end(), [](const std::string& p) {
            return p == "polarizability" || p == "hyperpolarizability"
                   || p == "raman";
          });
      const bool seeded = seen.contains("seed.dalton_dir");
      if (seen.contains("seed.freq_tol") && !(seeded && driven))
        throw std::runtime_error(
            b.name
            + ": seed.freq_tol does nothing here (it needs "
              "seed.dalton_dir and a driven-response property)");
      if (seeded && driven && !seen.contains("seed.freq_tol"))
        throw std::runtime_error(
            b.name
            + ": seed.freq_tol is required with seed.dalton_dir and a "
              "driven-response property (no default tolerance, B3)");
      if (seen.contains("seed.freq_tol")
          && !(controls.get<double>("seed.freq_tol") > 0.0))
        throw std::runtime_error(b.name + ": seed.freq_tol must be positive");
      if (seeded) {
        const auto dir = controls.get<std::string>("seed.dalton_dir");
        // Relative to where madqc was launched, as the unnamed path resolves
        // dalton.dir; every stage runs from its own directory.
        stage.seed_dalton_dir
            = std::filesystem::absolute(dir).lexically_normal().string();
        stage.seed_freq_tol = controls.get<double>("seed.freq_tol");
        // The engine's spellings, which its seed import reads (S5: no
        // engine change).
        rp.set_user_defined_value<std::string>("dalton.dir",
                                               stage.seed_dalton_dir);
        if (driven)
          rp.set_user_defined_value("seed.freq_tol", stage.seed_freq_tol);
      }
    }
    stage.checkpoint_every_iterations
        = controls.get<int>("checkpoint.every_iterations");
    stage.checkpoint_keep_all = controls.get<bool>("checkpoint.keep_all");
    if (seen.contains("checkpoint.every_iterations")
        && stage.checkpoint_every_iterations < 1)
      throw std::runtime_error(b.name
                               + ": checkpoint.every_iterations must be at "
                                 "least 1");

    // Bind the molecule/reference and translate property requests.
    stage.properties =
        controls.get<std::vector<std::string>>("requested_properties");
    if (b.kind == "moldft") {
      const auto name = controls.get<std::string>("molecule");
      if (!molecules.count(name))
        throw std::runtime_error("unresolved molecule " + name + " for " +
                                 b.name);
      if (seen.count("reference"))
        throw std::runtime_error("moldft takes molecule, not reference");
      for (const auto &k : seen)
        for (const auto& p : named_response_properties())
          if (k.find(p + ".") == 0)
            throw std::runtime_error(p + " belongs to molresponse");
      stage.params.get<madness::Molecule>() = molecules.at(name);
      // Current named-path capability limit: HF energy/dipole only.
      // This is not a limitation of moldft's underlying XC implementation.
      if (cp.get<std::string>("xc") != "hf")
        throw std::runtime_error("named moldft currently supports xc hf");
      if (stage.properties.empty())
        throw std::runtime_error("requested_properties must be explicit for " +
                                 b.name);
      for (const auto &p : stage.properties)
        if (p != "energy" && p != "dipole")
          throw std::runtime_error("unsupported moldft property " + p);
      cp.set_user_defined_value(
          "dipole", std::find(stage.properties.begin(), stage.properties.end(),
                              "dipole") != stage.properties.end());
      cp.set_user_defined_value("restart", std::string("none"));
      cp.set_user_defined_value("save", true);
    } else {
      stage.reference = controls.get<std::string>("reference");
      const bool in_deck = deck_kinds.contains(stage.reference);
      if (stage.reference.empty()
          || (in_deck
              && (!indices.contains(stage.reference)
                  || stages.at(indices.at(stage.reference)).kind != "moldft")))
        throw std::runtime_error("unresolved earlier moldft reference " +
                                 stage.reference + " for " + b.name);
      if (seen.count("molecule"))
        throw std::runtime_error("molresponse takes reference, not molecule");
      // DESIGN A1: nothing is copied from the reference stage. Its molecule,
      // charge, L and final thresh/k come from its published generation
      // when this stage runs; cp keeps only this stage's protocol.
      const auto own_protocol = cp.protocol();
      // D-017: the property blocks, validated and mapped onto the engine.
      translate_response_properties(b.name, controls, stage, rp);
      // D-061: econv on a molresponse stage is the excitation energy change
      // target of its excited-state bundle, per root, reported under the
      // engine check (D8, B8). A driven response has no local energy, so
      // without excited_states it is refused as accuracy.energy_change is.
      if (seen.contains("econv")) {
        if (std::find(stage.properties.begin(), stage.properties.end(),
                      "excited_states")
            == stage.properties.end())
          throw std::runtime_error(
              b.name
              + ": econv is not applicable: a driven response has no local "
                "energy (see the reference stage's energy change); on a "
                "stage that requests excited_states it is the excitation "
                "energy change target");
        stage.excited_states.omega_target = cp.econv();
      }
      // D-099 (B3): stall.window and stall.ratio drive the bundle's plateau
      // detector (the engine's, in both modes); the dipole legs' detector
      // takes progress.* instead, so on a stage without excited_states the
      // keys would do nothing and are refused.
      for (const char* key : {"stall.window", "stall.ratio"})
        if (seen.contains(key)
            && std::find(stage.properties.begin(), stage.properties.end(),
                         "excited_states")
                   == stage.properties.end())
          throw std::runtime_error(
              b.name + ": " + key
              + " drives the excited-state bundle's plateau detector; on a "
                "stage without excited_states it does nothing (the dipole "
                "legs' stall rules are progress.window and "
                "progress.minimum_relative_improvement)");
      // Validate both controls before executing any stage, including reuse.
      // B8 (Q3, Q5), D-092: as for moldft, letting the written targets decide
      // needs a target for every loop it would otherwise leave without one.
      // A driven-response property needs one of the three FD targets (its
      // legs' loop stops on them); a stage that requests only excited_states
      // needs econv (the bundle's loop stops on the engine check in both
      // modes, D-092, and the written econv is the stage verdict's target);
      // a stage with both needs the FD target and may add econv. progress.*
      // needs the targets to decide.
      const bool targets_decide
          = !rp.get<bool>(madness::convergence::use_default_criteria_key);
      if (targets_decide) {
        const bool driven
            = std::any_of(stage.properties.begin(), stage.properties.end(),
                          [](const std::string& property) {
                            return property != "excited_states";
                          });
        if (driven && molresponse_v3::fd_requested_accuracy_targets(rp).empty())
          throw std::runtime_error(
              b.name + ": " + madness::convergence::use_default_criteria_key
              + " false needs at least one written target for the "
                "driven-response legs (dconv, accuracy.bsh_residual, "
                "accuracy.kinetic_trace_change)");
        if (!driven && !stage.excited_states.omega_target)
          throw std::runtime_error(
              b.name + ": " + madness::convergence::use_default_criteria_key
              + " false on a stage that requests only excited_states needs "
                "econv, the excitation energy change target");
      }
      if (!targets_decide && rp.is_user_defined("progress.window"))
        throw std::runtime_error(
            b.name
            + ": progress.* stall rules act on your written targets; set "
            + madness::convergence::use_default_criteria_key + " false");
      molresponse_v3::fd_requested_progress(rp);
      // Defaults must come from this stage, not the upstream SCF maxiter.
      rp.set_user_defined_value("maxiter", rp.get<std::size_t>("maxiter"));
      valid_protocol(own_protocol);
      std::set<std::string> protocol_keys;
      for (double t : own_protocol)
        if (!protocol_keys
                 .insert(molresponse_v3::protocol_key(t, protocol_k(t)))
                 .second)
          throw std::runtime_error("response protocol rungs collide at engine "
                                   "protocol-key precision");
      // Reject downstream tightening instead of silently refining the source.
      // An in-deck reference is checked here, before anything runs; every
      // reference is checked again against its published generation.
      if (in_deck) {
        const auto& reference_cp = stages.at(indices.at(stage.reference))
                                       .params.get<CalculationParameters>();
        const int source_k = reference_cp.k() == -1
                                 ? protocol_k(reference_cp.protocol().back())
                                 : reference_cp.k();
        if (own_protocol.back() < reference_cp.protocol().back()
            || protocol_k(own_protocol.back()) > source_k) {
          // D-102: the deck is the contract; say what each side wrote.
          std::ostringstream text;
          namespace cv = madness::convergence;
          text << "response stage " << b.name << ": protocol final step thresh "
               << cv::detail::scientific(own_protocol.back()) << " (k "
               << protocol_k(own_protocol.back())
               << ") is tighter than reference " << stage.reference
               << ", written at final step thresh "
               << cv::detail::scientific(reference_cp.protocol().back())
               << " (k " << source_k << ")"
               << saved_reference_rung(world, output_root, stage.reference)
               << "; write the protocol you mean";
          throw std::runtime_error(text.str());
        }
      }
    }

    if (b.kind == "moldft") {
      // B8 (Q3, Q5): letting the written targets decide needs at least one,
      // and the progress.* stall rules act on written targets, so they need
      // the targets to decide.
      const bool targets_decide
          = !cp.get<bool>(madness::convergence::use_default_criteria_key);
      if (targets_decide && madness::SCF::requested_targets(cp).empty())
        throw std::runtime_error(
            b.name + ": " + madness::convergence::use_default_criteria_key
            + " false needs at least one written target (dconv, econv, "
              "accuracy.bsh_residual, accuracy.kinetic_trace_change)");
      if (!targets_decide && cp.is_user_defined("progress.window"))
        throw std::runtime_error(
            b.name
            + ": progress.* stall rules act on your written targets; set "
            + madness::convergence::use_default_criteria_key + " false");
    }

    // Derive electron counts/defaults, then apply the shared stage checks.
    valid_protocol(cp.protocol());
    const auto per_rung = cp.get<std::vector<int>>("maxiter_per_rung");
    if (b.kind == "moldft" && !per_rung.empty()
        && (per_rung.size() != 1 && per_rung.size() != cp.protocol().size()))
      throw std::runtime_error("maxiter list for " + b.name
                               + " needs one number or one per protocol rung");
    for (int cap : per_rung)
      if (cap <= 0)
        throw std::runtime_error("maxiter values must be positive in "
                                 + b.name);
    if (!std::isfinite(cp.dconv()) || cp.dconv() <= 0 || cp.maxiter() <= 0)
      throw std::runtime_error(
          "positive finite dconv and positive maxiter required");
    stage.params.set_derived_values();
    // A response stage has no molecule of its own; its reference's shell is
    // checked against the published generation when it runs.
    if (b.kind == "moldft"
        && (!cp.spin_restricted() || cp.nalpha() <= 0
            || cp.nalpha() != cp.nbeta()))
      throw std::runtime_error(
          "named stages currently require a closed-shell reference");
    if (!std::isfinite(cp.econv()) || cp.econv() <= 0 ||
        !std::isfinite(cp.L()) || cp.L() <= 0)
      throw std::runtime_error("positive finite econv and L required");
    if (!std::isfinite(rp.dconv()) || rp.dconv() <= 0 || rp.maxiter() == 0)
      throw std::runtime_error(
          "positive finite response dconv and positive maxiter required");
    // Only a fully validated stage becomes available to later references.
    stage.written = seen;
    indices.emplace(b.name, stages.size());
    stages.push_back(std::move(stage));
  }
  if (stages.empty())
    throw std::runtime_error("named input has no application stages");
  return stages;
}
// A3, step 6 (S3c2): one generation for the skeleton to publish
// (publish_generation). The application gathers its publish sample, from
// which the AppSpec's manifest field assembles the manifest (S4c3, D-054),
// all but the fields every manifest shares (S4c1, set_shared_fields) and the
// archives; it also says what to check before anything is written and how
// its archives are written.
struct GenerationWrite {
  // What AppSpec::manifest reads, as JSON (D-045's pattern).
  nlohmann::json publish_sample;
  // position.iterations_this_run: the engine iterations of this run so far.
  std::size_t iterations_this_run = 0;
  // D-077 (S4d2): this run's final generation of the stage. The skeleton
  // forms the shared verdict from its state records before it is published
  // and sets the publish sample's `converged` to it.
  bool final = false;
  // Rank 0, before the generation directory exists; throws why not.
  std::function<void(const madness::checkpoint::Manifest&)> check_before_write;
  // Every rank: write the archives into `directory` and add them to the
  // manifest (rank 0's is the one published); returns this rank's error.
  std::function<std::string(const std::filesystem::path& directory,
                            madness::checkpoint::Manifest&)>
      write_archives;
};
// A3, D-054 (S4c2): the AppSpec's publish_now at one save boundary: whether
// it publishes, and, its second output, the rung record of the sample (the
// published manifest's position.rung): null for an application whose
// position has no rung (response).
struct PublishDecision {
  bool publish = false;
  nlohmann::json rung = nullptr;
};

// D-058 (S4c1): the names one stage's run gives its generations,
// g-<stage>-<run>-<n>, n counting from 1 in publish order. Generations
// written before S4c1 are named g-<run>-<n> (moldft) and are still read: no
// reader parses a name (readers follow `latest` or a manifest).
struct GenerationSequence {
  std::string stage, run_id;
  std::size_t n = 0;
  std::string next() {
    return "g-" + stage + "-" + run_id + "-" + std::to_string(++n);
  }
};

// S4c1: the manifest fields every application's generation shares, filled
// by the skeleton before anything checks or writes the manifest: the stage,
// the generation's name, the iterations of this run, the generation each
// property record belongs to, and the working orbitals' state id, which is
// the generation (A1: a consumer requires state_id == generation; moldft
// only, response has no working orbitals).
inline void set_shared_fields(madness::checkpoint::Manifest& m,
                              const std::string& stage,
                              const std::string& generation,
                              std::size_t iterations_this_run) {
  m.stage = stage;
  m.generation = generation;
  m.position["iterations_this_run"] = iterations_this_run;
  for (auto& property : m.properties)
    property["generation"] = generation;
  if (m.working_orbitals.is_object())
    m.working_orbitals["state_id"] = generation;
}

// The skeleton's publish, as the run binding hands it to a checkpoint class:
// returns the generation as published and read back.
using GenerationPublisher
    = std::function<madness::checkpoint::Manifest(GenerationWrite)>;

// S4a (D-056): the per-state record every ladder application writes in
// numerical.states, one per state (B6, skeleton spec 3.4). D-106: v2 names
// the record's keys after the science (protocol_step, stop_reason, the
// frequency-dependent kind, the excitation energy and kinetic trace
// changes). v1 records, and response's records written before S4a under
// its former name, are still read, through upgrade_state_record.
inline constexpr const char* state_schema = "madqc.state.v2";
inline constexpr const char* state_schema_v1 = "madqc.state.v1";
inline constexpr const char* response_state_schema = "madqc.response.state.v1";

// D-101, D-106: a frequency as a state identity prints it: the shortest
// decimal that reads back to the same double (a written 0.04 is "0.04"), so
// the identity a generation saved and the one a request asks for are the
// same string exactly when the doubles are equal.
[[nodiscard]] inline std::string frequency_text(double frequency) {
  std::array<char, 32> buffer{};
  const auto [end, error]
      = std::to_chars(buffer.data(), buffer.data() + buffer.size(), frequency);
  if (error != std::errc())
    throw std::runtime_error("frequency_text: cannot print the frequency");
  return {buffer.data(), end};
}

// D-101: a frequency-dependent response state's identity, the one form the
// saved records, the request and upgrade_state_record use.
[[nodiscard]] inline std::string
frequency_dependent_identity(const std::string& perturbation,
                             double frequency) {
  return perturbation + "|" + frequency_text(frequency);
}

// D-110: an excited-state bundle in the deck's words. The saved id keeps the
// engine's spelling, es:<tda|full>_n<roots>, "full" being rpa (es_node_id in
// calc_manager.hpp); a line madqc prints says "rpa with 2 roots". Any other
// identity is returned as it is.
[[nodiscard]] inline std::string
excited_state_words(const std::string& identity) {
  const std::string prefix = "es:";
  const auto n = identity.rfind("_n");
  if (identity.rfind(prefix, 0) != 0 || n == std::string::npos
      || n < prefix.size())
    return identity;
  const std::string type = identity.substr(prefix.size(), n - prefix.size());
  const std::string roots = identity.substr(n + 2);
  if ((type != "tda" && type != "full") || roots.empty()
      || roots.find_first_not_of("0123456789") != std::string::npos)
    return identity;
  return (type == "full" ? "rpa" : "tda") + std::string(" with ") + roots
         + (roots == "1" ? " root" : " roots");
}

// D-106: a state record in the current schema's keys. A v1 (or pre-S4a
// response) record gets the renamed keys: rung -> protocol_step,
// iterations_at_rung -> iterations_at_protocol_step, stop -> stop_reason,
// kind fd -> frequency_dependent_response, the bundle's omega_residual ->
// excitation_energy_change (measurement and root rows), the ground state's
// delta_kinetic_energy -> absolute_kinetic_trace_change; and a
// frequency-dependent state's identity in the canonical frequency form. A
// current record is returned as it is. The generation on disk is never
// rewritten (D-008): this is how it is read.
[[nodiscard]] inline nlohmann::json
upgrade_state_record(const nlohmann::json& record) {
  const auto schema = record.value("schema", "");
  if (schema == state_schema)
    return record;
  if (schema != state_schema_v1 && schema != response_state_schema)
    return record; // validate_state_records names an unknown schema
  auto upgraded = record;
  const auto rename
      = [](nlohmann::json& object, const char* from, const char* to) {
          if (object.is_object() && object.contains(from)) {
            object[to] = object.at(from);
            object.erase(from);
          }
        };
  rename(upgraded, "rung", "protocol_step");
  rename(upgraded, "iterations_at_rung", "iterations_at_protocol_step");
  rename(upgraded, "stop", "stop_reason");
  if (upgraded.value("kind", "") == "fd")
    upgraded["kind"] = "frequency_dependent_response";
  // Rows madqc names itself rename; a frequency-dependent state's rows are
  // the engine's own (convergence and saved-state measurements, raw
  // evidence) and keep its spelling.
  if (upgraded.value("kind", "") == "excited_states"
      && upgraded.contains("engine_measurements"))
    rename(upgraded["engine_measurements"], "omega_residual",
           "excitation_energy_change");
  if (upgraded.value("kind", "") == "ground_state")
    for (const char* part : {"engine_measurements", "saved_measurements"})
      if (upgraded.contains(part))
        rename(upgraded[part], "delta_kinetic_energy",
               "absolute_kinetic_trace_change");
  if (upgraded.contains("roots"))
    for (auto& root : upgraded["roots"])
      rename(root, "omega_residual", "excitation_energy_change");
  // A pre-S4a response record has no kind; its driven states carry a
  // perturbation and a numeric frequency (a bundle's frequency is null).
  const bool driven
      = upgraded.value("kind", "") == "frequency_dependent_response"
        || (!upgraded.contains("kind")
            && !upgraded.value("perturbation", "").empty());
  if (driven && upgraded.value("frequency", nlohmann::json()).is_number())
    upgraded["identity"] = frequency_dependent_identity(
        upgraded.at("perturbation").get<std::string>(),
        upgraded.at("frequency").get<double>());
  upgraded["schema"] = state_schema;
  return upgraded;
}

// D-106: every record of a generation in the current schema's keys.
[[nodiscard]] inline nlohmann::json
upgrade_state_records(const nlohmann::json& records) {
  auto upgraded = nlohmann::json::array();
  for (const auto& record : records)
    upgraded.push_back(upgrade_state_record(record));
  return upgraded;
}

// A2, D-052: the checks on a generation's state records, the same on the
// write side and the read side. A record of the current schema has a kind
// and a phase on every measurement row (D-059; null for a row with no
// sample, D-063); one written before S4a may lack both.
inline void validate_state_records(const nlohmann::json& records) {
  const auto require = [](bool ok, const std::string& why) {
    if (!ok)
      throw std::runtime_error("state record: " + why);
  };
  require(records.is_array(), "numerical.states is not a list");
  for (const auto& written : records) {
    require(written.is_object(), "not an object");
    const auto schema = written.value("schema", "");
    require(schema == state_schema || schema == state_schema_v1
                || schema == response_state_schema,
            "unknown schema '" + schema + "'");
    // A record written before S4a may lack a kind and phases (D-059).
    const bool current = schema != response_state_schema;
    const auto record = upgrade_state_record(written);
    require(record.contains("key") && record.at("key").is_string()
                && record.contains("identity")
                && record.at("identity").is_string(),
            "missing key or identity");
    const auto where = record.at("key").get<std::string>();
    require(record.contains("protocol_step")
                && record.at("protocol_step").contains("thresh")
                && record.at("protocol_step").contains("k"),
            where + ": no protocol_step");
    require(record.contains("iterations_at_protocol_step")
                && record.at("iterations_at_protocol_step").is_number_integer(),
            where + ": no iterations_at_protocol_step");
    require(record.contains("stop_reason")
                && record.at("stop_reason").is_string()
                && madness::convergence::to_string(
                       madness::convergence::stop_reason_from_string(
                           record.at("stop_reason").get<std::string>()))
                       == record.at("stop_reason").get<std::string>(),
            where + ": unknown stop_reason");
    require(record.contains("diverged") && record.at("diverged").is_boolean()
                && record.contains("decided_by")
                && record.at("decided_by").is_string(),
            where + ": no diverged or decided_by");
    for (const char* part : {"engine_measurements", "saved_measurements"})
      require(record.contains(part) && record.at(part).is_object(),
              where + ": no " + part);
    if (!current)
      continue;
    require(record.contains("kind") && record.at("kind").is_string(),
            where + ": no kind");
    for (const char* part : {"engine_measurements", "saved_measurements"})
      for (const auto& row : record.at(part).items())
        require(row.value().is_object() && row.value().contains("phase"),
                where + ": " + part + "/" + row.key() + " has no phase");
  }
}

// B6, S4a: a state record as the planner reads it, for every application.
// The response-only fields are read when present; `saved_measurements`
// renames the saved state's measurements to the names the written targets
// use (response: fd_accuracy_measurements), and is the identity when empty.
[[nodiscard]] inline madness::planning::SavedState saved_state_from_record(
    const nlohmann::json& written,
    const std::function<nlohmann::json(const nlohmann::json&)>&
        saved_measurements
    = {}) {
  const auto schema = written.value("schema", "");
  if (schema != state_schema && schema != state_schema_v1
      && schema != response_state_schema)
    throw std::runtime_error("state record: unknown schema '" + schema + "'");
  const auto record = upgrade_state_record(written);
  madness::planning::SavedState state;
  state.identity = record.at("identity").get<std::string>();
  state.reference_generation = record.value("reference_generation", "");
  state.perturbation = record.value("perturbation", "");
  // An excited-state bundle has no frequency (null); moldft none at all.
  state.frequency
      = record.contains("frequency") && record.at("frequency").is_number()
            ? record.at("frequency").get<double>()
            : 0.0;
  state.rung = {record.at("protocol_step").at("thresh").get<double>(),
                record.at("protocol_step").at("k").get<int>()};
  state.iterations_at_rung
      = record.at("iterations_at_protocol_step").get<int>();
  state.stop = madness::convergence::stop_reason_from_string(
      record.at("stop_reason").get<std::string>());
  state.diverged = record.at("diverged").get<bool>();
  state.engine_measurements = record.at("engine_measurements");
  state.saved_measurements
      = saved_measurements ? saved_measurements(record.at("saved_measurements"))
                           : record.at("saved_measurements");
  return state;
}

// B6, S4b (D-048): what the plan needs that the stage does not hold, for
// every application: moldft the physical block of the live engine, response
// its reference's generation.
struct PlanContext {
  std::string identity;
  std::string reference_generation;
};

// One state's plan (D-064): the planner's decision, the requested rung it
// starts at, the saved state's rung and a seed's frequency, which the run
// binding turns into its engine's input, and the engine-side names it needs
// (response: the engine's node id and the state key that resumes; moldft:
// none).
struct PlannedState {
  std::string identity;
  madness::planning::LadderPlan plan;
  madness::planning::Rung start;
  std::optional<madness::planning::Rung> own_rung;
  std::optional<double> seed_frequency;
  nlohmann::json fields = nlohmann::json::object();
};

// B6 (D-107): the part of a state's continuation record both applications
// write: the planner's action in words, the protocol step it starts at, the
// iterations already done there, the reason, the unmet targets, and the
// summary. Each application adds its own fields (moldft: planner; response:
// state and its node fields).
[[nodiscard]] inline nlohmann::json
continuation_record(const PlannedState& state,
                    const std::vector<madness::planning::Rung>& ladder) {
  namespace planning = madness::planning;
  static const std::map<planning::Action, const char*> words
      = {{planning::Action::load, "loaded"},
         {planning::Action::resume, "resumed"},
         {planning::Action::seed, "seeded"},
         {planning::Action::fresh, "fresh"}};
  const auto& plan = state.plan;
  return {{"action", words.at(plan.action)},
          {"protocol_step",
           {{"index", plan.start_rung + 1},
            {"thresh", state.start.thresh},
            {"k", state.start.k}}},
          {"iterations_done_at_protocol_step", plan.iterations_done},
          {"why", plan.why},
          {"unmet", plan.unmet},
          {"unmet_not_enforceable", plan.unmet_not_enforceable},
          {"summary", planning::describe(plan, ladder)}};
}

// The plans of a stage's states, in request order, and what --continue did
// and why (B6), the results.json record in the application's form.
struct StagePlans {
  nlohmann::json continuation = nlohmann::json::object();
  std::vector<PlannedState> states;
};

} // namespace named_stages
