#pragma once

// Load MRA function declarations before ResponseProtocol includes its operators.
#include <madness/chem/ParameterManager.hpp>
#include <madness/chem/SCF.h>

#include <apps/molresponse/ResponseProtocol.hpp>
#include <apps/molresponse/solvers/fd_explicit_accuracy.hpp>
#include <apps/molresponse/solvers/response_metadata.hpp>
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
  };
  ExcitedStatesRequest excited_states;
  // A3: checkpoint cadence and retention (moldft); 0 is "decision points only".
  int checkpoint_every_iterations = 0;
  bool checkpoint_keep_all = false;
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

// Build and validate the whole sequence before any calculation starts.
// Names enter these lookup tables in file order, so references must point back.
inline std::vector<Stage> parse(madness::World &world,
                                const std::string &text) {
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
         "print_level"};
  // B8: dconv is the induced-density-change target of a molresponse stage; a
  // driven response has no local energy, so no energy target.
  const std::map<std::string, std::string> replaced_response_keys
      = {{"accuracy.density_change", "dconv (the density-change target)"},
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
      else if (b.kind == "molresponse" && response_keys.count(k) &&
               rp.parameter_exists(k))
        read_value(world, rp, line);
      else
        fail(line, "unsupported parameter " + k + " in " + b.kind);
    }

    if (b.kind == "molresponse") {
      // Validate both controls before executing any stage, including reuse.
      // B8 (Q3, Q5): as for moldft, letting the written targets decide needs
      // at least one, and progress.* needs the targets to decide.
      const bool targets_decide
          = !rp.get<bool>(madness::convergence::use_default_criteria_key);
      if (targets_decide
          && molresponse_v3::fd_requested_accuracy_targets(rp).empty())
        throw std::runtime_error(
            b.name + ": " + madness::convergence::use_default_criteria_key
            + " false needs at least one written target (dconv, "
              "accuracy.bsh_residual, accuracy.kinetic_trace_change)");
      if (!targets_decide && rp.is_user_defined("progress.window"))
        throw std::runtime_error(
            b.name
            + ": progress.* stall rules act on your written targets; set "
            + madness::convergence::use_default_criteria_key + " false");
      molresponse_v3::fd_requested_progress(rp);
    }

    // A3: a stage publishes inside a solve only where it has a save boundary
    // there; molresponse publishes at solve ends only (DESIGN §1.4).
    // checkpoint.keep_all is every application's (D-070).
    if (seen.contains("checkpoint.every_iterations") && b.kind != "moldft")
      throw std::runtime_error(
          b.name
          + ": checkpoint.every_iterations is not available for molresponse, "
            "which publishes at solve ends only");
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
            || protocol_k(own_protocol.back()) > source_k)
          throw std::runtime_error(
              "response protocol is tighter than reference " + stage.reference);
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
// numerical.states, one per state (B6, skeleton spec 3.4). Response's records
// written before S4a carry its former name and are still read.
inline constexpr const char* state_schema = "madqc.state.v1";
inline constexpr const char* response_state_schema = "madqc.response.state.v1";

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
  for (const auto& record : records) {
    require(record.is_object(), "not an object");
    const auto schema = record.value("schema", "");
    const bool current = schema == state_schema;
    require(current || schema == response_state_schema,
            "unknown schema '" + schema + "'");
    require(record.contains("key") && record.at("key").is_string()
                && record.contains("identity")
                && record.at("identity").is_string(),
            "missing key or identity");
    const auto where = record.at("key").get<std::string>();
    require(record.contains("rung") && record.at("rung").contains("thresh")
                && record.at("rung").contains("k"),
            where + ": no rung");
    require(record.contains("iterations_at_rung")
                && record.at("iterations_at_rung").is_number_integer(),
            where + ": no iterations_at_rung");
    require(record.contains("stop") && record.at("stop").is_string()
                && madness::convergence::to_string(
                       madness::convergence::stop_reason_from_string(
                           record.at("stop").get<std::string>()))
                       == record.at("stop").get<std::string>(),
            where + ": unknown stop");
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
    const nlohmann::json& record,
    const std::function<nlohmann::json(const nlohmann::json&)>&
        saved_measurements
    = {}) {
  const auto schema = record.value("schema", "");
  if (schema != state_schema && schema != response_state_schema)
    throw std::runtime_error("state record: unknown schema '" + schema + "'");
  madness::planning::SavedState state;
  state.identity = record.at("identity").get<std::string>();
  state.reference_generation = record.value("reference_generation", "");
  state.perturbation = record.value("perturbation", "");
  // An excited-state bundle has no frequency (null); moldft none at all.
  state.frequency
      = record.contains("frequency") && record.at("frequency").is_number()
            ? record.at("frequency").get<double>()
            : 0.0;
  state.rung = {record.at("rung").at("thresh").get<double>(),
                record.at("rung").at("k").get<int>()};
  state.iterations_at_rung = record.at("iterations_at_rung").get<int>();
  state.stop = madness::convergence::stop_reason_from_string(
      record.at("stop").get<std::string>());
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

// The plans of a stage's states, in request order, and what --continue did
// and why (B6), the results.json record in the application's form.
struct StagePlans {
  nlohmann::json continuation = nlohmann::json::object();
  std::vector<PlannedState> states;
};

} // namespace named_stages
