#pragma once
#include <madness/chem/CCParameters.h>
#include <madness/chem/CalculationParameters.h>
#include <madness/chem/SCFInput.h>
#include <algorithm>
#include <madness/chem/ResponseParameters.hpp>
#include <madness/chem/TDHF.h>
#include <madness/chem/oep.h>
#include <madness/chem/pcm.h>
#include <madness/mra/QCCalculationParametersBase.h>

#include <type_traits>

using namespace madness;
using path = std::filesystem::path;

struct OptimizationParameters : public QCCalculationParametersBase {
  static constexpr const char *tag = "optimization";
  OptimizationParameters(const OptimizationParameters &other) = default;

  OptimizationParameters(World &world, const commandlineparser &parser)
      : OptimizationParameters() {
    read_input_and_commandline_options(world, parser, tag);
  }
  OptimizationParameters() {
    // No `method` here: the reference method comes from --wf=<scf|nemo>, and
    // --optimize says to optimize its geometry. Naming it twice invites the two
    // to disagree.
    initialize<int>("maxiter", 20, "maximum number of geometry steps");
    initialize<bool>("initial_hessian", false,
                     "compute inital hessian for optimization");
    initialize<std::string>("algopt", "bfgs", "hessian update used by MolOpt",
                            {"bfgs", "sr1"});
    initialize<double>("maxstep", 0.1,
                       "maximum step in any cartesian coordinate (a.u.)");
    // Gradient tolerances are independent of electronic density convergence.
    // The driver derives the energy/step defaults from the protocol; explicit
    // input values always take precedence.
    initialize<double>("etol", 1.e-5, "convergence: energy change");
    initialize<double>("gtol", 1.e-4, "convergence: maximum gradient element");
    initialize<double>("xtol", 1.e-4, "convergence: maximum cartesian step");
    initialize<double>("value_precision", 1.e-5, "assumed precision of the energy");
    initialize<double>("gradient_precision", 1.e-5,
                       "assumed precision of the gradient");
  }

  void set_protocol_defaults(double threshold) {
    set_derived_value("etol", std::max(1.0e-6, 2.0 * threshold));
    set_derived_value("xtol", 1.0e-3);
    set_derived_value("value_precision", std::max(1.0e-6, threshold));
  }

  std::string get_tag() const override { return std::string(tag); }

  using QCCalculationParametersBase::read_input_and_commandline_options;
  using QCCalculationParametersBase::print;

  [[nodiscard]] int get_maxiter() const { return get<int>("maxiter"); }
  [[nodiscard]] bool get_initial_hessian() const {
    return get<bool>("initial_hessian");
  }
  [[nodiscard]] std::string get_algopt() const {
    return get<std::string>("algopt");
  }
  [[nodiscard]] double get_maxstep() const { return get<double>("maxstep"); }
  [[nodiscard]] double get_etol() const { return get<double>("etol"); }
  [[nodiscard]] double get_gtol() const { return get<double>("gtol"); }
  [[nodiscard]] double get_xtol() const { return get<double>("xtol"); }
  [[nodiscard]] double get_value_precision() const {
    return get<double>("value_precision");
  }
  [[nodiscard]] double get_gradient_precision() const {
    return get<double>("gradient_precision");
  }
};

/// Deck-level `io` block (roadmap change 5): run-wide I/O configuration.
/// `backend` selects the restart-archive family every task that persists
/// MRA state writes (today: response; moldft & friends as they adopt it) —
/// it generalizes the response-only `response.hdf5` knob, which is kept as
/// a working alias. Future parallel-HDF5 knobs (collective I/O, per-dataset
/// compression) belong in this block, not per-task.
struct IOParameters : public QCCalculationParametersBase {
  static constexpr const char *tag = "io";
  IOParameters(const IOParameters &other) = default;

  IOParameters(World &world, const commandlineparser &parser)
      : IOParameters() {
    read_input_and_commandline_options(world, parser, tag);
  }
  IOParameters() {
    initialize<std::string>(
        "backend", "native",
        "restart-archive backend for tasks that support it: native "
        "(MADNESS parallel archives) or hdf5 (single .h5 blobs; requires a "
        "-DMADNESS_ENABLE_HDF5=ON build)",
        {"native", "hdf5"});
    initialize<std::string>(
        "dalton.dir", "",
        "run-wide seed directory: a DALTON run (loose RSPVEC + molden.inp, or a "
        "unique *.tar.gz) that seeds every stage of this calculation — the SCF "
        "(molden orbitals), the frequency-dependent response legs (XDIPLEN/... "
        "records) and the excited states (EXCITLAB records). The response-block "
        "`dalton.dir` is kept as a working alias; this one wins when both are set.");
  }

  std::string get_tag() const override { return std::string(tag); }
  // D-128: a directory is a path: its case is kept (the deck reader and the
  // stored text lowercase every other value).
  std::set<std::string> case_sensitive_keys() const override { return {"dalton.dir"}; }

  using QCCalculationParametersBase::read_input_and_commandline_options;

  void print() const {
    madness::print("------------IO Parameters---------------");
    madness::print("Backend: ", get<std::string>("backend"));
    madness::print("-------------------------------------------");
  }

  [[nodiscard]] std::string backend() const {
    return get<std::string>("backend");
  }
  [[nodiscard]] bool hdf5() const { return backend() == "hdf5"; }
  [[nodiscard]] std::string dalton_dir() const {
    return get<std::string>("dalton.dir");
  }
};

template <typename... Groups> class ParameterManager {
  std::tuple<Groups...> groups_;
  commandlineparser parser_;
  nlohmann::json all_input_json_;

  World &world_;

  // helper to invoke each group’s JSON export:
  template <typename G> void addGroupJson() {
    const auto &g = std::get<G>(groups_);
    auto j = g.to_json_if_precedence("defined");
    // if (world_.rank() == 0) {
    //   madness::print("Group: ", G::tag, " JSON: ", j.dump(4));
    // }
    if (!j.is_null())
      all_input_json_[G::tag] = j;
  }

public:
  ParameterManager() : world_(World::get_default()) {}

  /// Empty independent groups, populated by a scoped input reader.
  explicit ParameterManager(World &world) : world_(world) {}

  /// "Master" ctor: takes any single intput file, JSON or plain-text
  // ParameterManager(World &w, const path &filename) : world_(w) {
  ParameterManager(World &w, const commandlineparser &parser)
      : parser_(parser), world_(w) {
    // parser_.set_keyval("input", filename);
    //
    const bool user_defined_prefix = parser_.key_exists("user_defined_prefix");
    std::string inputfile = parser.value("input");
    std::string prefix_from_input = commandlineparser::remove_extension(
        commandlineparser::base_name(inputfile));

    const path &filename = parser_.value("input");

    std::string json_text;
    if (world_.rank() == 0) {
      std::ifstream input(filename);
      std::ostringstream contents; contents << input.rdbuf();
      if (json::accept(contents.str())) json_text = contents.str();
    }
    world_.gop.broadcast_serializable(json_text, 0);
    if (!json_text.empty()) {
      initFromJson(json::parse(json_text));
      // invoke each group’s JSON parser:
    } else {
      // plain-text file
      initFromText(filename);
    }

    // Prefix precedence, highest first: --prefix on the command line, then a
    // `prefix` in the deck's dft block, then the input file's base name, then
    // the "mad" default. This must be layered in AFTER the deck is parsed:
    // set_derived_value does not override an already-defined value, so a
    // deck-level prefix survives while the file-name fallback still applies
    // when the deck is silent. CalculationParameters::prefix() is the single
    // source of truth (see prefix() below) — computing it up front into a
    // separate member shadowed the deck value.
    auto &cparam = get<CalculationParameters>();
    if (user_defined_prefix) {
      cparam.set_user_defined_value("prefix", parser_.value_raw("prefix"));
    } else if (prefix_from_input != "input") {
      cparam.set_derived_value("prefix", prefix_from_input);
    }

    set_derived_values();
  }

  /// here comes some logic for the calculation, e.g. the number of electrons
  /// derived from the molecule
  void set_derived_values() {
    auto &cparam = this->get<CalculationParameters>();
    cparam.set_derived_values(this->get<Molecule>());
    // Derived blocks inherit the numerical knobs of the dft block unless the
    // deck sets them in the derived block itself (set_derived_value never
    // overrides a user-defined value). The response `protocol` ladder already
    // comes from the dft block (there is no response-level protocol key; one in
    // the response block is reported as an unknown parameter and ignored);
    // dconv and maxiter used to be independent defaults (1e-6 / 25) that
    // silently disagreed with the dft block.
    auto &rp = this->get<ResponseParameters>();
    // Response dconv: when the deck does not set it, derive it from the finest
    // protocol threshold as 100 x thresh. The solvers gate at 5 x max(thresh,
    // dconv) (convergence_policy.hpp). The factor is empirical: the BSH
    // amplitude residual of a converged FD leg plateaus at ~2-6e-5 for
    // thresh 1e-6/k8 and at ~4-8e-6 for thresh 1e-8/k10 (closeout series,
    // 2026-09-11; the density residual keeps falling to ~1e-7), so dconv =
    // thresh (gate 5e-8) is unattainable and burned 60 iterations per leg,
    // while 100 x thresh gives the validated 1e-4 at 1e-6 and a gate of 5e-6
    // at 1e-8, just above the floor. A tighter series sets dconv explicitly.
    const auto proto = cparam.protocol();
    const double finest = proto.empty() ? 1.0e-6 : *std::min_element(proto.begin(), proto.end());
    rp.set_derived_value("dconv", 100.0 * finest);
    rp.set_derived_value("maxiter", static_cast<size_t>(cparam.maxiter()));
  }

  /// dump out the merged JSON
  [[nodiscard]] const nlohmann::json &getAllInputJson() const {
    return all_input_json_;
  }

  /// access a particular group by type:
  template <typename G> const G &get() const { return std::get<G>(groups_); }
  template <typename G> G &get() { return std::get<G>(groups_); }
  template <typename G> void set(const G &g) { std::get<G>(groups_) = g; }
  /// The run-level output prefix. Delegates to CalculationParameters so a
  /// `prefix` set in the deck (or later via set_derived_value) is honored by
  /// everything that names files or directories from it.
  std::string prefix() const { return get<CalculationParameters>().prefix(); }

  /// pretty-print everything
  void print_all() const { (print_group_if_defined<Groups>(), ...); }

private:
  [[nodiscard]] bool shared_scf_owner() const {
    const auto workflow = parser_.key_exists("workflow") ? parser_.value("workflow") :
        parser_.key_exists("wf") ? parser_.value("wf") :
        parser_.key_exists("w") ? parser_.value("w") : "scf";
    return workflow == "scf" || workflow == "response";
  }
  template <class G> G read_group() {
    if constexpr (std::is_same_v<G, CalculationParameters>) {
      if (shared_scf_owner()) return read_scf_parameters(world_, parser_);
    }
    return G(world_, parser_);
  }
  void initFromJson(const nlohmann::json &j) {
    (
        [&] {
          std::string input_error;
          try {
          if constexpr (std::is_same_v<Groups, CalculationParameters> ||
                        std::is_same_v<Groups, ResponseParameters>) {
            auto &group = std::get<Groups>(groups_);
            if constexpr (std::is_same_v<Groups, CalculationParameters>)
              if (shared_scf_owner()) group.enable_scf_inputs();
            if (j.contains(Groups::tag))
              for (const auto &[key, value] : j.at(Groups::tag).items())
                group.validate_input_value(key, value.is_string() ? value.template get<std::string>() : value.dump());
          }
          if (j.contains(Groups::tag)) {
            std::get<Groups>(groups_).from_json(j.at(Groups::tag));
          }
          } catch (const std::exception& error) { input_error = error.what(); }
          int failing_rank = input_error.empty() ? world_.size() : world_.rank();
          world_.gop.min(failing_rank);
          if (failing_rank != world_.size()) {
            world_.gop.broadcast_serializable(input_error, failing_rank);
            throw std::runtime_error(input_error);
          }
          // Layer command-line overrides (--group=key=val) ON TOP of the JSON,
          // exactly as initFromText gets them through the file+CLI parser.
          // Without this, JSON-format decks silently ignored CLI overrides
          // (review MED). Applies to defaults too when the tag is absent.
          // Only the QCCalculationParametersBase-derived groups have the
          // override machinery — Molecule (also in the tuple) does not, so
          // guard with if constexpr.
          if constexpr (std::is_base_of_v<QCCalculationParametersBase, Groups>) {
            std::get<Groups>(groups_).read_commandline_options(world_, parser_,
                                                               Groups::tag);
          } else if constexpr (std::is_same_v<Groups, Molecule>) {
            // Geometry convenience options override JSON just as they override
            // text decks. Molecule is not a parameter group, so it needs its
            // own collective reader rather than read_commandline_options.
            if (parser_.key_exists("geometry") || parser_.key_exists("molecule"))
              std::get<Groups>(groups_) = Molecule(world_, parser_);
            else if (!j.contains(Groups::tag))
              throw std::runtime_error("missing molecular geometry in JSON input");
          }
        }(),
        ...);

    all_input_json_ = j;
  }
  // 1) read from a plain-text “.inp” file
  void initFromText(const path &filename) {
    // parser_.set_keyval("input", filename);
    // invoke each group’s file+CLI parser:
    ((void)(std::get<Groups>(groups_) = read_group<Groups>()), ...);

    // collect JSON for any defined keys:
    ((void)addGroupJson<Groups>(), ...);
  }

  template <typename G> void print_group_if_defined() const {
    const auto &g = std::get<G>(groups_);
    // grab only the user-defined values:
    auto j = g.to_json_if_precedence("defined");
    // json.empty() is true if no user-defined values
    if (!j.empty()) {
      g.print();
    }
  }
  static bool is_json_file(const path &f) {
    std::ifstream input_file_stream(f);
    bool is_json = json::accept(input_file_stream);
    input_file_stream.close();
    return is_json;
  }
  static json read_json_file(const path &input_file) {
    std::ifstream input_file_stream(input_file);
    auto j = json::parse(input_file_stream);
    input_file_stream.close();
    return j;
  }
};

// Define a concrete aliased ParameterManager type
using Params =
    ParameterManager<CalculationParameters, ResponseParameters,
                     Nemo::NemoCalculationParameters, OptimizationParameters,
                     OEP_Parameters, CCParameters, TDHFParameters, PCMParameters,
                     Molecule, IOParameters>;
