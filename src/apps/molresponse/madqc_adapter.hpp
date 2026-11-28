#ifndef MOLRESPONSE_V3_MADQC_ADAPTER_HPP
#define MOLRESPONSE_V3_MADQC_ADAPTER_HPP

#include "solvers/shared_response_input.hpp"
#include "solvers/fd_explicit_accuracy.hpp"

// -----------------------------------------------------------------------------
// molresponse_v3_lib — madqc adapter (doc 16 R3 / doc 05 seam).
//
// Satisfies the duck-typed interface ResponseApplication<Library> expects
// (Applications.hpp): `Library::label()` + `Library::run_response(world, params,
// scf, outdir) -> Results{metadata, properties, vibrational_analysis,
// raman_spectra}`. So `ResponseApplication<molresponse_v3_lib>` runs the v3
// pipeline through the same madqc workflow path as v2's molresponse_lib —
// enabling a SAME-INPUT calc_info.json parity check (engine = v2 vs v3).
//
// It builds a v3 GroundState from the moldft restart ARCHIVE (resolved from
// scf_calc->work_dir, like v2's make_ground_context — NOT from the in-memory
// SCF, whose MOs may be unloaded on a restart-in-place run), maps the
// ResponseParameters input deck to a v3 ResponsePlan, and calls
// run_response_with_ground.
//
// SCOPE: multi-property mapping (R3b) — polarizability + hyperpolarizability +
// single-component raman + resonant/excited, deduped via merge_plans.
// This header lives in the v3 app and is included by madqc.cpp (the app links
// both MADchem and v3's GroundState.cpp); MADchem/WorkflowBuilders must NOT
// include it (that would be a circular library dependency — the engine is
// selected in madqc.cpp instead).
// -----------------------------------------------------------------------------

#include <apps/molresponse/deck_settings.hpp>
#include <apps/molresponse/orchestrator/response_workflow.hpp>
#include <apps/molresponse/solvers/dalton_import.hpp>
#include <apps/molresponse/solvers/dalton_gs_seed.hpp>   // GS seed from molden (dalton.dir)

#include <madness/chem/CalculationParameters.h>
#include <madness/chem/ParameterManager.hpp>   // Params
#include <madness/chem/ResponseParameters.hpp>
#include <madness/chem/SCF.h>
#include <madness/chem/SCFRun.h>
#include <nlohmann/json.hpp>
#include <madness/world/worldprofile.h>

#include <algorithm>
#include <filesystem>
#include <memory>
#include <vector>

namespace molresponse_v3 {

/// Run-wide seed directory: `io.dalton.dir` (ParameterManager IOParameters)
/// wins; the response-block `dalton.dir` is the alias. Same rule for HDF5:
/// `io.backend hdf5` or the alias `response.hdf5`.
template <typename Params>
inline std::string effective_dalton_dir(const Params &params) {
  const std::string io_dir = params.template get<IOParameters>().dalton_dir();
  if (!io_dir.empty()) return io_dir;
  return params.template get<ResponseParameters>().dalton_dir();
}
template <typename Params>
inline bool effective_hdf5(const Params &params) {
  return params.template get<IOParameters>().hdf5() ||
         params.template get<ResponseParameters>().hdf5();
}

// ---------------------------------------------------------------------------
// GROUND-STATE SEED HOOK (2026-09-09, seeding showcase). Installed by madqc on
// the SCF application (SCFApplication::set_pre_run_hook) when the deck carries
// `dalton.dir`. Runs collectively inside the SCF work dir right before the
// engine plans its restart:
//   * no-op when an archive <prefix>.restartdata* already exists (a real
//     restart always wins over a seed) or when `restart` is not auto/iterate;
//   * locates molden.inp in dalton.dir (same resolver as the FD/ES import),
//     fingerprints its geometry against the deck's Molecule (hard error);
//   * n_occ = (sum Z - charge)/2 (closed shell), L = deck `l`,
//     thresh = first protocol rung; writes <prefix>.restartdata via
//     write_gs_seed_from_molden, so RestartPlan(auto) -> restartdata/iterate.
// The projection happens at the deck's box; SCF re-sets its own
// FunctionDefaults afterwards, and load_mos re-projects k/thresh as needed.
// ---------------------------------------------------------------------------

/// D-096: the deck's print_level is moldft's scale (0 none, 1 final energy,
/// 2 iterations, 3 timings, the default, 10 debug); the engine's levels are
/// Silent, Normal, Verbose, Debug. Mapped as Drivers.hpp maps the scale for
/// the optimizer (3 -> 1): an unwritten named print_level runs the engine at
/// Normal, one row per iteration; the unnamed decks writing 1 are unchanged;
/// 4 to 9 add the channel detail (Verbose); 10 and above are Debug. Before
/// this the value was cast onto the enum, so the default ran at Debug.
/// Nothing but printing reads the level.
[[nodiscard]] inline PrintLevel engine_print_level(int deck_print_level) {
  if (deck_print_level <= 0)
    return PrintLevel::Silent;
  if (deck_print_level <= 3)
    return PrintLevel::Normal;
  if (deck_print_level <= 9)
    return PrintLevel::Verbose;
  return PrintLevel::Debug;
}

inline void seed_gs_from_dalton_dir(World &world, const Params &params,
                                    const std::filesystem::path &workdir,
                                    const std::string &dalton_dir) {
  namespace fs = std::filesystem;
  const auto &cp  = params.get<CalculationParameters>();
  const auto &mol = params.get<Molecule>();
  const std::string prefix = cp.prefix();
  const std::string mode   = cp.restart();
  if (mode != "auto" && mode != "iterate") {
    if (world.rank() == 0)
      print("[DALTON-SEED] GS: restart =", mode, "-> not seeding the ground state");
    return;
  }
  // Existing archive wins (rank 0 decides, collective broadcast).
  int have = 0;
  if (world.rank() == 0) {
    for (const auto &e : fs::directory_iterator(workdir)) {
      const std::string n = e.path().filename().string();
      if (n.rfind(prefix + ".restartdata", 0) == 0) { have = 1; break; }
    }
  }
  world.gop.broadcast(have, 0);
  if (have) {
    if (world.rank() == 0)
      print("[DALTON-SEED] GS:", prefix + ".restartdata*", "already present in",
            workdir.string(), "-> restart from it, not from the DALTON seed");
    return;
  }
  // Locate + fingerprint (rank 0), broadcast the molden path or the error.
  std::string molden, err, report;
  if (world.rank() == 0) {
    try {
      auto m = locate_dalton_dir(dalton_dir, (workdir / "dalton_import").string(),
                                 "", "", "");
      auto check = fingerprint_dalton_geometry(m, mol, 1e-4);
      report = check.report;
      if (!check.ok)
        throw std::runtime_error(
            "dalton import (GS seed): GEOMETRY FINGERPRINT MISMATCH\n" + check.report);
      molden = m.molden_path;
    } catch (const std::exception &ex) { err = ex.what(); }
  }
  world.gop.broadcast_serializable(err, 0);
  if (!err.empty()) throw std::runtime_error(err);
  world.gop.broadcast_serializable(molden, 0);
  world.gop.broadcast_serializable(report, 0);

  const double Z     = mol.total_nuclear_charge();
  const double nelec = Z - cp.charge();
  const long   ne    = std::lround(nelec);
  if (std::abs(nelec - static_cast<double>(ne)) > 1e-6 || ne % 2 != 0)
    throw std::runtime_error("dalton import (GS seed): closed-shell seed needs an even "
                             "electron count, got " + std::to_string(nelec));
  GsSeedOptions opt;
  opt.L        = cp.L();
  opt.thresh   = cp.protocol().empty() ? 1e-4 : cp.protocol().front();
  opt.xc       = cp.get<std::string>("xc");
  opt.localize = cp.get<std::string>("localize");
  opt.extra_prefixes = {prefix + ".gs_seed"};   // preserved copy (save_mos overwrites <prefix>.restartdata)
  opt.active_molecule = &mol;   // RestartPlan matches geometry at 1e-8 and eprec exactly
  if (world.rank() == 0) {
    print("[DALTON-SEED] GS: metadata molecule (stamped) eprec =", mol.parameters.eprec());
    for (std::size_t i = 0; i < mol.natom(); ++i) {
      const auto at = mol.get_atom(i);
      printf("[DALTON-SEED] GS:   atom %zu  Z=%d  %.10f %.10f %.10f\n", i, at.atomic_number, at.x, at.y, at.z);
    }
  }
  if (world.rank() == 0) {
    print("[DALTON-SEED] GS: seeding", prefix + ".restartdata", "from", molden,
          " n_occ =", ne / 2, " L =", opt.L, " thresh =", opt.thresh);
    print(report);
  }
  auto rep = write_gs_seed_from_molden(world, molden, static_cast<int>(ne / 2), prefix, opt);
  if (world.rank() == 0) {
    nlohmann::json j;
    j["seed"] = "dalton_import"; j["stage"] = "ground_state";
    j["molden"] = molden; j["dalton_dir"] = dalton_dir;
    j["n_occ"] = rep.n_occ; j["n_ao"] = rep.n_ao; j["n_mo"] = rep.n_mo;
    j["L"] = opt.L; j["thresh"] = opt.thresh;
    j["max_offdiag_pre_loewdin"] = rep.max_offdiag_pre;
    j["archive"] = rep.archive;
    j["preserved_copy"] = prefix + ".gs_seed.restartdata";
    std::ofstream out((workdir / (prefix + ".gs_seed.json")).string());
    out << j.dump(2) << "\n";
  }
}
} // namespace molresponse_v3

/// Global namespace (mirrors `molresponse_lib`) so madqc can write
/// `ResponseApplication<molresponse_v3_lib>`.

struct molresponse_v3_lib {
  /// Output subdir name + the interface ResponseApplication reads.
  static const char *label() { return "molresponse"; }

  /// Structured result returned to ResponseApplication (→ calc_info.json).
  /// What the response task adds to calc_info.json: the stage times
  /// ({"v3_timing": {plan_build, solve, assemble}}), for the wall-time lines at
  /// the end of the log. Every value is in results.json.
  struct Results {
    nlohmann::json metadata;
  };

  /// The reference is the live SCF of the same process, and the response
  /// runs in a session on that SCF's checkpoint (run_response_on_reference).
  /// A session needs the workflow to own the ground-state input and no
  /// separate reference prefix; any other case is refused.
  inline static Results
  run_response(madness::World &world, const Params &params,
               const std::shared_ptr<madness::SCF> &scf_calc,
               const std::filesystem::path &outdir,
               const std::filesystem::path &reference_prefix = {}) {
    const bool shared = reference_prefix.empty() && scf_calc->param.scf_input_owner;
    MADNESS_CHECK_THROW(shared,
        "response requires the workflow's own SCF reference: a separate reference "
        "prefix, or a workflow that does not own the ground-state input, is not supported");
    MADNESS_CHECK_THROW(scf_calc->last_calculation_outcome && scf_calc->last_calculation_outcome->complete,
                       "response requires a qualified SCF reference");
    const auto reference = scf_ground_reference(world, *scf_calc);
    const molresponse_v3::response_shared::GroundSource ground{reference.archive, reference.protocol, reference.identity};
    auto source_prefix = std::filesystem::path(reference.archive);
    source_prefix.replace_extension(); // the native loader appends .restartdata
    return run_response_on_reference(world, params, scf_calc->molecule, outdir,
                                     source_prefix, {},
                                     &scf_calc->param, &scf_calc->pcm_param,
                                     ground, scf_checkpoint_directory(*scf_calc));
  }

  /// DESIGN A1: the response on a saved reference alone: its orbital archive
  /// (`reference_prefix` + ".restartdata") and its molecule. Named stages call
  /// this with both read from the reference's published generation, and on
  /// --continue with the planner's decision per FD state (`state_plans`).
  inline static Results run_response_on_reference(
      madness::World& world, const Params& params,
      const madness::Molecule& reference_molecule,
      const std::filesystem::path& outdir,
      const std::filesystem::path& reference_prefix,
      const molresponse_v3::StatePlans& state_plans = {},
      const madness::CalculationParameters* verified_reference = nullptr,
      const madness::PCMParameters* verified_pcm = nullptr,
      const std::optional<molresponse_v3::response_shared::GroundSource>& shared_source = {},
      const std::filesystem::path& shared_root = {}) {
    using namespace madness;
    using namespace molresponse_v3;

    const auto &cp = params.get<CalculationParameters>();
    const auto &rp = params.get<ResponseParameters>();
    madness::iteration::TimingCollector<> response_timing;
    response_timing.begin("run",{{"clock_scope","driver_rank"}});
    response_timing.begin("response",{{"calculation_id","response"}});
    response_timing.begin("reference_and_input_setup");
    const auto shared_request = shared_response_request(rp);
    std::vector<double> protocol = cp.protocol();
    MADNESS_CHECK(!protocol.empty());

    // Explicit imported guesses use the user-requested protocol ladder.
    const std::string dalton_dir = effective_dalton_dir(params);
    if (shared_request) {
      // Check the written request before the planner can filter unknown values
      // and before loading/preparing any response reference or operators.
      const auto directions=rp.dipole_directions();
      const auto frequencies=rp.dipole_frequencies();
      const auto properties=rp.requested_properties();
      MADNESS_CHECK_THROW(std::all_of(properties.begin(),properties.end(),[](const auto& property){
            return property=="polarizability" || property=="hyperpolarizability" || property=="raman";}) &&
          !directions.empty() && directions.find_first_not_of("xyzXYZ")==std::string::npos &&
          !frequencies.empty() &&
          shared_source && !shared_root.empty() &&
          verified_reference,
          "shared response requires supported properties, valid dipole directions and frequencies");
      MADNESS_CHECK_THROW(!rp.excited_enable() || (verified_reference->spin_restricted() && rp.excited_num_states()>0),
          "shared excited response requires a closed-shell reference and positive root count");
      MADNESS_CHECK_THROW(!effective_hdf5(params),
          "shared response restart uses native checkpoints; HDF5 restart is not supported");
      std::set<std::string> frequency_keys,protocol_keys;
      for (double frequency:frequencies)
        MADNESS_CHECK_THROW(std::isfinite(frequency) &&
            frequency_keys.insert(freq_key(frequency)).second,
            "shared response frequencies must be finite and have distinct checkpoint identities");
      if (rp.nuclear()) {
        const auto nuclear_directions=rp.nuclear_directions();
        const auto nuclear_frequencies=rp.nuclear_frequencies();
        MADNESS_CHECK_THROW(verified_reference->spin_restricted(),
            "open-shell nuclear response source is not implemented");
        MADNESS_CHECK_THROW(!nuclear_directions.empty() &&
            nuclear_directions.find_first_not_of("xyzXYZ")==std::string::npos && !nuclear_frequencies.empty(),
            "nuclear response requires Cartesian directions and frequencies");
        std::set<std::string> nuclear_keys;
        for (double frequency:nuclear_frequencies)
          MADNESS_CHECK_THROW(std::isfinite(frequency) &&
              nuclear_keys.insert(freq_key(frequency)).second,
              "nuclear response frequencies must be finite and have distinct checkpoint identities");
      }
      double previous=std::numeric_limits<double>::infinity();
      for (double threshold:protocol) {
        const int order=default_k_for_thresh(threshold);
        MADNESS_CHECK_THROW(std::isfinite(threshold) && threshold>0 && threshold<previous &&
            protocol_keys.insert(protocol_key(threshold,order)).second,
            "shared response protocols must be positive, strictly decreasing and have distinct checkpoint identities");
        MADNESS_CHECK_THROW(cp.k()==-1 || cp.k()==order,
            "shared response fixed order differs from its scheduler representation");
        previous=threshold;
      }
      (void)madness::iteration::resolve(*shared_request,response_shared::model(),
          {protocol.back(),default_k_for_thresh(protocol.back())},true,"driven response");
    }

    // 1. Ground state from the moldft ARCHIVE (not the in-memory SCF). On a
    //    restart-in-place run madqc validates the SCF as "Ok" and never loads the
    //    MOs into memory (lib_.calc() just constructs the SCF), so scf_calc->amo is
    //    empty — building from the live SCF then segfaults in build_fock_matrices.
    //    Loading from the archive (exactly like v2's make_ground_context) avoids
    //    that AND lets prepare() reproject pristine MOs on each protocol climb.
    //
    //    Resolve the moldft work dir RELATIVE to the response `outdir`, exactly
    //    like v2's make_ground_context (MolresponseLib.hpp ~1154) and the
    //    CC2/TDHF/OEP applications. scf_calc->work_dir is stored relative to the
    //    top calc dir, but ResponseApplication::run has already chdir'd (ScopedCWD)
    //    into `outdir` (the response task dir), so using work_dir raw makes the
    //    archive lookup resolve against the wrong cwd → "could not find file:
    //    <work_dir>/<prefix>.restartdata" on multi-node madqc runs.
    namespace fs = std::filesystem;
    const double L = cp.L();
    set_response_protocol(world, L, protocol.front());
    // Named stages supply absolute paths and never change process cwd.
    const fs::path& source_prefix = reference_prefix;
    const std::string archive = source_prefix.string() + ".restartdata";
    const std::string fock_json = source_prefix.string() + ".fock.json"; // in the node archive
    GroundState gs
        = GroundState::from_archive(world, archive, reference_molecule, verified_reference, verified_pcm);
    // The node's Fock file is read only at the level its description records.
    if (shared_source) gs.set_fock_level(shared_source->protocol);
    const double thresh = FunctionDefaults<3>::get_thresh();
    auto coulop = poperatorT(
        CoulombOperatorPtr(world, gs.params().lo(), 0.001 * thresh));
    gs.prepare(world, 0.001 * thresh, coulop, fock_json);

    // 2. Map the input deck → a Tier-A plan (R3b). requested_properties +
    //    beta.*/excited.* knobs select which ResponsePropertyRequests to build;
    //    merge_plans dedupes shared FD states across them.
    std::vector<char> axes;
    for (char c : rp.dipole_directions()) {
      const char l = static_cast<char>(std::tolower(c));
      if (l == 'x' || l == 'y' || l == 'z') axes.push_back(l);
    }
    if (axes.empty()) axes = {'x', 'y', 'z'};
    const std::vector<double> freqs = rp.dipole_frequencies();
    const auto props = rp.requested_properties();
    auto wants = [&](const char *p) {
      return std::find(props.begin(), props.end(), std::string(p)) != props.end();
    };

    std::vector<ResponsePlan> plans;
    auto add = [&](ResponsePropertyRequest r) {
      // The 2PA tensor needs all three half-frequency dipole legs, even when
      // the separately requested polarizability selects only one direction.
      r.axes = r.tpa ? std::vector<char>{'x', 'y', 'z'} : axes;
      r.protocol_thresholds = protocol;
      plans.push_back(plan_one(r));
    };
    if (wants("polarizability")) {
      ResponsePropertyRequest r;
      r.kind = ResponsePropertyKind::Polarizability;
      r.frequencies = freqs;
      add(r);
    }
    // `quadratic true` is the legacy (v2-era) spelling of "compute beta" —
    // honor it alongside requested_properties so old decks keep their
    // hyperpolarizability instead of silently dropping it (M3 golden regen
    // caught this: the alpha+beta deck produced alpha-only output).
    if (wants("hyperpolarizability") || rp.quadratic()) {
      // Optical rectification needs signed-frequency legs in the quadratic
      // source, which the planner does not build yet (beta_vbc_pairs returns
      // no pair for OR): the legs would solve and beta would silently be
      // absent (review finding C10).
      if (rp.beta_or())
        throw std::runtime_error(
            "response: beta.or (optical rectification) is not implemented — its "
            "quadratic source needs signed-frequency legs. Drop beta.or for SHG "
            "(or static, at frequency 0) beta.");
      ResponsePropertyRequest r;
      r.kind = ResponsePropertyKind::Hyperpolarizability;
      r.beta_process = rp.beta_or() ? BetaProcess::OR : BetaProcess::SHG;
      // Deck `beta.frequencies` restricts the driver frequencies (and thus the
      // 2w legs) to the ones actually wanted; empty = the whole dipole grid.
      const auto bf = rp.beta_frequencies();
      r.frequencies = bf.empty() ? freqs : bf;
      add(r);
    }
    if (wants("raman")) {
      // Two spellings of the nuclear coordinates (D-020). raman.atoms and
      // raman.directions (lists; the named path's) span atoms x axes, every
      // atom or x, y, z when one list is empty: concrete nuclear legs and the
      // dipole x nuclear sources for each coordinate, derivative rows only,
      // no normal-mode analysis. With both empty, raman.nuc_atom and
      // raman.nuc_axis (the unnamed path's, until it is ported) select ONE
      // coordinate, or raman.nuc_atom -1 the FULL set: 3N nuclear legs, the
      // sources for every coordinate, and the normal-mode Raman analysis
      // (Hessian from the same nuclear response states, frequencies,
      // activities, depolarization ratios) under
      // properties.vibrational_analysis / properties.raman_spectra. Writing
      // both spellings is an error.
      ResponsePropertyRequest r;
      r.kind = ResponsePropertyKind::PolarizabilityGradient;
      r.gradient_mode = GradientMode::Nuclear;
      // Deck raman.frequencies restricts the optical frequencies, like
      // beta.frequencies; empty = the whole dipole grid.
      const auto rf = rp.raman_frequencies();
      r.frequencies = rf.empty() ? freqs : rf;
      const int natom = static_cast<int>(reference_molecule.natom());
      const auto atoms_written = rp.raman_atoms();
      std::vector<int> axes_written;
      for (char c : rp.raman_directions()) {
        const auto axis = std::string("xyz").find(
            static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
        if (axis == std::string::npos
            || std::find(axes_written.begin(), axes_written.end(),
                         static_cast<int>(axis))
                   != axes_written.end())
          throw std::runtime_error("response: raman.directions must be "
                                   "distinct axes out of xyz");
        axes_written.push_back(static_cast<int>(axis));
      }
      const bool list_form = !atoms_written.empty() || !axes_written.empty();
      if (list_form
          && (rp.is_user_defined("raman.nuc_atom")
              || rp.is_user_defined("raman.nuc_axis")))
        throw std::runtime_error(
            "response: give the nuclear coordinates once, as raman.atoms and "
            "raman.directions (lists) or as raman.nuc_atom and raman.nuc_axis "
            "(one coordinate), not both");
      if (list_form) {
        std::vector<int> atoms;
        for (int atom : atoms_written) {
          if (atom < 0 || atom >= natom
              || std::find(atoms.begin(), atoms.end(), atom) != atoms.end())
            throw std::runtime_error(
                "response: raman.atoms must be distinct atom indices from 0 "
                "to "
                + std::to_string(natom - 1) + " (the molecule has "
                + std::to_string(natom) + " atom(s))");
          atoms.push_back(atom);
        }
        if (atoms.empty())
          for (int atom = 0; atom < natom; ++atom)
            atoms.push_back(atom);
        if (axes_written.empty())
          axes_written = {0, 1, 2};
        for (int atom : atoms)
          for (int axis : axes_written)
            r.raman_coordinates.push_back({atom, axis});
        if (world.rank() == 0)
          print("response: polarizability gradient over",
                r.raman_coordinates.size(),
                "written nuclear coordinate(s); derivative rows only, no "
                "normal-mode analysis");
      } else if (rp.raman_nuc_atom() >= 0) {
        r.raman_nuc_atom = rp.raman_nuc_atom();
        r.raman_nuc_axis = rp.raman_nuc_axis();
        // Upper bound (review finding C13): an out-of-range atom otherwise
        // fails later, inside the nuclear-derivative functor, as an unnamed
        // exception.
        if (r.raman_nuc_atom >= natom)
          throw std::runtime_error(
              "response: raman.nuc_atom " + std::to_string(r.raman_nuc_atom)
              + " is out of range for a molecule of " + std::to_string(natom)
              + " atom(s) (atoms are numbered from 0)");
        if (r.raman_nuc_axis < 0 || r.raman_nuc_axis > 2)
          throw std::runtime_error(
              "response: raman.nuc_axis must be 0 (x), 1 (y) or 2 (z)");
      } else {
        r.raman_nuc_atom = -1;   // full set; the axis is ignored
        r.raman_nuc_axis = -1;
        if (world.rank() == 0)
          print("Raman: all atoms, directions x y z ("
                + std::to_string(3 * natom)
                + " nuclear coordinates, with normal-mode analysis)");
      }
      add(r);
    }
    if (rp.excited_enable()) {
      ResponsePropertyRequest r;
      r.kind = ResponsePropertyKind::PolarizabilityGradient;
      r.gradient_mode = GradientMode::Resonant;
      r.n_roots = static_cast<int>(rp.excited_num_states());
      // Roots only unless the deck also asks for the two-photon contraction
      // (excited.tpa) — the derived FD legs serve only the 2PA residue.
      r.tpa = rp.excited_tpa();
      add(r);
    }
    if (plans.empty()) {  // default: polarizability
      ResponsePropertyRequest r;
      r.kind = ResponsePropertyKind::Polarizability;
      r.frequencies = freqs;
      add(r);
    }

    // 3. Build the workflow input + settings; run the core with the ground
    //    state already loaded. archive_file is still passed through: the GS
    //    fingerprint gate hashes it (restart safety), and the F2 subworld
    //    fan-out needs it to load per-subworld ground states.
    ResponseWorkflowInput in;
    in.archive_file = archive;
    in.protocols = protocol;
    in.plan = merge_plans(plans);
    if (rp.nuclear()) {
      ResponsePlan nuclear_plan;
      for (char direction:rp.nuclear_directions()) {
        const char axis=static_cast<char>(std::tolower(static_cast<unsigned char>(direction)));
        MADNESS_CHECK_THROW(axis>='x' && axis<='z',"invalid nuclear response direction");
        for (double frequency:rp.nuclear_frequencies()) {
          MADNESS_CHECK_THROW(std::isfinite(frequency),"nuclear response frequency must be finite");
          // The all-atoms sentinel is the planner's full-Raman request. These
          // are independent nuclear states, without an implied mode analysis.
          for(int atom=0;atom<gs.molecule().natom();++atom)
            nuclear_plan.nuclear_fd.push_back({Perturbation::nuclear(atom,axis-'x'),frequency,protocol});
        }
      }
      in.plan=merge_plans({in.plan,nuclear_plan});
    }
    // excited.tda=false → Full (X,Y) ES bundle (default TDA).
    if (rp.excited_enable() && !rp.excited_tda())
      for (auto &e : in.plan.es) e.tda = false;
    // ResponseApplication::run has already chdir'd (ScopedCWD) into `outdir`, and
    // `outdir` is RELATIVE — so the calc dir is the cwd ("."). Using outdir here
    // would double the path (outdir/outdir) and the metadata would be written/read
    // in different places, leaving Output.properties empty. (outdir is still used
    // above to resolve the ground archive relative to this cwd.)
    in.settings.calc_dir = reference_prefix.empty() ? "." : outdir.string();
    in.settings.state_plans = state_plans;
    // excited.maxiter / guess_max_iter / maxsub (review C2), when the deck sets them.
    apply_deck_es_knobs(rp, in.settings);
    // Deck `subworlds N` -> the F2 state-parallel fan-out (same path as the
    // standalone --fd-subworlds flag; archive_file above makes it live).
    in.settings.fd_subworlds = std::max(0, rp.subworlds());
    // Deck `subworld_ranks R` -> universe-level split into subworlds of R ranks
    // that may span nodes (large-system regime); 0 keeps per-node packing.
    in.settings.fd_subworld_ranks = std::max(0, rp.subworld_ranks());
    // Deck `dalton.dir` + `seed.freq_tol` -> nearest-frequency DALTON guess for
    // the derived (two-photon) FD legs (calc_executor solve_fd seam).
    in.settings.dalton_dir    = dalton_dir;
    in.settings.seed_freq_tol = rp.seed_freq_tol();
    in.settings.diagnostic_properties = rp.diagnostic_properties();
    if (rp.diagnostic_properties())
      for (char axis : rp.dipole_directions())
        if (axis >= 'x' && axis <= 'z') in.settings.diagnostic_dipole_axes.push_back(axis - 'x');
    if (world.rank() == 0 && in.settings.fd_subworlds > 0) {
      print("response: deck subworlds =", in.settings.fd_subworlds,
            "(F2 state-parallel fan-out requested)");
      // Review io HIGH (early warning, not a hard gate): on a MULTI-NODE run
      // the subworlds write native archives at subworld-rank-count, but
      // property assembly reloads them at universe scale — the native np-guard
      // then aborts AFTER the full solve. Warn up front so the user isn't
      // surprised late; the HDF5 backend gathers to one client and reloads at
      // any np. (Single-node subworlds, where writer==reader np, are fine.)
#ifdef MADNESS_HAS_HDF5
      const bool hdf5_on = hdf5_io_enabled();
#else
      const bool hdf5_on = false;
#endif
      if (!hdf5_on)
        print("response: NOTE — subworlds with the NATIVE backend: per-subworld "
              "archives are treated as np-portable by default at property "
              "assembly (nio=1); set MADRESPONSE_STRICT_NP=1 to hard-block a "
              "cross-np reload instead. `response { io { backend hdf5 } }` "
              "sidesteps the question entirely on multi-node runs.");
    }
    in.settings.policy.dconv_user = rp.dconv();
    // Deck overrides for the iteration policy (2026-09-10), applied only when
    // the deck SETS them: the deck default kain=false predates the solvers'
    // KAIN-on default, so honouring it unconditionally would switch KAIN off
    // for every run. This policy is the FD solves' as is. The ES executors
    // build their own from it (es_iteration_policies): kain.min_residual and
    // the stall detector carry over, while kain / maxrotn reach ES through
    // apply_deck_es_knobs above.
    if (rp.is_user_defined("kain"))    in.settings.policy.kain    = rp.kain();
    if (rp.is_user_defined("maxrotn")) in.settings.policy.maxrotn = rp.maxrotn();
    if (rp.is_user_defined("kain.min_residual")) {
      in.settings.policy.kain_min_residual = rp.kain_min_residual();
      in.settings.es_kain_min_residual_written = true;  // D-122: it wins
    }
    in.settings.print_level = engine_print_level(rp.print_level());
    in.settings.verbosity = madness::iteration::verbosity(rp.print_level());
    // ES initial-guess knobs (deck: response { excited.guess virtual_ao,
    // excited.guess_basis aug-cc-pvtz }). Deck default (solid_harmonics /
    // aug-cc-pvdz) maps to the ExecutorSettings default — existing decks are
    // unchanged. virtual_ao is the energy-ordered AO-virtual guess; it is
    // required to reach totally-symmetric / radially-excited states on atoms,
    // which the angular-only solid-harmonic trials structurally cannot span.
    in.settings.es_guess       = parse_es_guess_mode(rp.excited_guess());
    in.settings.es_guess_basis = rp.excited_guess_basis();
    if (world.rank() == 0 && rp.excited_enable() &&
        in.settings.es_guess != ESGuessMode::SolidHarmonics)
      print("response: excited.guess =", to_string(in.settings.es_guess),
            (in.settings.es_guess == ESGuessMode::VirtualAO
                 ? std::string("(basis " + in.settings.es_guess_basis + ")")
                 : std::string()));
    // Deck-level HDF5 opt-in (response.hdf5 true) — the env var
    // MADRESPONSE_IO_HDF5 still works; the deck parameter wins when set.
    if (effective_hdf5(params)) {
#ifdef MADNESS_HAS_HDF5
      set_hdf5_io_enabled(true);
#else
      throw std::runtime_error(
          "response.hdf5 requested but this build has no HDF5 support — "
          "configure with -DMADNESS_ENABLE_HDF5=ON");
#endif
    }

#ifdef MADNESS_HAS_HDF5
    // Ground state in HDF5 (2026-09-09): with the HDF5 restart opt-in, mirror
    // moldft's native archive as <prefix>.restartdata.h5 (same stream, one
    // blob) so a calc dir can carry GS + response + ES in HDF5. Native stays
    // authoritative for moldft; GroundState::from_archive reads either.
    // Always (re)write: the DALTON GS seed hook writes <prefix>.restartdata.h5
    // BEFORE moldft runs (its HDF5 twin of the seed), so an existence check
    // here kept the unconverged seed as the .h5 GS (closeout attempt 9, 2026-09-10:
    // seed_h2o.restartdata.h5 == the 12:01:12 seed twin, native at 12:01:41).
    // Native was still preferred by from_archive, so nothing was mis-loaded,
    // but a native-less consumer would have read the seed. The preserved twin
    // lives under <prefix>.gs_seed.restartdata.h5.
    if (hdf5_io_enabled()) {
      // An explicit reference may be an immutable checkpoint generation.
      // The response-prepared mirror belongs to this response attempt.
      const auto mirror = reference_prefix.empty()
          ? archive + ".h5"
          : (outdir / "reference.restartdata.h5").string();
      gs.save_archive_hdf5(world, mirror);
    }
#endif

    if (shared_request) {
      MADNESS_CHECK_THROW(shared_source &&
          (!in.plan.fd.empty() || !in.plan.es.empty()),
          "shared response supports independent driven states");
      std::vector<madness::iteration::Protocol> ladder;
      for (double threshold:protocol) ladder.push_back({threshold,default_k_for_thresh(threshold)});
      std::vector<response_shared::TaskSpec> definitions;
      const auto add_definition=[&](const Perturbation& perturbation,double frequency) {
        MADNESS_CHECK_THROW(perturbation.kind==Perturbation::Kind::Dipole ||
            (perturbation.kind==Perturbation::Kind::NuclearDisplacement && gs.is_spin_restricted()),
            "response source is not implemented for this perturbation and spin layout");
        definitions.push_back({fd_node_id(perturbation,frequency),nlohmann::json{
          {"perturbation",perturbation.description()},{"frequency",frequency},
          {"layout",std::string(gs.is_spin_restricted() ? "closed_shell_" : "open_shell_") +
              (detail_exec::is_static_freq(frequency) ? "static" : "full")},
          {"exchange_tensor",in.settings.policy.exchange_tensor}}});
      };
      for (const auto& leg:in.plan.fd) add_definition(leg.pert,leg.freq);
      for (const auto& leg:in.plan.nuclear_fd) {
        if (leg.pert.is_all_atoms())
          for (int atom=0;atom<gs.molecule().natom();++atom)
            add_definition(Perturbation::nuclear(atom,leg.pert.axis),leg.freq);
        else add_definition(leg.pert,leg.freq);
      }
      for(const auto& bundle:in.plan.es) {
        const std::string model=bundle.tda?"tda":"full";
        definitions.push_back({es_node_id(bundle.tda,bundle.n_roots),
          {{"excited_model",model},{"layout","closed_shell_"+model}},
          response_shared::excited_task(bundle.n_roots),"es:"+model});
      }
      in.settings.shared_response=std::make_shared<response_shared::Session>(world,shared_root,*shared_source,
          *shared_request,ladder,definitions,std::move(response_timing),
          in.settings.verbosity);
      // Direct sources run once at the final protocol, after B/C qualify.
      // The binding then checks exact inputs instead of legacy convergence flags.
      for(const auto& request:expand_vbc_requests(in.plan.vbc,static_cast<int>(gs.molecule().natom()))) {
        StatePlan entry;entry.kind=StatePlan::Kind::fresh;
        entry.start_thresh=protocol.back();
        in.settings.state_plans[vbc_node_id(request.pert_b,request.pert_c,request.freq_b,request.freq_c)]=entry;
      }
      const auto plan_state=[&](const std::string& id,const auto& binding) {
        StatePlan entry;
        entry.kind=binding.reusable() ? StatePlan::Kind::load :
            binding.has_saved() ? StatePlan::Kind::resume : StatePlan::Kind::fresh;
        const auto start=ladder.at(binding.plan().start_protocol);
        entry.start_thresh=start.thresh;
        if (binding.reusable()) {entry.met=true;entry.loaded_key=protocol_key(start.thresh,start.k);}
        in.settings.state_plans[id]=entry;
      };
      for (const auto& [id,binding]:in.settings.shared_response->states) plan_state(id,*binding);
      for (const auto& [id,binding]:in.settings.shared_response->excited) plan_state(id,*binding);
    }

    // Deck `dalton.dir <path>` — the seed-from-directory import contract
    // (showcase W3; import-only, madness never invokes DALTON). Runs BEFORE
    // the workflow and writes a seed, <run>/seeds/response/<node name>/, for
    // each state with no save of its own; the state starts from it.
    // The GS is prepared at protocol.front() above (the documented
    // precondition); geometry-fingerprint / frequency mismatches throw here,
    // failing the response task loudly instead of silently solving cold.
    const bool need_response_import=
        std::any_of(in.settings.shared_response->states.begin(),in.settings.shared_response->states.end(),
          [](const auto& entry){return !entry.second->has_saved();}) ||
        std::any_of(in.settings.shared_response->excited.begin(),in.settings.shared_response->excited.end(),
          [](const auto& entry){return !entry.second->has_saved();});
    if (!dalton_dir.empty() && need_response_import) {
      in.settings.shared_response->timing().begin("import_guesses");
      // A seed belongs to this calculation when its settings (geometry and
      // model) match; the nodes it is for are those with no save of their own.
      DaltonSeedTargets targets{shared_root,in.settings.calc_dir,{},{}};
      for (const auto& [id,state]:in.settings.shared_response->states)
        if (!state->has_saved())
          targets.driven[state->name()]=response_seed::seed_settings(state->scientific_settings());
      for (const auto& [id,set]:in.settings.shared_response->excited)
        if (!set->has_saved())
          targets.excited[std::filesystem::path(set->folder()).filename().string()]=
              response_seed::seed_settings(set->scientific_settings());
      run_dalton_import(world, gs, reference_molecule, in.plan,
                        targets, dalton_dir, {}, {}, {}, 1e-4,
                        /*es_y_from_dalton=*/rp.seed_es_y() == "dalton");
      // D-118: the import's own time is the import_guesses node of
      // timings.json; what it seeded from is in each seed node's details.
      in.settings.shared_response->timing().end();
    }

    // Pass the RESOLVED fock path through (the ground-state node's
    // scf.fock.json): every per-rung re-prepare and subworld GS reload reads
    // it at the node's level and rebuilds the matrix at any other.
    ResponseWorkflowOutput out =
        run_response_with_ground(world, gs, L, fock_json, in);

    // 4. Map Output → Results: the stage times ride into calc_info.json for
    //    the wall-time lines; the values and the stop reason are in results.json.
    Results res;
    res.metadata = {{"v3_timing", std::move(out.timing)}};
    return res;
  }
};

#endif // MOLRESPONSE_V3_MADQC_ADAPTER_HPP
