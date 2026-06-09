#ifndef MOLRESPONSE_V3_MADQC_ADAPTER_HPP
#define MOLRESPONSE_V3_MADQC_ADAPTER_HPP

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
// It builds a v3 GroundState from madqc's IN-MEMORY SCF reference (the
// GroundState(world, shared_ptr<SCF>) ctor exists for exactly this), maps the
// ResponseParameters input deck to a v3 ResponsePlan, and calls
// run_response_with_ground (the archive-free core).
//
// R3a SCOPE: polarizability (alpha) only. beta / raman / excited mapping is R3b.
// This header lives in the v3 app and is included by madqc.cpp (the app links
// both MADchem and v3's GroundState.cpp); MADchem/WorkflowBuilders must NOT
// include it (that would be a circular library dependency — the engine is
// selected in madqc.cpp instead).
// -----------------------------------------------------------------------------

#include <apps/molresponse_v3/orchestrator/response_workflow.hpp>

#include <madness/chem/CalculationParameters.h>
#include <madness/chem/ParameterManager.hpp>   // Params
#include <madness/chem/ResponseParameters.hpp>
#include <madness/chem/SCF.h>
#include <madness/external/nlohmann_json/json.hpp>

#include <algorithm>
#include <filesystem>
#include <memory>
#include <vector>

/// Global namespace (mirrors `molresponse_lib`) so madqc can write
/// `ResponseApplication<molresponse_v3_lib>`.
struct molresponse_v3_lib {
  /// Output subdir name + the interface ResponseApplication reads.
  static const char *label() { return "molresponse_v3"; }

  /// Structured result returned to ResponseApplication (→ calc_info.json).
  struct Results {
    nlohmann::json metadata;
    nlohmann::json properties;
    nlohmann::json vibrational_analysis;  // empty for alpha (R3a)
    nlohmann::json raman_spectra;         // empty for alpha (R3a)
  };

  inline static Results
  run_response(madness::World &world, const Params &params,
               const std::shared_ptr<madness::SCF> &scf_calc,
               const std::filesystem::path &outdir) {
    using namespace madness;
    using namespace molresponse_v3;

    const auto &cp = params.get<CalculationParameters>();
    const auto &rp = params.get<ResponseParameters>();
    const std::vector<double> protocol = cp.protocol();
    MADNESS_CHECK(!protocol.empty());

    // 1. v3 GroundState from the live SCF; set the coarsest protocol + prepare.
    GroundState gs(world, scf_calc);
    const double L = cp.L();
    set_response_protocol(world, L, protocol.front());
    const double thresh = FunctionDefaults<3>::get_thresh();
    auto coulop = poperatorT(
        CoulombOperatorPtr(world, gs.params().lo(), 0.001 * thresh));
    gs.prepare(world, 0.001 * thresh, coulop, /*fock_json=*/"");

    // 2. Map the input deck → a Tier-A plan. R3a: polarizability only.
    ResponsePropertyRequest req;
    req.kind = ResponsePropertyKind::Polarizability;
    req.frequencies = rp.dipole_frequencies();
    req.protocol_thresholds = protocol;
    req.axes.clear();
    for (char c : rp.dipole_directions()) {
      const char l = static_cast<char>(std::tolower(c));
      if (l == 'x' || l == 'y' || l == 'z') req.axes.push_back(l);
    }
    if (req.axes.empty()) req.axes = {'x', 'y', 'z'};

    // 3. Build the workflow input + settings; run the archive-free core.
    ResponseWorkflowInput in;
    in.protocols = protocol;
    in.plan = plan_one(req);
    // ResponseApplication::run has already chdir'd (ScopedCWD) into `outdir`, and
    // `outdir` is RELATIVE — so the calc dir is the cwd ("."). Using outdir here
    // would double the path (outdir/outdir) and the metadata would be written/read
    // in different places, leaving Output.properties empty.
    (void)outdir;
    in.settings.calc_dir = ".";
    in.settings.max_iters = static_cast<int>(rp.maxiter());
    in.settings.policy.dconv_user = rp.dconv();
    in.settings.print_level =
        static_cast<PrintLevel>(std::max(0, std::min(3, rp.print_level())));

    ResponseWorkflowOutput out =
        run_response_with_ground(world, gs, L, /*fock_json=*/"", in);

    // 4. Map Output → Results. Stash v3 timing/diagnostics under metadata so
    //    they ride into the workflow's calc_info.json.
    Results res;
    res.metadata = std::move(out.metadata);
    res.properties = std::move(out.properties);
    if (world.rank() == 0) {
      res.metadata["engine"]         = "molresponse_v3";
      res.metadata["v3_timing"]      = std::move(out.timing);
      res.metadata["v3_diagnostics"] = std::move(out.diagnostics);
    }
    return res;
  }
};

#endif // MOLRESPONSE_V3_MADQC_ADAPTER_HPP
