// ===========================================================================
// test_tpa_pq_vs_vbc.cpp — gate for the PROMOTED state-free two-photon
// residue source (P,Q) (kernels/tpa_source_spec.hpp, source-unification
// step 3), and the measured comparison against beta's quadratic source V^{bc}.
//
// The (P,Q) assembly was derived and validated IN this test (see git history:
// the regrouping identities are pinned by test_tpa_regroup_identities, and the
// original inline build agreed with the c-grouped composition to 3e-7 on
// stored FD responses; that composition and its gate were removed in step 4c
// part 5b). It now lives in kernels/tpa_source_spec.hpp as data over the
// declarative engine (kernels/source_spec.hpp). The gates:
//
//   (2) ONE-ELECTRON gate: the full assemble_tpa_pq(B,C,vB,vC) adds exactly
//       tpa::tpa_moment_residue_1e's half(b,c)+half(c,b) content to the
//       two-electron assemble_tpa_pq(B,C): contraction(full) ==
//       contraction(2e) + S_1e.
//   (4) the symmetrized builder equals the sum of both orderings.
//   (5) the two-electron (P,Q) equals V^{bc}'s two-electron part after
//       Q-projection.
//
// Stand-in vectors are used by default (CI needs no external data beyond the
// ground-state fixture). They must NOT be Q-projected and must admix phi with
// irregular coefficients: Q-projected vectors are orthogonal to phi, which
// silently zeroes every <phi|...> factor, and pure Cartesian dipoles on a
// symmetric molecule zero further terms by selection rule. Either makes the
// test pass without testing anything.
//
//   test_tpa_pq_vs_vbc --archive=<moldft restartdata> [--thresh=X] [--k=N]
//       [--calc-dir=<madqc run folder> --freq=<omega_f/2> [--baxis=N --caxis=N]]
// With --calc-dir, b and c are the saved driven states dipole_<baxis> and
// dipole_<caxis> (default z) at --freq, found among the run folder's node
// folders (states/response/) and loaded at this test's level. A static state
// (X only) is used as (X, Y = X). A missing state stops the test.
// ===========================================================================

#include "../GroundState.hpp"
#include "../Perturbations.hpp"
#include "../ResponseProtocol.hpp"
#include "../kernels/common_ops.hpp"
#include "../kernels/tags.hpp"
#include "../kernels/tpa.hpp"
#include "../kernels/tpa_source_spec.hpp"
#include "../kernels/two_electron.hpp"
#include "../kernels/vbc.hpp"
#include "../solvers/build_response_ground_state.hpp"
#include "../solvers/response_state.hpp"
#include "../tools/saved_states.hpp"

#include <madness/misc/info.h>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

using namespace madness;
using namespace molresponse_v3;

using vecfuncT = std::vector<real_function_3d>;

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  int rc = 0;
  try {
    startup(world, argc, argv, true);
    {
      commandlineparser parser(argc, argv);
      if (!parser.key_exists("archive")) {
        if (world.rank() == 0)
          print("Usage: test_tpa_pq_vs_vbc --archive=<path> [--thresh=X] [--k=N] "
                "[--calc-dir=<run folder> --freq=<omega_f/2> [--baxis=N --caxis=N]]");
        finalize();
        return 1;
      }
      const std::string archive_path = parser.value_raw("archive");
      auto header = GroundState::read_archive_header(world, archive_path);
      const int override_k = parser.key_exists("k") ? std::stoi(parser.value("k")) : -1;
      const double thresh = parser.key_exists("thresh")
                                ? std::stod(parser.value("thresh"))
                                : default_thresh_for_k(header.k);
      set_response_protocol(world, header.L, thresh, override_k);

      Molecule molecule;
      auto dir = std::filesystem::path(archive_path).parent_path();
      for (const auto &name : {"moldft.calc_info.json", "mad.calc_info.json"}) {
        auto cand = dir / name;
        if (std::filesystem::exists(cand)) {
          std::ifstream ifs(cand);
          nlohmann::json j;
          ifs >> j;
          nlohmann::json mj;
          if (j.contains("tasks") && j["tasks"].is_array() && !j["tasks"].empty())
            mj = j["tasks"][0]["molecule"];
          else if (j.contains("molecule")) mj = j["molecule"];
          if (!mj.is_null()) molecule.from_json(mj);
          break;
        }
      }
      auto gs = GroundState::from_archive(world, archive_path, molecule);
      std::string fock_json;
      for (const auto &name : {"moldft.fock.json", "mad.fock.json"}) {
        auto cand = dir / name;
        if (std::filesystem::exists(cand)) { fock_json = cand.string(); break; }
      }
      const double t = FunctionDefaults<3>::get_thresh();
      auto coulop = poperatorT(CoulombOperatorPtr(world, gs.params().lo(), 0.001 * t));
      gs.prepare(world, 0.001 * t, coulop, fock_json);
      auto g0 = build_response_ground_state_closed_shell(
          world, gs, gs.hf_exchange_coefficient(), gs.params().lo());

      const vecfuncT &phi = g0.amo;
      const size_t n = phi.size();
      const double cxc = g0.c_xc;

      // ---- stand-in b, c, f vectors (see header note on why not Q-projected)
      auto d0 = dipole_perturbation(world, gs, 0);
      auto d1 = dipole_perturbation(world, gs, 1);
      auto d2 = dipole_perturbation(world, gs, 2);
      auto comb = [&](double a, const vecfuncT &u, double b, const vecfuncT &v,
                      double cp) {
        vecfuncT r = madness::copy(world, u);
        gaxpy(world, a, r, b, v);
        gaxpy(world, 1.0, r, cp, phi);
        truncate(world, r);
        return r;
      };
      vecfuncT xb = comb(1.0, d0,  0.23, d1,  0.31);
      vecfuncT yb = comb(1.0, d1, -0.41, d2,  0.17);
      vecfuncT xc = comb(1.0, d2,  0.13, d0, -0.29);
      vecfuncT yc = comb(1.0, d0,  0.53, d1,  0.11);
      const vecfuncT xf = comb(1.0, d1, -0.37, d2, -0.19);
      const vecfuncT yf = comb(1.0, d2,  0.71, d0,  0.43);

      // ---- optionally replace b,c with SAVED driven states ----------------
      // --calc-dir=<run folder> --freq=<omega_f/2> [--baxis=N --caxis=N]
      // The stand-ins above exercise the algebra; real responses give physical
      // magnitudes for the (p,q)-vs-V^{bc} comparison.
      std::string src = "stand-in (dipole combinations + phi admixture)";
      if (parser.key_exists("calc-dir")) {
        MADNESS_CHECK_THROW(parser.key_exists("freq"), "--calc-dir needs --freq");
        const std::string cdir = parser.value_raw("calc-dir");
        const double freq = std::stod(parser.value("freq"));
        const int ba = parser.key_exists("baxis") ? std::stoi(parser.value("baxis")) : 2;
        const int ca = parser.key_exists("caxis") ? std::stoi(parser.value("caxis")) : 2;
        auto load = [&](int axis) {
          const std::string pert = std::string("dipole_") + "xyz"[axis];
          const auto node = saved_states::find_driven(world, cdir, pert, freq);
          if (!node)
            throw std::runtime_error("no saved driven state " + pert + " at frequency " +
                                     parser.value("freq") + " in " + cdir + "/states/response");
          saved_states::print_loaded(world, *node);
          return saved_states::load_driven_xy(world, *node);
        };
        const auto rb = load(ba);
        const auto rcc = load(ca);
        xb = madness::copy(world, rb.x_alpha);
        yb = madness::copy(world, rb.y_alpha);
        xc = madness::copy(world, rcc.x_alpha);
        yc = madness::copy(world, rcc.y_alpha);
        src = "SAVED driven states  dipole_" + std::string(1, "xyz"[ba]) +
              " / dipole_" + std::string(1, "xyz"[ca]) + "  @ omega=" + parser.value("freq");
      }

      ResponseStateXY<ClosedShell> B, C, F;
      B.x_alpha = madness::copy(world, xb); B.y_alpha = madness::copy(world, yb);
      C.x_alpha = madness::copy(world, xc); C.y_alpha = madness::copy(world, yc);
      F.x_alpha = madness::copy(world, xf); F.y_alpha = madness::copy(world, yf);

      // ============ (P,Q) from the PROMOTED kernel (2e part) ==============
      auto PQ = tpa::assemble_tpa_pq(world, g0, B, C);

      const double sgrouped = inner(xf, PQ.x_alpha) + inner(yf, PQ.y_alpha);

      // ============ gate 2: full (1e + 2e) vs 2e (P,Q) + residue_1e =======
      const auto mu_b = dipole_operator(world, 0);
      const auto mu_c = dipole_operator(world, 2);
      auto PQ_full = tpa::assemble_tpa_pq(world, g0, B, C, mu_b, mu_c);
      // arrays laid out so S_1e(0,1) = half(b,c) + half(c,b) with (vB, vC)
      const std::array<ResponseStateXY<ClosedShell>, 3> resp{
          B, C, C};
      const std::array<real_function_3d, 3> ops{mu_b, mu_c, mu_c};
      const auto S1e = tpa::tpa_moment_residue_1e(world, g0, F, resp, ops, 1.0);
      const double ref_full = sgrouped + S1e(0L, 1L);
      const double got_full =
          inner(xf, PQ_full.x_alpha) + inner(yf, PQ_full.y_alpha);

      // ============ gate 4: symmetrized builder == two-ordering sum =======
      // tpa_pq_spec_sym exploits the family-D collapse (D(B,C)==B(C,B), the
      // adjoint identity Gbar=G^T); this gate keeps the CURRENT per-ordering
      // builder as the validation reference for the simplified one.
      {
        // Reference and sym both UNTRUNCATED raw engine output: the internal
        // truncations of assemble_tpa_pq differ per path (each ordering
        // truncated separately vs one final truncate) and inject O(thresh)
        // relative noise that has nothing to do with the identity.
        auto ref_bc = source_spec::assemble_source(
            world, g0, tpa::tpa_pq_spec(world, g0, B, C));
        auto ref_cb = source_spec::assemble_source(
            world, g0, tpa::tpa_pq_spec(world, g0, C, B));
        vecfuncT refx = std::move(ref_bc[0]);
        gaxpy(world, 1.0, refx, 1.0, ref_cb[0]);
        vecfuncT refy = std::move(ref_bc[1]);
        gaxpy(world, 1.0, refy, 1.0, ref_cb[1]);
        auto symr = source_spec::assemble_source(
            world, g0, tpa::tpa_pq_spec_sym(world, g0, B, C));
        auto nrm = [&](const vecfuncT &v) {
          return std::sqrt(std::abs(inner(world, v, v).sum()));
        };
        vecfuncT dx = madness::copy(world, symr[0]);
        gaxpy(world, 1.0, dx, -1.0, refx);
        vecfuncT dy = madness::copy(world, symr[1]);
        gaxpy(world, 1.0, dy, -1.0, refy);
        const double ex = nrm(dx) / std::max(1e-30, nrm(refx));
        const double ey = nrm(dy) / std::max(1e-30, nrm(refy));
        // Tolerance: the sym path applies the DAGGERED kernel where the
        // reference transposes the undaggered occupied block — the adjoint
        // identity holds to MRA operator-apply precision (measured ~1e-2 x
        // thresh; same scaling as the M1 identity in
        // test_tpa_regroup_identities), not to roundoff. A real bookkeeping
        // error shows at O(0.1-1) relative; gate at 0.1 x thresh.
        const double tol4 = 0.1 * t;
        const bool ok4 = ex < tol4 && ey < tol4;
        if (world.rank() == 0)
          print("gate 4 (sym == sum of orderings):  rel dx=", ex,
                " rel dy=", ey, ok4 ? " PASS" : " FAIL");
        if (!ok4) rc = 1;
      }

      // ============ (3) gate 5: (P,Q) two-electron == V^{bc} two-electron ==
      // Zero one-electron operators isolate the two-electron content. Since
      // 2026-09-09 (kernels/vbc.hpp in the equation's leg orientation) the two
      // builders encode the same source; they differ only by Q-projection
      // (P,Q leave the apply moves unprojected) and by the truncation points
      // of the pair-density family, so compare Q(P,Q) with V at thresh scale.
      // Before the fix this block printed rel 0.41 / 0.23 (HANDOFF §2a).
      // Both builders assume response vectors in the virtual space (V^{bc}'s
      // occupied-matrix term is Sum_k x^c_k F^b_kp, in the span of x^c), so the
      // comparison is made on Q-projected copies of b, c (stand-ins carry a
      // phi admixture; stored FD responses are already Q-space) and both
      // outputs are Q-projected. norm2 is collective: every norm is computed
      // here, on all ranks, before the rank-0 print block (job 2161695 hung
      // with the norms inside it).
      ResponseStateXY<ClosedShell> Bq, Cq;
      Bq.x_alpha = g0.Qa(B.x_alpha); Bq.y_alpha = g0.Qa(B.y_alpha);
      Cq.x_alpha = g0.Qa(C.x_alpha); Cq.y_alpha = g0.Qa(C.y_alpha);
      auto PQq = tpa::assemble_tpa_pq(world, g0, Bq, Cq);   // 2e only
      real_function_3d zop = madness::copy(phi[0]); zop.scale(0.0);
      auto V = vbc::compute_vbc<ClosedShell>(world, g0, Bq, Cq, zop, madness::copy(zop));
      vecfuncT PQx = g0.Qa(PQq.x_alpha);
      vecfuncT PQy = g0.Qa(PQq.y_alpha);
      vecfuncT Vx  = g0.Qa(V.x_alpha);
      vecfuncT Vy  = g0.Qa(V.y_alpha);

      auto nrm = [&](const vecfuncT &v) { return norm2(world, v); };
      vecfuncT dpx = madness::copy(world, PQx); gaxpy(world, 1.0, dpx, -1.0, Vx);
      vecfuncT dqy = madness::copy(world, PQy); gaxpy(world, 1.0, dqy, -1.0, Vy);
      const double nPQx = nrm(PQx), nPQy = nrm(PQy);
      const double nPQx_raw = nrm(PQq.x_alpha), nPQy_raw = nrm(PQq.y_alpha);
      const double nVx = nrm(Vx), nVy = nrm(Vy);
      const double ndpx = nrm(dpx), ndqy = nrm(dqy);
      const double relx5 = nVx > 1e-12 ? ndpx / nVx : ndpx;
      const double rely5 = nVy > 1e-12 ? ndqy / nVy : ndqy;
      // Tolerance: the two builders regroup the pair-density family at
      // different truncation points (four untruncated densities vs one
      // bundled, truncated one); on the two-electron part alone this is
      // 1.8e-5 / 1.6e-5 relative at thresh 1e-6 (h2o fixture, job 2161860),
      // i.e. ~18 x thresh. An orientation error shows at 0.1-1 relative
      // (0.41 / 0.23 before the fix). Gate at 50 x thresh, in line with
      // gates 1-2 (200 x thresh) on the same regrouping.
      const double tol5 = 50.0 * t;
      const bool ok5 = relx5 < tol5 && rely5 < tol5;
      if (!ok5) rc = 1;

      if (world.rank() == 0) {
        print("\n=== promoted (P,Q) source (tpa_source_spec)  vs  references ===");
        print("  thresh =", t, "  n_occ =", n, "  c_xc =", cxc);
        print("  b,c source:", src);
        print("\n--- gate 2: full (1e+2e) vs 2e (P,Q) + tpa_moment_residue_1e ---");
        printf("  two-electron <x^f|P>+<y^f|Q>    = %+.10f\n", sgrouped);
        printf("  reference  2e + S_1e            = %+.10f\n", ref_full);
        printf("  contraction of full (P,Q)       = %+.10f\n", got_full);
        const double relf = std::abs(ref_full) > 1e-12
                              ? std::abs(got_full - ref_full) / std::abs(ref_full)
                              : std::abs(got_full - ref_full);
        printf("  relative difference             = %.3e   %s\n", relf,
               relf < 200.0 * t ? "PASS" : "FAIL");
        print("\n--- gate 5: Q(P,Q) two-electron == Q(V^{bc}) two-electron (Q-space b,c) ---");
        printf("  %-28s ||P||=%10.6f  ||Q||=%10.6f   (unprojected: %10.6f %10.6f)\n",
               "Q(P,Q) total", nPQx, nPQy, nPQx_raw, nPQy_raw);
        printf("  %-28s ||x||=%10.6f  ||y||=%10.6f\n", "Q(V^{bc}) 2-electron", nVx, nVy);
        printf("  %-28s ||  ||=%10.6f  ||  ||=%10.6f\n", "Q(P,Q) - Q(V^{bc})", ndpx, ndqy);
        printf("  relative:  ||Q(P)-Vx||/||Vx|| = %.3e   ||Q(Q)-Vy||/||Vy|| = %.3e   tol=%.1e   %s\n",
               relx5, rely5, tol5, ok5 ? "PASS" : "FAIL");
        if (!(relf < 200.0 * t)) rc = 1;   // sticky: keep gate-4/5 verdicts
        print("\n", rc == 0 ? "PQ_SOURCE_SPEC VALIDATED" : "PQ_SOURCE_SPEC FAILED");
      }
      world.gop.broadcast(rc, 0);
    }
  } catch (const std::exception &e) {
    if (world.rank() == 0) print("EXCEPTION:", e.what());
    rc = 2;
  }
  finalize();
  return rc;
}
