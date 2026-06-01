// ===========================================================================
// test_calc_manager_run.cpp — end-to-end drive of the calc manager (doc 15,
// 15a). Sets up a ground state from a moldft archive (same recipe as the FD
// skeleton test), plans a polarizability request, then lets CalcManager +
// FdResponseExecutor schedule and solve every FD protocol step. Validates by reading
// back response_metadata.json and checking the expected FD states converged.
//
// This is an ALLOCATION test (it runs real MRA solves). Usage:
//   test_calc_manager_run --archive=<moldft restartdata> \
//       [--omega=0.0,0.057] [--axes=xyz] [--protocol=1e-4,1e-6] [--es-roots=N] \
//       [--maxiter=N] [--dconv=X] [--calc-dir=DIR] [--print-level=0..3]
// ===========================================================================

#include "../GroundState.hpp"
#include "../ResponseProtocol.hpp"
#include "../ResponsePropertyPlanner.hpp"
#include "../calc/calc_executor.hpp"
#include "../solvers/convergence_policy.hpp"
#include "../solvers/response_metadata.hpp"

#include <madness/external/nlohmann_json/json.hpp>
#include <madness/misc/info.h>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

using namespace madness;
using namespace molresponse_v3;

namespace {

std::vector<double> parse_csv_doubles(const std::string &s) {
  return parse_protocol_csv(s);  // reuses the protocol CSV parser
}

std::vector<char> parse_axes(const std::string &s) {
  std::vector<char> ax;
  for (char c : s) {
    char l = std::tolower(c);
    if (l == 'x' || l == 'y' || l == 'z') ax.push_back(l);
  }
  if (ax.empty()) ax = {'x', 'y', 'z'};
  return ax;
}

} // namespace

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  int rc = 0;
  try {
    startup(world, argc, argv, true);
    if (world.rank() == 0) {
      print("\n=========================================================");
      print(" molresponse_v3 CalcManager run test  (15a, FD/alpha)");
      print("=========================================================\n");
    }
    {
      commandlineparser parser(argc, argv);
      if (!parser.key_exists("archive")) {
        if (world.rank() == 0)
          print("Usage: test_calc_manager_run --archive=<path> "
                "[--omega=0.0,0.057] [--axes=xyz] [--protocol=1e-4,1e-6] "
                "[--maxiter=N] [--dconv=X] [--calc-dir=DIR] "
                "[--print-level=0..3]");
        finalize();
        return 1;
      }

      const std::string archive_path = parser.value_raw("archive");
      const std::vector<double> freqs =
          parser.key_exists("omega")
              ? parse_csv_doubles(parser.value("omega"))
              : std::vector<double>{0.0};
      const std::vector<char> axes =
          parse_axes(parser.key_exists("axes") ? parser.value("axes") : "xyz");
      const int max_iters = parser.key_exists("maxiter")
                                ? std::stoi(parser.value("maxiter")) : 25;
      const double dconv = parser.key_exists("dconv")
                               ? std::stod(parser.value("dconv")) : 1e-4;
      const int pl_int = parser.key_exists("print-level")
                             ? std::stoi(parser.value("print-level")) : 1;
      const PrintLevel print_level =
          static_cast<PrintLevel>(std::max(0, std::min(3, pl_int)));

      auto header = GroundState::read_archive_header(world, archive_path);

      std::vector<double> protocol;
      if (parser.key_exists("protocol"))
        protocol = parse_protocol_csv(parser.value("protocol"));
      else
        protocol.push_back(default_thresh_for_k(header.k));
      set_response_protocol(world, header.L, protocol.front());

      // Molecule (for natom -> nuclear expansion; harmless for pure alpha).
      Molecule molecule;
      auto archive_dir = std::filesystem::path(archive_path).parent_path();
      for (const auto &name : {"moldft.calc_info.json", "mad.calc_info.json"}) {
        auto candidate = archive_dir / name;
        if (std::filesystem::exists(candidate)) {
          std::ifstream ifs(candidate);
          nlohmann::json j; ifs >> j;
          nlohmann::json mol_json;
          if (j.contains("tasks") && j["tasks"].is_array() && !j["tasks"].empty())
            mol_json = j["tasks"][0]["molecule"];
          else if (j.contains("molecule"))
            mol_json = j["molecule"];
          if (!mol_json.is_null()) molecule.from_json(mol_json);
          break;
        }
      }
      auto gs = GroundState::from_archive(world, archive_path, molecule);
      if (world.rank() == 0) gs.print_info();

      std::string fock_json;
      for (const auto &name : {"moldft.fock.json", "mad.fock.json"}) {
        auto candidate = archive_dir / name;
        if (std::filesystem::exists(candidate)) {
          fock_json = candidate.string();
          break;
        }
      }
      const double cur_thresh = FunctionDefaults<3>::get_thresh();
      auto coulop = poperatorT(
          CoulombOperatorPtr(world, gs.params().lo(), 0.001 * cur_thresh));
      gs.prepare(world, 0.001 * cur_thresh, coulop, fock_json);

      ConvergencePolicy policy;
      policy.dconv_user = dconv;

      const std::string calc_dir =
          parser.key_exists("calc-dir") ? parser.value_raw("calc-dir")
                                        : std::string("calc_manager_run");

      // ---- Plan: polarizability, or (--es-roots=N) resonant gradient -----
      const int es_roots =
          parser.key_exists("es-roots") ? std::stoi(parser.value("es-roots")) : 0;
      ResponsePropertyRequest req;
      req.axes = axes;
      req.protocol_thresholds = protocol;
      if (es_roots > 0) {
        // ES(TDA) bundle + symbolic derived dipole FD; the calc manager solves
        // the bundle, then expands to FD at the converged excitation energies.
        req.kind          = ResponsePropertyKind::PolarizabilityGradient;
        req.gradient_mode = GradientMode::Resonant;
        req.n_roots       = es_roots;
      } else {
        req.kind        = ResponsePropertyKind::Polarizability;
        req.frequencies = freqs;
      }
      ResponsePlan plan = plan_one(req);

      if (world.rank() == 0) {
        print("\nRUN CONFIG:");
        print("  archive    =", archive_path);
        print("  shell      =",
              gs.is_spin_restricted() ? "ClosedShell" : "OpenShell");
        print("  omega      =", freqs);
        print("  protocol   =", protocol);
        print("  calc_dir   =", calc_dir);
        print("  mode       =", (es_roots > 0 ? "resonant-gradient (ES)" : "polarizability"));
        if (es_roots > 0) print("  es_roots   =", es_roots);
        print("  FD requests=", (int)plan.fd.size(),
              "  ES requests=", (int)plan.es.size());
      }

      // ---- Drive the calc manager ----------------------------------------
      CalcManager::Policy mgr_policy;
      mgr_policy.max_iters_per_step = max_iters;
      CalcManager mgr(plan, calc_dir, mgr_policy);
      mgr.build(molecule.natom());

      ExecutorContext ctx{world, gs, header.L, fock_json,
                          policy, print_level, calc_dir, max_iters};
      FdResponseExecutor exec(ctx);
      mgr.run(world, exec);

      // ---- Validate at the top protocol ----------------------------------
      const std::string top_key = protocol_key_at(protocol.back());
      if (world.rank() == 0) {
        auto meta = ResponseMetadata::load_or_create(
            calc_dir + "/response_metadata.json");
        const auto &j = meta.json();
        if (es_roots > 0) {
          // ES bundle converged with n_roots, plus the promoted derived FDs.
          bool es_ok = j["excited_states"].contains(top_key) &&
                       j["excited_states"][top_key].value("converged", false);
          int n_es = (es_ok && j["excited_states"][top_key].contains("roots"))
                         ? (int)j["excited_states"][top_key]["roots"].size() : 0;
          int derived = 0;
          if (j.contains("fd_states"))
            for (auto &kv : j["fd_states"].items())
              if (kv.value().contains(top_key))
                for (auto &fe : kv.value()[top_key].items())
                  if (fe.value().value("converged", false)) ++derived;
          const int want_derived = es_roots * (int)axes.size();
          print("  ES bundle  converged=", es_ok, "  roots=", n_es, "/", es_roots);
          print("  derived FD converged=", derived, "/", want_derived,
                " key=", top_key);
          bool ok = es_ok && n_es == es_roots && derived >= want_derived;
          print("\n", ok ? "PASSED" : "FAILED",
                " (ES + ", derived, "/", want_derived, " derived FD)");
          rc = ok ? 0 : 1;
        } else {
          int expected = 0, converged = 0;
          for (const auto &r : plan.fd) {
            ++expected;
            const std::string pert = r.pert.description();
            const std::string fk   = ResponseMetadata::freq_key(r.freq);
            bool ok = j["fd_states"].contains(pert) &&
                      j["fd_states"][pert].contains(top_key) &&
                      j["fd_states"][pert][top_key].contains(fk) &&
                      j["fd_states"][pert][top_key][fk].value("converged", false);
            if (ok) ++converged;
            print("  ", ok ? "[PASS]" : "[FAIL]", pert, "@", r.freq,
                  " key=", top_key);
          }
          print("\n", (converged == expected) ? "PASSED" : "FAILED", " (",
                converged, "/", expected, " FD states converged)");
          rc = (converged == expected) ? 0 : 1;
        }
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
