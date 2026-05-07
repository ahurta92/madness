// Validation harness for the v3 excited-state solver — Stage 1: TDA RHF.
//
// Reuses the protocol-setup pattern from test_solver.cpp:
//   1. Read archive header → set FunctionDefaults<3> via set_response_protocol
//   2. Load ground state → gs.prepare()
//   3. Run es_solve(world, ResponseType::TDA, num_roots, gs, ...)
//   4. Compare first N ω against the legacy H2 reference encoded in
//      MoleculeReference (matches the table in test_solver.cpp).

#include "ESSolver.hpp"
#include "ESSolverGuess.hpp"
#include "GroundState.hpp"
#include "ResponseFunctions.hpp"
#include "ResponseKernel.hpp"
#include "ResponseProtocol.hpp"

#include <madness/misc/info.h>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>
#include <madness/world/worldmem.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace madness;
using namespace molresponse_v3;

// -----------------------------------------------------------------------
// Reference data (mirrors `test_solver.cpp`).
// -----------------------------------------------------------------------
struct ESReference {
  std::string name;
  int num_roots;
  std::vector<double> omega_au;
};

// H2 (RHF). Same values as test_solver.cpp:42-43, sourced from legacy
// molresponse runs at k=8.
static const ESReference H2_ES_REF = {
    "H2", 4, {0.46535391, 0.47696774, 0.48065960, 0.48069615}};

// -----------------------------------------------------------------------
// Result reporter
// -----------------------------------------------------------------------
struct ESCheck {
  std::string name;
  bool passed;
  double expected;
  double actual;
  double tol;
  double error;
};

ESCheck check(const std::string &name, double expected, double actual,
              double tol) {
  double err = std::abs(actual - expected);
  return {name, err < tol, expected, actual, tol, err};
}

void print_check(const ESCheck &c) {
  print("  [", (c.passed ? "PASS" : "FAIL"), "]  ", c.name,
        ": expected=", c.expected, "  actual=", c.actual, "  err=", c.error,
        "  tol=", c.tol);
}

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  try {
    startup(world, argc, argv, true);

    if (world.rank() == 0) {
      print("\n============================================");
      print("  molresponse_v3 ES Solver — TDA RHF  (Stage 1)");
      print("============================================\n");
    }

    commandlineparser parser(argc, argv);
    if (!parser.key_exists("archive")) {
      if (world.rank() == 0) {
        print("Usage: test_es_solver --archive=<path> "
              "[--num-roots=N] [--maxiter=N] [--dconv=X] "
              "[--maxrotn=X] [--thresh=X] [--k=N]");
      }
      finalize();
      return 1;
    }

    // Use value_raw for paths so case is preserved (e.g. ~/Projects/H2_HF/...).
    std::string archive_path = parser.value_raw("archive");
    long num_roots = parser.key_exists("num-roots")
                         ? std::stol(parser.value("num-roots"))
                         : 4;
    long maxiter =
        parser.key_exists("maxiter") ? std::stol(parser.value("maxiter")) : 25;
    double dconv =
        parser.key_exists("dconv") ? std::stod(parser.value("dconv")) : 1e-4;
    double maxrotn =
        parser.key_exists("maxrotn") ? std::stod(parser.value("maxrotn")) : 0.5;
    double tol = parser.key_exists("tol") ? std::stod(parser.value("tol"))
                                          : 0.005; // 5 mhartree per root

    // ----- protocol from archive header (same pattern as test_solver) -----
    auto header = GroundState::read_archive_header(world, archive_path);
    double protocol_thresh = parser.key_exists("thresh")
                                 ? std::stod(parser.value("thresh"))
                                 : default_thresh_for_k(header.k);
    int override_k = parser.key_exists("k") ? std::stoi(parser.value("k")) : -1;

    if (world.rank() == 0) {
      print("Archive  k=", header.k, "  L=", header.L,
            "  protocol_thresh=", protocol_thresh);
      print("Solver   num_roots=", num_roots, "  maxiter=", maxiter,
            "  dconv=", dconv, "  maxrotn=", maxrotn);
    }
    set_response_protocol(world, header.L, protocol_thresh, override_k);

    // ----- molecule from calc_info.json -----
    Molecule molecule;
    auto archive_dir = std::filesystem::path(archive_path).parent_path();
    for (const auto &name : {"moldft.calc_info.json", "mad.calc_info.json"}) {
      auto candidate = archive_dir / name;
      if (std::filesystem::exists(candidate)) {
        std::ifstream ifs(candidate);
        nlohmann::json j;
        ifs >> j;
        nlohmann::json mol_json;
        if (j.contains("tasks") && j["tasks"].is_array() &&
            !j["tasks"].empty()) {
          mol_json = j["tasks"][0]["molecule"];
        } else if (j.contains("molecule")) {
          mol_json = j["molecule"];
        }
        if (!mol_json.is_null())
          molecule.from_json(mol_json);
        break;
      }
    }

    // ----- load ground state and prepare -----
    auto gs = GroundState::from_archive(world, archive_path, molecule);
    if (world.rank() == 0)
      gs.print_info();

    std::string fock_json;
    for (const auto &name : {"moldft.fock.json", "mad.fock.json"}) {
      auto candidate = archive_dir / name;
      if (std::filesystem::exists(candidate)) {
        fock_json = candidate.string();
        if (world.rank() == 0)
          print("Using Fock from:", fock_json);
        break;
      }
    }
    const double cur_thresh = FunctionDefaults<3>::get_thresh();
    auto coulop = poperatorT(
        CoulombOperatorPtr(world, gs.params().lo(), 0.001 * cur_thresh));
    gs.prepare(world, 0.001 * cur_thresh, coulop, fock_json);

    // ----- pick reference based on molecule -----
    const ESReference *ref = nullptr;
    if (gs.num_alpha() == 1 && gs.is_spin_restricted()) {
      ref = &H2_ES_REF;
    }
    if (!ref) {
      if (world.rank() == 0) {
        print("\nNo ES reference available for this molecule "
              "(num_alpha=",
              gs.num_alpha(), "  restricted=", gs.is_spin_restricted(), ")");
        print("Stage 1 supports only H2 RHF (num_alpha=1).");
      }
      finalize();
      return 1;
    }

    if (world.rank() == 0) {
      print("\nReference: ", ref->name, "  (legacy first ", ref->num_roots,
            " roots, k=8)");
      for (int i = 0; i < ref->num_roots; i++)
        print("  ω[", i, "]_ref =", ref->omega_au[i]);
      print("  per-root tolerance =", tol, " hartree\n");
    }

    // ----- solve -----
    if (world.rank() == 0)
      print("=== ES TDA RHF solve ===");
    auto solve = es_solve(world, ResponseType::TDA, num_roots, gs, maxiter,
                          dconv, maxrotn, PrintLevel::Verbose);

    // ----- compare first ref->num_roots vs reference -----
    std::vector<ESCheck> results;
    if (world.rank() == 0) {
      print("\n=== Final excitation energies ===");
      for (long s = 0; s < num_roots; s++) {
        print("  ω[", s, "] =", solve.omegas(s), "  res=", solve.residuals[s]);
      }
    }
    for (int i = 0; i < ref->num_roots && i < num_roots; i++) {
      results.push_back(check("omega_root_" + std::to_string(i),
                              ref->omega_au[i], solve.omegas(i), tol));
    }
    results.push_back(
        check("converged", 1.0, solve.converged ? 1.0 : 0.0, 0.5));

    if (world.rank() == 0) {
      print("\n=== Test Results vs Legacy Reference ===");
      int passed = 0, failed = 0;
      for (const auto &r : results) {
        print_check(r);
        (r.passed ? passed : failed)++;
      }
      print("\n  ", passed, " passed, ", failed, " failed\n");
    }

    bool all_passed = true;
    for (const auto &r : results)
      if (!r.passed)
        all_passed = false;

    world.gop.fence();
    finalize();
    return all_passed ? 0 : 1;
  } catch (const std::exception &e) {
    if (world.rank() == 0)
      print("ERROR:", e.what());
    finalize();
    return 2;
  }
}
