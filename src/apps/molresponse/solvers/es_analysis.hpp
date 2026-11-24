#ifndef MOLRESPONSE_V3_SOLVERS_ES_ANALYSIS_HPP
#define MOLRESPONSE_V3_SOLVERS_ES_ANALYSIS_HPP

// =========================================================================
// Post-convergence excited-state analysis (closed-shell).
//
// Ports the legacy TDDFT::analysis + TDDFT::analyze_vectors
// (src/apps/molresponse_legacy/iterate_excited.cc, TDDFT.cc) onto the v3
// ESSolver<Type, ClosedShell>::State. For each converged root it computes
// the transition properties expected of an excited-state report:
//
//   * transition dipole moments      d_a   = Σ_j <φ_j | r_a · (x_j [+ y_j])>
//   * dipole oscillator strength      f     = (2/3) |d|² ω
//   * transition quadrupole moments   Q_ab  = Σ_j <φ_j | r_a r_b · (x_j[+y_j])>
//   * dominant occupied contributions       (per-occupied |x_j| weights)
//
// Type selects the transition-moment assembly: TDA uses X only; Full adds
// the Y block (matching the legacy full-response sign convention — Y enters
// the moments with the same sign as X). Open-shell ES is out of scope, so
// these helpers are gated on ClosedShell.
//
// Collective discipline: every inner-product / norm here is collective on
// `world` and runs on all ranks; only the human-readable printing is
// rank-0-gated. Output is print-only for now — a dedicated chemistry
// property layer (JSON / HDF5 / cube) will consume the returned structs
// later, which is why analyze_es_bundle returns data rather than printing
// it directly.
// =========================================================================

#include "../GroundState.hpp"
#include "../kernels/tags.hpp"   // PrintLevel
#include "es_solver.hpp"

#include <madness/chem/SCF.h>            // MomentFunctor, DipoleFunctor, SCF::project_ao_basis_only
#include <madness/mra/mra.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <type_traits>
#include <vector>

namespace molresponse_v3 {

/// Per-root transition properties (closed-shell). Mirrors the data the
/// legacy TDDFT::analysis printed for one excited-state root.
struct ESRootAnalysis {
  bool dipole_available = false;
  bool quadrupole_available = false;
  double                                omega = 0.0;
  std::array<double, 3>                 dipole{{0.0, 0.0, 0.0}};
  double                                oscillator = 0.0;
  std::array<std::array<double, 3>, 3>  quadrupole{};
  /// |x_alpha[j]| per occupied orbital j (used for dominant-contribution
  /// reporting). Y is intentionally not folded in here — the legacy report
  /// listed the X and Y occupied weights side by side, so the Full printer
  /// recomputes the Y weights when it needs them.
  std::vector<double>                   orbital_norms;
};

namespace detail_es_analysis {

/// Cartesian coordinate functions x, y, z (degree-1 moments). Shared by the
/// transition-moment assembly and the per-orbital analyzer.
inline std::array<madness::real_function_3d, 3>
coordinate_functions(madness::World &world) {
  using namespace madness;
  const double tol = FunctionDefaults<3>::get_thresh();
  std::array<real_function_3d, 3> r;
  for (int a = 0; a < 3; ++a) {
    std::vector<int> dir(3, 0);
    dir[a] = 1;
    r[a] = real_factory_3d(world).functor(
        real_functor_3d(new MomentFunctor(dir)));
    r[a].truncate(tol);
  }
  return r;
}

}  // namespace detail_es_analysis

/// Compute transition properties for every root in a converged closed-shell
/// ES bundle. `ground_orbitals` are the occupied alpha orbitals
/// (gs.orbitals_alpha()). Pure computation — collective on all ranks, no
/// printing.
template <typename Type, typename Shell>
std::vector<ESRootAnalysis>
analyze_es_bundle(madness::World &world,
                  const std::vector<madness::real_function_3d> &ground_orbitals,
                  const typename ESSolver<Type, Shell>::State &state) {
  using namespace madness;
  static_assert(std::is_same_v<Shell, ClosedShell>,
                "analyze_es_bundle: closed-shell only");
  constexpr bool has_y = std::is_same_v<Type, Full>;

  const size_t m = state.roots.size();
  const size_t n = ground_orbitals.size();

  // ---- Bisection probe (TEMPORARY) ---------------------------------------
  // MADQC_ANALYZE_PROBE caps how far into the analysis we run, to localize the
  // teardown crash with one rebuild + a few quick runs:
  //   0  omega only (no MADNESS ops)
  //   1  + create coordinate functions (MomentFunctor over the box)
  //   2  + deep-copy + reconstruct the ground orbitals
  //   3  + norm2 + transition-dipole loop (inner/mul)
  //   4+ full (adds the quadrupole loop)            [default]
  int probe = 4;
  if (const char *p = std::getenv("MADQC_ANALYZE_PROBE")) probe = std::atoi(p);
  if (world.rank() == 0)
    madness::print("[ANALYZE] probe stage =", probe);

  std::vector<ESRootAnalysis> out(m);
  for (size_t i = 0; i < m; ++i)
    out[i].omega =
        (static_cast<long>(i) < state.omega.size()) ? state.omega(i) : 0.0;
  if (probe <= 0) return out;

  const auto r = detail_es_analysis::coordinate_functions(world);
  if (probe <= 1) return out;

  // Operate on a DEEP COPY of the ground orbitals — never gs's live handles.
  vecfuncT orbs = madness::copy(world, ground_orbitals);
  madness::reconstruct(world, orbs);
  if (probe <= 2) return out;

  // Stage 3: norms + transition dipole.
  for (size_t i = 0; i < m; ++i) {
    ESRootAnalysis &a = out[i];
    const auto &X = state.roots[i].x_alpha;

    a.orbital_norms.assign(n, 0.0);
    for (size_t j = 0; j < n; ++j) a.orbital_norms[j] = X[j].norm2();

    for (size_t j = 0; j < n; ++j) {
      for (int ax = 0; ax < 3; ++ax)
        a.dipole[ax] += inner(orbs[j], r[ax] * X[j]);
      if constexpr (has_y) {
        const auto &Y = state.roots[i].y_alpha;
        for (int ax = 0; ax < 3; ++ax)
          a.dipole[ax] += inner(orbs[j], r[ax] * Y[j]);
      }
    }
    for (auto &d : a.dipole) d *= -std::sqrt(2.0);

    a.oscillator = 2.0 / 3.0 *
                   (a.dipole[0] * a.dipole[0] + a.dipole[1] * a.dipole[1] +
                    a.dipole[2] * a.dipole[2]) *
                   a.omega;
    a.dipole_available = true;
  }
  if (probe <= 3) return out;

  // Stage 4: transition quadrupole.
  for (size_t i = 0; i < m; ++i) {
    ESRootAnalysis &a = out[i];
    const auto &X = state.roots[i].x_alpha;
    for (size_t j = 0; j < n; ++j) {
      for (int aa = 0; aa < 3; ++aa)
        for (int bb = 0; bb < 3; ++bb) {
          a.quadrupole[aa][bb] +=
              inner(orbs[j], r[aa] * r[bb] * X[j]);
          if constexpr (has_y) {
            const auto &Y = state.roots[i].y_alpha;
            a.quadrupole[aa][bb] +=
                inner(orbs[j], r[aa] * r[bb] * Y[j]);
          }
        }
    }
    for (auto &row : a.quadrupole)
      for (auto &q : row) q *= std::sqrt(2.0);
    a.quadrupole_available = true;
  }
  return out;
}

/// Print the transition-property block for a bundle (rank 0 only). `tda`
/// controls only the dominant-contribution layout (X-only vs X/Y columns);
/// when `tda` is false, `y_orbital_norms[i]` supplies the per-root Y weights.
inline void
print_es_analysis(madness::World &world,
                  const std::vector<ESRootAnalysis> &an, bool tda,
                  const std::vector<std::vector<double>> &y_orbital_norms = {}) {
  using namespace madness;
  if (world.rank() != 0) return;

  for (size_t i = 0; i < an.size(); ++i) {
    const ESRootAnalysis &a = an[i];
    const size_t n = a.orbital_norms.size();

    printf("   Response Function %d\t\t%7.8f a.u.", static_cast<int>(i), a.omega);
    print("\n   --------------------------------------------");
    printf("   Response Function %d\t\t%7.8f eV", static_cast<int>(i),
           a.omega * 27.2114);
    print("\n   --------------------------------------------");

    if (a.dipole_available) {
    print("\n   Transition Dipole Moments (a.u., e a0)");
    printf("   X: %7.8f   Y: %7.8f   Z: %7.8f\n", a.dipole[0], a.dipole[1],
           a.dipole[2]);

    printf("\n   Dipole Oscillator Strength: %7.8f\n", a.oscillator);

    } else print("   Transition dipole/oscillator strength: unavailable (analysis diagnostic)");

    if (a.quadrupole_available) {
    print("\n   Transition Quadrupole Moments (a.u., e a0^2; Cartesian r_a r_b, "
          "not traceless)");
    printf("   %16s %16s %16s\n", "X", "Y", "Z");
    printf("   X %16.8f %16.8f %16.8f\n", a.quadrupole[0][0], a.quadrupole[0][1],
           a.quadrupole[0][2]);
    printf("   Y %16.8f %16.8f %16.8f\n", a.quadrupole[1][0], a.quadrupole[1][1],
           a.quadrupole[1][2]);
    printf("   Z %16.8f %16.8f %16.8f\n", a.quadrupole[2][0], a.quadrupole[2][1],
           a.quadrupole[2][2]);

    } else print("   Transition quadrupole: unavailable (analysis diagnostic)");

    // Dominant contributions: top 5 occupied orbitals by |x_j| (legacy order).
    std::vector<size_t> order(n);
    for (size_t j = 0; j < n; ++j) order[j] = j;
    std::sort(order.begin(), order.end(), [&](size_t p, size_t q) {
      return a.orbital_norms[p] > a.orbital_norms[q];
    });

    if (tda) {
      print("\n   Dominant Contributions:");
      for (size_t j = 0; j < std::min(size_t(5), n); ++j)
        printf("   Occupied %d   %7.8f\n", static_cast<int>(order[j]),
               a.orbital_norms[order[j]]);
      print("\n");
    } else {
      const std::vector<double> &yn =
          (i < y_orbital_norms.size()) ? y_orbital_norms[i]
                                       : std::vector<double>(n, 0.0);
      print("\n   Dominant Contributions:");
      print("                  x          y");
      for (size_t j = 0; j < std::min(size_t(5), n); ++j) {
        const size_t o = order[j];
        printf("   Occupied %d   %7.8f %7.8f\n", static_cast<int>(o),
               a.orbital_norms[o], (o < yn.size()) ? yn[o] : 0.0);
      }
      print("\n");
    }
  }
}

/// Serialize a bundle's transition properties (omega, transition dipole,
/// oscillator strength, transition quadrupole, dominant occupied weights) to
/// JSON. This is the persistent, machine-readable record for later analysis;
/// the human-readable console report is print_es_analysis. `stable_index`
/// (if non-empty) attaches each root's stable identity / root_id.
inline nlohmann::json
es_analysis_to_json(const std::vector<ESRootAnalysis> &an, bool tda,
                    const std::vector<std::vector<double>> &y_norms,
                    const std::vector<int> &stable_index) {
  nlohmann::json roots = nlohmann::json::array();
  for (size_t i = 0; i < an.size(); ++i) {
    const ESRootAnalysis &a = an[i];
    nlohmann::json r;
    r["slot"]    = static_cast<int>(i);
    if (i < stable_index.size()) {
      r["stable_index"] = stable_index[i];
      r["root_id"]      = make_root_id(stable_index[i]);
    }
    r["omega"]                 = a.omega;
    r["omega_ev"]              = a.omega * 27.2114;
    if (a.dipole_available) {
    r["transition_dipole"]     = {a.dipole[0], a.dipole[1], a.dipole[2]};
    r["oscillator_strength"]   = a.oscillator;
    } else r["availability"]["transition_dipole_and_strength"] = "analysis diagnostic";
    if (a.quadrupole_available) {
    r["transition_quadrupole"] = {
        {a.quadrupole[0][0], a.quadrupole[0][1], a.quadrupole[0][2]},
        {a.quadrupole[1][0], a.quadrupole[1][1], a.quadrupole[1][2]},
        {a.quadrupole[2][0], a.quadrupole[2][1], a.quadrupole[2][2]}};
    } else r["availability"]["transition_quadrupole"] = "analysis diagnostic";
    r["occupied_weights_x"]    = a.orbital_norms;
    if (!tda && i < y_norms.size()) r["occupied_weights_y"] = y_norms[i];
    roots.push_back(r);
  }
  nlohmann::json j;
  j["analysis"]      = tda ? "tda" : "full";
  j["n_roots"]       = static_cast<int>(an.size());
  j["protocol_key"]  = protocol_key();
  j["units"]         = {{"omega", "hartree"}, {"omega_ev", "eV"},
                        {"transition_dipole", "a.u. (e a0)"},
                        {"transition_quadrupole",
                         "a.u. (e a0^2); Cartesian r_a r_b, not traceless"}};
  j["roots"]         = roots;
  return j;
}

}  // namespace molresponse_v3

#endif  // MOLRESPONSE_V3_SOLVERS_ES_ANALYSIS_HPP
