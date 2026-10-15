#ifndef MOLRESPONSE_V3_SOLVERS_RAMAN_ANALYSIS_HPP
#define MOLRESPONSE_V3_SOLVERS_RAMAN_ANALYSIS_HPP

// ===========================================================================
// Full vibrational Raman (post-solve, closed shell): the normal-mode analysis
// that turns the per-coordinate polarizability derivatives assemble_beta
// writes under properties/raman into a spectrum.
//
//   1. Hessian from the converged static nuclear-displacement responses
//      x^{(a,q)} (the same states the Raman sources used; no extra solve):
//        H_ij = < rho_i | dV/dQ_j >  (+ < d_axis rho_0 | dV/dQ_j > on the same
//        atom, the integration-by-parts term),   rho_i = 4 sum_k x^i_k phi_k.
//      The noisy diagonal is skipped and rebuilt from translational
//      invariance, the matrix is purified (heavier atom's element wins) and
//      symmetrized, and the nuclear-repulsion Hessian is added -- the recipe
//      of Nemo::hessian (chem/nemo.cc) and of the removed molresponse_v2
//      PropertyManager::compute_hessian, written for v3's ResponseStateXY.
//   2. Projected frequencies, normal modes and reduced masses from
//      chem/vibanal.h; Cartesian displacements per unit normal coordinate
//      L = M^{-1/2} * modes (Molecule::massweights()).
//   3. The 9 x 3N derivative table per optical frequency (from the raman rows
//      at this protocol), projected on L, reduced to alpha'^2, beta'^2,
//      activity and depolarization ratio (kernels/raman_spectrum.hpp).
//
// Records "vibrational_analysis" and "raman_spectra" in response_metadata.json
// (mirrored into the madqc task record), with units stated. Collective
// discipline: every MRA inner product runs on all ranks; JSON is rank 0.
// ===========================================================================

#include "../calc/calc_executor.hpp"      // ExecutorContext, detail_exec::load_fd_as_xy
#include "../kernels/raman_spectrum.hpp"
#include "../ResponsePropertyPlanner.hpp"
#include "response_metadata.hpp"

#include <madness/chem/molecular_functors.h>   // MolecularDerivativeFunctor
#include <madness/chem/vibanal.h>              // compute_frequencies, compute_reduced_mass
#include <madness/constants.h>
#include <madness/mra/mra.h>
#include <madness/tensor/tensor_json.hpp>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <optional>
#include <string>
#include <vector>

namespace molresponse_v3 {

/// True when the plan carries the all-atoms nuclear sentinel, i.e. the
/// request asked for every coordinate (raman.nuc_atom -1).
inline bool plan_requests_full_raman(const ResponsePlan &plan) {
  for (const auto &r : plan.nuclear_fd)
    if (r.pert.is_all_atoms()) return true;
  return false;
}

struct VibrationalAnalysis {
  madness::Tensor<double> hessian;           ///< 3N x 3N, a.u. (Eh / bohr^2)
  madness::Tensor<double> frequencies;       ///< 3N projected frequencies, a.u.
  madness::Tensor<double> normalmodes;       ///< 3N x 3N mass-weighted eigenvectors
  madness::Tensor<double> normalmodes_atomic;///< 3N x 3N Cartesian displacement per unit normal coordinate
  madness::Tensor<double> reducedmass;       ///< 3N, a.u. (m_e)
  double max_asymmetry = 0.0;                ///< largest |H_ij - H_ji|/2 before symmetrization
  std::vector<std::string> missing_legs;     ///< nuclear legs not on disk (empty on success)
};

/// Hessian + normal modes from the static nuclear-displacement responses in
/// `calc_dir` at the ACTIVE protocol (caller has set it and prepared gs).
/// Returns nullopt, with the missing leg names printed, if any leg is absent.
inline std::optional<VibrationalAnalysis>
compute_vibrational_analysis(ExecutorContext &ctx) {
  using namespace madness;
  World &world = ctx.world;
  GroundState &gs = ctx.gs;
  const Molecule &mol = gs.molecule();
  const long natom = static_cast<long>(mol.natom());
  const long dim = 3 * natom;

  VibrationalAnalysis out;
  const vecfuncT &phi0 = gs.orbitals_alpha();

  // Nuclear response densities rho_i = 4 sum_k x_k phi_k (closed shell, static).
  std::vector<real_function_3d> rho(static_cast<size_t>(dim));
  for (long a = 0; a < natom; ++a)
    for (int ax = 0; ax < 3; ++ax) {
      const Perturbation p = Perturbation::nuclear(static_cast<int>(a), ax);
      auto st = detail_exec::load_fd_as_xy<ClosedShell>(world, ctx.calc_dir, p, 0.0);
      if (!st) { out.missing_legs.push_back(p.description()); continue; }
      auto xphi = mul(world, st->x_alpha, phi0, true);
      real_function_3d r = sum(world, xphi, true);
      r.scale(4.0);
      r.truncate();
      rho[static_cast<size_t>(3 * a + ax)] = r;
    }
  if (!out.missing_legs.empty()) {
    if (world.rank() == 0) {
      print("[RAMAN] vibrational analysis skipped — nuclear legs missing at",
            protocol_key(), ":");
      for (const auto &m : out.missing_legs) print("[RAMAN]   ", m);
    }
    return std::nullopt;
  }

  // Ground-state density and its Cartesian derivatives (integration-by-parts term).
  real_function_3d rho0 =
      gs.scf().make_density(world, gs.occupations_alpha(), phi0);
  rho0.scale(2.0);
  std::vector<real_function_3d> drho0(3);
  for (int ax = 0; ax < 3; ++ax) {
    real_derivative_3d D = free_space_derivative<double, 3>(world, ax);
    real_function_3d rc = copy(rho0).refine();
    drho0[static_cast<size_t>(ax)] = D(rc);
  }

  Tensor<double> H(dim, dim);
  for (long iatom = 0; iatom < natom; ++iatom)
    for (long iaxis = 0; iaxis < 3; ++iaxis) {
      const long i = 3 * iatom + iaxis;
      for (long jatom = 0; jatom < natom; ++jatom)
        for (long jaxis = 0; jaxis < 3; ++jaxis) {
          const long j = 3 * jatom + jaxis;
          if (i == j) continue;   // noisy; rebuilt from translational invariance
          madchem::MolecularDerivativeFunctor mdf(mol, static_cast<int>(jatom),
                                                  static_cast<int>(jaxis));
          H(i, j) = inner(rho[static_cast<size_t>(i)], mdf);
          if (iatom == jatom) H(i, j) += inner(drho0[static_cast<size_t>(iaxis)], mdf);
        }
    }
  for (long i = 0; i < dim; ++i) H(i, i) = 0.0;

  // Purify: for unlike atoms keep the element whose response density belongs
  // to the heavier atom (the better-resolved one); for like atoms average.
  {
    Tensor<double> P = copy(H);
    for (long iatom = 0; iatom < natom; ++iatom)
      for (long iaxis = 0; iaxis < 3; ++iaxis) {
        const long i = 3 * iatom + iaxis;
        for (long jatom = 0; jatom < natom; ++jatom)
          for (long jaxis = 0; jaxis < 3; ++jaxis) {
            const long j = 3 * jatom + jaxis;
            const double mean = 0.5 * (P(i, j) + P(j, i));
            out.max_asymmetry = std::max(out.max_asymmetry, 0.5 * std::abs(P(i, j) - P(j, i)));
            const unsigned ZA = mol.get_atomic_number(static_cast<size_t>(iatom));
            const unsigned ZB = mol.get_atomic_number(static_cast<size_t>(jatom));
            if (ZA < ZB) P(i, j) = P(j, i);
            else if (ZA > ZB) P(j, i) = P(i, j);
            else { P(i, j) = mean; P(j, i) = mean; }
          }
      }
    H = P;
  }
  H += transpose(H);
  H.scale(0.5);
  // Translational invariance fixes the diagonal: sum over atoms of H(i, j+axis_i) = 0.
  for (long i = 0; i < dim; ++i) {
    double s = 0.0;
    for (long j = 0; j < dim; j += 3) s += H(i, j + (i % 3));
    H(i, i) = -s;
  }
  H += mol.nuclear_repulsion_hessian();

  Tensor<double> modes;
  out.frequencies = compute_frequencies(mol, H, modes, /*project_tr=*/true, /*print=*/false);
  out.normalmodes = modes;
  out.reducedmass = compute_reduced_mass(mol, modes);
  out.normalmodes_atomic = inner(mol.massweights(), modes);   // M^{-1/2} * eigenvectors
  out.hessian = H;

  if (world.rank() == 0) {
    print("[RAMAN] Hessian from", dim, "nuclear response states; max asymmetry",
          out.max_asymmetry);
    printf("[RAMAN] projected frequencies (cm^-1):");
    for (long i = 0; i < dim; ++i) printf(" %10.3f", out.frequencies(i) * constants::au2invcm);
    printf("\n");
    fflush(stdout);
  }
  return out;
}

/// Full Raman assembly at protocol `thresh`: Hessian + modes (collective), then
/// the per-mode spectrum from the raman rows recorded at this protocol (rank 0),
/// written as the "vibrational_analysis" and "raman_spectra" sections.
inline void assemble_raman_spectrum(ExecutorContext &ctx, const ResponsePlan &plan,
                                    double thresh) {
  using namespace madness;
  World &world = ctx.world;
  GroundState &gs = ctx.gs;
  if (!plan_requests_full_raman(plan)) return;
  if (!gs.is_spin_restricted()) {
    if (world.rank() == 0)
      print("[RAMAN] open-shell vibrational analysis is not implemented — SKIPPED.");
    return;
  }

  set_response_protocol(world, ctx.L, thresh);
  const double t0 = FunctionDefaults<3>::get_thresh();
  {
    auto coulop = poperatorT(CoulombOperatorPtr(world, gs.params().lo(), 0.001 * t0));
    gs.prepare(world, 0.001 * t0, coulop, ctx.fock_json);
  }
  const std::string key = protocol_key();
  if (world.rank() == 0)
    print("\n=== Raman vibrational analysis  protocol_key=", key, " ===");

  auto vib = compute_vibrational_analysis(ctx);
  world.gop.fence();
  if (!vib) return;
  if (world.rank() != 0) return;   // the rest is rank-0 tensor algebra + JSON

  const Molecule &mol = gs.molecule();
  const long natom = static_cast<long>(mol.natom());
  const long dim = 3 * natom;
  constexpr double tr_cutoff_cm1 = 1.0;
  const auto modes = raman::vibrational_modes(vib->frequencies, tr_cutoff_cm1);

  auto meta = ResponseMetadata::load_or_create(ctx.calc_dir + "/response_metadata.json");
  const auto &j = meta.json();

  // --- vibrational_analysis section --------------------------------------
  nlohmann::json va;
  va["protocol_key"] = key;
  va["n_atoms"] = natom;
  va["hessian"] = tensor_to_json(vib->hessian);
  va["frequencies"] = tensor_to_json(vib->frequencies);
  {
    Tensor<double> fcm = copy(vib->frequencies);
    fcm.scale(constants::au2invcm);
    va["frequencies_cm1"] = tensor_to_json(fcm);
  }
  va["normalmodes"] = tensor_to_json(vib->normalmodes);
  va["normalmodes_atomic"] = tensor_to_json(vib->normalmodes_atomic);
  va["reducedmass"] = tensor_to_json(vib->reducedmass);
  va["au2invcm"] = constants::au2invcm;
  va["max_hessian_asymmetry"] = vib->max_asymmetry;
  va["vibrational_mode_indices"] = modes;
  va["tr_cutoff_cm1"] = tr_cutoff_cm1;
  va["source"] = "Hessian from the static nuclear-displacement response states at this protocol; "
                 "diagonal from translational invariance; nuclear-repulsion Hessian added";
  va["units"] = {{"hessian", "Eh/bohr^2"}, {"frequencies", "hartree"},
                 {"frequencies_cm1", "cm^-1"},
                 {"normalmodes", "mass-weighted eigenvectors (dimensionless)"},
                 {"normalmodes_atomic", "bohr per unit mass-weighted normal coordinate (M^-1/2 * mode)"},
                 {"reducedmass", "m_e"}};

  // --- raman_spectra section ---------------------------------------------
  nlohmann::json rs;
  rs["protocol_key"] = key;
  rs["tr_cutoff_cm1"] = tr_cutoff_cm1;
  nlohmann::json vfreq = nlohmann::json::array();
  for (long m : modes) vfreq.push_back(vib->frequencies(m) * constants::au2invcm);
  rs["vibrational_frequencies_cm1"] = vfreq;
  rs["row_order"] = "row 3*A+B for dipole axes A,B in xyz; column 3*atom+axis";
  rs["units"] = {{"polarizability_derivatives", "a.u. (e^2 bohr^2 / Eh per bohr)"},
                 {"polarizability_derivatives_normal_modes", "a.u. per unit mass-weighted normal coordinate"},
                 {"alpha2_au", "a.u."}, {"beta2_au", "a.u."}, {"activity_au", "a.u. (45 alpha'^2 + 4 beta'^2)"},
                 {"alpha2", "Angstrom^4/amu"}, {"beta2", "Angstrom^4/amu"},
                 {"pol_int", "Angstrom^4/amu (activity, linear polarization)"},
                 {"depol_int", "Angstrom^4/amu (dep_ratio * pol_int)"},
                 {"dep_ratio", "dimensionless (3 beta'^2 / activity)"},
                 {"freq_cm1", "cm^-1"}};
  rs["au_to_A4_amu"] = raman::AU_TO_A4_AMU;

  std::vector<double> pol_freqs;
  for (const auto &r : plan.vbc)
    if (r.pert_c.is_all_atoms() &&
        std::find_if(pol_freqs.begin(), pol_freqs.end(),
                     [&](double f) { return std::abs(f - r.freq_b) < 1e-9; }) == pol_freqs.end())
      pol_freqs.push_back(r.freq_b);
  std::sort(pol_freqs.begin(), pol_freqs.end());
  rs["polarization_frequencies"] = pol_freqs;

  const bool have_rows = j.contains("properties") && j["properties"].contains("raman") &&
                         j["properties"]["raman"].contains(key);
  nlohmann::json derivs = nlohmann::json::object(), derivs_modes = nlohmann::json::object(),
                 spectra = nlohmann::json::object(), complete = nlohmann::json::object();
  static const char *axn = "xyz";
  for (double w : pol_freqs) {
    const std::string wkey = ResponseMetadata::freq_key(w);
    Tensor<double> D(9L, dim);
    long found = 0; bool all_conv = true;
    if (have_rows) {
      for (const auto &row : j["properties"]["raman"][key]) {
        if (std::abs(row.value("freq_b", -1.0) - w) > 1e-9) continue;
        if (std::abs(row.value("freq_c", -1.0)) > 1e-9) continue;
        const std::string A = row.value("A", ""), B = row.value("B", ""), C = row.value("C", "");
        if (A.size() != 1 || B.rfind("dipole_", 0) != 0 || C.rfind("nuc_", 0) != 0) continue;
        const int ia = static_cast<int>(std::string(axn).find(A[0]));
        const int ib = static_cast<int>(std::string(axn).find(B.back()));
        // C = nuc_<atom>_<axis>
        const auto us = C.find('_', 4);
        if (ia < 0 || ib < 0 || us == std::string::npos) continue;
        const int atom = std::stoi(C.substr(4, us - 4));
        const int ic = static_cast<int>(std::string(axn).find(C.back()));
        if (atom < 0 || atom >= natom || ic < 0) continue;
        D(3 * ia + ib, 3 * atom + ic) = row.value("beta", 0.0);
        ++found;
        all_conv = all_conv && row.value("converged", false);
      }
    }
    complete[wkey] = {{"rows_found", found}, {"rows_expected", 9 * dim},
                      {"all_rows_converged", all_conv}};
    if (found != 9 * dim) {
      print("[RAMAN] omega=", w, ": only", found, "of", 9 * dim,
            "derivative rows at", key, "— spectrum not assembled for this frequency");
      continue;
    }
    derivs[wkey] = tensor_to_json(D);
    if (modes.empty()) continue;
    Tensor<double> L(dim, static_cast<long>(modes.size()));
    for (size_t q = 0; q < modes.size(); ++q)
      for (long k = 0; k < dim; ++k) L(k, static_cast<long>(q)) = vib->normalmodes_atomic(k, modes[q]);
    const Tensor<double> Dq = raman::project_on_modes(D, L);
    derivs_modes[wkey] = tensor_to_json(Dq);
    nlohmann::json rows = nlohmann::json::array();
    print("[RAMAN] omega=", w, " mode   freq(cm^-1)   activity(A^4/amu)   depol   alpha'^2   beta'^2 (a.u.)");
    for (size_t q = 0; q < modes.size(); ++q) {
      const auto a = raman::mode_tensor(Dq, static_cast<long>(q));
      const auto inv = raman::invariants(a);
      const double fcm = vib->frequencies(modes[q]) * constants::au2invcm;
      rows.push_back({{"mode", static_cast<int>(q) + 1},
                      {"freq_cm1", fcm},
                      {"alpha2_au", inv.alpha2}, {"beta2_au", inv.beta2},
                      {"activity_au", inv.activity},
                      {"alpha2", inv.alpha2 * raman::AU_TO_A4_AMU},
                      {"beta2", inv.beta2 * raman::AU_TO_A4_AMU},
                      {"pol_int", inv.activity * raman::AU_TO_A4_AMU},
                      {"dep_ratio", inv.depolarization},
                      {"depol_int", inv.depolarization * inv.activity * raman::AU_TO_A4_AMU},
                      {"converged", all_conv}});
      printf("[RAMAN]   %3d  %12.3f  %16.4f  %7.4f  %10.5f %10.5f\n",
             static_cast<int>(q) + 1, fcm, inv.activity * raman::AU_TO_A4_AMU,
             inv.depolarization, inv.alpha2, inv.beta2);
    }
    fflush(stdout);
    spectra[wkey] = rows;
  }
  rs["polarizability_derivatives"] = derivs;
  rs["polarizability_derivatives_normal_modes"] = derivs_modes;
  rs["raman_spectra"] = spectra;
  rs["completeness"] = complete;

  meta.set_section("vibrational_analysis", va);
  meta.set_section("raman_spectra", rs);
  meta.save();
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_RAMAN_ANALYSIS_HPP
