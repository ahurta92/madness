#ifndef MADNESS_CHEM_DALTON_SEED_H
#define MADNESS_CHEM_DALTON_SEED_H

// Seed an SCF from a DALTON calculation: project the occupied molden orbitals
// onto the MRA basis at the run's box and first protocol threshold,
// Loewdin-orthonormalize them, and write them as <prefix>.restartdata with a
// RestartMetadata header, so `restart auto` resumes from them as from any other
// archive. The orbitals are a starting guess, never a result: the header records
// no dconv convergence, so the restart plan iterates.
//
// Why Loewdin: independently projected Gaussian MOs keep O(1e-4) mutual
// overlaps that SCF::load_mos does not remove, and a non-orthonormal first
// density sent water to -106.9 Ha.
//
//   write_gs_seed_from_molden: the projection and the archive;
//   seed_gs_from_dalton_dir:   the madqc pre-run hook (locate the DALTON files,
//                              check the geometry, count electrons, seed).

#include <madness/chem/ParameterManager.hpp>
#include <madness/chem/Restart.h>
#include <madness/chem/dalton_dir.h>
#include <madness/chem/molden_gto.h>
#include <madness/chem/molecule.h>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>
#include <madness/world/parallel_archive.h>
#include <madness/world/vector_archive.h>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace madness {

namespace dalton_seed_detail {

/// One DALTON MO (its AO-coefficient column) as a real MRA functor.
class DaltonMOFunctor : public FunctionFunctorInterface<double, 3> {
  const DaltonMoldenBasis &basis;
  std::vector<double> c;
  std::vector<coord_3d> centers;
public:
  DaltonMOFunctor(const DaltonMoldenBasis &b, std::vector<double> w)
      : basis(b), c(std::move(w)) {
    for (const auto &sh : basis.shells) centers.push_back({sh.cx, sh.cy, sh.cz});
  }
  double operator()(const coord_3d &r) const override {
    double val = 0.0, bf[9];
    for (size_t s = 0; s < basis.shells.size(); s++) {
      const auto &sh = basis.shells[s];
      sh.evaluate(r[0], r[1], r[2], bf);
      const int off = basis.ao_offsets[s];
      for (int k = 0; k < sh.n_ao; k++) val += c[static_cast<size_t>(off + k)] * bf[k];
    }
    return val;
  }
  std::vector<coord_3d> special_points() const override { return centers; }
};

/// Element symbol (molden labels carry index suffixes: H1, H2, O) -> Z.
inline int symbol_to_Z(const std::string &s) {
  std::string elem;
  for (char c : s) { if (std::isalpha(static_cast<unsigned char>(c))) elem += c; else break; }
  static const std::vector<std::pair<std::string, int>> tab = {
      {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5}, {"C", 6},
      {"N", 7}, {"O", 8}, {"F", 9}, {"Ne", 10}, {"Na", 11}, {"Mg", 12},
      {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16}, {"Cl", 17}, {"Ar", 18}};
  for (const auto &[sym, z] : tab) if (sym == elem) return z;
  throw std::runtime_error("dalton gs seed: unknown element symbol '" + s + "'");
}

/// The wavelet order for a projection threshold: the table SCF::set_protocol
/// uses, so the seed is projected at the k the SCF will run at.
inline int k_for_thresh(double thresh) {
  if (thresh >= 0.9e-2) return 4;
  if (thresh >= 0.9e-4) return 6;
  if (thresh >= 0.9e-6) return 8;
  if (thresh >= 0.9e-8) return 10;
  return 12;
}

} // namespace dalton_seed_detail

struct GsSeedOptions {
  double      L        = 200.0;   ///< box half-width; MUST match the run's `l`
  double      thresh   = 1e-4;    ///< projection thresh (k from the protocol table)
  double      energy   = 0.0;     ///< provenance only
  std::string xc       = "hf";    ///< provenance only (load_mos discards)
  std::string localize = "canon"; ///< provenance only
  int         nio      = 1;
  /// Extra archive prefixes to write the SAME seed to (e.g. "<prefix>.gs_seed"):
  /// moldft's save_mos overwrites <prefix>.restartdata on its first save, so a
  /// preserved copy is the only record of what the SCF started from.
  std::vector<std::string> extra_prefixes;
  /// Molecule to stamp into the RestartMetadata. RestartPlan matches the
  /// archive geometry against the deck at 1e-8 bohr and its eprec exactly;
  /// the molden coordinates differ from the deck at ~1e-8 (DALTON prints
  /// 10 digits) and carry no eprec, so a seed stamped with them was rejected
  /// as "geometry moved", and moldft fell back to the initial guess. The caller has already fingerprinted molden vs the
  /// active molecule at 1e-4; stamp the active molecule.
  const Molecule *active_molecule = nullptr;
  /// Called after each archive is written, with its name, an emitter that
  /// serializes the same content into another parallel archive, and the header
  /// and orbital count it holds. molresponse uses it to write an HDF5 twin of
  /// the seed; chem itself writes none.
  using Emitter  = std::function<void(archive::ParallelOutputArchive<archive::VectorOutputArchive> &)>;
  using AfterSave = std::function<void(World &, const std::string &archive_name, const Emitter &,
                                       const RestartMetadata &meta, unsigned int n_mo)>;
  AfterSave after_save;
};

struct GsSeedReport {
  int    n_occ = 0, n_ao = 0, n_mo = 0;
  double max_offdiag_pre = 0.0;   ///< max |S_ij|, i!=j, before Loewdin
  double max_norm_dev    = 0.0;   ///< max | ||phi_i|| - 1 | after
  std::string archive;            ///< <out_prefix>.restartdata
};

/// Write <out_prefix>.restartdata seeded from the molden's lowest n_occ MOs.
/// Collective. Sets FunctionDefaults<3> cell/k/thresh for the projection
/// (callers re-set their own protocol afterwards if it differs).
inline GsSeedReport
write_gs_seed_from_molden(World &world, const std::string &molden_path,
                          int n_occ, const std::string &out_prefix,
                          const GsSeedOptions &opt = {}) {
  using dalton_seed_detail::DaltonMOFunctor;
  using dalton_seed_detail::symbol_to_Z;
  GsSeedReport rep;

  DaltonMoldenResult molden = read_molden(molden_path);
  const int n_ao = molden.n_ao;
  if (n_occ <= 0 || n_occ > molden.n_mo)
    throw std::runtime_error("dalton gs seed: n_occ=" + std::to_string(n_occ) +
                             " out of range (molden n_mo=" + std::to_string(molden.n_mo) + ")");
  rep.n_occ = n_occ; rep.n_ao = n_ao; rep.n_mo = molden.n_mo;

  Molecule molecule;
  for (size_t a = 0; a < molden.atom_symbols.size(); ++a) {
    const int Z = symbol_to_Z(molden.atom_symbols[a]);
    molecule.add_atom(molden.coords[a][0], molden.coords[a][1], molden.coords[a][2],
                      static_cast<double>(Z), Z);
  }

  const int k = dalton_seed_detail::k_for_thresh(opt.thresh);
  Tensor<double> cell(3L, 2L);
  for (int i = 0; i < 3; i++) { cell(i, 0) = -opt.L; cell(i, 1) = opt.L; }
  FunctionDefaults<3>::set_cell(cell);
  FunctionDefaults<3>::set_k(k);
  FunctionDefaults<3>::set_thresh(opt.thresh);

  std::vector<real_function_3d> amo;
  for (int i = 0; i < n_occ; i++) {
    std::vector<double> col(molden.mo_coeffs.begin() + static_cast<ptrdiff_t>(i) * n_ao,
                            molden.mo_coeffs.begin() + static_cast<ptrdiff_t>(i + 1) * n_ao);
    std::shared_ptr<FunctionFunctorInterface<double, 3>> ff =
        std::make_shared<DaltonMOFunctor>(molden.basis, std::move(col));
    amo.push_back(FunctionFactory<double, 3>(world).functor(ff).thresh(opt.thresh)
                      .truncate_on_project());
  }
  {
    Tensor<double> S = matrix_inner(world, amo, amo, true);
    for (int i = 0; i < n_occ; i++)
      for (int j = 0; j < n_occ; j++)
        if (i != j) rep.max_offdiag_pre = std::max(rep.max_offdiag_pre, std::abs(S(i, j)));
  }
  amo = orthonormalize_symmetric(amo);

  Tensor<double> aeps(static_cast<long>(n_occ)), aocc(static_cast<long>(n_occ));
  for (int i = 0; i < n_occ; i++) {
    aeps(i) = (i < static_cast<int>(molden.mo_energies.size())) ? molden.mo_energies[i] : 0.0;
    aocc(i) = 1.0;
    rep.max_norm_dev = std::max(rep.max_norm_dev, std::abs(amo[i].norm2() - 1.0));
  }
  std::vector<int> aset(static_cast<size_t>(n_occ), 0);

  RestartMetadata meta;
  meta.current_energy       = opt.energy;
  meta.spin_restricted      = true;
  meta.L                    = opt.L;
  meta.k                    = k;
  const Molecule &stamp     = opt.active_molecule ? *opt.active_molecule : molecule;
  meta.molecule             = stamp;
  meta.xc                   = opt.xc;
  meta.localize             = opt.localize;
  meta.converged_for_thresh = opt.thresh;
  meta.representation       = Representation::mo;
  meta.eprec                = stamp.parameters.eprec();
  meta.madness_version      = MADNESS_PACKAGE_VERSION;
  // one id for the seed and its preserved copy: they are the same orbitals
  meta.archive_id           = new_archive_id(world);
  meta.origin               = "dalton-seed";
  meta.nalpha               = n_occ;
  meta.nbeta                = n_occ;
  meta.nmo_beta             = n_occ;
  auto emit = [&](auto &ar) {
    meta.write(ar);
    ar & static_cast<unsigned int>(amo.size());
    ar & aeps & aocc & aset;
    for (unsigned int i = 0; i < amo.size(); ++i) ar & amo[i];
  };
  auto write_to = [&](const std::string &prefix) {
    const std::string name = prefix + ".restartdata";
    {
      archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(
          world, name.c_str(), opt.nio);
      emit(ar);
    }
    world.gop.fence();
if (opt.after_save)
      opt.after_save(world, name, GsSeedOptions::Emitter(
          [&](archive::ParallelOutputArchive<archive::VectorOutputArchive> &ar) { emit(ar); }),
          meta, static_cast<unsigned int>(amo.size()));
    return name;
  };
  rep.archive = write_to(out_prefix);
  for (const auto &px : opt.extra_prefixes) write_to(px);

  if (world.rank() == 0)
    print("[DALTON-SEED] GS: wrote", n_occ, "occupied MOs ->", rep.archive,
          " (molden", molden_path, " n_ao =", n_ao, " n_mo =", molden.n_mo,
          " L =", opt.L, " k =", k, " thresh =", opt.thresh,
          " max|S_ij| pre-Loewdin =", rep.max_offdiag_pre, ")");
  return rep;
}

// ---------------------------------------------------------------------------
// Ground-state seed hook. madqc installs it on the SCF application
// (SCFApplication::set_pre_run_hook) when the deck carries `dalton.dir`. Runs collectively inside the SCF work dir right before the
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
inline void seed_gs_from_dalton_dir(World &world, const Params &params,
                                    const std::filesystem::path &workdir,
                                    const std::string &dalton_dir,
                                    GsSeedOptions::AfterSave after_save = {}) {
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
  opt.after_save = std::move(after_save);
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

} // namespace madness

#endif // MADNESS_CHEM_DALTON_SEED_H
