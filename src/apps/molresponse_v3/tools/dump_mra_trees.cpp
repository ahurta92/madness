// Offline MRA tree-skeleton dumper  (molresponse_v3 — R2/L4 export, structural target)
//
// Loads already-converged ground-state orbitals from a moldft restart archive
// and dumps each orbital's adaptive multiresolution octree as JSON, in TWO
// representations:
//
//   mo_<i>_recon.json  reconstructed: leaves hold scaling coefficients, so the
//                      per-node coeff norm = the function's L2 content in that
//                      box  ->  magnitude / spatial-support map.
//   mo_<i>_comp.json   compressed: internal nodes hold difference (wavelet)
//                      coefficients, so the per-node coeff norm = the local
//                      detail/error added at that scale  ->  error-at-resolution
//                      map (the essence of MRA).
//
// Output is consumed by the gecko `mra_tree` Python module, which turns the JSON
// into a flat box table + HDF5/HyperTreeGrid for ParaView and the analysis/viz
// layer (depth histogram, wavelet energy-per-scale spectrum, 2D slices, and the
// N x N product-overlap matrix used to study Coulomb/exchange partitioning).
//
// This tool LOADS ONLY; it never triggers an electronic solve (Tier-B sibling,
// doc 18). It reuses MADNESS's own tree dumper (print_tree_jsonfile), which
// already emits per-node {has_coeff, has_children, norm, norm_tree, snorm,
// dnorm, rank, dim} + owner, so no new tree-walking code is needed.
//
// All MADNESS Function ops below (reconstruct/compress, print_tree_jsonfile,
// norm2) are COLLECTIVE -- they run on every rank; only printing is gated on
// rank 0.
//
// Usage:
//   dump_mra_trees --archive=<prefix>.restartdata [--out=DIR]
//                  [--maxlevel=N] [--k=N] [--thresh=X] [--max-orbitals=N]

#include "../GroundState.hpp"
#include "../ResponseProtocol.hpp"

#include <madness/chem/atomutil.h>
#include <madness/external/nlohmann_json/json.hpp>
#include <madness/misc/info.h>
#include <madness/mra/funcplot.h>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace madness;
using namespace molresponse_v3;

namespace {

// Mirror the molecule sourcing the FD/ES skeletons use: read the geometry from
// the moldft calc_info JSON sitting next to the archive. from_archive needs a
// Molecule to derive electron counts; an empty one starts the closed-shell load
// path but matching the canonical path avoids surprises.
Molecule load_molecule_near(const std::string& archive_path) {
  Molecule molecule;
  auto archive_dir = std::filesystem::path(archive_path).parent_path();
  for (const auto& name : {"moldft.calc_info.json", "mad.calc_info.json"}) {
    auto candidate = archive_dir / name;
    if (std::filesystem::exists(candidate)) {
      std::ifstream ifs(candidate);
      nlohmann::json j;
      ifs >> j;
      nlohmann::json mol_json;
      if (j.contains("tasks") && j["tasks"].is_array() && !j["tasks"].empty())
        mol_json = j["tasks"][0]["molecule"];
      else if (j.contains("molecule"))
        mol_json = j["molecule"];
      if (!mol_json.is_null()) molecule.from_json(mol_json);
      break;
    }
  }
  return molecule;
}

// Write one Function's octree as JSON, MPI-safe. `Function::print_tree_json` is
// COLLECTIVE -- every rank must call it (rank 0 walks the tree and pulls remote
// nodes behind internal fences) -- but it emits content only on rank 0. So we
// buffer into a stringstream on every rank and let ONLY rank 0 open and write
// the file. (MADNESS's own print_tree_jsonfile opens the ofstream on every rank
// and writes the wrapper on every rank, so under -np>1 all ranks race on the
// same path and produce corrupt JSON.)
void write_tree_json(World& world, const Function<double, 3>& f,
                     const std::string& path) {
  std::ostringstream oss;
  const Tensor<double> cell = FunctionDefaults<3>::get_cell();
  oss << "{\"cell\":[";
  for (int d = 0; d < 3; ++d) {
    oss << "[" << cell(d, 0) << "," << cell(d, 1) << "]";
    if (d != 2) oss << ",";
  }
  oss << "],\"tree\":{";
  f.print_tree_json(oss);  // collective: call on every rank; fills oss on rank 0
  oss << "}}";
  if (world.rank() == 0) {
    std::ofstream os(path);
    os << oss.str();
  }
}

// Write the nuclear geometry as VTK POLYDATA points in BOHR -- the same frame as
// the box dumps -- with per-atom atomic number Z and covalent radius, so ParaView
// can glyph the atoms as element-sized spheres overlaid on the octree. The
// Molecule is a small replicated object (not a distributed Function), so this is
// pure local data: write on rank 0 only, no collective ops, no fence.
void write_geometry_vtk(World& world, const Molecule& mol,
                        const std::string& path) {
  if (world.rank() != 0) return;
  const std::size_t n = mol.natom();
  std::ostringstream pts, zs, rs;
  for (std::size_t i = 0; i < n; ++i) {
    const Atom& a = mol.get_atom(i);
    pts << a.x << " " << a.y << " " << a.z << "\n";
    zs << a.atomic_number << "\n";
    rs << get_atomic_data(a.atomic_number).covalent_radius << "\n";
  }
  std::ofstream os(path);
  os << "# vtk DataFile Version 3.0\nMRA molecule geometry (bohr)\nASCII\n"
     << "DATASET POLYDATA\nPOINTS " << n << " float\n"
     << pts.str() << "VERTICES " << n << " " << 2 * n << "\n";
  for (std::size_t i = 0; i < n; ++i) os << "1 " << i << "\n";
  os << "POINT_DATA " << n << "\n"
     << "SCALARS Z float 1\nLOOKUP_TABLE default\n" << zs.str()
     << "SCALARS radius float 1\nLOOKUP_TABLE default\n" << rs.str();
}

// Sample one orbital on a regular grid and write a Gaussian .cube (geometry
// baked in) so ParaView/VMD can render it as signed isosurfaces -- the actual
// orbital lobes/nodes, which the box dumps (one norm per box) cannot show.
//
// The grid is a cubic box snug around the molecule (bbox +/- `pad` bohr), NOT the
// full ~200-bohr cell, so the points land where the orbital has support. Sampling
// uses Function::eval_cube -- the same batched, distributed grid evaluator behind
// plotvtk/plotdx (each box fills its grid cells), which is COLLECTIVE (all ranks
// call it) and orders of magnitude faster than per-point evaluation. Only rank 0
// then writes the file.
void write_cube_file(World& world, Function<double, 3> f, const Molecule& mol,
                     const std::string& path, int npoints, double pad) {
  double lo[3] = {1e9, 1e9, 1e9}, hi[3] = {-1e9, -1e9, -1e9};
  for (std::size_t a = 0; a < mol.natom(); ++a) {
    const Atom& at = mol.get_atom(a);
    const double c[3] = {at.x, at.y, at.z};
    for (int d = 0; d < 3; ++d) {
      lo[d] = std::min(lo[d], c[d]);
      hi[d] = std::max(hi[d], c[d]);
    }
  }
  Tensor<double> cell(3, 2);
  for (int d = 0; d < 3; ++d) {
    lo[d] -= pad;
    hi[d] += pad;
    cell(d, 0) = lo[d];
    cell(d, 1) = hi[d];
  }

  // Collective, batched grid evaluation -> full grid (rank 0 holds it).
  const std::vector<long> npt(3, npoints);
  const Tensor<double> grid = f.eval_cube(cell, npt);

  if (world.rank() != 0) return;
  double delta[3];
  for (int d = 0; d < 3; ++d) delta[d] = (hi[d] - lo[d]) / (npoints - 1);
  FILE* file = std::fopen(path.c_str(), "w");
  if (!file) return;
  std::fprintf(file, "MADNESS MRA orbital (bohr)\n%s\n", path.c_str());
  std::fprintf(file, "%d %12.6f %12.6f %12.6f\n", int(mol.natom()),
               lo[0], lo[1], lo[2]);
  std::fprintf(file, "%d %12.6f %12.6f %12.6f\n", npoints, delta[0], 0.0, 0.0);
  std::fprintf(file, "%d %12.6f %12.6f %12.6f\n", npoints, 0.0, delta[1], 0.0);
  std::fprintf(file, "%d %12.6f %12.6f %12.6f\n", npoints, 0.0, 0.0, delta[2]);
  for (std::size_t a = 0; a < mol.natom(); ++a) {
    const Atom& at = mol.get_atom(a);
    std::fprintf(file, "%d %12.6f %12.6f %12.6f %12.6f\n", at.atomic_number,
                 double(at.atomic_number), at.x, at.y, at.z);
  }
  long per_line = 0;
  for (int i = 0; i < npoints; ++i)
    for (int j = 0; j < npoints; ++j)
      for (int k = 0; k < npoints; ++k) {
        std::fprintf(file, "%13.5e ", grid(i, j, k));
        if (++per_line == 6) {
          std::fprintf(file, "\n");
          per_line = 0;
        }
      }
  std::fprintf(file, "\n");
  std::fclose(file);
}

// Dump one orbital's octree in both representations. `stem` is the output path
// prefix (no extension). Collective: must be called on every rank.
void dump_function_trees(World& world, Function<double, 3> f,
                         const std::string& stem) {
  f.reconstruct();
  write_tree_json(world, f, stem + "_recon.json");
  f.compress();
  write_tree_json(world, f, stem + "_comp.json");
}

}  // namespace

int main(int argc, char** argv) {
  World& world = initialize(argc, argv);
  int rc = 0;
  try {
    startup(world, argc, argv, true);

    commandlineparser parser(argc, argv);
    if (!parser.key_exists("archive")) {
      if (world.rank() == 0) {
        print("Usage: dump_mra_trees --archive=<prefix>.restartdata "
              "[--out=DIR] [--maxlevel=N] [--k=N] [--thresh=X] "
              "[--max-orbitals=N] [--cube] [--cube-npoints=N] [--cube-pad=X]");
        print("  Loads ground-state orbitals at their native (k, thresh) and");
        print("  dumps per-orbital MRA octree JSON (reconstructed + compressed),");
        print("  geometry.vtk, and (with --cube) per-orbital .cube isosurfaces.");
      }
      finalize();
      return 1;
    }

    // Inner scope so MADNESS objects destruct before finalize().
    {
      const std::string archive_path = parser.value_raw("archive");
      const std::string out_dir =
          parser.key_exists("out") ? parser.value_raw("out")
                                   : std::string("mra_trees");
      const int max_orbitals = parser.key_exists("max-orbitals")
                                   ? std::stoi(parser.value("max-orbitals"))
                                   : -1;
      // Optional orbital isosurface export (.cube). Off by default -- it is the
      // expensive step (per-point grid evaluation).
      const bool cube = parser.key_exists("cube");
      const int cube_npoints = parser.key_exists("cube-npoints")
                                   ? std::stoi(parser.value("cube-npoints")) : 80;
      const double cube_pad = parser.key_exists("cube-pad")
                                  ? std::stod(parser.value("cube-pad")) : 6.0;

      auto header = GroundState::read_archive_header(world, archive_path);

      // Load at NATIVE resolution unless explicitly overridden. Matching the
      // archive's (k, thresh) means SCF::load_mos performs no reprojection, so
      // the tree we walk is exactly the one the solver produced. Reprojecting
      // would rebuild the tree and we'd be visualizing a different object.
      const int k =
          parser.key_exists("k") ? std::stoi(parser.value("k")) : header.k;
      const double thresh = parser.key_exists("thresh")
                                ? std::stod(parser.value("thresh"))
                                : header.converged_for_thresh;
      set_response_protocol(world, header.L, thresh, k);

      auto molecule = load_molecule_near(archive_path);
      auto gs = GroundState::from_archive(world, archive_path, molecule);
      if (world.rank() == 0) gs.print_info();

      if (world.rank() == 0) std::filesystem::create_directories(out_dir);
      world.gop.fence();

      write_geometry_vtk(world, gs.molecule(), out_dir + "/geometry.vtk");
      if (world.rank() == 0) {
        print("  wrote geometry.vtk  (", gs.molecule().natom(), "atoms )");
      }

      const vecfuncT& mos = gs.orbitals_alpha();
      const std::size_t n =
          (max_orbitals >= 0)
              ? std::min<std::size_t>(mos.size(),
                                      static_cast<std::size_t>(max_orbitals))
              : mos.size();

      if (world.rank() == 0) {
        print("DUMP_MRA  archive=", archive_path, "  out=", out_dir,
              "  k=", k, "  thresh=", thresh, "  n_orbitals=", n);
      }

      for (std::size_t i = 0; i < n; ++i) {
        const std::string stem = out_dir + "/mo_" + std::to_string(i);
        dump_function_trees(world, mos[i], stem);
        if (cube) {
          write_cube_file(world, mos[i], gs.molecule(), stem + ".cube",
                          cube_npoints, cube_pad);
        }
        const double nrm = mos[i].norm2();  // collective -- all ranks
        if (world.rank() == 0) print("  dumped mo_", i, "  norm2=", nrm);
      }
      world.gop.fence();
    }

    finalize();
    return rc;
  } catch (const std::exception& e) {
    if (world.rank() == 0) print("Error:", e.what());
    finalize();
    return 1;
  }
}
