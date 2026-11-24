// Offline MRA tree-skeleton dumper  (molresponse_v3 — R2/L4 export, structural target)
//
// Loads the saved ground-state orbitals of a madqc run folder and dumps each orbital's adaptive multiresolution octree as JSON, in TWO
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
// The run folder is --run-dir (default: the current directory). The ground
// state comes from its node, <run>/states/ground_state/, at the level its
// description records, with the molecule from the run's analysis files
// (saved_states.hpp:read_ground). With --fd and --es it also dumps the saved
// response states: it lists the node folders under <run>/states/response/,
// reads each description, and loads the driven states (--fd) and excited-state
// sets (--es) through the node module. It never writes into the run folder.
//
// Usage:
//   dump_mra_trees [--run-dir=<madqc run folder>] [--out=DIR]
//                  [--maxlevel=N] [--k=N] [--thresh=X] [--max-orbitals=N]
//                  [--fd] [--es]

#include "../GroundState.hpp"
#include "../ResponseProtocol.hpp"
#include "../kernels/tags.hpp"     // Static / Full / TDA, ClosedShell
#include "../kernels/static.hpp"   // Kernels<Static,ClosedShell>::compute_density
#include "../kernels/full.hpp"     // Kernels<Full,ClosedShell>::compute_density
#include "../kernels/tda.hpp"      // Kernels<TDA,ClosedShell>::compute_density
#include "../solvers/response_metadata.hpp"  // ResponseMetadata::freq_key (output names)
#include "saved_states.hpp"        // node folders under <run>/states/response/

#include <madness/chem/atomutil.h>
#include <nlohmann/json.hpp>
#include <madness/misc/info.h>
#include <madness/mra/funcplot.h>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

// Native HyperTreeGrid (.htg) export is OPTIONAL, behind MADNESS_ENABLE_VTK
// (CMake -> defines MADNESS_HAS_VTK). A no-VTK build compiles unchanged and
// --htg becomes a no-op with a one-line notice; the JSON/.vtk path is untouched.
#ifdef MADNESS_HAS_VTK
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkHyperTreeGrid.h>
#include <vtkHyperTreeGridNonOrientedGeometryCursor.h>
#include <vtkNew.h>
#include <vtkXMLHyperTreeGridWriter.h>
// Overlapping-AMR (function values on the adaptive octree -> .vthb)
#include <vtkAMRBox.h>
#include <vtkFloatArray.h>
#include <vtkOverlappingAMR.h>
#include <vtkUniformGrid.h>
#include <vtkXMLUniformGridAMRWriter.h>
#include <array>
#endif

// Native MRA coefficient export (.mad.h5) is OPTIONAL, behind MADNESS_ENABLE_HDF5
// (CMake -> defines MADNESS_HAS_HDF5). It reuses the io-hdf5 structured Function
// writer (solvers/function_hdf5_io.hpp); a no-HDF5 build compiles unchanged and
// --coeffs becomes a one-line notice (the header is itself a no-op without the macro).
#ifdef MADNESS_HAS_HDF5
#include "../solvers/function_hdf5_io.hpp"
#endif

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

using namespace madness;
using namespace molresponse_v3;

namespace {

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

// ---- native HyperTreeGrid (.htg) export -----------------------------------
//
// Build a vtkHyperTreeGrid from the LIVE octree (a single root tree over the
// whole simulation cell, branch factor 2) instead of the legacy one-hexahedron-
// per-leaf VTK. Geometry is implicit in the tree path, so the file is ~20-30x
// smaller and ParaView gets a native Maximum-Level depth control + per-level
// opacity. Mirrors the validated reference write_htg() in
// madness_studies/scripts/view_boxes.py (child order, rms formula, the cursor
// SetGlobalIndexStart gotcha, the bohr frame, the array names) -- but builds
// from the live tree rather than the dumped _boxes.vtk.
#ifdef MADNESS_HAS_VTK
struct HtgNode {
  double norm = 0.0;
  double dnorm = -1.0;
  int owner = 0;
  bool has_children = false;
};
// Keyed by (level, lx, ly, lz) so std::map ordering needs no Vector operator<.
using HtgKey = std::tuple<int, Translation, Translation, Translation>;
using HtgMap = std::map<HtgKey, HtgNode>;

// rank-0-only recursion mirroring FunctionImpl::do_print_tree_json (mraimpl.h):
// coeffs.find(key).get() pulls remote nodes while the OTHER ranks service the
// request inside the trailing gop.fence() in collect_live_tree(). The norm
// matches FunctionNode::print_json (funcimpl.h); owner = coeffs.owner(key), the
// live-container pmap owner -- the same source do_print_tree_json writes and the
// legacy _boxes.vtk uses (the recon-JSON owner can be unreliable; this is not).
template <typename ImplT>
void collect_live_tree_r(const ImplT& impl, const Key<3>& key, HtgMap& out) {
  const auto& coeffs = impl.get_coeffs();
  auto it = coeffs.find(key).get();
  if (it == coeffs.end()) return;  // absent child of a has_children node: skip
  const auto& node = it->second;
  HtgNode rec;
  rec.norm = node.has_coeff() ? node.coeff().normf() : 0.0;
  if (rec.norm < 1e-12) rec.norm = 0.0;
  rec.dnorm = node.get_dnorm();
  rec.owner = coeffs.owner(key);
  rec.has_children = node.has_children();
  const auto& l = key.translation();
  out[{int(key.level()), l[0], l[1], l[2]}] = rec;
  if (node.has_children())
    for (KeyChildIterator<3> kit(key); kit; ++kit)
      collect_live_tree_r(impl, kit.key(), out);
}

// Collective: rank 0 walks the tree; the fence lets its remote finds complete.
void collect_live_tree(World& world, const Function<double, 3>& f, HtgMap& out) {
  const auto impl = f.get_impl();
  if (world.rank() == 0)
    collect_live_tree_r(*impl, impl->get_cdata().key0, out);
  world.gop.fence();
}

// Build the single-root HTG from the node map and write it (rank 0 only). DFS
// tracks (level, translation) explicitly; child order c -> (c&1,(c>>1)&1,(c>>2)&1)
// with l_child = 2*l + (i,j,k); rms = norm / sqrt(box_volume),
// box_volume = cell_volume / 2^(3*level). Cell arrays match the legacy .vtk
// names so existing ParaView coloring (level/norm/rms/owner) works unchanged.
void write_htg(World& world, const HtgMap& nodes, const std::string& path) {
  if (world.rank() != 0 || nodes.empty()) return;
  const Tensor<double> cell = FunctionDefaults<3>::get_cell();
  double cell_volume = 1.0;
  for (int d = 0; d < 3; ++d) cell_volume *= (cell(d, 1) - cell(d, 0));

  vtkNew<vtkHyperTreeGrid> htg;
  htg->Initialize();
  htg->SetDimensions(2, 2, 2);  // 1x1x1 cells -> a single root tree
  htg->SetBranchFactor(2);
  for (int d = 0; d < 3; ++d) {
    vtkNew<vtkDoubleArray> coord;
    coord->SetNumberOfValues(2);
    coord->SetValue(0, cell(d, 0));
    coord->SetValue(1, cell(d, 1));
    if (d == 0) htg->SetXCoordinates(coord.Get());
    else if (d == 1) htg->SetYCoordinates(coord.Get());
    else htg->SetZCoordinates(coord.Get());
  }

  vtkNew<vtkDoubleArray> a_level, a_norm, a_rms, a_dnorm, a_owner;
  a_level->SetName("level");
  a_norm->SetName("norm");
  a_rms->SetName("rms");
  a_dnorm->SetName("dnorm");
  a_owner->SetName("owner");

  vtkNew<vtkHyperTreeGridNonOrientedGeometryCursor> cur;
  htg->InitializeNonOrientedGeometryCursor(cur.Get(), 0, /*create=*/true);
  cur->SetGlobalIndexStart(0);  // without this GetGlobalNodeIndex() == -1

  std::function<void(int, Translation, Translation, Translation)> rec =
      [&](int level, Translation lx, Translation ly, Translation lz) {
        const vtkIdType gid = cur->GetGlobalNodeIndex();
        auto it = nodes.find({level, lx, ly, lz});
        const HtgNode hn = (it != nodes.end()) ? it->second : HtgNode{};
        const double box_volume = cell_volume / std::pow(2.0, 3.0 * level);
        a_level->InsertValue(gid, double(level));
        a_norm->InsertValue(gid, hn.norm);
        a_rms->InsertValue(
            gid, box_volume > 0.0 ? hn.norm / std::sqrt(box_volume) : 0.0);
        a_dnorm->InsertValue(gid, hn.dnorm);
        a_owner->InsertValue(gid, double(hn.owner));
        if (it != nodes.end() && it->second.has_children) {
          cur->SubdivideLeaf();
          for (int c = 0; c < 8; ++c) {
            cur->ToChild(c);
            rec(level + 1, 2 * lx + (c & 1), 2 * ly + ((c >> 1) & 1),
                2 * lz + ((c >> 2) & 1));
            cur->ToParent();
          }
        }
      };
  rec(0, 0, 0, 0);  // root = (level 0, {0,0,0}) = FunctionImpl cdata.key0

  htg->GetCellData()->AddArray(a_level.Get());
  htg->GetCellData()->AddArray(a_norm.Get());
  htg->GetCellData()->AddArray(a_rms.Get());
  htg->GetCellData()->AddArray(a_dnorm.Get());
  htg->GetCellData()->AddArray(a_owner.Get());

  vtkNew<vtkXMLHyperTreeGridWriter> writer;
  writer->SetFileName(path.c_str());
  // Inline ASCII, NOT the default Appended mode: in VTK 9.1 the appended-data
  // offset-fixup pass mangles the header (writes a RangeMax/offset fragment over
  // byte 0, clobbering <VTKFile>), so ParaView/VTK 9.5 rejects it as not
  // well-formed. ASCII writes data inline (no seek-back). Switch to
  // SetDataModeToBinary() later for the size win once validated -- it is also
  // inline (base64) and avoids the buggy appended path.
  writer->SetDataModeToAscii();
  writer->SetInputData(htg.Get());
  writer->Write();
}

// Collect (collective) + write (rank 0). Safe to call on every rank.
void write_htg_from_function(World& world, const Function<double, 3>& f,
                             const std::string& path) {
  HtgMap nodes;
  collect_live_tree(world, f, nodes);
  write_htg(world, nodes, path);
}

// ---- adaptive function-value export: vtkOverlappingAMR (.vthb) -------------
//
// PROTOTYPE (see AMR_EXPORT_NOTES.md). HTG stores ONE scalar per box (skeleton);
// OverlappingAMR stores the FUNCTION sampled on a small m^3 uniform grid PER LEAF
// BOX across refinement levels, so ParaView's AMR Contour gives smooth isosurfaces
// at sub-box resolution while keeping adaptivity (fine blocks only where the tree
// refined) -- the adaptive replacement for the dense .cube. Octree => refinement
// ratio 2; origin/spacing in bohr (aligns with geometry.vtk / the .htg).
//
// The MADNESS-side mapping below (leaf walk, box geometry, level grouping,
// index-space AMR box, per-level spacing) is the load-bearing logic. The VTK
// OverlappingAMR API is version-sensitive -- the assembly calls are a best-effort
// first cut to verify at the first -DMADNESS_ENABLE_VTK build (NOTES §VTK-API).

// Leaf (level,lx,ly,lz) list replicated on ALL ranks: rank 0 walks the tree
// (collect_live_tree_r), keeps the leaves (no children), and broadcasts the flat
// list so the eval_cube sampling loop stays collective + in lockstep everywhere.
std::vector<std::array<long, 4>>
collect_leaf_keys(World& world, const Function<double, 3>& f) {
  std::vector<long> flat;
  if (world.rank() == 0) {
    HtgMap nodes;
    const auto impl = f.get_impl();
    collect_live_tree_r(*impl, impl->get_cdata().key0, nodes);
    for (const auto& [key, nd] : nodes)
      if (!nd.has_children) {
        flat.push_back(std::get<0>(key));
        flat.push_back(std::get<1>(key));
        flat.push_back(std::get<2>(key));
        flat.push_back(std::get<3>(key));
      }
  }
  world.gop.fence();  // pair rank 0's remote coeffs.find()s
  std::size_t nflat = flat.size();
  world.gop.broadcast(&nflat, 1, 0);
  flat.resize(nflat);
  if (nflat) world.gop.broadcast(flat.data(), nflat, 0);
  std::vector<std::array<long, 4>> leaves(nflat / 4);
  for (std::size_t i = 0; i < leaves.size(); ++i)
    leaves[i] = {flat[4 * i], flat[4 * i + 1], flat[4 * i + 2], flat[4 * i + 3]};
  return leaves;
}

// Sample f on every leaf box (m^3 uniform grid via the collective eval_cube; rank
// 0 receives each grid) and assemble a vtkOverlappingAMR; write <path> on rank 0.
void write_amr(World& world, Function<double, 3> f, const std::string& path,
               int m) {
  f.reconstruct();
  const auto leaves = collect_leaf_keys(world, f);
  const Tensor<double> cell = FunctionDefaults<3>::get_cell();
  double lo[3], width[3];
  for (int d = 0; d < 3; ++d) {
    lo[d] = cell(d, 0);
    width[d] = cell(d, 1) - cell(d, 0);
  }
  long Lmax = 0, Lmin = leaves.empty() ? 0 : leaves.front()[0];
  for (const auto& k : leaves) {
    Lmax = std::max(Lmax, k[0]);
    Lmin = std::min(Lmin, k[0]);
  }

  // collective sampling: ALL ranks call eval_cube per leaf; rank 0 keeps the grids
  struct Block { long n, lx, ly, lz; Tensor<double> g; };
  std::vector<Block> blocks;
  const std::vector<long> npt(3, m);
  for (const auto& k : leaves) {
    Tensor<double> bc(3, 2);
    for (int d = 0; d < 3; ++d) {
      const double h = width[d] / std::pow(2.0, double(k[0]));  // leaf box width
      const double dx = h / m;                                  // m-grid sub-cell width
      // Sample at the m CELL CENTERS, not box corners: the outer corner of a
      // boundary leaf sits exactly on the sim-cell edge, where eval_cube ->
      // plot_cube_kernel maps it to x=1.0(+fp eps) and trips its x in [0,1]
      // assertion. Centers are strictly interior (x in (0,1)) and also align the
      // sampled values with the m-cell block's cell-data positions (fixes the
      // half-cell offset noted in AMR_EXPORT_NOTES.md simplification #2).
      bc(d, 0) = lo[d] + double(k[d + 1]) * h + 0.5 * dx;       // first cell center
      bc(d, 1) = lo[d] + double(k[d + 1] + 1) * h - 0.5 * dx;   // last cell center
    }
    Tensor<double> g = f.eval_cube(bc, npt);  // collective; rank 0 holds the grid
    if (world.rank() == 0) blocks.push_back({k[0], k[1], k[2], k[3], copy(g)});
  }
  if (world.rank() != 0) return;

  // --- assemble on rank 0 (VTK OverlappingAMR) ---
  // RE-BASE the hierarchy to the populated level range [Lmin, Lmax]: MADNESS
  // leaves never sit at the top tree levels (the root box is always subdivided),
  // so AMR level 0 must map to Lmin, not the MADNESS root. An empty AMR level 0
  // segfaults vtkXMLUniformGridAMRReader::UpdateInformation (it cannot build the
  // metadata). AMR level index = n - Lmin; spacing/origin stay physical (bohr).
  const int nlev = int(Lmax - Lmin + 1);
  std::vector<int> per_level(nlev, 0);
  for (const auto& b : blocks) per_level[b.n - Lmin]++;
  // Review MED: the rebase guarantees a populated level 0 (Lmin), but an
  // INTERIOR level can still be empty if the leaf depths have a gap — that can
  // trip the same vtkXMLUniformGridAMRReader metadata path as an empty level 0.
  // Warn here; the full fix (compact to only populated levels with the correct
  // per-level refinement ratio) needs ParaView validation and belongs to the
  // viz/pymra thread.
  for (int lev = 1; lev + 1 < nlev; ++lev)
    if (per_level[lev] == 0)
      print("[--amr] WARNING: interior AMR level", lev, "(MADNESS n =",
            Lmin + lev, ") has NO leaf boxes; the .vthb may fail to load in "
            "some readers. Report to the viz thread if ParaView rejects it.");

  vtkNew<vtkOverlappingAMR> amr;
  amr->Initialize(nlev, per_level.data());
  double origin[3] = {lo[0], lo[1], lo[2]};
  amr->SetOrigin(origin);
  amr->SetGridDescription(VTK_XYZ_GRID);
  for (int lev = 0; lev < nlev; ++lev) {
    const long n = Lmin + lev;
    double sp[3];
    for (int d = 0; d < 3; ++d) sp[d] = (width[d] / std::pow(2.0, double(n))) / m;
    amr->SetSpacing(lev, sp);
    if (lev > 0) amr->SetRefinementRatio(lev - 1, 2);  // octree => 2x per level
  }
  std::vector<int> next(nlev, 0);
  for (const auto& b : blocks) {
    const long n = b.n;
    const int lev = int(n - Lmin), j = next[lev]++;
    const long l[3] = {b.lx, b.ly, b.lz};
    int boxlo[3], boxhi[3];
    for (int d = 0; d < 3; ++d) {
      boxlo[d] = int(l[d] * m);
      boxhi[d] = boxlo[d] + m - 1;
    }
    vtkAMRBox box(boxlo, boxhi);
    amr->SetAMRBox(lev, j, box);

    vtkNew<vtkUniformGrid> grid;
    double gorigin[3], sp[3];
    for (int d = 0; d < 3; ++d) {
      const double h = width[d] / std::pow(2.0, double(n));
      gorigin[d] = lo[d] + double(l[d]) * h;
      sp[d] = h / m;
    }
    grid->SetOrigin(gorigin);
    grid->SetSpacing(sp);
    grid->SetDimensions(m + 1, m + 1, m + 1);  // m^3 cells
    vtkNew<vtkFloatArray> vals;
    vals->SetName("psi");
    vals->SetNumberOfValues(vtkIdType(m) * m * m);
    for (int iz = 0; iz < m; ++iz)       // eval_cube g(ix,iy,iz); VTK cell x-fastest
      for (int iy = 0; iy < m; ++iy)
        for (int ix = 0; ix < m; ++ix)
          vals->SetValue((vtkIdType(iz) * m + iy) * m + ix, float(b.g(ix, iy, iz)));
    grid->GetCellData()->AddArray(vals.Get());
    amr->SetDataSet(lev, j, grid.Get());
  }
  vtkNew<vtkXMLUniformGridAMRWriter> w;
  w->SetFileName(path.c_str());
  w->SetDataModeToAscii();  // see write_htg: default Appended mode mangles VTK 9.1 output
  w->SetInputData(amr.Get());
  w->Write();
}

void write_amr_from_function(World& world, const Function<double, 3>& f,
                             const std::string& path, int m) {
  write_amr(world, f, path, m);
}
#else
// No-VTK build: collective-safe no-ops (nothing happens on any rank).
inline void write_htg_from_function(World&, const Function<double, 3>&,
                                    const std::string&) {}
inline void write_amr_from_function(World&, const Function<double, 3>&,
                                    const std::string&, int) {}
#endif  // MADNESS_HAS_VTK

// ---- native MRA coefficient export: .mad.h5 (io-hdf5 structured writer) -----
//
// Reuse the io-hdf5 structured Function writer (solvers/function_hdf5_io.hpp) to emit
// each function's RECONSTRUCTED octree as one <stem>.mad.h5: /meta (k, thresh, cell,
// tree_state, ...), /keys [n_nodes x (3+NDIM)] = (level, lx,ly,lz, has_children,
// has_coeff) per node, /coeffs [n_coeff_nodes x k^3] = the scaling-coefficient tensor
// per leaf. Unlike --amr (which point-SAMPLES each leaf onto an m^3 grid), this is the
// function's native multiresolution basis, so a reader can reconstruct the MADNESS
// expansion exactly at any point/zoom. "Don't roll a new serializer"
// we write through the same writer io-hdf5 built.
//
// The structured writer is double-only and MADNESS_CHECKs world.size()==1 (the local
// container == the whole tree only at NP=1) and t.size()==k^3 per coeff node (so the
// tree MUST be reconstructed -- a compressed (2k)^3 node would trip the check). We
// therefore reconstruct first and skip cleanly under NP>1. The skip is uniform across
// ranks (every rank evaluates world.size()>1 identically) -> no collective divergence.
#ifdef MADNESS_HAS_HDF5
void write_coeffs_hdf5(World& world, Function<double, 3> f,
                       const std::string& path) {
  if (world.size() > 1) {
    if (world.rank() == 0)
      print("  [HDF5] --coeffs needs NP=1 (structured writer is single-rank); "
            "skipping", path);
    return;
  }
  f.reconstruct();  // leaves -> k^3 scaling coefficients (the shape /coeffs expects)
  molresponse_v3::save_function_hdf5(f, path);
}
#else
// No-HDF5 build: collective-safe no-op (nothing happens on any rank).
inline void write_coeffs_hdf5(World&, Function<double, 3>, const std::string&) {}
#endif  // MADNESS_HAS_HDF5

// Dump one orbital's octree in both representations. `stem` is the output path
// prefix (no extension). With `htg`, also write the native HyperTreeGrid next to
// each JSON (<stem>_boxes.htg from the reconstructed tree, <stem>_error.htg from
// the compressed tree). Collective: must be called on every rank.
void dump_function_trees(World& world, Function<double, 3> f,
                         const std::string& stem, bool htg = false) {
  f.reconstruct();
  write_tree_json(world, f, stem + "_recon.json");
  if (htg) write_htg_from_function(world, f, stem + "_boxes.htg");
  f.compress();
  write_tree_json(world, f, stem + "_comp.json");
  if (htg) write_htg_from_function(world, f, stem + "_error.htg");
}

// What every dumped function gets: the trees, and the optional cube, AMR and
// native-coefficient files. Collective.
struct DumpOptions {
  std::string out_dir;
  int max_orbitals = -1;
  bool cube = false;
  int cube_npoints = 80;
  double cube_pad = 6.0;
  bool htg = false;
  bool amr = false;
  int amr_m = 8;
  bool coeffs = false;
};

void dump_one(World& world, const real_function_3d& f, const std::string& stem,
              const std::string& what, const Molecule& mol, const DumpOptions& o) {
  dump_function_trees(world, f, stem, o.htg);
  if (o.cube) write_cube_file(world, f, mol, stem + ".cube", o.cube_npoints, o.cube_pad);
  if (o.amr) write_amr_from_function(world, f, stem + ".vthb", o.amr_m);
  if (o.coeffs) write_coeffs_hdf5(world, f, stem + ".mad.h5");
  const double nrm = f.norm2();  // collective
  if (world.rank() == 0) print("  dumped", what, " norm2=", nrm);
}

std::size_t orbital_count(std::size_t n, int max_orbitals) {
  return max_orbitals >= 0 ? std::min<std::size_t>(n, static_cast<std::size_t>(max_orbitals)) : n;
}

// The ground orbitals at the active k, for the response and transition
// densities: a LOCAL copy reprojected when the ground state was loaded at
// another k (the dumped ground trees stay native).
ResponseGroundState ground_at_active_k(World& world, const vecfuncT& gs_amo) {
  const int k = FunctionDefaults<3>::get_k();
  const double thresh = FunctionDefaults<3>::get_thresh();
  ResponseGroundState rgs;
  rgs.amo.reserve(gs_amo.size());
  for (const auto& o : gs_amo)
    rgs.amo.push_back(o.k() == k ? o : project(o, k, thresh, false));
  world.gop.fence();
  return rgs;
}

// One driven state at its saved level: each response orbital (x, and y for a
// full layout) and the response density rho^(1). Kernels::compute_density is
// the single source of the closed-shell convention (static folds Y=X -> factor
// 4; full uses x+y -> factor 2). Output stems keep the earlier names,
// fd_<perturbation>__<protocol>__<frequency>_{x<i>,y<i>,rho1}. Collective.
template <typename Type>
void dump_driven(World& world, const saved_states::SavedNode& node, const vecfuncT& gs_amo,
                 const Molecule& mol, const DumpOptions& o) {
  using Storage = std::conditional_t<std::is_same_v<Type, Full>, ResponseStateXY<ClosedShell>,
                                     ResponseStateX<ClosedShell>>;
  const auto& d = node.description;
  const auto store = saved_states::load_driven<Storage>(world, node);
  const std::string label = "fd_" + d.settings.at("perturbation").get<std::string>() + "__" +
                            protocol_key(d.protocol.thresh, d.protocol.k) + "__" +
                            ResponseMetadata::freq_key(d.settings.at("frequency").get<double>());
  auto dump_block = [&](const vecfuncT& blk, const char* tag) {
    for (std::size_t i = 0; i < orbital_count(blk.size(), o.max_orbitals); ++i)
      dump_one(world, blk[i], o.out_dir + "/" + label + "_" + tag + std::to_string(i),
               label + " " + tag + std::to_string(i), mol, o);
  };
  dump_block(store.x_alpha, "x");
  if constexpr (std::is_same_v<Type, Full>) dump_block(store.y_alpha, "y");
  const auto rgs = ground_at_active_k(world, gs_amo);
  const real_function_3d rho1 = Kernels<Type, ClosedShell>::compute_density(world, rgs, store);
  dump_one(world, rho1, o.out_dir + "/" + label + "_rho1", label + " rho1", mol, o);
}

// One excited-state set at its saved level: x_i (and y_i for full) per root
// and the transition density (TDA: 2 sum phi_i x_i; full: 2 sum phi_i
// (x_i + y_i)). Output stems keep the earlier names, es_<slot>__<protocol>_*
// (pymra.web's field contract); with more than one set in the run the model is
// added to the protocol part, es_<slot>__<protocol>-<model>_*, so the sets do
// not overwrite each other. Returns the number of roots. Collective.
template <typename Type>
int dump_excited(World& world, const saved_states::SavedNode& node, bool name_model,
                 const vecfuncT& gs_amo, const Molecule& mol, const DumpOptions& o) {
  const auto& d = node.description;
  const auto state = saved_states::load_excited<Type>(world, node);
  std::string key = protocol_key(d.protocol.thresh, d.protocol.k);
  if (name_model) key += "-" + d.settings.at("excited_model").get<std::string>();
  const auto rgs = ground_at_active_k(world, gs_amo);
  const int nr = static_cast<int>(state.roots.size());
  for (int f = 0; f < nr; ++f) {
    const auto& root = state.roots[f];
    const std::string label = "es_" + std::to_string(f) + "__" + key;
    const std::size_t no = orbital_count(root.x_alpha.size(), o.max_orbitals);
    for (std::size_t i = 0; i < no; ++i)
      dump_one(world, root.x_alpha[i], o.out_dir + "/" + label + "_x" + std::to_string(i),
               label + " x" + std::to_string(i), mol, o);
    if constexpr (std::is_same_v<Type, Full>)
      for (std::size_t i = 0; i < no && i < root.y_alpha.size(); ++i)
        dump_one(world, root.y_alpha[i], o.out_dir + "/" + label + "_y" + std::to_string(i),
                 label + " y" + std::to_string(i), mol, o);
    const real_function_3d tdens = Kernels<Type, ClosedShell>::compute_density(world, rgs, root);
    dump_one(world, tdens, o.out_dir + "/" + label + "_tdens", label + " tdens", mol, o);
  }
  return nr;
}

// The saved driven states (fd) and excited-state sets (es) of the run folder
// `run`, found by listing its node folders and reading their descriptions.
// Stops when the run folder has nothing of an asked kind. Collective.
struct ResponseNodes {
  std::vector<saved_states::SavedNode> driven, excited;
};
ResponseNodes find_response_states(World& world, const std::filesystem::path& run, bool fd, bool es) {
  ResponseNodes found;
  for (auto& node : saved_states::response_nodes(world, run)) {
    if (fd && saved_states::is_driven(node)) found.driven.push_back(node);
    if (es && saved_states::is_excited(node)) found.excited.push_back(node);
  }
  if (fd && found.driven.empty())
    throw std::runtime_error("--fd: no saved driven state under " + (run / "states/response").string());
  if (es && found.excited.empty())
    throw std::runtime_error("--es: no saved excited-state set under " + (run / "states/response").string());
  return found;
}

// Dump the found states, each at its saved level. Open-shell driven states are
// listed and skipped, as before. Collective.
void dump_response_states(World& world, double L, const std::filesystem::path& run, const ResponseNodes& found,
                          const vecfuncT& gs_amo, const Molecule& mol, const DumpOptions& o) {
  const auto& driven = found.driven;
  const auto& excited = found.excited;

  auto at_saved_level = [&](const saved_states::SavedNode& node) {
    set_response_protocol(world, L, node.description.protocol.thresh, node.description.protocol.k, false);
    saved_states::print_loaded(world, node);
  };
  int npoints = 0;
  for (const auto& node : driven) {
    at_saved_level(node);
    if (!saved_states::closed_shell(node)) {
      if (world.rank() == 0) print("  [FD] skip, not closed shell:", saved_states::layout(node));
      continue;
    }
    if (saved_states::full_layout(node)) dump_driven<Full>(world, node, gs_amo, mol, o);
    else dump_driven<Static>(world, node, gs_amo, mol, o);
    ++npoints;
  }
  if (!driven.empty() && world.rank() == 0) print("DUMP_MRA_FD  run=", run.string(), "  states_dumped=", npoints);

  int nroots = 0;
  for (const auto& node : excited) {
    at_saved_level(node);
    const int n = static_cast<int>(node.description.details.at("n_roots").get<std::size_t>());
    if (world.rank() == 0) print("  [ES] set", node.folder.filename().string(), " layout=",
                                 saved_states::layout(node), " roots=", n);
    if (saved_states::full_layout(node))
      nroots += dump_excited<Full>(world, node, excited.size() > 1, gs_amo, mol, o);
    else
      nroots += dump_excited<TDA>(world, node, excited.size() > 1, gs_amo, mol, o);
  }
  if (!excited.empty() && world.rank() == 0)
    print("DUMP_MRA_ES  run=", run.string(), "  sets=", int(excited.size()), "  roots_dumped=", nroots);
}

}  // namespace

int main(int argc, char** argv) {
  World& world = initialize(argc, argv);
  int rc = 0;
  try {
    startup(world, argc, argv, true);

    commandlineparser parser(argc, argv);
    if (parser.key_exists("help") || parser.key_exists("archive")) {
      if (world.rank() == 0) {
        if (parser.key_exists("archive"))
          print("  --archive is replaced by --run-dir: the ground state is read from "
                "<run>/states/ground_state/.");
        print("Usage: dump_mra_trees [--run-dir=DIR] "
              "[--out=DIR] [--maxlevel=N] [--k=N] [--thresh=X] "
              "[--max-orbitals=N] [--cube] [--cube-npoints=N] [--cube-pad=X] "
              "[--htg] [--amr] [--amr-m=N] [--coeffs] [--fd] [--es]");
        print("  Loads the ground-state orbitals of the madqc run folder --run-dir");
        print("  (default: the current directory) from states/ground_state/ at the");
        print("  level its description records (--k/--thresh override), and");
        print("  dumps per-orbital MRA octree JSON (reconstructed + compressed),");
        print("  geometry.vtk, and (with --cube) per-orbital .cube isosurfaces.");
        print("  With --htg (needs a -DMADNESS_ENABLE_VTK=ON build), also writes");
        print("  each tree as a native VTK HyperTreeGrid (<stem>_boxes.htg /");
        print("  <stem>_error.htg) next to the legacy JSON/.vtk. With --amr,");
        print("  writes vtkOverlappingAMR (<stem>.vthb; function values on an");
        print("  --amr-m^3 grid/leaf box -> adaptive isosurfaces) for the ground");
        print("  MOs and (with --fd) the response orbitals + rho^(1).");
        print("  With --coeffs (needs -DMADNESS_ENABLE_HDF5=ON, NP=1), writes each");
        print("  function's NATIVE MRA coefficients as <stem>.mad.h5 (/meta+/keys+");
        print("  /coeffs: per-leaf n,l,s[k^3] + cell,k,thresh) -- the exact basis");
        print("  for offline reconstruction, not a point sample like .cube/.vthb.");
        print("  With --fd, also dumps the saved closed-shell driven states of the");
        print("  run folder (the older name --fd-calc-dir of --run-dir is still");
        print("  accepted) from its node folders");
        print("  under states/response/: x_i (+ y_i for a dynamic state) and rho^(1),");
        print("  stems fd_<perturbation>__<protocol>__<frequency>_*.");
        print("  With --es, also dumps every excited-state set there (TDA or full):");
        print("  x_i (+ y_i for full) per root and the transition density, stems");
        print("  es_<root>__<protocol>_*. Each state is dumped at its saved level.");
      }
      finalize();
      return 1;
    }

    // Inner scope so MADNESS objects destruct before finalize().
    {
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
      // Response states of a madqc run folder: driven states (--fd) and
      // excited-state sets (--es) from its node folders, into the same out dir,
      // so one mra_tree.py run yields a combined ground + response analysis.
      // --fd-calc-dir is the older name of --run-dir.
      const bool dump_fd = parser.key_exists("fd");
      const bool dump_es = parser.key_exists("es");
      const std::string run_dir =
          parser.key_exists("run-dir")       ? parser.value_raw("run-dir")
          : parser.key_exists("fd-calc-dir") ? parser.value_raw("fd-calc-dir")
                                             : std::string(".");
      // Find the ground state and the response states before any output, so a
      // run folder without them stops with nothing written.
      const auto ground = saved_states::read_ground(world, run_dir);
      ResponseNodes response_states;
      if (dump_fd || dump_es)
        response_states = find_response_states(world, run_dir, dump_fd, dump_es);
      // Native HyperTreeGrid (.htg) export of every tree, in addition to the
      // JSON/.vtk. Needs a -DMADNESS_ENABLE_VTK=ON build; otherwise a no-op.
      const bool htg = parser.key_exists("htg");
      // Adaptive function-value export (.vthb, vtkOverlappingAMR) -- samples each
      // leaf box on an m^3 grid (--amr-m, default 8). Covers the ground MOs and
      // (with --fd) the response orbitals + rho^(1), matching --cube/--htg.
      const bool amr = parser.key_exists("amr");
      const int amr_m = parser.key_exists("amr-m")
                            ? std::stoi(parser.value("amr-m")) : 8;
      // Native MRA coefficient export (.mad.h5) via the io-hdf5 structured writer:
      // per-leaf (n, l, s[k^3]) + header (cell, k, thresh). Needs a
      // -DMADNESS_ENABLE_HDF5=ON build and NP=1; otherwise skipped with a notice.
      // Covers the ground MOs and (with --fd) the response orbitals + rho^(1),
      // matching --cube/--amr.
      const bool coeffs = parser.key_exists("coeffs");
#ifndef MADNESS_HAS_VTK
      if ((htg || amr) && world.rank() == 0)
        print("  [VTK] --htg/--amr requested but this binary was built without "
              "VTK; reconfigure with -DMADNESS_ENABLE_VTK=ON. Skipping .htg/.vthb "
              "(JSON/.vtk/.cube output is unaffected).");
#endif
#ifndef MADNESS_HAS_HDF5
      if (coeffs && world.rank() == 0)
        print("  [HDF5] --coeffs requested but this binary was built without "
              "HDF5; reconfigure with -DMADNESS_ENABLE_HDF5=ON. Skipping .mad.h5 "
              "(JSON/.vtk/.cube/.vthb output is unaffected).");
#endif

      // Load at the node's saved level unless explicitly overridden. Matching
      // it means SCF::load_mos performs no reprojection, so the tree we walk
      // is exactly the one the solver produced. Reprojecting would rebuild the
      // tree and we'd be visualizing a different object.
      const int k = parser.key_exists("k") ? std::stoi(parser.value("k"))
                                           : ground.protocol().k;
      const double thresh = parser.key_exists("thresh")
                                ? std::stod(parser.value("thresh"))
                                : ground.protocol().thresh;
      auto gs = saved_states::load_ground(world, ground, thresh, k);
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
        print("DUMP_MRA  run=", run_dir, "  out=", out_dir,
              "  k=", k, "  thresh=", thresh, "  n_orbitals=", n);
      }

      for (std::size_t i = 0; i < n; ++i) {
        const std::string stem = out_dir + "/mo_" + std::to_string(i);
        dump_function_trees(world, mos[i], stem, htg);
        if (cube) {
          write_cube_file(world, mos[i], gs.molecule(), stem + ".cube",
                          cube_npoints, cube_pad);
        }
        if (amr) write_amr_from_function(world, mos[i], stem + ".vthb", amr_m);
        if (coeffs) write_coeffs_hdf5(world, mos[i], stem + ".mad.h5");
        const double nrm = mos[i].norm2();  // collective -- all ranks
        if (world.rank() == 0) print("  dumped mo_", i, "  norm2=", nrm);
      }
      world.gop.fence();

      // Response states, after the ground dump: each is loaded at its own saved
      // level (this changes FunctionDefaults). ground.L is the box length.
      if (dump_fd || dump_es) {
        const DumpOptions o{out_dir, max_orbitals, cube, cube_npoints, cube_pad,
                            htg, amr, amr_m, coeffs};
        dump_response_states(world, ground.L, run_dir, response_states, mos,
                             gs.molecule(), o);
        world.gop.fence();
      }
    }

    finalize();
    return rc;
  } catch (const std::exception& e) {
    if (world.rank() == 0) print("Error:", e.what());
    finalize();
    return 1;
  }
}
