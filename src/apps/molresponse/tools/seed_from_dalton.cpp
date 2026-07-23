// seed_from_dalton — build a molresponse TDA excited-state SEED bundle by
// reading a DALTON RSPVEC binary and Molden file directly (no Python/HDF5
// intermediate). Projects the per-occupied-orbital response function
//   x_i(r) = sum_mu D[mu,i] chi_mu(r),  D = C_vir @ X.T
// into MADNESS Function<double,3> via a DaltonResponseFunctor, wraps all N
// functions as a one-root ESSolver<TDA,ClosedShell>::State, and calls
// save_es_roots(..., converged=false). The resulting bundle is identical in
// schema to the HDF5-seeded one produced by seed_es_from_hdf5; the ES node
// Resumes from it on the next run.
//
// Constraints: NP=1 only (same as seed_es_from_hdf5). Cell is set to
// [-200,200]^3 here (not read from an archive); k and thresh must match
// the solver protocol you want to warm-start.
//
// Usage:
//   seed_from_dalton --rspvec=RSPVEC --molden=molden.inp --n-occ=N
//                    --root=0 --omega=<au> --calc-dir=DIR
//                    [--thresh=1e-4] [--scale=1.41421356]

#include "dalton_rspvec.hpp"
#include "dalton_gto.hpp"

#include "../ResponseProtocol.hpp"
#include "../kernels/tags.hpp"
#include "../solvers/response_state.hpp"
#include "../solvers/es_solver.hpp"
#include "../solvers/es_save_load.hpp"

#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <cmath>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace madness;
using namespace molresponse_v3;

namespace {

class DaltonResponseFunctor : public madness::FunctionFunctorInterface<double, 3> {
    const DaltonMoldenBasis& basis;
    std::vector<double> weights;    // D[:,i] for orbital i
    std::vector<madness::coord_3d> centers;

public:
    DaltonResponseFunctor(const DaltonMoldenBasis& b, std::vector<double> w)
        : basis(b), weights(std::move(w)) {
        for (const auto& sh : basis.shells)
            centers.push_back({sh.cx, sh.cy, sh.cz});
    }

    double operator()(const madness::coord_3d& r) const override {
        double val = 0.0;
        double bf[9];  // max 9 components for g shell
        for (size_t s = 0; s < basis.shells.size(); s++) {
            const auto& sh = basis.shells[s];
            sh.evaluate(r[0], r[1], r[2], bf);
            const int off = basis.ao_offsets[s];
            for (int k = 0; k < sh.n_ao; k++)
                val += weights[static_cast<size_t>(off + k)] * bf[k];
        }
        return val;
    }

    std::vector<madness::coord_3d> special_points() const override {
        return centers;
    }
};

} // namespace

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    startup(world, argc, argv, true);

    commandlineparser parser(argc, argv);

    if (world.size() != 1) {
        if (world.rank() == 0)
            print("ERROR: seed_from_dalton must run on a single rank (NP=1).");
        finalize();
        return 2;
    }

    if (!parser.key_exists("rspvec") || !parser.key_exists("molden") ||
        !parser.key_exists("n-occ")  || !parser.key_exists("calc-dir") ||
        !(parser.key_exists("root") || parser.key_exists("roots")) ||
        !(parser.key_exists("omega") || parser.key_exists("omegas"))) {
        if (world.rank() == 0) {
            print("Usage: seed_from_dalton --rspvec=RSPVEC --molden=molden.inp "
                  "--n-occ=N --root=IDX --omega=<au> --calc-dir=DIR");
            print("  multi-root: --roots=0,1,2 --omegas=w0,w1,w2 (one N-root bundle)");
            print("  [--thresh=1e-4] [--scale=<sqrt2>]");
        }
        finalize();
        return 2;
    }

    const std::string rspvec_path = parser.value_raw("rspvec");
    const std::string molden_path = parser.value_raw("molden");
    const std::string calc_dir   = parser.value_raw("calc-dir");
    const int         n_occ      = std::stoi(parser.value("n-occ"));
    // Multi-root: --roots=0,1,2 --omegas=w0,w1,w2 build ONE N-root seed
    // bundle (what a production --es-roots=N solve warm-starts from).
    // --root/--omega remain as the single-root spelling.
    std::vector<int>    root_list;
    std::vector<double> omega_list;
    {
        std::stringstream rs(parser.key_exists("roots") ? parser.value("roots")
                                                        : parser.value("root"));
        std::stringstream os(parser.key_exists("omegas")
                                 ? parser.value("omegas")
                                 : parser.value("omega"));
        std::string t;
        while (std::getline(rs, t, ','))
            if (!t.empty()) root_list.push_back(std::stoi(t));
        while (std::getline(os, t, ','))
            if (!t.empty()) omega_list.push_back(std::stod(t));
    }
    if (root_list.empty() || root_list.size() != omega_list.size()) {
        if (world.rank() == 0)
            print("ERROR: --roots and --omegas must be same-length CSV lists.");
        finalize();
        return 2;
    }
    const double      thresh     =
        parser.key_exists("thresh") ? std::stod(parser.value("thresh")) : 1e-4;
    const double      scale      =
        parser.key_exists("scale") ? std::stod(parser.value("scale"))
                                   : std::sqrt(2.0);

    const int k = default_k_for_thresh(thresh);
    const std::string key    = protocol_key(thresh, k);
    const std::string bundle = calc_dir + "/es__" + key;

    {
        // Read RSPVEC
        auto [info, entries] = read_rspvec(rspvec_path);
        for (int ri : root_list) {
            if (ri < 0 || ri >= static_cast<int>(entries.size())) {
                if (world.rank() == 0)
                    print("ERROR: root index", ri, "out of range (file has",
                          static_cast<int>(entries.size()), "entries).");
                finalize();
                return 2;
            }
        }

        // Read Molden
        DaltonMoldenResult molden = read_molden(molden_path);
        const int n_mo  = molden.n_mo;
        const int n_ao  = molden.n_ao;
        const int n_vir = n_mo - n_occ;

        if (n_vir <= 0) {
            if (world.rank() == 0)
                print("ERROR: n_occ=", n_occ, ">= n_mo=", n_mo);
            finalize();
            return 2;
        }

        // Set MADNESS cell and discretization parameters BEFORE projecting
        Tensor<double> cell(3L, 2L);
        for (int i = 0; i < 3; i++) { cell(i, 0) = -200.0; cell(i, 1) = 200.0; }
        FunctionDefaults<3>::set_cell(cell);
        FunctionDefaults<3>::set_k(k);
        FunctionDefaults<3>::set_thresh(thresh);

        // One N-root State: per root, split the RSPVEC entry, build
        // D = C_vir @ X.T and project one Function per occupied orbital.
        ESSolver<TDA, ClosedShell>::State s;
        const size_t NR = root_list.size();
        s.omega = Tensor<double>(static_cast<long>(NR));
        std::vector<double> xnorms;
        for (size_t rr = 0; rr < NR; ++rr) {
        const auto& entry = entries[static_cast<size_t>(root_list[rr])];

        // Split RSPVEC entry into X (and optional Y)
        auto [X_flat, Y_flat] = split_ov(entry.vec, n_occ, n_vir);
        // X_flat: n_occ*n_vir, row-major (occ outer, vir inner)

        // Check TDA norm: ||X||^2 should be ~0.5
        double xnorm2 = 0.0;
        for (double v : X_flat) xnorm2 += v * v;
        if (world.rank() == 0 && std::abs(xnorm2 - 0.5) > 0.05)
            print("WARNING: root", root_list[rr], "||X||^2 =", xnorm2,
                  "(expected ~0.5 for TDA)");
        xnorms.push_back(xnorm2);

        // Compute D = C_vir @ X.T: D[mu, i] = sum_a C[mu, n_occ+a] * X[i, a]
        // D is n_ao x n_occ; stored column-major (D_col_i = D[0..n_ao-1, i]).
        // mo_coeffs is column-major: C[mu + n_ao*mo_idx]
        std::vector<double> D(static_cast<size_t>(n_ao * n_occ), 0.0);
        for (int mu = 0; mu < n_ao; mu++) {
            for (int i = 0; i < n_occ; i++) {
                double val = 0.0;
                for (int a = 0; a < n_vir; a++) {
                    // C[mu, n_occ+a] = mo_coeffs[mu + n_ao*(n_occ+a)]
                    double c_vir = molden.mo_coeffs[
                        static_cast<size_t>(mu + n_ao * (n_occ + a))];
                    // X[i, a] = X_flat[i*n_vir + a]
                    double x_ia = X_flat[
                        static_cast<size_t>(i * n_vir + a)];
                    val += c_vir * x_ia;
                }
                // D[mu, i] = val -> column i at position mu
                D[static_cast<size_t>(mu + n_ao * i)] = val;
            }
        }

        // Project one Function<double,3> per occupied orbital
        ResponseStateX<ClosedShell> root;
        for (int i = 0; i < n_occ; i++) {
            std::vector<double> D_col_i(
                D.begin() + static_cast<ptrdiff_t>(i * n_ao),
                D.begin() + static_cast<ptrdiff_t>((i + 1) * n_ao));

            // Declare as base type so FunctionFactory::functor() picks the
            // shared_ptr<FunctionFunctorInterface> overload, not the generic
            // template callable overload (which would try op(coord) on the ptr).
            std::shared_ptr<madness::FunctionFunctorInterface<double, 3>> ffunctor =
                std::make_shared<DaltonResponseFunctor>(molden.basis, std::move(D_col_i));

            Function<double, 3> fn =
                FunctionFactory<double, 3>(world)
                    .functor(ffunctor)
                    .thresh(thresh)
                    .truncate_on_project();

            if (scale != 1.0) fn.scale(scale);
            root.x_alpha.push_back(fn);
        }

        s.roots.push_back(std::move(root));
        s.omega(static_cast<long>(rr)) = omega_list[rr];
        }  // per-root loop

        s.iter = 0;
        s.last_bsh_residual     = std::vector<double>(NR, 1.0);
        s.last_density_residual = std::vector<double>(NR, 1.0);
        s.last_omega_residual   = std::vector<double>(NR, 1.0);

        save_es_roots<TDA, ClosedShell>(world, s, bundle, /*converged=*/false);

        if (world.rank() == 0) {
            print("seed_from_dalton: wrote", static_cast<int>(NR),
                  "root(s) x", n_occ, "orbital functions");
            print("  rspvec    =", rspvec_path);
            print("  molden    =", molden_path);
            for (size_t rr = 0; rr < NR; ++rr)
                print("  root", root_list[rr], ": omega(au) =", omega_list[rr],
                      "  ||X||^2 =", xnorms[rr], "  scale =", scale);
            print("  n_occ =", n_occ, "  n_vir =", n_vir,
                  "  n_ao =", n_ao, "  n_mo =", n_mo);
            print("  protocol_key =", key, "  k =", k, "  thresh =", thresh);
            print("  bundle    =", bundle);
            print("  Run (NP=1): test_calc_manager_run --archive=<gs> --calc-dir=",
                  calc_dir, " --es-roots=", static_cast<int>(NR),
                  " --protocol=", thresh, " -> ES node Resumes from seed.");
        }
    }  // all Function/State objects destruct here, before finalize()

    finalize();
    return 0;
}
