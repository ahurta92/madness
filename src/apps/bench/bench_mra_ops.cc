/*
  bench_mra_ops.cc  --  microbenchmark of MADNESS MRA operations vs
  (number of functions n) x (protocol = thresh, k).

  Goal: emit one machine-readable (JSONL) record per
  (operation, n, k, thresh, P, threads) capturing wall time plus the
  quantities that drive the performance model in
  docs/parallel_runtime_guide/  --  tree size, leaf count, per-rank node
  imbalance, RMI message/byte traffic, and resident memory.

  Build:  add_subdirectory(bench) in src/apps/CMakeLists.txt, then
          ninja bench_mra_ops
  Run:    MAD_NUM_THREADS=8 mpiexec -np 4 ./bench_mra_ops \
              --k=6,8,10 --thresh=1e-4,1e-6 --n=1,4,16 --reps=3 \
              --out=bench.jsonl
  Analyze: python3 analyze.py bench.jsonl
*/

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/mra/vmra.h>
#include <madness/world/worldrmi.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace madness;

// ----------------------------------------------------------------------------
// Test-function generator: each function is a sum of a few Gaussians at
// deterministic (seeded) centers, so the tree for a given n is reproducible
// across protocols -- the only things that change N_leaf are k and thresh.
// special_points() forces refinement at the centers, mimicking orbital cusps.
// ----------------------------------------------------------------------------
class GaussianBlob : public FunctionFunctorInterface<double, 3> {
    std::vector<coord_3d> centers_;
    double inv2w2_;
public:
    GaussianBlob(std::vector<coord_3d> c, double width)
        : centers_(std::move(c)), inv2w2_(1.0 / (2.0 * width * width)) {}

    double operator()(const coord_3d& r) const override {
        double s = 0.0;
        for (const auto& c : centers_) {
            const double dx = r[0]-c[0], dy = r[1]-c[1], dz = r[2]-c[2];
            s += std::exp(-inv2w2_ * (dx*dx + dy*dy + dz*dz));
        }
        return s;
    }
    std::vector<coord_3d> special_points() const override { return centers_; }
};

static std::vector<real_function_3d>
make_functions(World& world, int n, int gpf, double width, double L, unsigned base_seed) {
    std::vector<real_function_3d> v;
    v.reserve(n);
    const double half = 0.4 * L;   // keep centers well inside the cell
    for (int i = 0; i < n; ++i) {
        std::mt19937 gen(base_seed + 7919u * (unsigned)i);
        std::uniform_real_distribution<double> U(-half, half);
        std::vector<coord_3d> centers(gpf);
        for (int g = 0; g < gpf; ++g) {
            coord_3d c; c[0] = U(gen); c[1] = U(gen); c[2] = U(gen);
            centers[g] = c;
        }
        auto functor = std::shared_ptr<FunctionFunctorInterface<double,3>>(
            new GaussianBlob(centers, width));
        v.push_back(real_factory_3d(world).functor(functor));
    }
    world.gop.fence();
    truncate(world, v);            // canonical reconstructed+truncated state
    return v;
}

// ----------------------------------------------------------------------------
// Metrics helpers
// ----------------------------------------------------------------------------
struct TreeMetrics {
    double nodes_total = 0, leaves_total = 0;   // global, summed over the vector
    double nodes_rank_min = 0, nodes_rank_max = 0;  // local-node count, across ranks
    double max_depth = 0;
};

static std::size_t count_leaves_local(const real_function_3d& f) {
    std::size_t c = 0;
    const auto& coeffs = f.get_impl()->get_coeffs();
    for (auto it = coeffs.begin(); it != coeffs.end(); ++it)
        if (!it->second.has_children()) ++c;   // leaf == no children
    return c;
}

static TreeMetrics tree_metrics(World& world, const std::vector<real_function_3d>& v) {
    TreeMetrics m;
    std::size_t local_nodes = 0, local_leaves = 0;
    double maxd = 0;
    for (const auto& f : v) {
        local_nodes  += f.get_impl()->get_coeffs().size();
        local_leaves += count_leaves_local(f);
        maxd = std::max<double>(maxd, (double)f.max_depth());
    }
    double nodes = (double)local_nodes, leaves = (double)local_leaves;
    double nmin = (double)local_nodes, nmax = (double)local_nodes;
    world.gop.sum(&nodes, 1);
    world.gop.sum(&leaves, 1);
    world.gop.min(&nmin, 1);
    world.gop.max(&nmax, 1);
    world.gop.max(&maxd, 1);
    m.nodes_total = nodes; m.leaves_total = leaves;
    m.nodes_rank_min = nmin; m.nodes_rank_max = nmax; m.max_depth = maxd;
    return m;
}

struct RmiDelta { double nmsg_sent=0, nbyte_sent=0, nmsg_recv=0, nbyte_recv=0; };

static RmiDelta rmi_snapshot() {
    const auto& s = RMI::get_stats();
    return { (double)s.nmsg_sent, (double)s.nbyte_sent,
             (double)s.nmsg_recv, (double)s.nbyte_recv };
}
static RmiDelta rmi_diff_global(World& world, const RmiDelta& a, const RmiDelta& b) {
    RmiDelta d{ b.nmsg_sent-a.nmsg_sent, b.nbyte_sent-a.nbyte_sent,
                b.nmsg_recv-a.nmsg_recv, b.nbyte_recv-a.nbyte_recv };
    world.gop.sum(&d.nmsg_sent, 1);  world.gop.sum(&d.nbyte_sent, 1);
    world.gop.sum(&d.nmsg_recv, 1);  world.gop.sum(&d.nbyte_recv, 1);
    return d;
}

static double proc_status_kb(const char* key) {
    std::ifstream f("/proc/self/status");
    std::string line;
    const std::size_t klen = std::strlen(key);
    while (std::getline(f, line)) {
        if (line.compare(0, klen, key) == 0) {
            std::istringstream iss(line.substr(klen));
            long kb = 0; iss >> kb; return (double)kb;
        }
    }
    return -1.0;
}

// ----------------------------------------------------------------------------
// Timing: warmup once, then `reps` timed runs; report median/min/max.
// `prepare()` builds fresh inputs OUTSIDE the timed region; `run()` is timed.
// ----------------------------------------------------------------------------
struct Timing { double median=0, tmin=0, tmax=0; int reps=0; RmiDelta rmi; };

template <class Prepare, class Run>
static Timing time_op(World& world, int reps, Prepare prepare, Run run) {
    Timing t; t.reps = reps;
    { auto in = prepare(); world.gop.fence(); run(in); world.gop.fence(); } // warmup
    std::vector<double> ts;
    RmiDelta last{};
    for (int r = 0; r < reps; ++r) {
        auto in = prepare();
        world.gop.fence();
        RmiDelta s0 = rmi_snapshot();
        double t0 = wall_time();
        run(in);
        world.gop.fence();
        double t1 = wall_time();
        RmiDelta s1 = rmi_snapshot();
        ts.push_back(t1 - t0);
        last = rmi_diff_global(world, s0, s1);
    }
    std::sort(ts.begin(), ts.end());
    t.tmin = ts.front(); t.tmax = ts.back();
    t.median = ts[ts.size()/2];
    t.rmi = last;
    return t;
}

// ----------------------------------------------------------------------------
// JSONL output (manual, no external json dependency)
// ----------------------------------------------------------------------------
static void emit(std::ostream& os, const std::string& op, int n, int k, double thresh,
                 int P, int threads, double L,
                 const Timing& t, const TreeMetrics& m) {
    const double imbalance = (m.nodes_total > 0)
        ? m.nodes_rank_max * (double)P / m.nodes_total : 1.0;
    const double bytes_est = m.leaves_total * std::pow((double)k, 3) * 8.0;
    os << "{"
       << "\"op\":\"" << op << "\","
       << "\"n\":" << n << ",\"k\":" << k << ",\"thresh\":" << thresh << ","
       << "\"P\":" << P << ",\"threads\":" << threads << ",\"L\":" << L << ","
       << "\"time_s\":{\"median\":" << t.median << ",\"min\":" << t.tmin
            << ",\"max\":" << t.tmax << ",\"reps\":" << t.reps << "},"
       << "\"tree\":{\"nodes\":" << m.nodes_total << ",\"leaves\":" << m.leaves_total
            << ",\"nodes_rank_min\":" << m.nodes_rank_min
            << ",\"nodes_rank_max\":" << m.nodes_rank_max
            << ",\"imbalance\":" << imbalance
            << ",\"max_depth\":" << m.max_depth << "},"
       << "\"coeff_bytes_est\":" << bytes_est << ","
       << "\"rmi\":{\"nmsg_sent\":" << t.rmi.nmsg_sent
            << ",\"nbyte_sent\":" << t.rmi.nbyte_sent
            << ",\"nmsg_recv\":" << t.rmi.nmsg_recv
            << ",\"nbyte_recv\":" << t.rmi.nbyte_recv << "}"
       << "}\n";
}

// ----------------------------------------------------------------------------
// Tiny CLI helpers
// ----------------------------------------------------------------------------
static std::string opt(int argc, char** argv, const std::string& key, const std::string& def) {
    const std::string pfx = "--" + key + "=";
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a.compare(0, pfx.size(), pfx) == 0) return a.substr(pfx.size());
    }
    return def;
}
static std::vector<double> parse_doubles(const std::string& s) {
    std::vector<double> out; std::stringstream ss(s); std::string tok;
    while (std::getline(ss, tok, ',')) if (!tok.empty()) out.push_back(std::stod(tok));
    return out;
}
static std::vector<int> parse_ints(const std::string& s) {
    std::vector<int> out; for (double d : parse_doubles(s)) out.push_back((int)d); return out;
}

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    startup(world, argc, argv);

    const std::vector<int>    ks      = parse_ints(opt(argc, argv, "k", "6,8,10"));
    const std::vector<double> threshs = parse_doubles(opt(argc, argv, "thresh", "1e-4,1e-6"));
    const std::vector<int>    ns      = parse_ints(opt(argc, argv, "n", "1,4,16"));
    const double L     = std::stod(opt(argc, argv, "L", "20.0"));
    const int    reps  = std::stoi(opt(argc, argv, "reps", "3"));
    const int    gpf   = std::stoi(opt(argc, argv, "gaussians", "3"));
    const double width = std::stod(opt(argc, argv, "width", "0.6"));
    const std::string outpath = opt(argc, argv, "out", "bench_results.jsonl");
    const std::string optype  = opt(argc, argv, "operator", "bsh"); // bsh | coulomb

    const int P = world.size();
    const int threads = ThreadPool::size() + 1; // workers + main (excludes RMI thread)

    std::ofstream os;
    if (world.rank() == 0) {
        os.open(outpath);
        os << "{\"record\":\"meta\",\"P\":" << P << ",\"threads\":" << threads
           << ",\"L\":" << L << ",\"gaussians\":" << gpf << ",\"width\":" << width
           << ",\"reps\":" << reps << ",\"operator\":\"" << optype << "\"}\n";
        print("bench_mra_ops: P=", P, " threads/rank=", threads,
              " ks=", ks.size(), " threshs=", threshs.size(), " ns=", ns.size());
    }

    for (int k : ks) {
        for (double thresh : threshs) {
            FunctionDefaults<3>::set_k(k);
            FunctionDefaults<3>::set_thresh(thresh);
            FunctionDefaults<3>::set_cubic_cell(-L/2, L/2);
            FunctionDefaults<3>::set_refine(true);
            FunctionDefaults<3>::set_initial_level(2);
            FunctionDefaults<3>::set_truncate_mode(1);

            // Operator built once per protocol (excluded from per-op timing).
            const double lo = 1e-4;
            real_convolution_3d Op = (optype == "coulomb")
                ? CoulombOperator(world, lo, thresh)
                : BSHOperator3D(world, 1.0, lo, thresh);

            for (int n : ns) {
                // Two independent reconstructed+truncated operand sets.
                std::vector<real_function_3d> A = make_functions(world, n, gpf, width, L, 11u);
                std::vector<real_function_3d> B = make_functions(world, n, gpf, width, L, 101u);
                TreeMetrics m = tree_metrics(world, A);

                auto record = [&](const std::string& name, const Timing& t) {
                    if (world.rank() == 0) emit(os, name, n, k, thresh, P, threads, L, t, m);
                };

                // project: cost of constructing the n functions from functors.
                record("project", time_op(world, reps,
                    [&]{ return 0; },
                    [&](int){ auto v = make_functions(world, n, gpf, width, L, 11u); }));

                // compress: reconstructed -> compressed (bottom-up sweep).
                record("compress", time_op(world, reps,
                    [&]{ reconstruct(world, A); return copy(world, A); },
                    [&](std::vector<real_function_3d>& v){ compress(world, v); }));

                // reconstruct: compressed -> reconstructed (top-down sweep).
                record("reconstruct", time_op(world, reps,
                    [&]{ reconstruct(world, A); auto v = copy(world, A); compress(world, v); return v; },
                    [&](std::vector<real_function_3d>& v){ reconstruct(world, v); }));

                // truncate: norm-based pruning (two sweeps, two fences).
                record("truncate", time_op(world, reps,
                    [&]{ reconstruct(world, A); return copy(world, A); },
                    [&](std::vector<real_function_3d>& v){ truncate(world, v); }));

                // gaxpy: A += B  (compressed form; copy of A so A is untouched).
                record("gaxpy", time_op(world, reps,
                    [&]{ auto v = copy(world, A); compress(world, v); compress(world, B); return v; },
                    [&](std::vector<real_function_3d>& v){ gaxpy(world, 1.0, v, 1.0, B, true); }));

                // multiply: elementwise A[i]*B[i]  (reconstructed inputs).
                record("multiply", time_op(world, reps,
                    [&]{ reconstruct(world, A); reconstruct(world, B); return 0; },
                    [&](int){ auto v = mul(world, A, B); }));

                // matrix_inner: n x n Gram matrix <A|B> (one all-reduce).
                record("matrix_inner", time_op(world, reps,
                    [&]{ compress(world, A); compress(world, B); return 0; },
                    [&](int){ Tensor<double> g = matrix_inner(world, A, B); }));

                // apply: BSH/Coulomb convolution on each of n functions (the heavy op).
                record("apply", time_op(world, reps,
                    [&]{ reconstruct(world, A); return 0; },
                    [&](int){ auto v = apply(world, Op, A); }));

                if (world.rank() == 0)
                    print("  done n=", n, " k=", k, " thresh=", thresh,
                          " leaves=", (long)m.leaves_total, " imbalance=", m.nodes_rank_max*(double)P/std::max(1.0,m.nodes_total));
            }
        }
    }

    if (world.rank() == 0) { os.close(); print("wrote ", outpath); }
    world.gop.fence();
    finalize();
    return 0;
}
