// ===========================================================================
// CalcManager scheduling-core tests (doc 15, 15a). Pure C++, no MPI / no MRA
// runtime — exercises build_dag, reconcile_rung, deps_ready, and schedule
// against hand-written response_metadata.json blobs.
// ===========================================================================

#include "../calc/calc_manager.hpp"
#include "../ResponsePropertyPlanner.hpp"

#include <madness/external/nlohmann_json/json.hpp>

#include <algorithm>
#include <cstdio>
#include <string>

using namespace molresponse_v3;
using nlohmann::json;

namespace {

int failed = 0;

#define EXPECT(cond, label)                                                \
  do {                                                                     \
    if (cond) { std::printf("  [PASS]  %s\n", label); }                    \
    else      { std::printf("  [FAIL]  %s\n", label); ++failed; }          \
  } while (0)

const CalcNode *find_id(const std::vector<CalcNode> &dag, const std::string &id) {
  for (const auto &n : dag)
    if (n.id == id) return &n;
  return nullptr;
}

bool wave_has(const std::vector<WorkItem> &w, const std::string &id) {
  for (const auto &it : w)
    if (it.node->id == id) return true;
  return false;
}

NodeAction action_in(const std::vector<WorkItem> &w, const std::string &id) {
  for (const auto &it : w)
    if (it.node->id == id) return it.action;
  return NodeAction::Skip; // sentinel "not present"
}

// Minimal empty aggregate.
json empty_meta() {
  return json{{"schema_version", 1},
              {"protocols", json::object()},
              {"fd_states", json::object()},
              {"excited_states", json::object()},
              {"properties", json::object()}};
}

// Register a protocol key in the registry.
void put_protocol(json &m, double thresh) {
  const int k = default_k_for_thresh(thresh);
  m["protocols"][protocol_key(thresh, k)] = {{"thresh", thresh}, {"k", k},
                                             {"index", -1}};
}

// Write an fd_states entry.
void put_fd(json &m, const std::string &pert, double thresh, double freq,
            bool converged, bool diverged = false) {
  const std::string key = rung_key(thresh);
  const std::string fk  = ResponseMetadata::freq_key(freq);
  put_protocol(m, thresh);
  m["fd_states"][pert][key][fk] = {{"freq", freq},
                                   {"converged", converged},
                                   {"diverged", diverged}};
}

void put_es(json &m, double thresh, bool converged, bool diverged = false) {
  const std::string key = rung_key(thresh);
  put_protocol(m, thresh);
  m["excited_states"][key] = {{"converged", converged},
                              {"diverged", diverged}};
}

} // namespace

int main() {
  const std::vector<double> P = {1e-4, 1e-6};

  // ====== build_dag: polarizability on x,y,z ===============================
  std::printf("=== build_dag: alpha(0.057) x/y/z ===\n");
  {
    ResponsePropertyRequest r;
    r.kind = ResponsePropertyKind::Polarizability;
    r.frequencies = {0.057};
    r.protocol_thresholds = P;
    auto dag = build_dag(plan_one(r), /*n_atoms=*/3);
    EXPECT(dag.size() == 3, "3 FD nodes");
    EXPECT(find_id(dag, "fd:dipole_x@f0.05700") != nullptr, "dipole_x node id");
    const CalcNode *zx = find_id(dag, "fd:dipole_z@f0.05700");
    EXPECT(zx && zx->kind == CalcKind::FD && std::abs(zx->freq - 0.057) < 1e-9,
           "dipole_z node fields");
    EXPECT(zx && zx->seed_from.empty(),
           "no static seed hint when no ω=0 in plan");
  }

  // ====== build_dag: SHG seeds dynamics from static of same channel ========
  std::printf("=== build_dag: static-seed hint ===\n");
  {
    // Two dipole_x FD: ω=0 and ω=0.1 -> dynamic should hint at the static.
    ResponsePlan plan;
    plan.fd.push_back({Perturbation::dipole(0), 0.0, P});
    plan.fd.push_back({Perturbation::dipole(0), 0.1, P});
    auto dag = build_dag(plan, 0);
    const CalcNode *dyn = find_id(dag, "fd:dipole_x@f0.10000");
    EXPECT(dyn && dyn->seed_from == "fd:dipole_x@f0.00000",
           "dynamic seeds from same-channel static");
    const CalcNode *stat = find_id(dag, "fd:dipole_x@f0.00000");
    EXPECT(stat && stat->seed_from.empty(), "static has no seed hint");
  }

  // ====== build_dag: nuclear_all expands to 3N ============================
  std::printf("=== build_dag: nuclear_all expansion ===\n");
  {
    ResponsePropertyRequest r;
    r.kind = ResponsePropertyKind::PolarizabilityGradient;
    r.gradient_mode = GradientMode::Nuclear;
    r.frequencies = {0.0};
    r.protocol_thresholds = P;
    auto plan = plan_one(r);
    auto dag = build_dag(plan, /*n_atoms=*/2);  // 2 atoms × 3 axes = 6 nuc nodes
    int nuc = 0;
    for (const auto &n : dag)
      if (n.kind == CalcKind::NuclearFD) ++nuc;
    EXPECT(nuc == 6, "2 atoms × 3 axes = 6 NuclearFD nodes");
    EXPECT(find_id(dag, "fd:nuc_0_x@f0.00000") != nullptr, "nuc atom0 x");
    EXPECT(find_id(dag, "fd:nuc_1_z@f0.00000") != nullptr, "nuc atom1 z");
    for (const auto &n : dag)
      EXPECT(!(n.kind == CalcKind::NuclearFD && n.pert.is_all_atoms()),
             "no unexpanded all-atoms sentinel remains");
  }

  // ====== build_dag: resonance Raman -> ES + symbolic derived FD ===========
  std::printf("=== build_dag: resonant gradient (ES + derived) ===\n");
  {
    ResponsePropertyRequest r;
    r.kind = ResponsePropertyKind::PolarizabilityGradient;
    r.gradient_mode = GradientMode::Resonant;
    r.n_roots = 4;
    r.protocol_thresholds = P;
    auto dag = build_dag(plan_one(r), 0);
    const CalcNode *es = find_id(dag, "es:tda_n4");
    EXPECT(es && es->kind == CalcKind::ES && es->n_roots == 4, "ES node");
    const CalcNode *d = find_id(dag, "dfd:dipole_x@*");
    EXPECT(d && d->is_symbolic(), "derived FD is symbolic");
    EXPECT(d && d->hard_deps.size() == 1 && d->hard_deps[0] == "es:tda_n4",
           "derived FD hard-deps the ES bundle");
  }

  // ====== reconcile_rung: the 5-row table =================================
  std::printf("=== reconcile_rung table ===\n");
  {
    CalcNode fd;
    fd.kind = CalcKind::FD;
    fd.pert = Perturbation::dipole(0);
    fd.freq = 0.057;
    fd.protocols = P;
    fd.id = fd_node_id(fd.pert, fd.freq);

    EXPECT(reconcile_rung(fd, empty_meta(), 1e-6) == NodeAction::Fresh,
           "absent -> Fresh");

    json m1 = empty_meta();
    put_fd(m1, "dipole_x", 1e-6, 0.057, /*converged=*/true);
    EXPECT(reconcile_rung(fd, m1, 1e-6) == NodeAction::Skip,
           "converged at rung -> Skip");

    json m2 = empty_meta();
    put_fd(m2, "dipole_x", 1e-4, 0.057, /*converged=*/true);  // coarser only
    EXPECT(reconcile_rung(fd, m2, 1e-6) == NodeAction::Restart,
           "converged at coarser rung -> Restart");

    json m3 = empty_meta();
    put_fd(m3, "dipole_x", 1e-6, 0.057, /*converged=*/false, /*diverged=*/false);
    EXPECT(reconcile_rung(fd, m3, 1e-6) == NodeAction::Resume,
           "present, not converged, not diverged -> Resume");

    json m4 = empty_meta();
    put_fd(m4, "dipole_x", 1e-6, 0.057, /*converged=*/false, /*diverged=*/true);
    EXPECT(reconcile_rung(fd, m4, 1e-6) == NodeAction::Fresh,
           "present, diverged -> Fresh");
  }

  // ====== reconcile_rung: ES path =========================================
  std::printf("=== reconcile_rung ES ===\n");
  {
    CalcNode es;
    es.kind = CalcKind::ES; es.tda = true; es.n_roots = 3; es.protocols = P;
    es.id = es_node_id(true, 3);
    json m = empty_meta();
    put_es(m, 1e-6, /*converged=*/true);
    EXPECT(reconcile_rung(es, m, 1e-6) == NodeAction::Skip, "ES converged -> Skip");
    json m2 = empty_meta();
    put_es(m2, 1e-4, /*converged=*/true);
    EXPECT(reconcile_rung(es, m2, 1e-6) == NodeAction::Restart,
           "ES coarser converged -> Restart");
  }

  // ====== deps_ready: VBC-like synthetic dep ==============================
  std::printf("=== deps_ready (synthetic VBC dep) ===\n");
  {
    std::vector<CalcNode> dag;
    CalcNode b; b.kind = CalcKind::FD; b.pert = Perturbation::dipole(0);
    b.freq = 0.05; b.protocols = P; b.id = fd_node_id(b.pert, b.freq);
    CalcNode c; c.kind = CalcKind::FD; c.pert = Perturbation::dipole(1);
    c.freq = 0.07; c.protocols = P; c.id = fd_node_id(c.pert, c.freq);
    CalcNode vbc; vbc.kind = CalcKind::VBC; vbc.id = "vbc:test";
    vbc.protocols = P; vbc.hard_deps = {b.id, c.id};
    dag = {b, c, vbc};

    json m = empty_meta();
    EXPECT(!deps_ready(vbc, dag, m, 1e-6), "deps absent -> not ready");
    put_fd(m, "dipole_x", 1e-6, 0.05, true);
    EXPECT(!deps_ready(vbc, dag, m, 1e-6), "one dep ready -> still not ready");
    put_fd(m, "dipole_y", 1e-6, 0.07, true);
    EXPECT(deps_ready(vbc, dag, m, 1e-6), "both deps converged -> ready");
    // Same deps but at a coarser rung must NOT satisfy a finer rung.
    json m2 = empty_meta();
    put_fd(m2, "dipole_x", 1e-4, 0.05, true);
    put_fd(m2, "dipole_y", 1e-4, 0.07, true);
    EXPECT(!deps_ready(vbc, dag, m2, 1e-6),
           "deps converged only at coarser rung -> not ready at finer");
  }

  // ====== schedule: two-phase, channel-parallel ===========================
  std::printf("=== schedule: alpha x/y/z over [1e-4,1e-6] ===\n");
  {
    ResponsePropertyRequest r;
    r.kind = ResponsePropertyKind::Polarizability;
    r.frequencies = {0.0, 0.1};            // static + dynamic per axis
    r.protocol_thresholds = P;
    auto dag = build_dag(plan_one(r), 0);  // 6 FD nodes (3 axes × 2 freqs)
    auto waves = schedule(dag, P, empty_meta());

    // Expect: A1 (3 anchors = the ω=0 of each axis), A2 (3 dynamics @1e-4),
    //         then B (all 6 @1e-6).
    EXPECT(waves.size() == 3, "3 waves: A1, A2, B");
    if (waves.size() == 3) {
      EXPECT(waves[0].size() == 3, "A1 = 3 channel anchors");
      EXPECT(wave_has(waves[0], "fd:dipole_x@f0.00000") &&
             wave_has(waves[0], "fd:dipole_y@f0.00000") &&
             wave_has(waves[0], "fd:dipole_z@f0.00000"),
             "A1 holds the static (anchor) of each channel");
      EXPECT(waves[1].size() == 3, "A2 = 3 dynamics @ P0");
      EXPECT(wave_has(waves[1], "fd:dipole_x@f0.10000"),
             "A2 holds the dynamic frequencies");
      EXPECT(waves[2].size() == 6, "B = all 6 nodes @ 1e-6");
      for (const auto &it : waves[0])
        EXPECT(std::abs(it.thresh - 1e-4) < 1e-12 &&
               it.action == NodeAction::Fresh,
               "A1 items: thresh=1e-4, Fresh");
    }
  }

  // ====== schedule: skips converged, symbolic excluded ====================
  std::printf("=== schedule: skip converged + exclude symbolic ===\n");
  {
    // Resonant: ES bundle + symbolic derived. Symbolic node must never appear.
    ResponsePropertyRequest r;
    r.kind = ResponsePropertyKind::PolarizabilityGradient;
    r.gradient_mode = GradientMode::Resonant;
    r.n_roots = 2;
    r.protocol_thresholds = P;
    auto dag = build_dag(plan_one(r), 0);
    auto waves = schedule(dag, P, empty_meta());
    for (const auto &w : waves)
      EXPECT(!wave_has(w, "dfd:dipole_x@*"),
             "symbolic derived FD never scheduled");
    // Only the ES node is schedulable -> at least one wave with es:tda_n2.
    bool es_seen = false;
    for (const auto &w : waves) es_seen |= wave_has(w, "es:tda_n2");
    EXPECT(es_seen, "ES node is scheduled");

    // Now mark the ES bundle converged at both rungs -> everything Skips.
    json m = empty_meta();
    put_es(m, 1e-4, true);
    put_es(m, 1e-6, true);
    auto waves2 = schedule(dag, P, m);
    EXPECT(waves2.empty(), "fully-converged plan schedules no work");
  }

  std::printf("\n%s  (%d failures)\n", failed ? "FAILED" : "PASSED", failed);
  return failed ? 1 : 0;
}
