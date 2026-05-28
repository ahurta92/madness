// ===========================================================================
// PropertyPlanner round-trip and dedupe tests. Pure C++, no MPI / no World.
// ===========================================================================

#include "../PropertyPlanner.hpp"

#include <algorithm>
#include <cstdio>
#include <string>

using namespace molresponse_v3;

namespace {

int failed = 0;

#define EXPECT(cond, label)                                                \
  do {                                                                     \
    if (cond) { std::printf("  [PASS]  %s\n", label); }                    \
    else      { std::printf("  [FAIL]  %s\n", label); ++failed; }          \
  } while (0)

bool has_fd(const ResponsePlan &p, const std::string &pert, double freq,
            double tol = 1e-7) {
  for (const auto &r : p.fd) {
    if (r.pert.description() == pert && std::abs(r.freq - freq) < tol)
      return true;
  }
  return false;
}

bool has_dfd(const ResponsePlan &p, const std::string &pert,
             const std::string &es_root_id) {
  for (const auto &r : p.derived_fd) {
    if (r.pert.description() == pert && r.es_root_id == es_root_id)
      return true;
  }
  return false;
}

bool fd_protocols_eq(const ResponsePlan &p, const std::string &pert,
                     double freq, const std::vector<double> &expected) {
  for (const auto &r : p.fd) {
    if (r.pert.description() == pert && std::abs(r.freq - freq) < 1e-7) {
      return r.protocols == expected;
    }
  }
  return false;
}

} // namespace

int main() {
  const std::vector<double> P = {1e-4, 1e-6};

  // ----- α(ω=0.057) on {x,y,z} -> 3 FD, no ES, no derived ------------------
  std::printf("=== alpha(0.057) on {x,y,z} ===\n");
  {
    PropertyRequest r;
    r.kind = PropertyKind::Alpha;
    r.frequencies = {0.057};
    r.protocol_thresholds = P;
    auto plan = plan_one(r);
    EXPECT(plan.fd.size() == 3,        "3 FD requests");
    EXPECT(plan.es.empty(),            "no ES requests");
    EXPECT(plan.derived_fd.empty(),    "no derived-FD requests");
    EXPECT(has_fd(plan, "dipole_x", 0.057), "dipole_x @ 0.057");
    EXPECT(has_fd(plan, "dipole_y", 0.057), "dipole_y @ 0.057");
    EXPECT(has_fd(plan, "dipole_z", 0.057), "dipole_z @ 0.057");
    EXPECT(plan.fd[0].protocols == P,  "protocol ramp propagated to FD");
  }

  // ----- α at 3 freqs on {x,z} -> 6 FD --------------------------------------
  std::printf("=== alpha({0, 0.04, 0.057}) on {x,z} ===\n");
  {
    PropertyRequest r;
    r.kind = PropertyKind::Alpha;
    r.frequencies = {0.0, 0.04, 0.057};
    r.axes = {'x','z'};
    r.protocol_thresholds = P;
    auto plan = plan_one(r);
    EXPECT(plan.fd.size() == 6,                       "6 FD requests");
    EXPECT(has_fd(plan, "dipole_x", 0.0)   && has_fd(plan, "dipole_z", 0.0),   "static pair");
    EXPECT(has_fd(plan, "dipole_x", 0.04)  && has_fd(plan, "dipole_z", 0.04),  "ω=0.04 pair");
    EXPECT(has_fd(plan, "dipole_x", 0.057) && has_fd(plan, "dipole_z", 0.057), "ω=0.057 pair");
  }

  // ----- β-SHG(ω=0.057) on {x,y,z} -> 6 FD: 3@ω + 3@2ω ----------------------
  std::printf("=== beta-SHG(0.057) on {x,y,z} ===\n");
  {
    PropertyRequest r;
    r.kind = PropertyKind::Beta;
    r.beta_process = BetaProcess::SHG;
    r.frequencies = {0.057};
    r.protocol_thresholds = P;
    auto plan = plan_one(r);
    EXPECT(plan.fd.size() == 6,                       "6 FD (3 @ ω + 3 @ 2ω)");
    EXPECT(has_fd(plan, "dipole_x", 0.057) && has_fd(plan, "dipole_x", 0.114), "x at ω and 2ω");
    EXPECT(has_fd(plan, "dipole_y", 0.057) && has_fd(plan, "dipole_y", 0.114), "y at ω and 2ω");
    EXPECT(has_fd(plan, "dipole_z", 0.057) && has_fd(plan, "dipole_z", 0.114), "z at ω and 2ω");
  }

  // ----- β-OR(ω=0.057) -> 6 FD: 3@ω + 3@0 -----------------------------------
  std::printf("=== beta-OR(0.057) on {x,y,z} ===\n");
  {
    PropertyRequest r;
    r.kind = PropertyKind::Beta;
    r.beta_process = BetaProcess::OR;
    r.frequencies = {0.057};
    r.protocol_thresholds = P;
    auto plan = plan_one(r);
    EXPECT(plan.fd.size() == 6,                       "6 FD (3 @ 0 + 3 @ ω)");
    EXPECT(has_fd(plan, "dipole_x", 0.0)   && has_fd(plan, "dipole_x", 0.057), "x at 0 and ω");
  }

  // ----- α(ω) + β-SHG(ω) merged -> 6 unique FD (α's ω-set ⊂ β's ω-set) ------
  std::printf("=== merge: alpha(0.057) + beta-SHG(0.057) ===\n");
  {
    PropertyRequest a;
    a.kind = PropertyKind::Alpha;
    a.frequencies = {0.057};
    a.protocol_thresholds = P;
    PropertyRequest b;
    b.kind = PropertyKind::Beta;
    b.beta_process = BetaProcess::SHG;
    b.frequencies = {0.057};
    b.protocol_thresholds = P;
    auto plan = merge_plans({plan_one(a), plan_one(b)});
    EXPECT(plan.fd.size() == 6, "dedupe: 6 unique FD (not 9 = 3 alpha + 6 beta)");
  }

  // ----- Protocol union on dedupe ------------------------------------------
  std::printf("=== merge: protocol-set union on same (pert,freq) ===\n");
  {
    PropertyRequest a;
    a.kind = PropertyKind::Alpha;
    a.frequencies = {0.057};
    a.axes = {'x'};
    a.protocol_thresholds = {1e-4};
    PropertyRequest b = a;
    b.protocol_thresholds = {1e-6};
    auto plan = merge_plans({plan_one(a), plan_one(b)});
    EXPECT(plan.fd.size() == 1,                                   "still 1 FD after merge");
    EXPECT(fd_protocols_eq(plan, "dipole_x", 0.057, {1e-4, 1e-6}), "union(coarse→fine) = [1e-4, 1e-6]");
  }

  // ----- ResonantRaman(n_roots=4) -> 1 ES + 3 derived FD with "*" -----------
  std::printf("=== resonant Raman (n_roots=4) ===\n");
  {
    PropertyRequest r;
    r.kind = PropertyKind::ResonantRaman;
    r.n_roots = 4;
    r.protocol_thresholds = P;
    auto plan = plan_one(r);
    EXPECT(plan.fd.empty(),                            "no direct FD");
    EXPECT(plan.es.size() == 1 && plan.es[0].n_roots == 4 && plan.es[0].tda,
                                                       "1 ES request (TDA, 4 roots)");
    EXPECT(plan.derived_fd.size() == 3,                "3 derived FD (one per axis)");
    EXPECT(has_dfd(plan, "dipole_x", "*") &&
           has_dfd(plan, "dipole_y", "*") &&
           has_dfd(plan, "dipole_z", "*"),             "all three axes, root_id='*'");
  }

  // ----- α + Raman together: no spurious dedupe between FD and derived-FD ---
  std::printf("=== merge: alpha + Raman keep separate ===\n");
  {
    PropertyRequest a;
    a.kind = PropertyKind::Alpha;
    a.frequencies = {0.057};
    a.protocol_thresholds = P;
    PropertyRequest r;
    r.kind = PropertyKind::ResonantRaman;
    r.n_roots = 2;
    r.protocol_thresholds = P;
    auto plan = merge_plans({plan_one(a), plan_one(r)});
    EXPECT(plan.fd.size() == 3,         "3 alpha FD");
    EXPECT(plan.es.size() == 1,         "1 ES bundle");
    EXPECT(plan.derived_fd.size() == 3, "3 derived FD");
  }

  std::printf("\n%s: %d failure(s)\n",
              failed == 0 ? "ALL PASS" : "FAILED", failed);
  return failed == 0 ? 0 : 1;
}
