// ===========================================================================
// ResponsePropertyPlanner round-trip and dedupe tests. Pure C++, no MPI.
// ===========================================================================

#include "../ResponsePropertyPlanner.hpp"

#include <algorithm>
#include <cmath>
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

bool has_fd(const std::vector<FDRequest> &v, const std::string &pert,
            double freq, double tol = 1e-7) {
  for (const auto &r : v)
    if (r.pert.description() == pert && std::abs(r.freq - freq) < tol)
      return true;
  return false;
}

bool has_dfd(const ResponsePlan &p, const std::string &pert,
             const std::string &es_root_id) {
  for (const auto &r : p.derived_fd)
    if (r.pert.description() == pert && r.es_root_id == es_root_id)
      return true;
  return false;
}

bool fd_protocols_eq(const ResponsePlan &p, const std::string &pert,
                     double freq, const std::vector<double> &expected) {
  for (const auto &r : p.fd)
    if (r.pert.description() == pert && std::abs(r.freq - freq) < 1e-7)
      return r.protocols == expected;
  return false;
}

} // namespace

int main() {
  const std::vector<double> P = {1e-4, 1e-6};

  // ----- Polarizability(ω=0.057) on {x,y,z} -> 3 FD ------------------------
  std::printf("=== Polarizability(0.057) on {x,y,z} ===\n");
  {
    ResponsePropertyRequest r;
    r.kind = ResponsePropertyKind::Polarizability;
    r.frequencies = {0.057};
    r.protocol_thresholds = P;
    auto plan = plan_one(r);
    EXPECT(plan.fd.size() == 3,     "3 FD requests");
    EXPECT(plan.es.empty(),         "no ES");
    EXPECT(plan.derived_fd.empty(), "no derived FD");
    EXPECT(plan.nuclear_fd.empty(), "no nuclear FD");
    EXPECT(has_fd(plan.fd, "dipole_x", 0.057) &&
           has_fd(plan.fd, "dipole_y", 0.057) &&
           has_fd(plan.fd, "dipole_z", 0.057), "x/y/z @ 0.057");
    EXPECT(plan.fd[0].protocols == P, "protocol ramp propagated");
  }

  // ----- Hyperpolarizability SHG(ω=0.057) -> 6 FD: 3@ω + 3@2ω --------------
  std::printf("=== Hyperpolarizability-SHG(0.057) ===\n");
  {
    ResponsePropertyRequest r;
    r.kind = ResponsePropertyKind::Hyperpolarizability;
    r.beta_process = BetaProcess::SHG;
    r.frequencies = {0.057};
    r.protocol_thresholds = P;
    auto plan = plan_one(r);
    EXPECT(plan.fd.size() == 6, "6 FD (3 @ ω + 3 @ 2ω)");
    EXPECT(has_fd(plan.fd, "dipole_z", 0.057) && has_fd(plan.fd, "dipole_z", 0.114),
           "z at ω and 2ω");
  }

  // ----- Hyperpolarizability OR(ω=0.057) -> 6 FD: 3@0 + 3@ω ----------------
  std::printf("=== Hyperpolarizability-OR(0.057) ===\n");
  {
    ResponsePropertyRequest r;
    r.kind = ResponsePropertyKind::Hyperpolarizability;
    r.beta_process = BetaProcess::OR;
    r.frequencies = {0.057};
    r.protocol_thresholds = P;
    auto plan = plan_one(r);
    EXPECT(plan.fd.size() == 6, "6 FD (3 @ 0 + 3 @ ω)");
    EXPECT(has_fd(plan.fd, "dipole_x", 0.0) && has_fd(plan.fd, "dipole_x", 0.057),
           "x at 0 and ω");
  }

  // ----- PolarizabilityGradient Nuclear (vibrational Raman) ----------------
  std::printf("=== PolarizabilityGradient Nuclear(ω=0.057) ===\n");
  {
    ResponsePropertyRequest r;
    r.kind = ResponsePropertyKind::PolarizabilityGradient;
    r.gradient_mode = GradientMode::Nuclear;
    r.frequencies = {0.057};
    r.protocol_thresholds = P;
    auto plan = plan_one(r);
    EXPECT(plan.fd.size() == 3,           "3 dipole FD @ ω");
    EXPECT(has_fd(plan.fd, "dipole_x", 0.057), "dipole_x @ ω present");
    EXPECT(plan.nuclear_fd.size() == 3,   "3 nuclear-disp FD (all-atoms, xyz)");
    EXPECT(has_fd(plan.nuclear_fd, "nuc_*_x", 0.0) &&
           has_fd(plan.nuclear_fd, "nuc_*_y", 0.0) &&
           has_fd(plan.nuclear_fd, "nuc_*_z", 0.0), "nuc_*_{x,y,z} @ 0");
    EXPECT(plan.es.empty() && plan.derived_fd.empty(), "no ES / no derived FD");
    // Full set: one dipole x nuclear source per (optical axis, displacement
    // axis), carrying the all-atoms sentinel for build_dag to expand.
    EXPECT(plan.vbc.size() == 9, "9 sentinel VBC (3 optical x 3 displacement axes)");
    bool all_sentinel = !plan.vbc.empty();
    for (const auto &v : plan.vbc)
      all_sentinel = all_sentinel && v.pert_c.is_all_atoms() &&
                     v.pert_b.kind == Perturbation::Kind::Dipole &&
                     std::abs(v.freq_b - 0.057) < 1e-9 && v.freq_c == 0.0;
    EXPECT(all_sentinel, "every VBC is dipole@ω x nuc_*@0");
  }

  // ----- Single coordinate: concrete nuclear leg, 3 VBC, no sentinel -------
  std::printf("=== PolarizabilityGradient Nuclear single (atom 1, z) ===\n");
  {
    ResponsePropertyRequest r;
    r.kind = ResponsePropertyKind::PolarizabilityGradient;
    r.gradient_mode = GradientMode::Nuclear;
    r.frequencies = {0.0};
    r.raman_nuc_atom = 1;
    r.raman_nuc_axis = 2;
    r.protocol_thresholds = P;
    auto plan = plan_one(r);
    EXPECT(plan.nuclear_fd.size() == 1 && has_fd(plan.nuclear_fd, "nuc_1_z", 0.0),
           "1 concrete nuclear leg nuc_1_z");
    EXPECT(plan.vbc.size() == 3, "3 VBC (one per optical axis)");
    bool none_sentinel = true;
    for (const auto &v : plan.vbc) none_sentinel = none_sentinel && !v.pert_c.is_all_atoms();
    EXPECT(none_sentinel, "no sentinel in the single-coordinate plan");
  }

  // ----- Nuclear axes are full 3N even if optical axes subset --------------
  std::printf("=== Nuclear gradient: full 3N displacement even w/ axes={z} ===\n");
  {
    ResponsePropertyRequest r;
    r.kind = ResponsePropertyKind::PolarizabilityGradient;
    r.gradient_mode = GradientMode::Nuclear;
    r.axes = {'z'};
    r.frequencies = {0.057};
    r.protocol_thresholds = P;
    auto plan = plan_one(r);
    EXPECT(plan.fd.size() == 1,         "1 dipole FD (axes={z})");
    EXPECT(plan.nuclear_fd.size() == 3, "still 3 nuclear-disp axes (need full 3N)");
    EXPECT(plan.vbc.size() == 3,        "3 sentinel VBC (1 optical x 3 displacement)");
  }

  // ----- PolarizabilityGradient Resonant: roots only (default) -------------
  std::printf("=== PolarizabilityGradient Resonant(n_roots=4) roots-only ===\n");
  {
    ResponsePropertyRequest r;
    r.kind = ResponsePropertyKind::PolarizabilityGradient;
    r.gradient_mode = GradientMode::Resonant;
    r.n_roots = 4;
    r.protocol_thresholds = P;
    auto plan = plan_one(r);
    EXPECT(plan.fd.empty() && plan.nuclear_fd.empty(),  "no direct FD");
    EXPECT(plan.es.size() == 1 && plan.es[0].n_roots == 4 && plan.es[0].tda,
                                                        "1 ES (TDA, 4 roots)");
    EXPECT(plan.derived_fd.empty(),  "roots only: no derived FD without tpa");
  }

  // ----- PolarizabilityGradient Resonant + tpa: adds the 2PA FD legs -------
  std::printf("=== PolarizabilityGradient Resonant(n_roots=4) + tpa ===\n");
  {
    ResponsePropertyRequest r;
    r.kind = ResponsePropertyKind::PolarizabilityGradient;
    r.gradient_mode = GradientMode::Resonant;
    r.n_roots = 4;
    r.tpa = true;
    r.protocol_thresholds = P;
    auto plan = plan_one(r);
    EXPECT(plan.es.size() == 1 && plan.es[0].n_roots == 4 && plan.es[0].tda,
                                                        "1 ES (TDA, 4 roots)");
    EXPECT(plan.derived_fd.size() == 3,                 "3 derived FD (per axis)");
    EXPECT(has_dfd(plan, "dipole_x", "*") &&
           has_dfd(plan, "dipole_y", "*") &&
           has_dfd(plan, "dipole_z", "*"),              "all axes, root_id='*'");
  }

  // ----- dedupe: Polarizability(ω) + Hyperpol-SHG(ω) -> 6 unique FD --------
  std::printf("=== merge: Polarizability + Hyperpol-SHG (dedupe) ===\n");
  {
    ResponsePropertyRequest a;
    a.kind = ResponsePropertyKind::Polarizability;
    a.frequencies = {0.057};
    a.protocol_thresholds = P;
    ResponsePropertyRequest b;
    b.kind = ResponsePropertyKind::Hyperpolarizability;
    b.beta_process = BetaProcess::SHG;
    b.frequencies = {0.057};
    b.protocol_thresholds = P;
    auto plan = merge_plans({plan_one(a), plan_one(b)});
    EXPECT(plan.fd.size() == 6, "6 unique FD (not 9)");
  }

  // ----- protocol-set union on dedupe --------------------------------------
  std::printf("=== merge: protocol-set union ===\n");
  {
    ResponsePropertyRequest a;
    a.kind = ResponsePropertyKind::Polarizability;
    a.frequencies = {0.057};
    a.axes = {'x'};
    a.protocol_thresholds = {1e-4};
    ResponsePropertyRequest b = a;
    b.protocol_thresholds = {1e-6};
    auto plan = merge_plans({plan_one(a), plan_one(b)});
    EXPECT(plan.fd.size() == 1, "still 1 FD after merge");
    EXPECT(fd_protocols_eq(plan, "dipole_x", 0.057, {1e-4, 1e-6}),
           "union(coarse→fine) = [1e-4, 1e-6]");
  }

  // ----- vibrational Raman + α merge: nuclear_fd stays separate ------------
  std::printf("=== merge: alpha + vibrational Raman ===\n");
  {
    ResponsePropertyRequest a;
    a.kind = ResponsePropertyKind::Polarizability;
    a.frequencies = {0.057};
    a.protocol_thresholds = P;
    ResponsePropertyRequest g;
    g.kind = ResponsePropertyKind::PolarizabilityGradient;
    g.gradient_mode = GradientMode::Nuclear;
    g.frequencies = {0.057};
    g.protocol_thresholds = P;
    auto plan = merge_plans({plan_one(a), plan_one(g)});
    EXPECT(plan.fd.size() == 3,         "dipole FD deduped to 3 (shared @ 0.057)");
    EXPECT(plan.nuclear_fd.size() == 3, "3 nuclear-disp FD");
  }

  std::printf("\n%s: %d failure(s)\n",
              failed == 0 ? "ALL PASS" : "FAILED", failed);
  return failed == 0 ? 0 : 1;
}
