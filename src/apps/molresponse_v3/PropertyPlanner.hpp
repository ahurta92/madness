#ifndef MOLRESPONSE_V3_PROPERTY_PLANNER_HPP
#define MOLRESPONSE_V3_PROPERTY_PLANNER_HPP

// ===========================================================================
// PropertyPlanner — turn a user-facing PropertyRequest into the list of
// response states (FD + ES) the calc manager needs to solve. Pure planning:
// no World, no I/O, no solver dispatch. Output is a ResponsePlan; the calc
// manager (separate concern) executes it and the results flow back through
// response_metadata.json (doc 13).
//
// Three properties supported:
//
//   PropertyKind::Alpha           α(ω) — needs FD states (dipole_a, ω) for
//                                 each requested axis a and frequency ω.
//   PropertyKind::Beta            β(ω₁,ω₂,ω₃) — `BetaProcess` picks the
//                                 triplet from the user-provided driver ω;
//                                 unique |ω|s become FD states for every axis.
//                                   Static : (0,0,0)
//                                   SHG    : (ω, ω, -2ω)   → {ω, 2ω}
//                                   OR     : (-ω, ω, 0)   → {ω, 0}
//                                   EOPE   : (-ω, ω, 0)   → {ω, 0}
//   PropertyKind::ResonantRaman   Needs an ES bundle (n_roots) AND FD states
//                                 at each excitation energy ω_root. The
//                                 latter are symbolic in the plan — see
//                                 DerivedFDRequest below.
//
// Dependency model: `derived_fd` entries are SYMBOLIC FD requests whose
// frequency comes from a future ES root. The es_root_id == "*" sentinel
// means "one per converged root in the matching ES bundle". The calc
// manager expands "*" after the ES bundle converges; the planner stays
// pure (no future-eigenvalue dependency in the planning step itself).
// ===========================================================================

#include "Perturbations.hpp"            // Perturbation, Perturbation::dipole

#include <algorithm>
#include <cstdio>
#include <functional>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace molresponse_v3 {

enum class PropertyKind { Alpha, Beta, ResonantRaman };

/// β-process frequency convention (driver ω → triplet).
enum class BetaProcess { Static, SHG, OR, EOPE };

/// User-facing request — one logical property at one set of frequencies.
struct PropertyRequest {
  PropertyKind        kind;
  /// α/β: driver ω(s). For α every ω becomes the FD freq; for β the
  /// `beta_process` triplet is built from each ω. Raman: ignored — the ES
  /// eigenvalues become the FD frequencies after the bundle converges.
  std::vector<double> frequencies;
  /// Cartesian directions: 'x','y','z' (case-insensitive). For α each
  /// becomes a single FD state per ω; for β the same axes are emitted
  /// at every triplet frequency; for Raman a derived FD is emitted per
  /// axis per future root.
  std::vector<char>   axes         = {'x','y','z'};
  BetaProcess         beta_process = BetaProcess::SHG;
  int                 n_roots      = 0;           // Raman only
  std::vector<double> protocol_thresholds;        // ramp; all states get this
};

/// One concrete FD state the calc manager must solve.
struct FDRequest {
  Perturbation        pert;
  double              freq;
  std::vector<double> protocols;
};

/// One ES bundle the calc manager must solve.
struct ESRequest {
  bool                tda = true;   // false → Full/RPA
  int                 n_roots;
  std::vector<double> protocols;
};

/// Symbolic FD request — `es_root_id` references a future ES output.
/// "*" expands to one record per converged root post-ES.
struct DerivedFDRequest {
  Perturbation        pert;
  std::string         es_root_id;   // "es_root_0001" or "*"
  std::vector<double> protocols;
};

struct ResponsePlan {
  std::vector<FDRequest>        fd;
  std::vector<ESRequest>        es;
  std::vector<DerivedFDRequest> derived_fd;
};

// ---------------------------------------------------------------------------
namespace detail_planner {

inline int axis_index(char a) {
  switch (a) {
    case 'x': case 'X': return 0;
    case 'y': case 'Y': return 1;
    case 'z': case 'Z': return 2;
    default: return -1;
  }
}

/// Unique |ω| set the calc manager needs to solve to assemble β at one
/// driver ω. β components mix these FD states (and possibly their
/// complex conjugates) — that's the assembler's concern.
inline std::vector<double> beta_freqs(BetaProcess proc, double omega) {
  std::set<double> uniq;
  switch (proc) {
    case BetaProcess::Static: uniq.insert(0.0); break;
    case BetaProcess::SHG:    uniq.insert(omega); uniq.insert(2.0 * omega); break;
    case BetaProcess::OR:
    case BetaProcess::EOPE:   uniq.insert(omega); uniq.insert(0.0); break;
  }
  return {uniq.begin(), uniq.end()};
}

/// Sorted (coarse → fine) deduplicated union of two protocol lists.
/// Used by merge_plans to combine ramps when the same state is requested
/// at different precisions across multiple PropertyRequests.
inline std::vector<double> union_protocols(std::vector<double> a,
                                            const std::vector<double> &b) {
  for (double t : b) a.push_back(t);
  std::sort(a.begin(), a.end(), std::greater<double>{});
  a.erase(std::unique(a.begin(), a.end()), a.end());
  return a;
}

/// Stable freq key (matches the doc-13 archive freq_key precision)
/// used for FD dedupe across nominally-equal floats.
inline std::string fd_freq_key(double f) {
  char buf[16];
  std::snprintf(buf, sizeof buf, "%.5f", f);
  return buf;
}

} // namespace detail_planner

// ---------------------------------------------------------------------------
/// Merge & dedupe across plans. Identical FD requests collapse with their
/// protocol sets unioned (so "α at ω with ramp [1e-4]" + "β-SHG at ω with
/// ramp [1e-6]" yields one FD at ω with ramp [1e-4, 1e-6]). Same for ES on
/// (tda, n_roots) and derived FD on (pert.description(), es_root_id).
inline ResponsePlan merge_plans(const std::vector<ResponsePlan> &plans) {
  ResponsePlan out;

  auto fd_key = [](const FDRequest &r) {
    return r.pert.description() + "@" + detail_planner::fd_freq_key(r.freq);
  };
  auto es_key = [](const ESRequest &r) {
    return std::string(r.tda ? "tda" : "full") + "_" + std::to_string(r.n_roots);
  };
  auto dfd_key = [](const DerivedFDRequest &r) {
    return r.pert.description() + "@" + r.es_root_id;
  };

  // Simple O(N) lookup in a parallel keys vector — counts are small
  // (handful of states per plan).
  auto upsert_fd = [&](const FDRequest &r) {
    static thread_local std::vector<std::string> keys;
    (void)keys;  // suppress unused warning if header included w/o use
    const auto k = fd_key(r);
    for (size_t i = 0; i < out.fd.size(); ++i) {
      if (fd_key(out.fd[i]) == k) {
        out.fd[i].protocols =
            detail_planner::union_protocols(out.fd[i].protocols, r.protocols);
        return;
      }
    }
    out.fd.push_back(r);
  };
  auto upsert_es = [&](const ESRequest &r) {
    const auto k = es_key(r);
    for (size_t i = 0; i < out.es.size(); ++i) {
      if (es_key(out.es[i]) == k) {
        out.es[i].protocols =
            detail_planner::union_protocols(out.es[i].protocols, r.protocols);
        return;
      }
    }
    out.es.push_back(r);
  };
  auto upsert_dfd = [&](const DerivedFDRequest &r) {
    const auto k = dfd_key(r);
    for (size_t i = 0; i < out.derived_fd.size(); ++i) {
      if (dfd_key(out.derived_fd[i]) == k) {
        out.derived_fd[i].protocols = detail_planner::union_protocols(
            out.derived_fd[i].protocols, r.protocols);
        return;
      }
    }
    out.derived_fd.push_back(r);
  };

  for (const auto &p : plans) {
    for (const auto &r : p.fd)         upsert_fd(r);
    for (const auto &r : p.es)         upsert_es(r);
    for (const auto &r : p.derived_fd) upsert_dfd(r);
  }
  return out;
}

inline ResponsePlan merge_plans(ResponsePlan a, ResponsePlan b) {
  return merge_plans(std::vector<ResponsePlan>{std::move(a), std::move(b)});
}

// ---------------------------------------------------------------------------
/// Plan one PropertyRequest. Output is self-deduped (merge_plans pass).
inline ResponsePlan plan_one(const PropertyRequest &req) {
  ResponsePlan plan;
  using detail_planner::axis_index;
  using detail_planner::beta_freqs;

  switch (req.kind) {
    case PropertyKind::Alpha: {
      for (double w : req.frequencies) {
        for (char ax : req.axes) {
          int i = axis_index(ax);
          if (i < 0) continue;
          plan.fd.push_back({Perturbation::dipole(i), w,
                             req.protocol_thresholds});
        }
      }
      break;
    }
    case PropertyKind::Beta: {
      for (double w : req.frequencies) {
        for (double bw : beta_freqs(req.beta_process, w)) {
          for (char ax : req.axes) {
            int i = axis_index(ax);
            if (i < 0) continue;
            plan.fd.push_back({Perturbation::dipole(i), bw,
                               req.protocol_thresholds});
          }
        }
      }
      break;
    }
    case PropertyKind::ResonantRaman: {
      plan.es.push_back({/*tda=*/true, req.n_roots, req.protocol_thresholds});
      for (char ax : req.axes) {
        int i = axis_index(ax);
        if (i < 0) continue;
        plan.derived_fd.push_back({Perturbation::dipole(i),
                                   /*es_root_id=*/"*",
                                   req.protocol_thresholds});
      }
      break;
    }
  }
  return merge_plans(std::vector<ResponsePlan>{std::move(plan)});
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_PROPERTY_PLANNER_HPP
