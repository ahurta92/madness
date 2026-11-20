#pragma once

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace molresponse_v3 {
struct RootComparison {
  std::optional<std::size_t> previous;
  double phase = 1.0;
  std::string reason = "comparison_unavailable";
};

// Conservative measurement eligibility, not a solver or root-ordering policy.
// A dominant match must beat every competitor in both directions. These
// internal cutoffs are provisional until the final numerical validation; they
// are independent of DALTON import diagnostics and never rotate a function.
struct RootComparisonLimits {
  double minimum_overlap = 0.8;
  double minimum_gap = 0.1;
};

/// Rows are current roots, columns are previous prepared roots. Supply the
/// positive Euclidean inner product and L2 norms over ALL stored X/Y/spin
/// functions. The signed RPA metric remains exclusively in the eigensolver.
/// Only mutual, unambiguous matches qualify comparison history; incomplete
/// matching is intentional and gives unavailable changes, never invented zeros.
inline std::vector<RootComparison> match_prepared_roots(
    const std::vector<std::vector<double>>& overlaps,
    const std::vector<double>& current_norms,
    const std::vector<double>& previous_norms,
    RootComparisonLimits limits = {}) {
  const auto n=current_norms.size(), m=previous_norms.size();
  if (overlaps.size()!=n || !std::isfinite(limits.minimum_overlap) ||
      limits.minimum_overlap<=0.0 || limits.minimum_overlap>1.0 ||
      !std::isfinite(limits.minimum_gap) || limits.minimum_gap<=0.0 || limits.minimum_gap>1.0)
    throw std::invalid_argument("invalid root comparison inputs");
  struct Best {
    std::optional<std::size_t> index;
    double value=-1.0, runner_up=0.0;
    void add(std::size_t candidate,double score) {
      if (score>value) {runner_up=std::max(0.0,value);value=score;index=candidate;}
      else runner_up=std::max(runner_up,score);
    }
  };
  std::vector<Best> rows(n),columns(m);
  std::vector<bool> valid_rows(n,true),valid_columns(m,true);
  std::vector<RootComparison> result(n);
  const auto valid_norm=[](double norm){return std::isfinite(norm) && norm>0.0;};
  for (std::size_t j=0;j<m;++j) valid_columns[j]=valid_norm(previous_norms[j]);
  for (std::size_t i=0;i<n;++i) {
    if (overlaps[i].size()!=m) throw std::invalid_argument("ragged root overlap matrix");
    valid_rows[i]=valid_norm(current_norms[i]);
    for (std::size_t j=0;j<m;++j) {
      if (!valid_norm(current_norms[i]) || !valid_norm(previous_norms[j])) continue;
      const double score=std::abs((overlaps[i][j]/current_norms[i])/previous_norms[j]);
      if (!std::isfinite(score) || score>1.0+1e-8) {
        valid_rows[i]=false;valid_columns[j]=false;continue;
      }
      rows[i].add(j,std::min(1.0,score));
      columns[j].add(i,std::min(1.0,score));
    }
  }
  for (std::size_t i=0;i<n;++i) {
    if (!valid_rows[i]) {result[i].reason="invalid_or_zero_root_norm";continue;}
    if (!rows[i].index) continue;
    const auto j=*rows[i].index;
    if (!valid_columns[j]) {result[i].reason="invalid_previous_root";continue;}
    if (rows[i].value<limits.minimum_overlap) {result[i].reason="weak_root_overlap";continue;}
    if (columns[j].index!=i || rows[i].value-rows[i].runner_up<limits.minimum_gap ||
        columns[j].value-columns[j].runner_up<limits.minimum_gap) {
      result[i].reason="ambiguous_root_correspondence";continue;
    }
    result[i].previous=j;
    result[i].phase=overlaps[i][j]<0.0 ? -1.0 : 1.0;
    result[i].reason.clear();
  }
  return result;
}
} // namespace molresponse_v3
