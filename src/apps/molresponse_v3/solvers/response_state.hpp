#ifndef MOLRESPONSE_V3_SOLVERS_RESPONSE_STATE_HPP
#define MOLRESPONSE_V3_SOLVERS_RESPONSE_STATE_HPP

// =========================================================================
// Storage shapes — two templates, each specialized on Shell.
//
//   ResponseStateX<ClosedShell>   : { x_alpha }
//   ResponseStateX<OpenShell>     : { x_alpha, x_beta }
//   ResponseStateXY<ClosedShell>  : { x_alpha, y_alpha }
//   ResponseStateXY<OpenShell>    : { x_alpha, y_alpha, x_beta, y_beta }
//
// Four distinct types via two structures × two shells. The TYPE
// determines exactly which fields are accessible — no empty-vector
// guard, no "is this populated" check; a closed-shell kernel cannot
// even spell `state.x_beta` because that member does not exist on
// `ResponseStateX<ClosedShell>`.
//
// Save / load: there are two LOGICAL bodies (X-only vs X-and-Y),
// each instantiated for both shells. Four save() / load() function
// definitions in this file, but they share the same algorithm
// pattern per shape.
// =========================================================================

#include "../kernels/tags.hpp"

#include <madness/mra/mra.h>
#include <madness/world/parallel_archive.h>

#include <string>
#include <vector>

namespace molresponse_v3 {

using namespace madness;

// ---- ResponseStateX<ClosedShell> ----------------------------------------
template <>
struct ResponseStateX<ClosedShell> {
  std::vector<real_function_3d> x_alpha;

  std::size_t num_alpha() const { return x_alpha.size(); }

  // ---- flatten / from_flat / flat_size ----
  // Pack all component vecfuncTs into a single flat vector (KAIN
  // operates on flat vectors). The flatten/from_flat round-trip
  // preserves all Storage data. Order is stable per Storage type so
  // KAIN's stored iterates index the same components iter-to-iter.
  std::size_t flat_size() const { return x_alpha.size(); }
  std::vector<real_function_3d> flatten() const { return x_alpha; }
  void from_flat(const std::vector<real_function_3d>& v) {
    x_alpha = v;
  }

  void save(World &world, const std::string &filename) const {
    archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(
        world, filename.c_str(), 1);
    const std::size_t na = x_alpha.size();
    ar & na;
    for (const auto &f : x_alpha) ar & f;
  }

  static ResponseStateX load(World &world, const std::string &filename) {
    archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(
        world, filename.c_str(), 1);
    ResponseStateX s;
    std::size_t na;
    ar & na;
    s.x_alpha.resize(na);
    for (auto &f : s.x_alpha) ar & f;
    return s;
  }
};

// ---- ResponseStateX<OpenShell> ------------------------------------------
template <>
struct ResponseStateX<OpenShell> {
  std::vector<real_function_3d> x_alpha;
  std::vector<real_function_3d> x_beta;

  std::size_t num_alpha() const { return x_alpha.size(); }
  std::size_t num_beta()  const { return x_beta.size(); }

  std::size_t flat_size() const { return x_alpha.size() + x_beta.size(); }
  std::vector<real_function_3d> flatten() const {
    auto v = x_alpha;
    v.insert(v.end(), x_beta.begin(), x_beta.end());
    return v;
  }
  void from_flat(const std::vector<real_function_3d>& v) {
    const std::size_t na = x_alpha.size();
    const std::size_t nb = x_beta.size();
    x_alpha = std::vector<real_function_3d>(
        v.begin(), v.begin() + na);
    x_beta  = std::vector<real_function_3d>(
        v.begin() + na, v.begin() + na + nb);
  }

  void save(World &world, const std::string &filename) const {
    archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(
        world, filename.c_str(), 1);
    const std::size_t na = x_alpha.size();
    const std::size_t nb = x_beta.size();
    ar & na & nb;
    for (const auto &f : x_alpha) ar & f;
    for (const auto &f : x_beta)  ar & f;
  }

  static ResponseStateX load(World &world, const std::string &filename) {
    archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(
        world, filename.c_str(), 1);
    ResponseStateX s;
    std::size_t na, nb;
    ar & na & nb;
    s.x_alpha.resize(na);
    s.x_beta.resize(nb);
    for (auto &f : s.x_alpha) ar & f;
    for (auto &f : s.x_beta)  ar & f;
    return s;
  }
};

// ---- ResponseStateXY<ClosedShell> ---------------------------------------
template <>
struct ResponseStateXY<ClosedShell> {
  std::vector<real_function_3d> x_alpha;
  std::vector<real_function_3d> y_alpha;

  std::size_t num_alpha() const { return x_alpha.size(); }

  std::size_t flat_size() const { return x_alpha.size() + y_alpha.size(); }
  std::vector<real_function_3d> flatten() const {
    auto v = x_alpha;
    v.insert(v.end(), y_alpha.begin(), y_alpha.end());
    return v;
  }
  void from_flat(const std::vector<real_function_3d>& v) {
    const std::size_t na = x_alpha.size();
    x_alpha = std::vector<real_function_3d>(
        v.begin(), v.begin() + na);
    y_alpha = std::vector<real_function_3d>(
        v.begin() + na, v.begin() + 2 * na);
  }

  void save(World &world, const std::string &filename) const {
    archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(
        world, filename.c_str(), 1);
    const std::size_t na = x_alpha.size();
    ar & na;
    for (const auto &f : x_alpha) ar & f;
    for (const auto &f : y_alpha) ar & f;
  }

  static ResponseStateXY load(World &world, const std::string &filename) {
    archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(
        world, filename.c_str(), 1);
    ResponseStateXY s;
    std::size_t na;
    ar & na;
    s.x_alpha.resize(na);
    s.y_alpha.resize(na);
    for (auto &f : s.x_alpha) ar & f;
    for (auto &f : s.y_alpha) ar & f;
    return s;
  }
};

// ---- ResponseStateXY<OpenShell> -----------------------------------------
template <>
struct ResponseStateXY<OpenShell> {
  std::vector<real_function_3d> x_alpha;
  std::vector<real_function_3d> y_alpha;
  std::vector<real_function_3d> x_beta;
  std::vector<real_function_3d> y_beta;

  std::size_t num_alpha() const { return x_alpha.size(); }
  std::size_t num_beta()  const { return x_beta.size(); }

  std::size_t flat_size() const {
    return 2 * (x_alpha.size() + x_beta.size());
  }
  std::vector<real_function_3d> flatten() const {
    auto v = x_alpha;
    v.insert(v.end(), y_alpha.begin(), y_alpha.end());
    v.insert(v.end(), x_beta.begin(),  x_beta.end());
    v.insert(v.end(), y_beta.begin(),  y_beta.end());
    return v;
  }
  void from_flat(const std::vector<real_function_3d>& v) {
    const std::size_t na = x_alpha.size();
    const std::size_t nb = x_beta.size();
    x_alpha = std::vector<real_function_3d>(
        v.begin(),                  v.begin() + na);
    y_alpha = std::vector<real_function_3d>(
        v.begin() + na,             v.begin() + 2 * na);
    x_beta  = std::vector<real_function_3d>(
        v.begin() + 2 * na,         v.begin() + 2 * na + nb);
    y_beta  = std::vector<real_function_3d>(
        v.begin() + 2 * na + nb,    v.begin() + 2 * na + 2 * nb);
  }

  void save(World &world, const std::string &filename) const {
    archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(
        world, filename.c_str(), 1);
    const std::size_t na = x_alpha.size();
    const std::size_t nb = x_beta.size();
    ar & na & nb;
    for (const auto &f : x_alpha) ar & f;
    for (const auto &f : y_alpha) ar & f;
    for (const auto &f : x_beta)  ar & f;
    for (const auto &f : y_beta)  ar & f;
  }

  static ResponseStateXY load(World &world, const std::string &filename) {
    archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(
        world, filename.c_str(), 1);
    ResponseStateXY s;
    std::size_t na, nb;
    ar & na & nb;
    s.x_alpha.resize(na);
    s.y_alpha.resize(na);
    s.x_beta.resize(nb);
    s.y_beta.resize(nb);
    for (auto &f : s.x_alpha) ar & f;
    for (auto &f : s.y_alpha) ar & f;
    for (auto &f : s.x_beta)  ar & f;
    for (auto &f : s.y_beta)  ar & f;
    return s;
  }
};

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_RESPONSE_STATE_HPP
