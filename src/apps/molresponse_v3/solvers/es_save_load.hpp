#ifndef MOLRESPONSE_V3_SOLVERS_ES_SAVE_LOAD_HPP
#define MOLRESPONSE_V3_SOLVERS_ES_SAVE_LOAD_HPP

// =========================================================================
// Solve / save / restart for ESSolver<Type, Shell>.
//
// On disk, an "ES roots directory" is laid out as:
//
//   my_dir/
//   ├── roots.json    ← sidecar: type, shell, n_roots, k, thresh, per-root
//   │                   metadata (slot, omega, residuals, file name)
//   ├── root_0        ← MADNESS parallel archive (one per slot)
//   ├── root_1
//   └── ...
//
// Per-root binary content uses the existing `ResponseStateX<Shell>::save`
// / `::load` machinery in response_state.hpp — this header only adds the
// sidecar + a directory-level wrapper.
//
// `roots.json` is the *index*: it tells the loader (a) what type/shell to
// expect, (b) which file holds which slot, (c) the ω / residual recorded
// at save time, (d) the per-root STABLE identity (stable_index / root_id /
// display_name) plus the bundle's slot_permutation (slot -> stable_index).
// ω is recorded for sanity (caller compares against a reference) — load
// does not permute by ω. The slot_permutation is what lets a restart match
// roots by stable identity rather than by raw array slot (doc 03, Inc 9).
//
// On load: `iter` is reset to 0 (fresh count for the resumed solver),
// `rho_alpha_prev` is left empty (ESSolver::step's first-iter check is
// already gated by `iter <= 1`). The KAIN buffer of the resumed solver
// starts empty too — that's the right call for a restart test.
//
// Cross-protocol restart: the bundle's native k/thresh is recorded but
// NOT enforced. The caller's protocol-step `prepare` hook re-projects
// each root to the active k/thresh via `madness::project`, same path the
// per-protocol stepper already uses.
//
// Collective discipline: parallel archive I/O is collective on `world`
// and must be called by every rank. The JSON sidecar read/write is
// rank-0-only (filesystem op); rank 0 broadcasts the parsed metadata.
// =========================================================================

#include "es_root_identity.hpp"
#include "es_solver.hpp"
#include "response_state.hpp"
#include "../ResponseProtocol.hpp"   // protocol_key(thresh, k)
#include "../kernels/tags.hpp"

#include <madness/external/nlohmann_json/json.hpp>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace molresponse_v3 {

namespace detail_save_load {

template <typename Type>
inline const char *type_tag() {
  if constexpr (std::is_same_v<Type, TDA>)  return "tda";
  else if constexpr (std::is_same_v<Type, Full>) return "full";
  else return "unknown";
}

template <typename Shell>
inline const char *shell_tag() {
  if constexpr (std::is_same_v<Shell, ClosedShell>) return "closed_shell";
  else if constexpr (std::is_same_v<Shell, OpenShell>) return "open_shell";
  else return "unknown";
}

inline std::string root_file(int slot) {
  return "root_" + std::to_string(slot);
}

} // namespace detail_save_load

/// Save an ESSolver state to a directory. Writes one parallel archive
/// per root + a `roots.json` sidecar with per-slot metadata.
template <typename Type, typename Shell>
void save_es_roots(madness::World &world,
                   const typename ESSolver<Type, Shell>::State &state,
                   const std::string &dir,
                   bool converged) {
  const int n_roots = static_cast<int>(state.roots.size());
  MADNESS_CHECK(n_roots > 0);

  if (world.rank() == 0) {
    std::filesystem::create_directories(dir);
  }
  world.gop.fence();

  // Per-root binary archives — collective.
  for (int s = 0; s < n_roots; ++s) {
    const std::string path = dir + "/" + detail_save_load::root_file(s);
    state.roots[s].save(world, path);
  }

  if (world.rank() == 0) {
    const double thresh = madness::FunctionDefaults<3>::get_thresh();

    // Stable identity: use what the state carries, else assign 0..N-1 so a
    // save always records identity. Display names are grouped at this
    // protocol boundary (tol = 10 * thresh; matches doc 03).
    std::vector<int> stable_index = state.stable_index;
    assign_initial_stable_index(stable_index, n_roots);
    const std::vector<std::string> display_names =
        assign_display_names(state.omega, 10.0 * thresh);

    nlohmann::json j;
    j["type"]      = detail_save_load::type_tag<Type>();
    j["shell"]     = detail_save_load::shell_tag<Shell>();
    j["n_roots"]   = n_roots;
    const int k_now = madness::FunctionDefaults<3>::get_k();
    j["k"]            = k_now;
    j["thresh"]       = thresh;
    j["protocol_key"] = protocol_key(thresh, k_now);  // doc 13 join key
    j["iter"]         = state.iter;
    j["converged"]    = converged;

    // slot_permutation[slot] = stable_index — the cross-protocol root map.
    j["slot_permutation"] = stable_index;

    nlohmann::json roots_arr = nlohmann::json::array();
    const long n_omega = state.omega.size();
    for (int s = 0; s < n_roots; ++s) {
      nlohmann::json entry;
      entry["slot"]         = s;
      entry["stable_index"] = stable_index[s];
      entry["root_id"]      = make_root_id(stable_index[s]);
      entry["display_name"] =
          (static_cast<size_t>(s) < display_names.size()) ? display_names[s]
                                                          : std::string();
      entry["omega"] = (s < n_omega) ? state.omega(s)
                                     : std::numeric_limits<double>::quiet_NaN();
      entry["bsh_residual"] =
          (static_cast<size_t>(s) < state.last_bsh_residual.size())
              ? state.last_bsh_residual[s]
              : std::numeric_limits<double>::quiet_NaN();
      entry["density_residual"] =
          (static_cast<size_t>(s) < state.last_density_residual.size())
              ? state.last_density_residual[s]
              : std::numeric_limits<double>::quiet_NaN();
      entry["file"] = detail_save_load::root_file(s);
      roots_arr.push_back(entry);
    }
    j["roots"] = roots_arr;

    std::ofstream out(dir + "/roots.json");
    out << j.dump(2) << "\n";
  }
  world.gop.fence();
}

/// Load an ESSolver state previously written by `save_es_roots`.
///
/// Validates that the sidecar's `type` / `shell` / `n_roots` match the
/// requested template instantiation. Per-slot ω + residuals are
/// restored; `iter` resets to 0 and `rho_alpha_prev` stays empty so the
/// resumed solver picks up from a clean count.
template <typename Type, typename Shell>
typename ESSolver<Type, Shell>::State
load_es_roots(madness::World &world, const std::string &dir) {
  using State = typename ESSolver<Type, Shell>::State;

  // Rank 0 reads + parses; broadcast the n_roots scalar and the per-slot
  // ω / residuals via Tensor<double> + vector<double> hops. The
  // expected type/shell strings are baked in at compile time, so we
  // only need rank 0 to assert them — no broadcast needed for those.
  int n_roots = 0;
  int iter_at_save = 0;
  madness::Tensor<double> omega_recorded;
  std::vector<double> bsh_residual_recorded;
  std::vector<double> drho_residual_recorded;
  std::vector<int> stable_index_recorded;

  if (world.rank() == 0) {
    std::ifstream in(dir + "/roots.json");
    if (!in) {
      throw std::runtime_error("load_es_roots: cannot open " + dir +
                               "/roots.json");
    }
    nlohmann::json j;
    in >> j;

    const std::string saved_type  = j.value("type",  std::string{});
    const std::string saved_shell = j.value("shell", std::string{});
    const std::string want_type   = detail_save_load::type_tag<Type>();
    const std::string want_shell  = detail_save_load::shell_tag<Shell>();
    if (saved_type != want_type || saved_shell != want_shell) {
      throw std::runtime_error(
          "load_es_roots: saved type/shell = " + saved_type + "/" +
          saved_shell + ", requested = " + want_type + "/" + want_shell);
    }

    n_roots      = j.value("n_roots", 0);
    iter_at_save = j.value("iter", 0);
    MADNESS_CHECK(n_roots > 0);
    if (!j.contains("roots") || !j["roots"].is_array() ||
        static_cast<int>(j["roots"].size()) != n_roots) {
      throw std::runtime_error(
          "load_es_roots: roots[] array missing or wrong length in " +
          dir + "/roots.json");
    }

    omega_recorded         = madness::Tensor<double>(n_roots);
    bsh_residual_recorded  .assign(n_roots, 0.0);
    drho_residual_recorded .assign(n_roots, 0.0);
    stable_index_recorded  .assign(n_roots, 0);
    for (int s = 0; s < n_roots; ++s) {
      const auto &e = j["roots"][s];
      const int slot = e.value("slot", -1);
      MADNESS_CHECK(slot == s);
      omega_recorded(s)         = e.value("omega", 0.0);
      bsh_residual_recorded[s]  = e.value("bsh_residual", 0.0);
      drho_residual_recorded[s] = e.value("density_residual", 0.0);
      // Stable identity: legacy files predate it — fall back to slot so
      // the loaded bundle still carries a well-formed identity vector.
      stable_index_recorded[s]  = e.value("stable_index", s);
    }
  }

  // Broadcast the metadata to every rank.
  world.gop.broadcast(n_roots, 0);
  world.gop.broadcast(iter_at_save, 0);
  if (world.rank() != 0) {
    omega_recorded         = madness::Tensor<double>(n_roots);
    bsh_residual_recorded  .assign(n_roots, 0.0);
    drho_residual_recorded .assign(n_roots, 0.0);
    stable_index_recorded  .assign(n_roots, 0);
  }
  world.gop.broadcast(omega_recorded.ptr(), n_roots, 0);
  world.gop.broadcast(bsh_residual_recorded.data(),  n_roots, 0);
  world.gop.broadcast(drho_residual_recorded.data(), n_roots, 0);
  world.gop.broadcast(stable_index_recorded.data(),  n_roots, 0);

  // Per-root binary archives — collective.
  State s;
  s.roots.reserve(n_roots);
  for (int b = 0; b < n_roots; ++b) {
    const std::string path = dir + "/" + detail_save_load::root_file(b);
    s.roots.push_back(
        std::remove_reference_t<decltype(s.roots[0])>::load(world, path));
  }

  s.omega                 = std::move(omega_recorded);
  s.last_bsh_residual     = std::move(bsh_residual_recorded);
  s.last_density_residual = std::move(drho_residual_recorded);
  s.stable_index          = std::move(stable_index_recorded);
  s.iter                  = 0;       // fresh count for the resumed solver
  s.diverged              = false;
  // rho_alpha_prev left empty — first-iter density-residual check in
  // ESSolver::step() is already gated by iter <= 1.

  if (world.rank() == 0) {
    madness::print("[LOAD] es_roots: n_roots=", n_roots,
                   " iter_at_save=", iter_at_save,
                   " dir=", dir);
    for (int b = 0; b < n_roots; ++b) {
      madness::print("       slot=", b,
                     "  root_id=", make_root_id(s.stable_index[b]),
                     "  omega=", s.omega(b),
                     "  bsh_res=", s.last_bsh_residual[b],
                     "  drho_res=", s.last_density_residual[b]);
    }
  }

  return s;
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_ES_SAVE_LOAD_HPP
