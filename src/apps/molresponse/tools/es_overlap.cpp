// es_overlap — response-metric overlap of a SEED excited-state vector with each
// root of a converged excited-state set. Identifies which converged MADNESS
// root corresponds to the operator-selected seed by CHARACTER (max
// |<seed|root>|), not by energy — robust when the target sits in a
// dense/Rydberg manifold where energy-ordering is ambiguous. This is the "(b)
// overlap-based pinning" step of the dalton-cis-seed targeting thread: pick the
// state that best matches the cheap-calc (operator-selected) eigenvector.
//
// Overlap uses the response metric rs::metric_inner (TDA: sum_i <x_i|x_i>;
// full: <X|X> - <Y|Y>), normalized:
//   ov_i = |<seed|root_i>| / sqrt(<seed|seed> <root_i|root_i>).
//
// Usage:
//   es_overlap --seed-dir=<folder> --bundle-dir=<folder>
// Each folder is an excited-state set node (<run>/states/response/excited-*/)
// or an excited-state seed (<run>/seeds/response/excited-*/). The model (TDA or
// full) comes from the descriptions and must be the same for both; --full is
// accepted and not needed. Both are loaded at the --bundle-dir set's saved
// level. Prints the full |<seed_r|root_i>| matrix (all seed roots x all
// converged roots) plus the best match per seed root. Nothing is written.

#include "../kernels/response_space_ops.hpp"  // rs::metric_inner
#include "../kernels/tags.hpp"                // TDA, Full, ClosedShell
#include "saved_states.hpp"

#include <madness/misc/info.h>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace madness;
using namespace molresponse_v3;

namespace {

template <typename Type>
void run_overlap(World &world, const saved_states::SavedNode &seed_node,
                 const saved_states::SavedNode &bundle_node) {
  auto seed = saved_states::load_excited<Type>(world, seed_node);
  auto conv = saved_states::load_excited<Type>(world, bundle_node);
  MADNESS_CHECK(!seed.roots.empty());
  const std::size_t R = seed.roots.size(), N = conv.roots.size();

  std::vector<double> snn(R), rnn(N);
  for (std::size_t r = 0; r < R; ++r) snn[r] = rs::metric_inner(seed.roots[r], seed.roots[r]);
  for (std::size_t i = 0; i < N; ++i) rnn[i] = rs::metric_inner(conv.roots[i], conv.roots[i]);
  std::vector<std::vector<double>> ov(R, std::vector<double>(N, 0.0));
  for (std::size_t r = 0; r < R; ++r)
    for (std::size_t i = 0; i < N; ++i) {
      const double sr = rs::metric_inner(seed.roots[r], conv.roots[i]);  // collective
      ov[r][i] = (snn[r] > 0.0 && rnn[i] > 0.0) ? std::abs(sr) / std::sqrt(snn[r] * rnn[i]) : 0.0;
    }

  if (world.rank() == 0) {
    print("\n=== es_overlap ===");
    print("  seed   :", seed_node.folder.string(), " n_roots=", (int)R);
    print("  bundle :", bundle_node.folder.string(), " n_roots=", (int)N);
    std::printf("  converged roots: ");
    for (std::size_t i = 0; i < N; ++i)
      std::printf("  %2zu:%8.4f", i, i < (std::size_t)conv.omega.size() ? conv.omega(long(i)) : 0.0);
    std::printf("  (au)\n  |<seed_r|root_i>| (rows = seed roots, omega_seed in au):\n");
    for (std::size_t r = 0; r < R; ++r) {
      const double os = r < (std::size_t)seed.omega.size() ? seed.omega(long(r)) : 0.0;
      std::printf("  seed %2zu %8.4f :", r, os);
      std::size_t best = 0;
      for (std::size_t i = 1; i < N; ++i)
        if (ov[r][i] > ov[r][best]) best = i;
      double sumsq = 0.0;
      for (double v : ov[r]) sumsq += v * v;
      for (std::size_t i = 0; i < N; ++i) std::printf(" %9.6f", ov[r][i]);
      std::printf("   best=%zu (%.6f)  ||P_conv seed||=%.6f\n", best, ov[r][best], std::sqrt(sumsq));
    }
    std::printf("  (||P_conv seed|| = norm of the seed's projection onto the span of the "
                "converged roots; 1 = the seed lies in that span)\n");
  }
}

}  // namespace

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  int rc = 0;
  try {
    startup(world, argc, argv, true);
    commandlineparser parser(argc, argv);
    if (!parser.key_exists("seed-dir") || !parser.key_exists("bundle-dir")) {
      if (world.rank() == 0)
        print("Usage: es_overlap --seed-dir=<excited-state set or seed folder> "
              "--bundle-dir=<excited-state set or seed folder>");
      finalize();
      return 2;
    }
    {
      const auto seed = saved_states::require_node(world, parser.value_raw("seed-dir"));
      const auto bundle = saved_states::require_node(world, parser.value_raw("bundle-dir"));
      for (const auto *node : {&seed, &bundle}) {
        MADNESS_CHECK_THROW(saved_states::is_excited(*node),
                            ("not an excited-state set or seed: " + node->folder.string()).c_str());
        MADNESS_CHECK_THROW(saved_states::closed_shell(*node),
                            ("not a closed-shell excited-state set: " + node->folder.string()).c_str());
      }
      MADNESS_CHECK_THROW(saved_states::full_layout(seed) == saved_states::full_layout(bundle),
                          "the seed and the set are of different models (TDA and full)");
      // Both at the converged set's saved level; the seed is projected to it.
      const auto &level = bundle.description.protocol;
      set_response_protocol(world, saved_states::cell_half_width(bundle), level.thresh, level.k, false);
      saved_states::print_loaded(world, seed);
      saved_states::print_loaded(world, bundle);
      if (saved_states::full_layout(bundle)) run_overlap<Full>(world, seed, bundle);
      else run_overlap<TDA>(world, seed, bundle);
    }
  } catch (const std::exception &e) {
    if (world.rank() == 0) print("Error:", e.what());
    rc = 1;
  }
  finalize();
  return rc;
}
