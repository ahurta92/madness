// es_overlap — response-metric overlap of a SEED excited-state vector with each
// root of a converged ES bundle. Identifies which converged MADNESS root
// corresponds to the operator-selected seed by CHARACTER (max |<seed|root>|),
// not by energy — robust when the target sits in a dense/Rydberg manifold where
// energy-ordering is ambiguous. This is the "(b) overlap-based pinning" step of
// the dalton-cis-seed targeting thread: pick the state that best matches the
// cheap-calc (operator-selected) eigenvector.
//
// Overlap uses the response metric rs::metric_inner (TDA: sum_i <x_i|x_i>),
// normalized: ov_i = |<seed|root_i>| / sqrt(<seed|seed> <root_i|root_i>).
//
// Usage:
//   es_overlap --seed-dir=<es__KEY seed bundle> --bundle-dir=<converged es bundle> [--full]
// Prints the full |<seed_r|root_i>| matrix (all seed roots x all converged
// roots) plus the best match per seed root. --full loads Full (X,Y) bundles and
// uses the RPA metric <X|X> - <Y|Y>; default TDA (X only).

// The reader of the older excited-state bundle layout (es__<key>/roots.json)
// is gone. Reading excited-state node folders (states/response/excited-*/,
// seeds/response/excited-*/) arrives in step 4c part 5d; until then this
// program only says so.

#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

using namespace madness;

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  if (world.rank() == 0)
    print("es_overlap: not available in this version. It read the older "
          "excited-state bundle layout, which is gone; reading node folders "
          "(states/response/, seeds/response/) arrives in the next part.");
  finalize();
  return 2;
}
