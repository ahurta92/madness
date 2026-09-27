#pragma once
// The HDF5 twin of the DALTON ground-state seed. The seed itself (projection,
// archive, madqc hook) lives in chem (madness/chem/dalton_seed.h); only this
// callback stays in molresponse, because the HDF5 I/O is molresponse's.
// Pass it as GsSeedOptions::after_save or to seed_gs_from_dalton_dir.

#include <madness/chem/dalton_seed.h>
#include "function_hdf5_io.hpp"   // save_parallel_archive_hdf5 (no-op without MADNESS_HAS_HDF5)
#include "restart_hdf5.hpp"       // its /restart attributes

namespace molresponse_v3 {

/// A no-op unless MADNESS was built with HDF5 and the HDF5 backend is on
/// (io backend hdf5): then each seed archive gets a <name>.h5 twin, with the
/// header duplicated as /restart attributes, so the DALTON starting point can
/// be inspected next to the MADNESS result.
inline madness::GsSeedOptions::AfterSave gs_seed_hdf5_twin() {
#ifdef MADNESS_HAS_HDF5
  return [](madness::World &world, const std::string &name,
            const madness::GsSeedOptions::Emitter &emit,
            const madness::RestartMetadata &meta, unsigned int n_mo) {
    if (hdf5_io_enabled())
      save_parallel_archive_hdf5(world, name + ".h5", /*deflate=*/0,
                                 [&](auto &ar) { emit(ar); },
                                 [&](hid_t file) { write_restart_attributes(file, meta, n_mo); });
  };
#else
  return {};
#endif
}

} // namespace molresponse_v3
