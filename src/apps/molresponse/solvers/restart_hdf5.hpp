// =========================================================================
// restart_hdf5.hpp — the restartdata header as HDF5 attributes.
//
// A ground-state `<prefix>.restartdata.h5` holds the same stream as the native
// archive (RestartMetadata, then the orbitals) inside one opaque byte dataset
// (function_hdf5_io.hpp, save_parallel_archive_hdf5). The stream stays the
// source of truth, so one loader serves both formats; this file duplicates the
// header as attributes of a `/restart` group, written into the same tmp file
// before the same rename, so that the archive's identity (archive_id), its
// Hamiltonian and its electron counts can be read without the blob:
//
//   /restart  schema, header_version, archive_id (16 hex digits), origin,
//             representation, hamiltonian_key (json), localize, molecule
//             (json), nalpha, nbeta, nmo_alpha, nmo_beta, spin_restricted,
//             k, L, converged_for_thresh, converged_for_dconv, current_energy,
//             madness_version
//
// A reader that finds both must find them agreeing; see GroundState.
// =========================================================================
#ifndef MOLRESPONSE_V3_SOLVERS_RESTART_HDF5_HPP
#define MOLRESPONSE_V3_SOLVERS_RESTART_HDF5_HPP

#ifdef MADNESS_HAS_HDF5

#include "function_hdf5_io.hpp"

#include <madness/chem/Restart.h>
#include <nlohmann/json.hpp>

#include <filesystem>
#include <optional>
#include <string>

namespace molresponse_v3 {

inline constexpr int kRestartAttrSchema = 1;

/// what the `/restart` group records
struct RestartAttributes {
  int schema = 0;
  unsigned int header_version = 0;
  madness::ArchiveId archive_id = 0;
  std::string origin;
  int representation = 0;
  nlohmann::json hamiltonian_key;
  std::string localize;
  nlohmann::json molecule;
  int nalpha = -1, nbeta = -1, nmo_alpha = 0, nmo_beta = -1;
  bool spin_restricted = true;
  int k = 0;
  double L = 0.0;
  double converged_for_thresh = 1.e10, converged_for_dconv = 1.e10;
  double current_energy = 1.e10;
  std::string madness_version;
};

/// write \p meta as the `/restart` group of the open file \p file (rank 0)
inline void write_restart_attributes(hid_t file, const madness::RestartMetadata &meta,
                                     const unsigned int nmo_alpha) {
  namespace h = detail_function_hdf5;
  hid_t g = H5Gcreate2(file, "restart", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  MADNESS_CHECK_THROW(g >= 0, "write_restart_attributes: H5Gcreate2(/restart) failed");
  h::write_attr_i(g, "schema", kRestartAttrSchema);
  h::write_attr_i(g, "header_version", madness::RestartMetadata::CURRENT_VERSION);
  h::write_attr_s(g, "archive_id", madness::archive_id_to_string(meta.archive_id));
  h::write_attr_s(g, "origin", meta.origin);
  h::write_attr_i(g, "representation", static_cast<int>(meta.representation));
  h::write_attr_s(g, "hamiltonian_key", meta.hamiltonian_key().to_json().dump());
  h::write_attr_s(g, "localize", meta.localize);
  h::write_attr_s(g, "molecule", meta.molecule.to_json().dump());
  h::write_attr_i(g, "nalpha", meta.nalpha);
  h::write_attr_i(g, "nbeta", meta.nbeta);
  h::write_attr_i(g, "nmo_alpha", nmo_alpha);
  h::write_attr_i(g, "nmo_beta", meta.nmo_beta);
  h::write_attr_i(g, "spin_restricted", meta.spin_restricted ? 1 : 0);
  h::write_attr_i(g, "k", meta.k);
  h::write_attr_d(g, "L", meta.L);
  h::write_attr_d(g, "converged_for_thresh", meta.converged_for_thresh);
  h::write_attr_d(g, "converged_for_dconv", meta.converged_for_dconv);
  h::write_attr_d(g, "current_energy", meta.current_energy);
  h::write_attr_s(g, "madness_version", meta.madness_version);
  H5Gclose(g);
}

/// the `/restart` group of \p path, without reading the blob; rank-local
///
/// nullopt if the file is absent or has no such group (written before it
/// existed), or cannot be read.
inline std::optional<RestartAttributes> peek_restartdata_hdf5(const std::string &path) {
  namespace h = detail_function_hdf5;
  if (!std::filesystem::exists(path)) return std::nullopt;
  h::H5ErrorScope quiet;
  hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) return std::nullopt;
  std::optional<RestartAttributes> out;
  if (H5Lexists(file, "restart", H5P_DEFAULT) > 0) {
    hid_t g = H5Gopen2(file, "restart", H5P_DEFAULT);
    if (g >= 0) {
      try {
        RestartAttributes a;
        a.schema = int(h::read_attr_i(g, "schema"));
        a.header_version = unsigned(h::read_attr_i(g, "header_version"));
        a.archive_id = madness::archive_id_from_string(h::read_attr_s(g, "archive_id"));
        a.origin = h::read_attr_s(g, "origin");
        a.representation = int(h::read_attr_i(g, "representation"));
        a.hamiltonian_key = nlohmann::json::parse(h::read_attr_s(g, "hamiltonian_key"));
        a.localize = h::read_attr_s(g, "localize");
        a.molecule = nlohmann::json::parse(h::read_attr_s(g, "molecule"));
        a.nalpha = int(h::read_attr_i(g, "nalpha"));
        a.nbeta = int(h::read_attr_i(g, "nbeta"));
        a.nmo_alpha = int(h::read_attr_i(g, "nmo_alpha"));
        a.nmo_beta = int(h::read_attr_i(g, "nmo_beta"));
        a.spin_restricted = h::read_attr_i(g, "spin_restricted") != 0;
        a.k = int(h::read_attr_i(g, "k"));
        a.L = h::read_attr_d(g, "L");
        a.converged_for_thresh = h::read_attr_d(g, "converged_for_thresh");
        a.converged_for_dconv = h::read_attr_d(g, "converged_for_dconv");
        a.current_energy = h::read_attr_d(g, "current_energy");
        a.madness_version = h::read_attr_s(g, "madness_version");
        out = a;
      } catch (...) {
        out.reset();
      }
      H5Gclose(g);
    }
  }
  H5Fclose(file);
  return out;
}

}  // namespace molresponse_v3

#endif  // MADNESS_HAS_HDF5
#endif  // MOLRESPONSE_V3_SOLVERS_RESTART_HDF5_HPP
