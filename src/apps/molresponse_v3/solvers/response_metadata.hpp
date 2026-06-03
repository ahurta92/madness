#ifndef MOLRESPONSE_V3_SOLVERS_RESPONSE_METADATA_HPP
#define MOLRESPONSE_V3_SOLVERS_RESPONSE_METADATA_HPP

// =========================================================================
// ResponseMetadata — the unified response_metadata.json writer (doc 13).
//
// Single source of truth for the aggregate metadata file FD and ES
// persistence share. Schema:
//
//   {
//     "schema_version": 1,
//     "protocols":       { "<protocol_key>": {thresh, k, index} },
//     "fd_states":       { "<pert>": { "<protocol_key>": { "<freq_key>":
//                          {freq, type, shell, converged, iter,
//                           bsh_residual, archive} } } },
//     "excited_states":  { "<protocol_key>": {type, shell, n_roots,
//                          bundle_dir, converged, slot_permutation,
//                          roots[]} },
//     "properties":      { "<name>": { "<protocol_key>": [ {...}, ... ] } }
//   }
//
// The writer is pure JSON — no World, no MPI dependency. CALLERS guard with
//   if (world.rank() == 0) { ResponseMetadata::load_or_create(...) ... save(); }
//   world.gop.fence();
// so other ranks observe the file after rank 0 has finished writing.
//
// Atomic-write contract: save() writes to `<path>.tmp` then renames over
// `<path>`. On POSIX the rename is atomic within a filesystem, so a crash
// mid-write leaves either the old file intact or the new one — never a
// half-written aggregate. Two concurrent writers on the same file are not
// supported; serialize at the orchestration layer.
//
// The writer takes payloads as raw nlohmann::json blobs rather than typed
// structs: keeps the FD/ES save paths in charge of what to record, and the
// writer in charge of where it lands. Helper builders (typed -> json) can
// live next to the solver code that owns the types.
// =========================================================================

#include <madness/external/nlohmann_json/json.hpp>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

namespace molresponse_v3 {

class ResponseMetadata {
public:
  static constexpr int kSchemaVersion = 1;

  /// Load the file if it exists, otherwise seed a fresh schema. Rank-0 only
  /// (filesystem op); the caller is responsible for guarding with rank().
  static ResponseMetadata load_or_create(const std::string &path) {
    ResponseMetadata m;
    m.path_ = path;
    if (std::filesystem::exists(path)) {
      std::ifstream in(path);
      if (!in) {
        throw std::runtime_error("ResponseMetadata: cannot read " + path);
      }
      in >> m.j_;
      const int v = m.j_.value("schema_version", 0);
      if (v != kSchemaVersion) {
        throw std::runtime_error(
            "ResponseMetadata: schema_version " + std::to_string(v) +
            " in " + path + " not recognized (expected " +
            std::to_string(kSchemaVersion) + ")");
      }
      // Ensure required top-level keys exist for upsert paths below.
      for (const char *k : {"protocols", "fd_states",
                             "excited_states", "properties", "vbc_states"}) {
        if (!m.j_.contains(k)) m.j_[k] = nlohmann::json::object();
      }
    } else {
      m.j_ = nlohmann::json{
          {"schema_version", kSchemaVersion},
          {"protocols",      nlohmann::json::object()},
          {"fd_states",      nlohmann::json::object()},
          {"excited_states", nlohmann::json::object()},
          {"properties",     nlohmann::json::object()},
          {"vbc_states",     nlohmann::json::object()},
      };
    }
    return m;
  }

  /// Upsert a protocol registry entry. Overwrites any existing entry at
  /// `key` (the (thresh, k) physical identity is stable — only `index` can
  /// legitimately change across runs with different ramps).
  void set_protocol(const std::string &key, double thresh, int k, int index) {
    j_["protocols"][key] = nlohmann::json{
        {"thresh", thresh}, {"k", k}, {"index", index}};
  }

  /// Upsert one FD point: fd_states/<pert>/<protocol_key>/<freq_key> = entry.
  void set_fd_state(const std::string &pert,
                    const std::string &protocol_key,
                    const std::string &freq_key,
                    const nlohmann::json &entry) {
    j_["fd_states"][pert][protocol_key][freq_key] = entry;
  }

  /// Replace the ES bundle entry at this protocol. ES is bundle-at-a-time —
  /// the whole entry is rewritten on save, not patched root-by-root.
  void set_es_bundle(const std::string &protocol_key,
                     const nlohmann::json &entry) {
    j_["excited_states"][protocol_key] = entry;
  }

  /// Upsert one VBC quadratic source: vbc_states/<vbc_id>/<protocol_key> = entry.
  void set_vbc_state(const std::string &vbc_id,
                     const std::string &protocol_key,
                     const nlohmann::json &entry) {
    j_["vbc_states"][vbc_id][protocol_key] = entry;
  }

  /// Append a property record to properties/<name>/<protocol_key>[]. The
  /// caller stamps it with whatever provenance is meaningful (es_root_id,
  /// fd_freq, the value, etc. — see doc 13's matching contract).
  void add_property(const std::string &name,
                    const std::string &protocol_key,
                    const nlohmann::json &record) {
    auto &arr = j_["properties"][name][protocol_key];
    if (!arr.is_array()) arr = nlohmann::json::array();
    arr.push_back(record);
  }

  /// Atomic write. Throws on filesystem error.
  void save() const {
    const std::string tmp = path_ + ".tmp";
    {
      std::ofstream out(tmp);
      if (!out) throw std::runtime_error(
          "ResponseMetadata: cannot open for write: " + tmp);
      out << j_.dump(2) << "\n";
    }
    std::filesystem::rename(tmp, path_);
  }

  /// Direct read access — for tests and the (rare) case a caller needs to
  /// query something the typed setters don't cover.
  const nlohmann::json &json() const { return j_; }
  const std::string    &path() const { return path_; }

  /// Canonical frequency key for fd_states / archive naming. f%.5f matches
  /// the dimensionless precision used in v2 and is enough to distinguish
  /// any physically relevant frequency we'd solve at. Doc 13.
  static std::string freq_key(double freq) {
    char buf[32];
    std::snprintf(buf, sizeof buf, "f%.5f", freq);
    return {buf};
  }

private:
  nlohmann::json j_;
  std::string    path_;
};

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_RESPONSE_METADATA_HPP
