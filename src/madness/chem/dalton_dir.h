#ifndef MADNESS_CHEM_DALTON_DIR_H
#define MADNESS_CHEM_DALTON_DIR_H

// Locate the output of a DALTON run and check that it is for the molecule at
// hand. Import-only: MADNESS never runs DALTON.
//
//   locate_dalton_dir: finds molden, RSPVEC and .out in a directory by
//     convention (loose files, or a unique *.tar.gz extracted on rank 0), and
//     reads the basis, method and molden geometry into a DaltonManifest;
//   fingerprint_dalton_geometry: compares that geometry with a Molecule, atom
//     by atom, within a tolerance in bohr.
//
// Consumers: the SCF ground-state seed (chem) and the molresponse response
// seeds (apps/molresponse/solvers/dalton_import.hpp).

#include <madness/chem/fnv1a64.h>
#include <madness/chem/molden_gto.h>
#include <madness/chem/molecule.h>
#include <madness/misc/info.h>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace madness {

// -------------------------------------------------------------------------
// Manifest — everything the import knows about the DALTON directory.
// -------------------------------------------------------------------------
struct DaltonManifest {
  std::string dir;           // the dalton.dir as given
  std::string molden_path;   // resolved molden file
  std::string rspvec_path;   // resolved RSPVEC binary
  std::string out_path;      // resolved .out ("" = none found; provenance-degraded)
  std::string basis  = "unknown";   // from the .out
  std::string method = "unknown";   // "HF" | "unknown" (non-HF is rejected earlier)
  std::string geometry_hash;        // FNV-1a-64 over the canonical geometry string
  std::vector<int> atomic_numbers;                // molden [Atoms]
  std::vector<std::array<double, 3>> coords;      // bohr

  /// The provenance a consumer records for a seeded run: the DALTON files,
  /// basis, method and geometry hash, plus the MADNESS build that wrote the
  /// MRA artifacts (molresponse stamps it into response_metadata.json).
  nlohmann::json provenance() const {
    return nlohmann::json{
        {"dalton_dir",    dir},
        {"molden",        molden_path},
        {"rspvec",        rspvec_path},
        {"out",           out_path},
        {"basis",         basis},
        {"method",        method},
        {"geometry_hash", geometry_hash},
        // Which MADNESS build wrote the MRA seed artifacts. Mandatory since
        // the archive format is build-locked ("incompatible MADNESS version"
        // on load) — without this a rejected seed dir is undiagnosable.
        {"mra_writer", nlohmann::json{
            {"version",    madness::info::version()},
            {"git_commit", madness::info::git_commit()},
        }},
    };
  }
};


namespace detail_dalton_import {

namespace fs = std::filesystem;

/// Candidates in `dir` whose filename matches: exact `name`, or suffix `ext`.
inline std::vector<std::string> find_candidates(const std::string &dir,
                                                const std::string &exact,
                                                const std::string &suffix) {
  std::vector<std::string> hits;
  std::error_code ec;
  // Exact conventional name wins outright.
  if (!exact.empty() && fs::exists(dir + "/" + exact, ec))
    return {dir + "/" + exact};
  for (const auto &e : fs::directory_iterator(dir, ec)) {
    if (ec || !e.is_regular_file()) continue;
    const std::string name = e.path().filename().string();
    if (!suffix.empty() && name.size() > suffix.size() &&
        name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0)
      hits.push_back(e.path().string());
  }
  std::sort(hits.begin(), hits.end());
  return hits;
}

inline std::string list_join(const std::vector<std::string> &v) {
  std::string s;
  for (const auto &x : v) { if (!s.empty()) s += ", "; s += x; }
  return s;
}

/// Resolve one artifact slot: explicit override > exact conventional name >
/// unique suffix match. Zero or multiple suffix matches -> "" (caller decides
/// whether that is fatal); an override that does not exist IS fatal here.
inline std::string resolve_slot(const std::string &dir,
                                const std::string &override_path,
                                const std::string &exact,
                                const std::string &suffix,
                                std::string &err, const char *what) {
  if (!override_path.empty()) {
    if (!fs::exists(override_path)) {
      err = std::string("dalton import: explicit ") + what + " override '" +
            override_path + "' does not exist";
    }
    return override_path;
  }
  auto c = find_candidates(dir, exact, suffix);
  if (c.size() == 1) return c[0];
  if (c.size() > 1) {
    err = std::string("dalton import: multiple ") + what + " candidates in " +
          dir + " (" + list_join(c) + ") — pass the explicit override";
  }
  return {};
}

} // namespace detail_dalton_import

// -------------------------------------------------------------------------
// locate — by-convention artifact discovery (+ tarball extraction).
// Pure rank-0 function: performs filesystem reads and (only when molden/
// RSPVEC are not loose and a unique *.tar.gz exists) an extraction into
// `extract_dir` via a shell-out to `tar`. Throws on any ambiguity.
//
// `allow_extract = false` makes it read-only: the tarball branch is skipped,
// so it is safe to call from inside a subworld, where several rank-0s would
// otherwise extract into the same directory at once (B24). Prefer the
// prepared-seed registry above; this is the fallback for a directory whose
// artifacts are already loose.
// -------------------------------------------------------------------------
inline DaltonManifest locate_dalton_dir(const std::string &dir,
                                        const std::string &extract_dir,
                                        const std::string &molden_override = {},
                                        const std::string &rspvec_override = {},
                                        const std::string &out_override = {},
                                        bool allow_extract = true) {
  namespace fs = std::filesystem;
  using namespace detail_dalton_import;

  if (!fs::exists(dir) || !fs::is_directory(dir))
    throw std::runtime_error("dalton import: dalton.dir '" + dir +
                             "' is not a directory");

  DaltonManifest m;
  m.dir = dir;
  std::string err;

  m.molden_path = resolve_slot(dir, molden_override, "molden.inp",
                               ".molden.inp", err, "molden");
  if (err.empty() && m.molden_path.empty())
    m.molden_path = resolve_slot(dir, "", "", ".molden", err, "molden");
  if (!err.empty()) throw std::runtime_error(err);

  m.rspvec_path = resolve_slot(dir, rspvec_override, "RSPVEC",
                               ".RSPVEC", err, "RSPVEC");
  if (!err.empty()) throw std::runtime_error(err);

  // Not loose -> try the DALTON output tarball (unique *.tar.gz). Extraction
  // is documented behavior: `tar -xzf <tar> -C <extract_dir> RSPVEC molden.inp`
  // (only the two members), so pre-extraction is never required but always
  // honored. tar's exit code is ignored on purpose (one member may be absent);
  // what matters is which files exist afterwards.
  if ((m.molden_path.empty() || m.rspvec_path.empty()) && allow_extract) {
    auto tars = find_candidates(dir, "", ".tar.gz");
    if (tars.size() > 1)
      throw std::runtime_error(
          "dalton import: molden/RSPVEC not loose in " + dir +
          " and multiple *.tar.gz candidates (" + list_join(tars) +
          ") — extract manually or pass explicit overrides");
    if (tars.size() == 1) {
      fs::create_directories(extract_dir);
      const std::string cmd = "tar -xzf '" + tars[0] + "' -C '" + extract_dir +
                              "' RSPVEC molden.inp 2>/dev/null";
      // glibc marks system() warn_unused_result and (void) does not silence
      // it under gcc. tar exits nonzero if EITHER member is absent, and a
      // GS-only directory legitimately lacks RSPVEC — so partial extraction
      // is decided by the existence checks below, not by tar's rc. Only a
      // total failure (nothing extracted) is an error here.
      const int tar_rc = std::system(cmd.c_str());
      if (tar_rc != 0 && !fs::exists(extract_dir + "/molden.inp") &&
          !fs::exists(extract_dir + "/RSPVEC"))
        throw std::runtime_error(
            "dalton import: extraction of '" + tars[0] +
            "' produced neither RSPVEC nor molden.inp (tar exit " +
            std::to_string(tar_rc) + ")");
      if (m.molden_path.empty() && fs::exists(extract_dir + "/molden.inp"))
        m.molden_path = extract_dir + "/molden.inp";
      if (m.rspvec_path.empty() && fs::exists(extract_dir + "/RSPVEC"))
        m.rspvec_path = extract_dir + "/RSPVEC";
    }
  }
  if (m.molden_path.empty())
    throw std::runtime_error(
        "dalton import: no molden file in " + dir +
        " (looked for molden.inp, *.molden.inp, *.molden, and inside a unique "
        "*.tar.gz) — pass --dalton-molden=PATH");
  if (m.rspvec_path.empty())
    throw std::runtime_error(
        "dalton import: no RSPVEC in " + dir +
        " (looked for RSPVEC, *.RSPVEC, and inside a unique *.tar.gz) — pass "
        "--dalton-rspvec=PATH");

  {
    auto outs = find_candidates(dir, "", ".out");
    if (!out_override.empty()) {
      if (!fs::exists(out_override))
        throw std::runtime_error("dalton import: explicit out override '" +
                                 out_override + "' does not exist");
      m.out_path = out_override;
    } else if (outs.size() == 1) {
      m.out_path = outs[0];
    } else if (outs.size() > 1) {
      // Several *.out: keep the ones that are DALTON outputs (banner in the
      // first lines). A SLURM stdout named <job>.out next to <dal>_<mol>.out is
      // the common case (closeout jobs 2162014/15/25 aborted here).
      std::vector<std::string> dalton_outs;
      for (const auto &o : outs) {
        std::ifstream in(o);
        std::string line;
        for (int i = 0; i < 80 && std::getline(in, line); ++i) {
          // banner lines of a real DALTON output (a SLURM log that merely
          // says "Running DALTON" must not match)
          if (line.find("This is output from DALTON") != std::string::npos ||
              line.find("Dalton - An Electronic Structure Program") != std::string::npos) {
            dalton_outs.push_back(o);
            break;
          }
        }
      }
      if (dalton_outs.size() == 1)
        m.out_path = dalton_outs[0];
      else
        throw std::runtime_error(
            "dalton import: multiple *.out candidates in " + dir + " (" +
            list_join(outs) + "), " + std::to_string(dalton_outs.size()) +
            " of them DALTON outputs — pass --dalton-out=PATH");
    }
    // zero .out files: allowed (provenance-degraded; caller warns).
  }

  // Geometry from the molden (authoritative: it is what the MOs live on).
  {
    DaltonMoldenResult mol = read_molden(m.molden_path);
    m.atomic_numbers = mol.atomic_numbers;
    m.coords         = mol.coords;
    std::uint64_t h = kFnv1a64Basis;
    char line[128];
    for (size_t i = 0; i < m.coords.size(); ++i) {
      std::snprintf(line, sizeof line, "%d %.6f %.6f %.6f\n",
                    m.atomic_numbers[i], m.coords[i][0], m.coords[i][1],
                    m.coords[i][2]);
      h = fnv1a64_update(h, line, std::strlen(line));
    }
    char hex[17];
    std::snprintf(hex, sizeof hex, "%016llx",
                  static_cast<unsigned long long>(h));
    m.geometry_hash = hex;
  }

  // Method + basis from the .out (when present).
  if (!m.out_path.empty()) {
    std::ifstream in(m.out_path);
    std::string line;
    const std::regex basis_re("Basis set used is \"([^\"]+)\"");
    const std::regex wf_re("Wave function type\\s+-+\\s*(.+?)\\s*-+\\s*$");
    std::smatch sm;
    while (std::getline(in, line)) {
      if (m.basis == "unknown" && std::regex_search(line, sm, basis_re))
        m.basis = sm[1];
      if (m.method == "unknown" && std::regex_search(line, sm, wf_re))
        m.method = sm[1];
      if (m.basis != "unknown" && m.method != "unknown") break;
    }
  }
  return m;
}

// -------------------------------------------------------------------------
// Geometry fingerprint — DALTON geometry vs the active Molecule.
// -------------------------------------------------------------------------
struct DaltonGeometryCheck {
  bool        ok = false;
  double      max_dev = 0.0;   // bohr
  std::string report;          // both geometries side by side + verdict
};

inline DaltonGeometryCheck
fingerprint_dalton_geometry(const DaltonManifest &m,
                            const madness::Molecule &mol,
                            double tol_bohr = 1e-4) {
  DaltonGeometryCheck c;
  std::ostringstream r;
  r << "dalton import geometry fingerprint (bohr; tol=" << tol_bohr << "):\n";
  r << "  DALTON (" << m.molden_path << ")  vs  active Molecule\n";
  const size_t nd = m.coords.size(), nm = mol.natom();
  if (nd != nm) {
    r << "  ATOM-COUNT MISMATCH: dalton natom=" << nd
      << "  molecule natom=" << nm << "\n";
    for (size_t i = 0; i < nd; ++i)
      r << "    dalton  Z=" << m.atomic_numbers[i] << "  " << m.coords[i][0]
        << " " << m.coords[i][1] << " " << m.coords[i][2] << "\n";
    for (size_t i = 0; i < nm; ++i) {
      const auto &a = mol.get_atom(static_cast<unsigned int>(i));
      r << "    active  Z=" << a.atomic_number << "  " << a.x << " " << a.y
        << " " << a.z << "\n";
    }
    c.ok = false;
    c.max_dev = std::numeric_limits<double>::infinity();
    c.report = r.str();
    return c;
  }
  bool z_ok = true;
  for (size_t i = 0; i < nd; ++i) {
    const auto &a = mol.get_atom(static_cast<unsigned int>(i));
    const double dx = m.coords[i][0] - a.x;
    const double dy = m.coords[i][1] - a.y;
    const double dz = m.coords[i][2] - a.z;
    const double dev = std::max({std::abs(dx), std::abs(dy), std::abs(dz)});
    c.max_dev = std::max(c.max_dev, dev);
    const bool zmatch =
        (m.atomic_numbers[i] == static_cast<int>(a.atomic_number));
    z_ok = z_ok && zmatch;
    char buf[256];
    std::snprintf(buf, sizeof buf,
                  "  atom %2zu  Z %2d|%2u  dalton % .8f % .8f % .8f   active "
                  "% .8f % .8f % .8f   |dev| %.3e%s\n",
                  i, m.atomic_numbers[i], a.atomic_number, m.coords[i][0],
                  m.coords[i][1], m.coords[i][2], a.x, a.y, a.z, dev,
                  zmatch ? "" : "  <-- ELEMENT MISMATCH");
    r << buf;
  }
  c.ok = z_ok && (c.max_dev <= tol_bohr);
  r << "  max deviation = " << c.max_dev << " bohr  ->  "
    << (c.ok ? "MATCH" : "MISMATCH") << "\n";
  c.report = r.str();
  return c;
}


} // namespace madness

#endif // MADNESS_CHEM_DALTON_DIR_H
