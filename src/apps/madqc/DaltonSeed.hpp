#pragma once
// S5, D-112: the Dalton seed's member kind, held by the applications (C5: the
// storage core does not know it). A seeded stage's generations each hold a
// copy of the files the engine's Dalton reader reads, under seed/, listed in
// the manifest's seed block by name, size and SHA-256. One implementation
// for the molresponse stage (S5) and the moldft stage (D-112).
#include "Sha256.hpp"
#include <apps/molresponse/solvers/dalton_import.hpp>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>

namespace named_stages {

// What a seed reads: a molresponse stage's excited states and driven legs
// read molden.inp and RSPVEC; a moldft stage's orbitals read molden.inp only.
enum class SeedUse { response, ground_state };

// The seed block of a stage seeded from `dir`: each file the engine's reader
// locates there that this use reads (molden.inp; RSPVEC for a response seed;
// a Dalton .out when present), with its source, size and SHA-256. The reader
// is the first validator: it refuses an ambiguous directory; loose files only
// (no tarball extraction on the named path); a .out naming another method
// than HF is refused here. A ground-state seed also records the method the
// .out names, and that it could not be checked when there is none. Rank 0.
[[nodiscard]] inline nlohmann::json
dalton_seed_block(const std::string& dir, SeedUse use = SeedUse::response) {
  const auto found
      = molresponse_v3::locate_dalton_dir(dir, "", {}, {}, {}, false);
  if (found.method != "unknown" && found.method != "HF")
    throw std::runtime_error("Dalton seed: method " + found.method
                             + " (only HF seeds are validated)");
  nlohmann::json files = nlohmann::json::array();
  for (const auto& path :
       {found.molden_path,
        use == SeedUse::response ? found.rspvec_path : std::string(),
        found.out_path})
    if (!path.empty())
      files.push_back(
          {{"path", "seed/" + std::filesystem::path(path).filename().string()},
           {"source", path},
           {"bytes", std::filesystem::file_size(path)},
           {"sha256", sha256_file(path)}});
  nlohmann::json block = {{"dalton_dir", dir}, {"files", std::move(files)}};
  if (use == SeedUse::ground_state) {
    block["method"] = found.method;
    if (found.out_path.empty())
      block["method_check"]
          = "no Dalton .out in the directory: the method could not be "
            "checked (only HF seeds are validated)";
  }
  return block;
}

// The copies a generation holds against its seed block: present, size and
// SHA-256 as listed; then, for a response seed, the reader on the copies (a
// ground-state seed's copies hold molden.inp alone, which the reader does not
// accept without RSPVEC; their sums match the block made from the validated
// source). Returns what is wrong (empty when nothing is). Rank 0.
[[nodiscard]] inline std::string
check_seed_files(const std::filesystem::path& generation,
                 const nlohmann::json& seed, SeedUse use = SeedUse::response) {
  namespace fs = std::filesystem;
  std::string problems;
  const auto note = [&](const std::string& what) {
    problems += (problems.empty() ? "" : "; ") + what;
  };
  for (const auto& f : seed.at("files")) {
    const auto copy = generation / f.at("path").get<std::string>();
    if (!fs::is_regular_file(fs::symlink_status(copy)))
      note(f.at("path").get<std::string>() + " missing");
    else if (fs::file_size(copy) != f.at("bytes").get<std::uintmax_t>()
             || sha256_file(copy.string()) != f.at("sha256").get<std::string>())
      note(f.at("path").get<std::string>() + " altered");
  }
  if (problems.empty() && use == SeedUse::response) {
    try {
      (void)dalton_seed_block((generation / "seed").string());
    } catch (const std::exception& e) {
      note(e.what());
    }
  }
  return problems.empty() ? std::string{} : "Dalton seed files: " + problems;
}

// A generation directory's copies of its seed files, from their sources, then
// checked against the block; returns what check_seed_files returns. Rank 0,
// inside a publish's write_archives, so a generation without them stays
// unpublished (D-052).
[[nodiscard]] inline std::string
copy_seed_files(const std::filesystem::path& generation,
                const nlohmann::json& seed, SeedUse use = SeedUse::response) {
  std::filesystem::create_directories(generation / "seed");
  for (const auto& f : seed.at("files"))
    std::filesystem::copy_file(f.at("source").get<std::string>(),
                               generation / f.at("path").get<std::string>());
  return check_seed_files(generation, seed, use);
}

// The manifest's seed block from a publish sample's: the files by their path
// under seed/ (their sources move to source_files). The caller adds what the
// engine recorded seeding from (seeded_from).
[[nodiscard]] inline nlohmann::json
manifest_seed_block(const nlohmann::json& sample_seed) {
  auto seed = sample_seed;
  for (auto& f : seed.at("files"))
    f.erase("source");
  seed["source_files"] = sample_seed.at("files");
  return seed;
}

// A manifest's seed block lists the copied seed files by path under seed/,
// size and SHA-256 (the files themselves are checked where the generation's
// directory is read: at publish and before a load). Returns what is wrong
// with its form (empty when nothing is).
[[nodiscard]] inline std::string
seed_block_form_problem(const nlohmann::json& seed) {
  if (!seed.contains("files") || !seed.at("files").is_array()
      || seed.at("files").empty())
    return "Dalton seed block without files";
  for (const auto& f : seed.at("files"))
    if (!(f.at("path").get<std::string>().starts_with("seed/")
          && f.at("bytes").is_number_unsigned()
          && f.at("sha256").get<std::string>().size() == 64))
      return "Dalton seed block: malformed file entry";
  return {};
}

} // namespace named_stages
