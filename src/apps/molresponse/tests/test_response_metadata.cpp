// =========================================================================
// Round-trip test for ResponseMetadata (Inc 13b).
//
// Pure C++ — no MPI, no MADNESS World. Covers:
//   1. fresh load_or_create seeds the schema
//   2. upsert protocol / FD state / ES bundle / property
//   3. save -> reload preserves everything bit-for-bit
//   4. atomic write leaves no *.tmp lingering
//   5. set_protocol upsert: last write wins
//   6. add_property appends to an array (multiple records per protocol)
//   7. freq_key formatting matches doc 13 ("f%.5f" -> "f0.05700")
// =========================================================================

#include "../solvers/gs_fingerprint.hpp"
#include "../solvers/response_metadata.hpp"

#include <unistd.h>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>

namespace {

int failed = 0;

#define EXPECT(cond, label)                                                \
  do {                                                                     \
    if (cond) { std::printf("  [PASS]  %s\n", label); }                    \
    else      { std::printf("  [FAIL]  %s\n", label); ++failed; }          \
  } while (0)

} // namespace

int main() {
  using molresponse_v3::ResponseMetadata;

  // Use a process-unique scratch path so parallel CI runs don't collide.
  const std::filesystem::path tmp =
      std::filesystem::temp_directory_path() /
      ("rmeta_test_" + std::to_string(::getpid()));
  std::filesystem::create_directories(tmp);
  const std::string path = (tmp / "response_metadata.json").string();
  std::filesystem::remove(path);

  // ---- freq_key formatting -------------------------------------------------
  std::printf("=== freq_key ===\n");
  EXPECT(ResponseMetadata::freq_key(0.057)  == "f0.05700", "freq_key(0.057)  == f0.05700");
  EXPECT(ResponseMetadata::freq_key(0.0)    == "f0.00000", "freq_key(0.0)    == f0.00000");
  EXPECT(ResponseMetadata::freq_key(0.4811) == "f0.48110", "freq_key(0.4811) == f0.48110");

  // ---- fresh load: schema is seeded ----------------------------------------
  std::printf("=== fresh load_or_create ===\n");
  {
    auto m = ResponseMetadata::load_or_create(path);
    const auto &j = m.json();
    EXPECT(j["schema_version"] == 1,             "schema_version == 1");
    EXPECT(j.contains("protocols"),              "protocols/ present");
    EXPECT(j.contains("fd_states"),              "fd_states/ present");
    EXPECT(j.contains("excited_states"),         "excited_states/ present");
    EXPECT(j.contains("properties"),             "properties/ present");
  }

  // ---- upsert + save -------------------------------------------------------
  std::printf("=== upsert + save ===\n");
  {
    auto m = ResponseMetadata::load_or_create(path);
    m.set_protocol("1e-04_k6", 1e-4, 6, 0);
    m.set_protocol("1e-06_k8", 1e-6, 8, 1);

    nlohmann::json fd_entry = {
        {"freq",         0.057},
        {"type",         "full"},
        {"shell",        "closed_shell"},
        {"converged",    true},
        {"iter",         8},
        {"bsh_residual", 1.2e-7},
        {"archive",      "dipole_x__1e-06_k8__f0.05700"}};
    m.set_fd_state("dipole_x", "1e-06_k8",
                   ResponseMetadata::freq_key(0.057), fd_entry);

    nlohmann::json es_entry = {
        {"type",             "tda"},
        {"shell",            "closed_shell"},
        {"n_roots",          2},
        {"bundle_dir",       "es_bundle__1e-06_k8"},
        {"converged",        true},
        {"slot_permutation", nlohmann::json::array({0, 1})},
        {"roots", nlohmann::json::array({
            {{"stable_index", 0}, {"root_id", "es_root_0000"}, {"slot", 0},
             {"omega", 0.468}, {"display_name", "S1"}},
            {{"stable_index", 1}, {"root_id", "es_root_0001"}, {"slot", 1},
             {"omega", 0.481}, {"display_name", "S2"}}})}};
    m.set_es_bundle("1e-06_k8", es_entry);

    m.add_property("resonant_raman", "1e-06_k8",
                   {{"es_root_id", "es_root_0001"},
                    {"fd_freq",    0.481},
                    {"value",      0.0}});

    m.save();
    EXPECT(std::filesystem::exists(path),                  "file written");
    EXPECT(!std::filesystem::exists(path + ".tmp"),        "no .tmp lingers");
  }

  // ---- reload: every field round-trips -------------------------------------
  std::printf("=== reload round-trip ===\n");
  {
    auto m = ResponseMetadata::load_or_create(path);
    const auto &j = m.json();

    EXPECT(j["protocols"]["1e-06_k8"]["thresh"] == 1e-6,           "protocol thresh");
    EXPECT(j["protocols"]["1e-06_k8"]["k"]      == 8,              "protocol k");
    EXPECT(j["protocols"]["1e-06_k8"]["index"]  == 1,              "protocol index");

    const auto &fd = j["fd_states"]["dipole_x"]["1e-06_k8"]["f0.05700"];
    EXPECT(fd["freq"]      == 0.057,                               "fd freq");
    EXPECT(fd["converged"] == true,                                "fd converged");
    EXPECT(fd["archive"]   == "dipole_x__1e-06_k8__f0.05700",      "fd archive");

    const auto &es = j["excited_states"]["1e-06_k8"];
    EXPECT(es["n_roots"]                       == 2,               "es n_roots");
    EXPECT(es["bundle_dir"]                    == "es_bundle__1e-06_k8", "es bundle_dir");
    EXPECT(es["slot_permutation"]              == nlohmann::json::array({0, 1}), "es slot_permutation");
    EXPECT(es["roots"][1]["root_id"]           == "es_root_0001", "es root[1] id");
    EXPECT(es["roots"][1]["display_name"]      == "S2",           "es root[1] display");

    const auto &props = j["properties"]["resonant_raman"]["1e-06_k8"];
    EXPECT(props.is_array() && props.size() == 1,                 "1 property record");
    EXPECT(props[0]["es_root_id"]              == "es_root_0001", "prop es_root_id");
    EXPECT(props[0]["fd_freq"]                 == 0.481,          "prop fd_freq");
  }

  // ---- upsert overwrites; add_property appends -----------------------------
  std::printf("=== upsert overwrites; properties append ===\n");
  {
    auto m = ResponseMetadata::load_or_create(path);
    // bump the index on an existing protocol
    m.set_protocol("1e-06_k8", 1e-6, 8, 99);
    // add a second property record at the same protocol
    m.add_property("resonant_raman", "1e-06_k8",
                   {{"es_root_id", "es_root_0000"},
                    {"fd_freq",    0.468},
                    {"value",      0.0}});
    m.save();

    auto m2 = ResponseMetadata::load_or_create(path);
    EXPECT(m2.json()["protocols"]["1e-06_k8"]["index"] == 99,
           "set_protocol last-write-wins");
    EXPECT(m2.json()["properties"]["resonant_raman"]["1e-06_k8"].size() == 2,
           "add_property appended (now 2 records)");
  }

  // ---- property-row identity: same identity replaces, not appends ----------
  std::printf("=== property identity upsert ===\n");
  {
    auto m = ResponseMetadata::load_or_create(path);
    // Re-assemble the es_root_0000 row (restart flow): value changes, identity
    // (es_root_id, fd_freq) does not -> REPLACE in place, still 2 records.
    m.add_property("resonant_raman", "1e-06_k8",
                   {{"es_root_id", "es_root_0000"},
                    {"fd_freq",    0.468},
                    {"value",      42.0}});
    // Alpha identity is (omega, directions); beta/raman is (A,B,C,freq_b,freq_c).
    m.add_property("alpha", "1e-06_k8",
                   {{"omega", 0.0}, {"directions", "z"}, {"converged", false}});
    m.add_property("alpha", "1e-06_k8",
                   {{"omega", 0.0}, {"directions", "z"}, {"converged", true}});
    m.add_property("alpha", "1e-06_k8",
                   {{"omega", 0.057}, {"directions", "z"}, {"converged", true}});
    // Same identity at a DIFFERENT protocol appends (history across protocols).
    m.add_property("alpha", "1e-04_k6",
                   {{"omega", 0.0}, {"directions", "z"}, {"converged", true}});
    m.save();

    auto m2 = ResponseMetadata::load_or_create(path);
    const auto &rr = m2.json()["properties"]["resonant_raman"]["1e-06_k8"];
    EXPECT(rr.size() == 2, "same-identity raman row replaced (still 2)");
    bool found = false;
    for (const auto &r : rr)
      if (r.value("es_root_id", std::string()) == "es_root_0000") {
        EXPECT(r.value("value", -1.0) == 42.0, "replaced row carries new value");
        found = true;
      }
    EXPECT(found, "es_root_0000 row present after upsert");
    const auto &al = m2.json()["properties"]["alpha"]["1e-06_k8"];
    EXPECT(al.size() == 2, "alpha: 2 identities (omega 0.0 and 0.057)");
    EXPECT(al[0].value("converged", false) == true,
           "alpha omega=0 row replaced by the converged re-run");
    EXPECT(m2.json()["properties"]["alpha"]["1e-04_k6"].size() == 1,
           "different protocol_key keeps its own history");
  }

  // ---- run_summary/dropped_work (full-replace upsert) -----------------------
  std::printf("=== dropped_work ===\n");
  {
    auto m = ResponseMetadata::load_or_create(path);
    m.set_dropped_work(nlohmann::json::array(
        {{{"id", "vbc:dipole_z__dipole_z@0.057_0.057"},
          {"reason", "prerequisites never converged"}}}));
    m.save();
    auto m2 = ResponseMetadata::load_or_create(path);
    EXPECT(m2.json()["run_summary"]["dropped_work"].size() == 1,
           "dropped_work recorded");
    m2.set_dropped_work(nlohmann::json::array());
    m2.save();
    auto m3 = ResponseMetadata::load_or_create(path);
    EXPECT(m3.json()["run_summary"]["dropped_work"].empty(),
           "a completing run clears dropped_work");
  }

  // ---- schema_version mismatch is rejected ---------------------------------
  std::printf("=== schema_version guard ===\n");
  {
    const std::string bad = (tmp / "bad.json").string();
    {
      std::ofstream o(bad);
      o << "{\"schema_version\":999,\"protocols\":{},\"fd_states\":{},"
           "\"excited_states\":{},\"properties\":{}}";
    }
    bool threw = false;
    try { (void)ResponseMetadata::load_or_create(bad); }
    catch (const std::exception &) { threw = true; }
    EXPECT(threw, "unrecognized schema_version throws");
  }

  // ---- GS-archive fingerprint (restart-safety gate) -------------------------
  std::printf("=== gs_fingerprint ===\n");
  {
    using molresponse_v3::gs_archive_fingerprint;
    using molresponse_v3::gs_archive_parts;
    using molresponse_v3::gs_fingerprint_verdict;
    using molresponse_v3::GsGateVerdict;
    using molresponse_v3::metadata_has_response_states;

    // Fake a 2-part parallel archive: <base>.00000 + <base>.00001.
    const std::string base = (tmp / "gs.restartdata").string();
    { std::ofstream a(base + ".00000", std::ios::binary); a << "orbitals-A"; }
    { std::ofstream b(base + ".00001", std::ios::binary); b << "part-two"; }

    EXPECT(gs_archive_parts(base).size() == 2, "2 archive parts discovered");
    const auto fp1 = gs_archive_fingerprint(base);
    EXPECT(fp1.nparts == 2,               "fingerprint nparts == 2");
    EXPECT(fp1.bytes == 18,               "fingerprint bytes == 18");
    EXPECT(fp1.hex.size() == 16,          "fingerprint hex is 16 chars");
    EXPECT(gs_archive_fingerprint(base).hex == fp1.hex,
           "fingerprint is deterministic");

    // A single flipped byte (phase-flip stand-in) changes the fingerprint.
    { std::ofstream a(base + ".00000", std::ios::binary); a << "orbitals-B"; }
    const auto fp2 = gs_archive_fingerprint(base);
    EXPECT(fp2.hex != fp1.hex, "changed archive bytes -> changed fingerprint");

    // Missing archive throws.
    bool threw = false;
    try { (void)gs_archive_fingerprint((tmp / "nope").string()); }
    catch (const std::exception &) { threw = true; }
    EXPECT(threw, "missing archive throws");

    // Verdicts against metadata: fresh dir -> FreshDir; states without a
    // stamp -> MissingStamp; stamped -> Match / Mismatch.
    const std::string gpath = (tmp / "gs_meta.json").string();
    auto m = ResponseMetadata::load_or_create(gpath);
    EXPECT(!metadata_has_response_states(m.json()), "fresh dir has no states");
    EXPECT(gs_fingerprint_verdict(m.json(), fp1.hex) == GsGateVerdict::FreshDir,
           "fresh dir -> FreshDir");

    m.set_fd_state("dipole_z", "1e-06_k8", "f0.00000", {{"status", "Converged"}});
    EXPECT(metadata_has_response_states(m.json()), "fd_state counts as states");
    EXPECT(gs_fingerprint_verdict(m.json(), fp1.hex) == GsGateVerdict::MissingStamp,
           "states without stamp -> MissingStamp");

    m.set_ground_state(base, fp1.hex, fp1.bytes, fp1.nparts);
    EXPECT(gs_fingerprint_verdict(m.json(), fp1.hex) == GsGateVerdict::Match,
           "matching stamp -> Match");
    EXPECT(gs_fingerprint_verdict(m.json(), fp2.hex) == GsGateVerdict::Mismatch,
           "different archive -> Mismatch");

    // The stamp round-trips through save/load like every other block.
    m.save();
    auto m2 = ResponseMetadata::load_or_create(gpath);
    EXPECT(m2.ground_state_fingerprint() == fp1.hex,
           "ground_state stamp survives save/reload");
    EXPECT(m2.json()["ground_state"]["nparts"] == 2,
           "ground_state block carries nparts");
  }

  // ---- stale shard sweep (F2 restart safety) --------------------------------
  std::printf("=== merge_stale_state_shards ===\n");
  {
    const auto cdir = tmp / "shard_calc";
    std::filesystem::create_directories(cdir);

    // Canonical with one FD state; two stranded shards (as an interrupted
    // 3-subworld run would leave them — note g2 wrote nothing and left none).
    {
      auto c = ResponseMetadata::load_or_create(
          (cdir / "response_metadata.json").string());
      c.set_fd_state("dipole_x", "1e-04_k6", "f0.00000", {{"status", "Converged"}});
      c.save();
      auto s0 = ResponseMetadata::load_or_create(
          (cdir / "response_metadata.group0.json").string());
      s0.set_fd_state("dipole_y", "1e-04_k6", "f0.00000", {{"status", "Converged"}});
      s0.set_protocol("1e-04_k6", 1e-4, 6, 0);
      s0.save();
      auto s1 = ResponseMetadata::load_or_create(
          (cdir / "response_metadata.group1.json").string());
      s1.set_fd_state("dipole_z", "1e-04_k6", "f0.00000", {{"status", "Running"}});
      s1.set_vbc_state("vbc_zz", "1e-04_k6", {{"status", "Converged"}});
      s1.save();
      // Decoys the sweep must NOT eat: non-numeric gid, and a foreign json.
      std::ofstream((cdir / "response_metadata.groupX.json").string())
          << "{\"schema_version\":1}";
      std::ofstream((cdir / "other.json").string()) << "{}";
    }

    const int merged =
        ResponseMetadata::merge_stale_state_shards(cdir.string());
    EXPECT(merged == 2, "two stranded shards merged");
    auto c = ResponseMetadata::load_or_create(
        (cdir / "response_metadata.json").string());
    const auto &fd = c.json()["fd_states"];
    EXPECT(fd.contains("dipole_x") && fd.contains("dipole_y") &&
           fd.contains("dipole_z"), "canonical holds pre-existing + both shards");
    EXPECT(c.json()["vbc_states"].contains("vbc_zz"), "shard vbc_states merged");
    EXPECT(c.json()["protocols"].contains("1e-04_k6"), "shard protocol unioned");
    EXPECT(!std::filesystem::exists(cdir / "response_metadata.group0.json") &&
           !std::filesystem::exists(cdir / "response_metadata.group1.json"),
           "merged shards removed");
    EXPECT(std::filesystem::exists(cdir / "response_metadata.groupX.json"),
           "non-shard decoy untouched");
    EXPECT(ResponseMetadata::merge_stale_state_shards(cdir.string()) == 0,
           "second sweep is a no-op (idempotent)");
  }

  // ---- kill-mid-save simulation: stranded truncated tmp -------------------
  // A writer killed between opening <path>.tmp and the rename leaves a
  // truncated tmp behind. It must not affect reads (the good file is intact),
  // and the next save must replace it and install atomically.
  std::printf("=== stranded truncated tmp ===\n");
  {
    const auto cdir = tmp / "atomic";
    std::filesystem::create_directories(cdir);
    const std::string mpath = (cdir / "response_metadata.json").string();
    {
      auto m = ResponseMetadata::load_or_create(mpath);
      m.set_protocol("1e-06_k8", 1e-6, 8, 0);
      m.save();
    }
    {  // simulate the kill: garbage where the tmp lives
      std::ofstream junk(mpath + ".tmp");
      junk << "{\"schema_version\": 1, \"proto";  // truncated JSON
    }
    {
      auto m = ResponseMetadata::load_or_create(mpath);   // reads the GOOD file
      EXPECT(m.json()["protocols"].contains("1e-06_k8"),
             "stranded tmp does not shadow the good index");
      m.set_protocol("1e-07_k10", 1e-7, 10, 1);
      m.save();                                            // replaces the tmp
    }
    EXPECT(!std::filesystem::exists(mpath + ".tmp"),
           "re-save consumed/replaced the stranded tmp");
    auto m = ResponseMetadata::load_or_create(mpath);
    EXPECT(m.json()["protocols"].contains("1e-06_k8") &&
           m.json()["protocols"].contains("1e-07_k10"),
           "index intact after kill-simulated save cycle");
  }

  // ---- drop_fd_rung: a redefined protocol key drops its FD entries ----------
  std::printf("=== drop_fd_rung ===\n");
  {
    auto m = ResponseMetadata::load_or_create((tmp / "drop.json").string());
    for (const char *pert : {"dipole_x", "dipole_y"})
      for (const char *key : {"1e-04_k6", "1e-06_k8"})
        m.set_fd_state(pert, key, ResponseMetadata::freq_key(0.0),
                       {{"freq", 0.0}, {"converged", true}});
    m.drop_fd_rung("1e-06_k8");
    const auto &fd = m.json()["fd_states"];
    EXPECT(!fd["dipole_x"].contains("1e-06_k8") && !fd["dipole_y"].contains("1e-06_k8"),
           "every channel's entries at the redefined key are gone");
    EXPECT(fd["dipole_x"].contains("1e-04_k6") && fd["dipole_y"].contains("1e-04_k6"),
           "other keys are kept");
  }

  // ---- a leg saved again drops the quadratic sources built from it ---------
  std::printf("=== drop_vbc_built_from ===\n");
  {
    auto m = ResponseMetadata::load_or_create((tmp / "vbc.json").string());
    const std::string f0 = ResponseMetadata::freq_key(0.0);
    const std::string f2 = ResponseMetadata::freq_key(0.02);
    const nlohmann::json leg = {{"freq", 0.02}, {"converged", true}, {"iter", 5}};
    for (const char *key : {"1e-04_k6", "1e-06_k8"}) {
      m.set_fd_state("dipole_z", key, f2, leg);
      m.set_fd_state("nuc_1_z", key, f0, leg);
      for (const char *id : {"vbc:dipole_z__dipole_z@f0.02000_f0.02000",
                             "vbc:dipole_x__dipole_z@f0.02000_f0.02000",
                             "vbc:dipole_z__nuc_1_z@f0.02000_f0.00000",
                             "vbc:dipole_z__dipole_z@f0.00000_f0.00000",
                             "vbc:dipole_x__nuc_1_z@f0.02000_f0.00000"})
        m.set_vbc_state(id, key, {{"converged", true}});
    }
    m.set_fd_state("dipole_z", "1e-06_k8", f2, leg);  // identical: a re-merge
    EXPECT(m.json()["vbc_states"]["vbc:dipole_z__dipole_z@f0.02000_f0.02000"]
               .contains("1e-06_k8"),
           "an identical entry keeps the sources");
    auto resolved = leg;
    resolved["iter"] = 3;
    m.set_fd_state("dipole_z", "1e-06_k8", f2, resolved);  // solved again
    const auto &vbc = m.json()["vbc_states"];
    for (const char *id : {"vbc:dipole_z__dipole_z@f0.02000_f0.02000",
                           "vbc:dipole_x__dipole_z@f0.02000_f0.02000",
                           "vbc:dipole_z__nuc_1_z@f0.02000_f0.00000"})
      EXPECT(!vbc[id].contains("1e-06_k8") && vbc[id].contains("1e-04_k6"),
             "a source built from the leg is dropped at that rung only");
    EXPECT(vbc["vbc:dipole_z__dipole_z@f0.00000_f0.00000"].contains("1e-06_k8"),
           "the same channel at another frequency is not an input");
    EXPECT(vbc["vbc:dipole_x__nuc_1_z@f0.02000_f0.00000"].contains("1e-06_k8"),
           "a source without the leg is kept");
    m.drop_fd_rung("1e-04_k6");
    EXPECT(!vbc["vbc:dipole_x__nuc_1_z@f0.02000_f0.00000"].contains("1e-04_k6"),
           "a redefined key drops its sources too");
  }

  // ---- a removed driven node is no longer named ----------------------------
  std::printf("=== drop_fd_node ===\n");
  {
    auto m = ResponseMetadata::load_or_create((tmp / "prune.json").string());
    const std::string old_f = ResponseMetadata::freq_key(0.16205);
    const std::string new_f = ResponseMetadata::freq_key(0.16139);
    const auto entry = [](const std::string &node) {
      return nlohmann::json{{"converged", true}, {"node", node}};
    };
    m.set_fd_state("dipole_x", "1e-04_k6", old_f, entry("states/response/dipole_x-f0.16205"));
    m.set_fd_state("dipole_x", "1e-04_k6", new_f, entry("states/response/dipole_x-f0.16139"));
    m.set_fd_state("dipole_y", "1e-04_k6", old_f, entry("states/response/dipole_y-f0.16205"));
    m.set_vbc_state("vbc:dipole_x__dipole_x@f0.16205_f0.16205", "1e-04_k6", {{"converged", true}});
    m.set_vbc_state("vbc:dipole_x__dipole_x@f0.16139_f0.16139", "1e-04_k6", {{"converged", true}});
    m.drop_fd_node("states/response/dipole_x-f0.16205");
    m.drop_fd_node("states/response/dipole_y-f0.16205");
    m.drop_fd_node("states/response/dipole_z-f0.16205");  // never recorded
    const auto &fd = m.json()["fd_states"];
    EXPECT(!fd["dipole_x"]["1e-04_k6"].contains(old_f), "the removed node's entry is gone");
    EXPECT(fd["dipole_x"]["1e-04_k6"].contains(new_f), "another node's entry is kept");
    EXPECT(!fd.contains("dipole_y"), "a channel left empty is dropped");
    const auto &vbc = m.json()["vbc_states"];
    EXPECT(!vbc["vbc:dipole_x__dipole_x@f0.16205_f0.16205"].contains("1e-04_k6"),
           "a source built from the removed node is dropped");
    EXPECT(vbc["vbc:dipole_x__dipole_x@f0.16139_f0.16139"].contains("1e-04_k6"),
           "a source built from the kept node stays");
  }

  // ---- best_usable_fd_source_key: an aliased entry is not a source --------
  std::printf("=== best_usable_fd_source_key, aliased key ===\n");
  {
    const std::string f0 = ResponseMetadata::freq_key(0.0);
    nlohmann::json meta = {
        {"protocols", {{"1e-04_k6", {{"thresh", 1e-4}, {"k", 6}}},
                       {"1e-06_k8", {{"thresh", 1.01e-6}, {"k", 8}}}}},
        {"fd_states", {{"dipole_z", {{"1e-04_k6", {{f0, {{"converged", true}}}}},
                                     {"1e-06_k8", {{f0, {{"converged", true}}}}}}}}}};
    // The saved 1.01e-6 entry is coarser than a 1e-6 target, so before this
    // rule it was the exact source and the state resumed from it.
    EXPECT(molresponse_v3::best_usable_fd_source_key(meta, "dipole_z", f0, 1e-6, 8, "1e-06_k8") ==
               "1e-04_k6",
           "an entry saved under another threshold of the same key is skipped");
    EXPECT(molresponse_v3::best_usable_fd_source_key(meta, "dipole_z", f0, 1.01e-6, 8, "1e-06_k8") ==
               "1e-06_k8",
           "the same threshold is still the exact source");
  }

  // ---- best_usable_es_source_key: an aliased bundle is not a source -------
  std::printf("=== best_usable_es_source_key, aliased key ===\n");
  {
    const nlohmann::json bundle = {{"type", "tda"}, {"shell", "closed_shell"},
                                   {"n_roots", 1}, {"converged", true}};
    nlohmann::json meta = {
        {"protocols", {{"1e-04_k6", {{"thresh", 1e-4}, {"k", 6}}},
                       {"1e-06_k8", {{"thresh", 1.01e-6}, {"k", 8}}}}},
        {"excited_states", {{"1e-04_k6", bundle}, {"1e-06_k8", bundle}}}};
    EXPECT(molresponse_v3::best_usable_es_source_key(meta, "tda", "closed_shell", 1,
                                                     1e-6, 8, "1e-06_k8") == "1e-04_k6",
           "a bundle saved under another threshold of the same key is skipped");
    EXPECT(molresponse_v3::best_usable_es_source_key(meta, "tda", "closed_shell", 1,
                                                     1.01e-6, 8, "1e-06_k8") == "1e-06_k8",
           "the same threshold is still the exact source");
  }

  std::filesystem::remove_all(tmp);
  std::printf("\n%s: %d failure(s)\n",
              failed == 0 ? "ALL PASS" : "FAILED", failed);
  return failed == 0 ? 0 : 1;
}
