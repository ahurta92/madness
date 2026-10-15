/// Storage checks deliberately use real archives and filesystem failures, no chemistry solve.
#include <madness/chem/StageCheckpoint.h>
#include <madness/world/MADworld.h>
#include <madness/world/parallel_archive.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <csignal>
#include <sys/resource.h>

using namespace madness;
namespace cp = madness::checkpoint;
namespace fs = std::filesystem;
namespace {
int failures = 0;
int checks = 0;
void check(bool ok, const std::string& message) {
    ++checks;
    if (!ok) { ++failures; std::cerr << "FAILED: " << message << '\n'; }
}
template <class F> void rejects(F&& f, const std::string& message) {
    bool threw = false;
    try { f(); } catch (const std::exception&) { threw = true; }
    check(threw, message);
}
void write_text(const fs::path& path, const std::string& text) {
    std::ofstream out(path, std::ios::binary);
    out.exceptions(std::ios::failbit | std::ios::badbit);
    out << text;
    out.close();
}
cp::Manifest manifest(const std::string& generation) {
    cp::Manifest m;
    m.stage = "hf"; m.method = "nemo"; m.generation = generation;
    m.physical = {{"geometry", "He"}, {"occupations", {2.0}}};
    m.numerical = {{"k", 6}, {"thresh", 1e-5}};
    m.reference = {{"stage", "reference"}, {"generation", "ref-1"}};
    m.working_orbitals = {{"ordering", {0}}, {"frozen", 0}};
    m.seed = {{"stage", "seed"}, {"generation", "seed-2"}};
    m.position = {{"iteration", 3}, {"rung", 0}, {"measurement_phase", "post-update"}};
    m.metrics = {{"energy_change", 0.02}};
    m.targets = {{"econv", 1e-6}};
    m.stop_reason = "stalled";
    m.properties = {{{"generation", generation}, {"energy", -2.8}, {"units", "hartree"}}};
    return m;
}
// Use the actual runtime layout: cookie on every shard, nio only on shard zero.
cp::Manifest save(World& world, const fs::path& stage, const std::string& generation, int writers = 0) {
    if (writers == 0) writers = world.size();
    auto m = manifest(generation);
    const auto dir = stage / generation;
    if (world.rank() == 0) fs::create_directories(dir);
    world.gop.fence();
    archive::ParallelOutputArchive out(world, (dir / "orbitals").c_str(), writers);
    int value = 42 + world.rank();
    if (out.is_io_node()) out.local_archive() & value;
    check(cp::close_archive(out).empty(), "successful archive closes cleanly on writers and nonwriters");
    world.gop.fence();
    if (world.rank() == 0)
        m.archives.push_back(cp::inspect_archive(dir, "orbitals", "orbitals", {}, writers));
    return m;
}
void read_payload(World& world, const fs::path& stage, const cp::Manifest& m) {
    archive::ParallelInputArchive in(world, (stage / m.generation / "orbitals").c_str());
    if (in.is_io_node()) {
        int value = 0;
        in.local_archive() & value;
        check(value == 42 + world.rank(), "selected generation loads the original payload");
    }
    in.close();
}
void storage(World& world, const fs::path& stage) {
    auto first = save(world, stage, "g1");
    cp::publish(world, stage, first);
    auto selected = cp::read_published(world, stage);
    check(selected.generation == "g1", "select published generation");
    check(selected.stop_reason == "stalled", "storage validity does not require convergence");
    check(!selected.metrics.contains("density_change"), "missing metrics stay missing");
    if (world.rank() == 0)
        check(cp::to_json_record(cp::read_manifest(stage / "g1")) == cp::to_json_record(first),
              "all manifest metadata round trips");
    read_payload(world, stage, selected);

    auto second = save(world, stage, "g2");
    check(cp::read_published(world, stage).generation == "g1", "ignore unpublished newer files");
    // A failure to create the manifest temporary file leaves latest untouched.
    if (world.rank() == 0) fs::create_directory(stage / "g2" / ".manifest.tmp");
    world.gop.fence();
    rejects([&] { cp::publish(world, stage, second); }, "collective failure before manifest publication");
    check(cp::read_published(world, stage).generation == "g1", "old pointer survives manifest failure");
    if (world.rank() == 0) fs::remove(stage / "g2" / ".manifest.tmp");
    world.gop.fence();
    // Fail after manifest publication, before the atomic latest replacement.
    if (world.rank() == 0) fs::create_directory(stage / ".latest-g2.tmp");
    world.gop.fence();
    rejects([&] { cp::publish(world, stage, second); }, "collective failure before pointer replacement");
    check(cp::read_published(world, stage).generation == "g1", "old pointer survives pointer failure");
    read_payload(world, stage, cp::read_published(world, stage));
    if (world.rank() == 0) check(fs::is_regular_file(stage / "g2" / "manifest.json"),
                               "pointer failure occurs after manifest publication");

    auto third = save(world, stage, "g3");
    rejects([&] { cp::publish(world, stage, third,
                  world.rank() == world.size() - 1 ? "archive close failed" : ""); },
            "a writer error on any rank prevents publication on every rank");
    check(cp::read_published(world, stage).generation == "g1", "writer failure retains old generation");
    cp::publish(world, stage, third);
    check(cp::read_published(world, stage).generation == "g3", "advance pointer after successful publication");
    if (world.rank() == 0) check(cp::validate_files(stage / "g1", first).valid(), "previous generation retained");
    rejects([&] { cp::publish(world, stage, third); }, "never overwrite a published manifest");

    // Rank-zero-only filesystem/schema tests; no collectives in the callable.
    if (world.rank() == 0) {
        const auto dir = stage / "g3";
        auto empty = manifest("empty");
        check(!cp::validate_files(stage / "empty", empty).valid(),
              "empty inventory still requires a real generation directory");
        fs::create_directory(stage / "empty");
        check(cp::validate_files(stage / "empty", empty).valid(),
              "existing empty generation has a valid empty inventory");
        auto bad = third;
        bad.archives[0].indices = {-1};
        check(!cp::validate_files(dir, bad).valid(), "negative component index rejected before publication");
        bad = third;
        bad.version = 99;
        check(!cp::validate_files(dir, bad).valid(), "unsupported format rejected");
        bad = third;
        bad.archives[0].shards[0].path = "../g1/orbitals.00000";
        check(!cp::validate_files(dir, bad).valid(), "escaping member rejected");
        bad.archives[0].shards[0].path = "/etc/passwd";
        check(!cp::validate_files(dir, bad).valid(), "absolute member rejected");
        bad = third;
        bad.properties[0]["generation"] = "g1";
        check(!cp::validate_files(dir, bad).valid(), "properties bound to generation");
        bad = third;
        bad.archives.push_back(bad.archives.front());
        check(!cp::validate_files(dir, bad).valid(), "duplicate inventory rejected");
        auto record = cp::to_json_record(third);
        record["version"] = 99;
        write_text(dir / "manifest.json", record.dump());
        rejects([&] { cp::read_manifest(dir); }, "unsupported on-disk manifest rejected");
        record["version"] = 1;
        record["archives"][0]["writers"] = -1;
        write_text(dir / "manifest.json", record.dump());
        rejects([&] { cp::read_manifest(dir); }, "negative writer count rejected");
        write_text(dir / "manifest.json", cp::to_json_record(third).dump());

        // Independent two-shard fixture also exercises nonzero shards in a serial run.
        const auto fixture = stage / "fixture";
        fs::create_directory(fixture);
        { archive::BinaryFstreamOutputArchive out((fixture / "pair.00000").string());
          int writers = 2; out & writers; out.close(); }
        { archive::BinaryFstreamOutputArchive out((fixture / "pair.00001").string());
          int payload = 123; out & payload; out.close(); }
        auto multi = manifest("fixture");
        multi.archives.push_back(cp::inspect_archive(fixture, "pair", "pair", {0, 1}, 2));
        check(cp::validate_files(fixture, multi, 2).valid(), "two-shard inventory valid with two readers");
        check(!cp::validate_files(fixture, multi, 1).valid(), "fewer readers than writers rejected");
        bad = multi;
        bad.archives[0].writers = 1;
        bad.archives[0].shards.pop_back();
        check(!cp::validate_files(fixture, bad, 2).valid(), "header/manifest writer mismatch rejected");
        const auto shard = fixture / "pair.00001";
        const auto original = fixture / "backup";
        fs::rename(shard, original);
        auto missing = cp::validate_files(fixture, multi, 2);
        check(!missing.valid() && missing.message().find("pair.00001") != std::string::npos,
              "missing nonzero shard diagnosed by member");
        fs::create_symlink(original, shard);
        check(!cp::validate_files(fixture, multi, 2).valid(), "symlink member rejected");
        fs::remove(shard);
        fs::rename(original, shard);
        { std::fstream bytes(shard, std::ios::binary | std::ios::in | std::ios::out);
          bytes.seekp(-1, std::ios::end); bytes.put('\x7f'); }
        auto corrupt = cp::validate_files(fixture, multi, 2);
        check(!corrupt.valid() && corrupt.message().find("pair.00001") != std::string::npos,
              "same-size corruption of nonzero shard diagnosed");
        // Recompute the checksum to isolate header checking from integrity checks.
        const auto zero = fixture / "pair.00000";
        write_text(zero, std::string("archive\0", 8));
        bad = multi;
        bad.archives[0].shards[0].bytes = fs::file_size(zero);
        bad.archives[0].shards[0].checksum = cp::checksum_file(zero);
        auto truncated = cp::validate_files(fixture, bad, 2);
        check(!truncated.valid() && truncated.message().find("pair.00000") != std::string::npos,
              "truncated writer header rejected even with matching integrity metadata");
        write_text(stage / "crc-vector", "123456789");
        check(cp::checksum_file(stage / "crc-vector") == "crc64-ecma182:6c40df5f0b497347",
              "CRC-64 standard check vector");
        write_text(stage / "crc-empty", "");
        check(cp::checksum_file(stage / "crc-empty") == "crc64-ecma182:0000000000000000",
              "CRC-64 empty file");
    }
    world.gop.fence();
    // Damage the nonzero writer's member in MPI runs: all ranks must reject it
    // before constructing a collective archive, even though shard zero survives.
    const auto published = cp::read_published(world, stage);
    const auto member = stage / "g3" / published.archives[0].shards.back().path;
    const auto backup = stage / "member-backup";
    if (world.rank() == 0) {
        fs::copy_file(member, backup);
        std::fstream bytes(member, std::ios::binary | std::ios::in | std::ios::out);
        bytes.seekg(-1, std::ios::end);
        const auto value = static_cast<char>(bytes.get() ^ 1);
        bytes.seekp(-1, std::ios::end);
        bytes.put(value);
    }
    world.gop.fence();
    rejects([&] { cp::read_published(world, stage); }, "corrupt last-writer shard rejected collectively");
    if (world.rank() == 0) fs::rename(backup, member);
    world.gop.fence();
    check(cp::read_published(world, stage).generation == "g3", "restored generation remains readable");
    if (world.rank() == 0) fs::remove(member);
    world.gop.fence();
    rejects([&] { cp::read_published(world, stage); }, "missing last-writer shard rejected collectively");

    const auto single_stage = stage / "single-writer";
    auto single = save(world, single_stage, "g1", 1);
    cp::publish(world, single_stage, single);
    read_payload(world, single_stage, cp::read_published(world, single_stage));
}
// A real buffered stream failure, including a nonzero writer in the MPI run.
// Limit only the local write/close section; restore process settings before any
// further filesystem work or collective operation.
void failed_write(World& world, const fs::path& stage) {
    auto old = save(world, stage, "old");
    cp::publish(world, stage, old);
    const auto dir = stage / "short-write";
    if (world.rank() == 0) fs::create_directory(dir);
    world.gop.fence();
    archive::ParallelOutputArchive output(world, (dir / "orbitals").c_str(), world.size());
    output.local_archive().flush();
    std::string error;
    struct rlimit original{};
    struct sigaction previous{}, ignored{};
    const bool failing_rank = world.rank() == world.size() - 1;
    if (failing_rank) {
        MADNESS_CHECK(getrlimit(RLIMIT_FSIZE, &original) == 0);
        ignored.sa_handler = SIG_IGN;
        sigemptyset(&ignored.sa_mask);
        MADNESS_CHECK(sigaction(SIGXFSZ, &ignored, &previous) == 0);
        auto limited = original;
        limited.rlim_cur = 64;
        MADNESS_CHECK(setrlimit(RLIMIT_FSIZE, &limited) == 0);
    }
    std::array<char, 4096> payload{};
    output.local_archive().store(payload.data(), payload.size());
    error = cp::close_archive(output);
    if (failing_rank) {
        MADNESS_CHECK(setrlimit(RLIMIT_FSIZE, &original) == 0);
        MADNESS_CHECK(sigaction(SIGXFSZ, &previous, nullptr) == 0);
        check(!error.empty(), "checked close reports an actual short write");
    }
    world.gop.fence();
    auto next = manifest("short-write");
    if (world.rank() == 0) {
        next.archives.push_back(cp::inspect_archive(dir, "orbitals", "orbitals", {}, world.size()));
        check(next.archives[0].shards.back().bytes == 64, "failure injection actually truncated the archive");
    }
    rejects([&] { cp::publish(world, stage, next, error); }, "actual write failure rejects publication collectively");
    const auto selected = cp::read_published(world, stage);
    check(selected.generation == "old", "actual write failure preserves latest");
    read_payload(world, stage, selected);
}
} // namespace
int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    std::string root;
    if (world.rank() == 0) {
        char pattern[] = "/tmp/madness-checkpoint-XXXXXX";
        const auto* result = mkdtemp(pattern);
        if (result) root = result;
    }
    world.gop.broadcast_serializable(root, 0);
    if (root.empty()) { finalize(); return 1; }
    try { storage(world, root); failed_write(world, fs::path(root) / "write-failure"); }
    catch (const std::exception& error) { check(false, error.what()); }
    world.gop.sum(failures);
    world.gop.sum(checks);
    if (world.rank() == 0) {
        std::cout << "stage checkpoint: " << checks - failures << '/' << checks << " checks passed\n";
        fs::remove_all(root);
    }
    finalize();
    return failures ? 1 : 0;
}
