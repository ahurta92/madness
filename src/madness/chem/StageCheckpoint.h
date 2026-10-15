/// \file StageCheckpoint.h
/// \brief Immutable stage generations over MADNESS binary parallel archives.
///
/// A stage directory contains generation directories and a small `latest` JSON
/// pointer. Callers own unique generation names and the complete required-member
/// list. All ranks must close their archives with close_archive() before publish(),
/// passing its result (and any other local error) through local_error. Never modify a saved
/// generation. Readers call read_published() collectively BEFORE constructing a
/// ParallelInputArchive. The shared filesystem must be visible to rank zero.
///
/// This layer checks storage, not scientific compatibility or convergence. The
/// application owns the meaning of the JSON metadata, including metric phases,
/// reference/orbital identities and resetting acceleration history on restart.
/// CRC-64 detects accidental corruption; it does not authenticate untrusted data.
/// Concurrent mutation by outside processes, MPI process loss and power-loss
/// durability are not supported. A crashed publisher can leave `.publish-lock`;
/// readers still work, but a subsequent publisher requires manual lock recovery.
#ifndef MADNESS_CHEM_STAGECHECKPOINT_H_INCLUDED
#define MADNESS_CHEM_STAGECHECKPOINT_H_INCLUDED

#include <madness/world/world.h>
#include <madness/world/worldgop.h>
#include <madness/world/parallel_archive.h>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <fcntl.h>
#include <unistd.h>

namespace madness::checkpoint {
namespace fs = std::filesystem;
using json = nlohmann::json;

struct Shard {
    std::string path;
    std::uint64_t bytes = 0;
    std::string checksum;
};
struct Archive {
    std::string role;
    std::string basename;
    std::vector<int> indices;
    int writers = 0;
    std::vector<Shard> shards;
};
struct Manifest {
    int version = 1;
    std::string stage, method, generation;
    json physical = json::object(), numerical = json::object();
    json reference = nullptr, working_orbitals = nullptr, seed = nullptr;
    json position = json::object(), metrics = json::object(), targets = json::object();
    std::string stop_reason;
    json properties = json::array();
    std::vector<Archive> archives;
};
struct Validation {
    std::vector<std::string> diagnostics;
    [[nodiscard]] bool valid() const { return diagnostics.empty(); }
    [[nodiscard]] std::string message() const {
        std::string result;
        for (const auto& item : diagnostics) result += item + "\n";
        return result;
    }
};

namespace detail {
inline void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error("checkpoint: " + message);
}
inline bool safe_relative(const std::string& name) {
    if (name.empty() || name.find('\0') != std::string::npos) return false;
    const fs::path path(name);
    if (path.is_absolute() || path.generic_string() != name) return false;
    return std::all_of(path.begin(), path.end(), [](const auto& part) {
        return !part.empty() && part != "." && part != "..";
    });
}
inline bool generation_name(const std::string& name) {
    return safe_relative(name) && name.find('/') == std::string::npos && name.front() != '.'
           && name != "latest";
}
inline std::string shard_name(const std::string& basename, int writer) {
    std::ostringstream out;
    out << basename << '.' << std::setfill('0') << std::setw(5) << writer;
    return out.str();
}
// Check every component below the caller-owned root, not just the final file.
inline fs::path member_path(const fs::path& dir, const std::string& name) {
    require(safe_relative(name), "unsafe member path: " + name);
    require(fs::is_directory(fs::symlink_status(dir)), "generation is not a real directory: " + dir.string());
    auto path = dir;
    for (const auto& part : fs::path(name)) {
        path /= part;
        require(!fs::is_symlink(fs::symlink_status(path)), "symlink member: " + name);
    }
    return path;
}
inline std::uint64_t unsigned_integer(const json& value) {
    require(value.is_number_integer(), "expected an integer");
    if (!value.is_number_unsigned()) require(value.get<std::int64_t>() >= 0, "negative integer");
    return value.get<std::uint64_t>();
}
inline int integer(const json& value) {
    const auto n = unsigned_integer(value);
    require(n <= static_cast<std::uint64_t>(std::numeric_limits<int>::max()), "integer out of range");
    return static_cast<int>(n);
}
inline void structure(const Manifest& m) {
    require(m.version == 1, "unsupported manifest version " + std::to_string(m.version));
    require(!m.stage.empty() && !m.method.empty(), "stage and method are required");
    require(generation_name(m.generation), "unsafe generation: " + m.generation);
    require(m.physical.is_object() && m.numerical.is_object() && m.position.is_object()
            && m.metrics.is_object() && m.targets.is_object(), "settings and evidence must be objects");
    for (const auto* field : {&m.reference, &m.working_orbitals, &m.seed})
        require(field->is_null() || field->is_object(), "identity/provenance must be null or an object");
    require(m.properties.is_array(), "properties must be an array");
    for (const auto& property : m.properties)
        require(property.is_object() && property.at("generation") == m.generation,
                "property belongs to a different generation");
    std::set<std::string> basenames, members;
    for (const auto& a : m.archives) {
        require(!a.role.empty(), "archive role is required");
        require(std::all_of(a.indices.begin(), a.indices.end(), [](int index) { return index >= 0; }),
                a.basename + ": negative component index");
        require(safe_relative(a.basename), "unsafe archive basename: " + a.basename);
        require(basenames.insert(a.basename).second, "duplicate archive: " + a.basename);
        require(a.writers > 0 && a.shards.size() == static_cast<std::size_t>(a.writers),
                a.basename + ": writer count does not match shard inventory");
        for (int i = 0; i < a.writers; ++i) {
            const auto& s = a.shards[static_cast<std::size_t>(i)];
            require(safe_relative(s.path) && s.path == shard_name(a.basename, i),
                    "unsafe or inconsistent shard path: " + s.path);
            require(members.insert(s.path).second, "duplicate shard: " + s.path);
            require(s.checksum.size() == 30 && s.checksum.substr(0, 14) == "crc64-ecma182:"
                    && s.checksum.find_first_not_of("0123456789abcdef", 14) == std::string::npos,
                    s.path + ": unsupported checksum");
        }
    }
}
inline json read_json(const fs::path& path) {
    require(fs::is_regular_file(fs::symlink_status(path)), "not a regular file: " + path.string());
    std::ifstream input(path);
    require(input.good(), "cannot read " + path.string());
    return json::parse(input);
}
// Strict, bounded header read. The existing binary archive stream does not turn
// short reads into exceptions; using it here could accept an uninitialized nio.
inline int header(const fs::path& path, bool rank_zero) {
    std::ifstream input(path, std::ios::binary);
    input.exceptions(std::ios::failbit | std::ios::badbit);
    std::array<char, sizeof(ARCHIVE_COOKIE)> cookie{};
    input.read(cookie.data(), cookie.size());
    require(std::memcmp(cookie.data(), ARCHIVE_COOKIE, cookie.size()) == 0,
            path.string() + ": invalid archive cookie");
    int writers = 0;
    if (rank_zero) {
        unsigned char tag = 0;
        input.read(reinterpret_cast<char*>(&tag), sizeof(tag));
        require(tag == archive::archive_typeinfo<int>::cookie, path.string() + ": invalid writer-count tag");
        input.read(reinterpret_cast<char*>(&writers), sizeof(writers));
    }
    return writers;
}
inline void write_exclusive(const fs::path& path, const std::string& content) {
    const int fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0600);
    require(fd >= 0, "cannot create " + path.string() + ": " + std::strerror(errno));
    std::size_t offset = 0;
    while (offset < content.size()) {
        const auto written = ::write(fd, content.data() + offset, content.size() - offset);
        if (written < 0 && errno == EINTR) continue;
        if (written <= 0) {
            const std::string error = std::strerror(errno);
            ::close(fd);
            throw std::runtime_error("checkpoint: cannot write " + path.string() + ": " + error);
        }
        offset += static_cast<std::size_t>(written);
    }
    require(::close(fd) == 0, "cannot close " + path.string());
}
inline void collective_error(World& world, std::string error) {
    int source = error.empty() ? world.size() : world.rank();
    world.gop.min(source);
    if (source < world.size()) {
        world.gop.broadcast_serializable(error, source);
        throw std::runtime_error(error);
    }
}
template <class F> inline std::string on_root(World& world, F&& function) {
    std::string error, result;
    if (world.rank() == 0) {
        try { result = function(); }
        catch (const std::exception& e) { error = e.what(); }
        catch (...) { error = "checkpoint: unknown coordinator error"; }
    }
    collective_error(world, error);
    world.gop.broadcast_serializable(result, 0);
    return result;
}
} // namespace detail

/// Streaming CRC-64/ECMA-182: polynomial 0x42f0e1eba9ea3693, init/xorout 0,
/// no reflection. Fixed width and byte order independent; no external dependency.
inline std::string checksum_file(const fs::path& path) {
    static const auto table = [] {
        std::array<std::uint64_t, 256> values{};
        for (std::size_t i = 0; i < values.size(); ++i) {
            std::uint64_t v = static_cast<std::uint64_t>(i) << 56;
            for (int bit = 0; bit < 8; ++bit)
                v = (v & (std::uint64_t{1} << 63)) ? (v << 1) ^ UINT64_C(0x42f0e1eba9ea3693) : v << 1;
            values[i] = v;
        }
        return values;
    }();
    std::ifstream input(path, std::ios::binary);
    detail::require(input.good(), "cannot read " + path.string());
    std::array<char, 65536> buffer{};
    std::uint64_t crc = 0;
    while (input) {
        input.read(buffer.data(), buffer.size());
        for (std::streamsize i = 0; i < input.gcount(); ++i)
            crc = table[(crc >> 56) ^ static_cast<unsigned char>(buffer[static_cast<std::size_t>(i)])] ^ (crc << 8);
    }
    detail::require(input.eof() && !input.bad(), "read failed: " + path.string());
    std::ostringstream out;
    out << "crc64-ecma182:" << std::hex << std::setfill('0') << std::setw(16) << crc;
    return out.str();
}

inline json to_json_record(const Manifest& m) {
    json archives = json::array();
    for (const auto& a : m.archives) {
        json shards = json::array();
        for (const auto& s : a.shards)
            shards.push_back({{"path", s.path}, {"bytes", s.bytes}, {"checksum", s.checksum}});
        archives.push_back({{"role", a.role}, {"basename", a.basename}, {"indices", a.indices},
                            {"writers", a.writers}, {"shards", shards}});
    }
    return {{"version", m.version}, {"stage", m.stage}, {"method", m.method}, {"generation", m.generation},
            {"physical", m.physical}, {"numerical", m.numerical}, {"reference", m.reference},
            {"working_orbitals", m.working_orbitals}, {"seed", m.seed}, {"position", m.position},
            {"metrics", m.metrics}, {"targets", m.targets}, {"stop_reason", m.stop_reason},
            {"properties", m.properties}, {"archives", archives}};
}
inline Manifest from_json_record(const json& j) {
    Manifest m;
    m.version = detail::integer(j.at("version"));
    detail::require(m.version == 1, "unsupported manifest version");
    j.at("stage").get_to(m.stage); j.at("method").get_to(m.method); j.at("generation").get_to(m.generation);
    m.physical = j.at("physical"); m.numerical = j.at("numerical"); m.reference = j.at("reference");
    m.working_orbitals = j.at("working_orbitals"); m.seed = j.at("seed"); m.position = j.at("position");
    m.metrics = j.at("metrics"); m.targets = j.at("targets"); j.at("stop_reason").get_to(m.stop_reason);
    m.properties = j.at("properties");
    detail::require(j.at("archives").is_array(), "archives must be an array");
    for (const auto& entry : j.at("archives")) {
        Archive a;
        entry.at("role").get_to(a.role); entry.at("basename").get_to(a.basename);
        detail::require(entry.at("indices").is_array() && entry.at("shards").is_array(),
                        "indices and shards must be arrays");
        for (const auto& index : entry.at("indices")) a.indices.push_back(detail::integer(index));
        a.writers = detail::integer(entry.at("writers"));
        for (const auto& s : entry.at("shards"))
            a.shards.push_back({s.at("path").get<std::string>(), detail::unsigned_integer(s.at("bytes")),
                                s.at("checksum").get<std::string>()});
        m.archives.push_back(std::move(a));
    }
    detail::structure(m);
    return m;
}

/// Read metadata only. Use validate_files or collective read_published before loading.
inline Manifest read_manifest(const fs::path& generation_dir) {
    auto m = from_json_record(detail::read_json(detail::member_path(generation_dir, "manifest.json")));
    detail::require(generation_dir.filename() == m.generation, "manifest generation/directory mismatch");
    return m;
}

/// Validate every listed shard. reader_count is the consuming World size; omit it
/// for a storage-only inspection. No verdict about convergence is produced.
inline Validation validate_files(const fs::path& generation_dir, const Manifest& m,
                                 int reader_count = std::numeric_limits<int>::max()) {
    Validation result;
    try {
        detail::structure(m);
        detail::require(generation_dir.filename() == m.generation, "manifest generation/directory mismatch");
        detail::require(fs::is_directory(fs::symlink_status(generation_dir)), "generation is not a real directory");
        detail::require(reader_count > 0, "reader count must be positive");
    } catch (const std::exception& error) {
        result.diagnostics.emplace_back(error.what());
        return result;
    }
    for (const auto& a : m.archives) {
        if (a.writers > reader_count)
            result.diagnostics.push_back(a.basename + ": reader count is smaller than archive writer count");
        for (std::size_t i = 0; i < a.shards.size(); ++i) {
            const auto& s = a.shards[i];
            try {
                const auto path = detail::member_path(generation_dir, s.path);
                detail::require(fs::is_regular_file(fs::symlink_status(path)), "missing or non-regular shard");
                detail::require(fs::file_size(path) == s.bytes, "byte count mismatch");
                detail::require(checksum_file(path) == s.checksum, "checksum mismatch");
                const int writers = detail::header(path, i == 0);
                if (i == 0) detail::require(writers == a.writers, "manifest/header writer-count mismatch");
            } catch (const std::exception& error) {
                result.diagnostics.push_back(s.path + ": " + error.what());
            }
        }
    }
    return result;
}

/// Noncollective inventory builder for already CLOSED archives. Pass the actual
/// I/O writer count (ParallelOutputArchive may clamp the requested count).
inline Archive inspect_archive(const fs::path& dir, const std::string& role, const std::string& basename,
                               const std::vector<int>& indices, int writers) {
    detail::require(writers > 0, "writer count must be positive");
    Archive a{role, basename, indices, writers, {}};
    for (int i = 0; i < writers; ++i) {
        const auto name = detail::shard_name(basename, i);
        const auto path = detail::member_path(dir, name);
        detail::require(fs::is_regular_file(fs::symlink_status(path)), "missing shard: " + name);
        const int stored = detail::header(path, i == 0);
        if (i == 0) detail::require(stored == writers, "manifest/header writer-count mismatch: " + name);
        a.shards.push_back({name, fs::file_size(path), checksum_file(path)});
    }
    return a;
}

/// Noncollective checked close. Return an error for publish(local_error), without
/// throwing on one writer while the others proceed to a collective. Call on all
/// ranks, including nonwriters. Inventory construction errors must likewise be
/// caught by the caller and passed to publish; no rank may skip that collective.
[[nodiscard]] inline std::string close_archive(archive::ParallelOutputArchive<>& output) {
    try {
        output.close();
        if (output.is_io_node() && !output.local_archive().good())
            return "checkpoint: archive write/flush/close failed";
    } catch (const std::exception& error) {
        return error.what();
    } catch (...) {
        return "checkpoint: unknown archive close error";
    }
    return {};
}

/// Collective publication; only rank zero's stage_dir and manifest are used.
/// Every rank supplies its own checked-close error, if any. Adapters remain
/// responsible for providing the complete intended payload/member inventory.
inline void publish(World& world, const fs::path& stage_dir, const Manifest& m,
                    const std::string& local_error = {}) {
    detail::collective_error(world, local_error);
    world.gop.fence();
    detail::on_root(world, [&] {
        detail::structure(m);
        detail::require(fs::is_directory(fs::symlink_status(stage_dir)), "stage is not a real directory");
        const auto lock_path = stage_dir / ".publish-lock";
        detail::require(fs::create_directory(lock_path), "another publisher or stale .publish-lock");
        struct Lock {
            fs::path path;
            ~Lock() { std::error_code error; fs::remove(path, error); }
        } lock{lock_path};
        const auto dir = detail::member_path(stage_dir, m.generation);
        const auto manifest_path = dir / "manifest.json";
        detail::require(!fs::exists(fs::symlink_status(manifest_path)), "generation already has a manifest");
        const auto validation = validate_files(dir, m, world.size());
        detail::require(validation.valid(), validation.message());
        detail::write_exclusive(dir / ".manifest.tmp", to_json_record(m).dump(2) + "\n");
        fs::rename(dir / ".manifest.tmp", manifest_path);
        const auto pointer_tmp = stage_dir / (".latest-" + m.generation + ".tmp");
        detail::write_exclusive(pointer_tmp, json{{"generation", m.generation}}.dump() + "\n");
        fs::rename(pointer_tmp, stage_dir / "latest");
        return std::string{};
    });
}

/// Collective preflight: select ONLY the published generation, validate all its
/// shards and writer counts on rank zero, then broadcast either metadata or the
/// failure. Never construct a collective archive inside rank zero's preflight.
inline Manifest read_published(World& world, const fs::path& stage_dir) {
    const auto record = detail::on_root(world, [&] {
        const auto pointer = detail::read_json(detail::member_path(stage_dir, "latest"));
        const auto generation = pointer.at("generation").get<std::string>();
        detail::require(detail::generation_name(generation), "unsafe published generation");
        const auto dir = detail::member_path(stage_dir, generation);
        const auto m = read_manifest(dir);
        const auto validation = validate_files(dir, m, world.size());
        detail::require(validation.valid(), validation.message());
        return to_json_record(m).dump();
    });
    return from_json_record(json::parse(record));
}
} // namespace madness::checkpoint
#endif
