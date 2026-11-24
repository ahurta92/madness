#pragma once

// Collective helpers for work that only rank zero does: run a function there,
// make every rank throw its error, read and write small files, close a native
// archive without throwing on one writer alone.
#include <madness/chem/MetricRecords.h>
#include <madness/world/world.h>
#include <madness/world/worldgop.h>
#include <madness/world/parallel_archive.h>
#include <nlohmann/json.hpp>

#include <cerrno>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <fcntl.h>
#include <unistd.h>

namespace madness::checkpoint {
using json = nlohmann::json;

namespace detail {
inline void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error("checkpoint: " + message);
}
inline json read_json(const std::filesystem::path& path) {
    require(std::filesystem::is_regular_file(std::filesystem::symlink_status(path)),
            "not a regular file: " + path.string());
    std::ifstream input(path);
    require(input.good(), "cannot read " + path.string());
    return json::parse(input);
}
inline void write_exclusive(const std::filesystem::path& path, const std::string& content) {
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
// Replace a small JSON file: write a temporary name, then rename over it.
inline void replace_json(const std::filesystem::path& path, const json& record) {
    records::finite_numbers(record);
    const auto temp = path.parent_path() / ("." + path.filename().string() + ".tmp");
    if (std::filesystem::is_regular_file(std::filesystem::symlink_status(temp)))
        std::filesystem::remove(temp); // abandoned, never a recovery source
    write_exclusive(temp, record.dump(2) + "\n");
    require(!std::filesystem::is_symlink(std::filesystem::symlink_status(path)), "symlink record destination");
    std::filesystem::rename(temp, path);
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

/// Noncollective checked close. Returns the error instead of throwing on one
/// writer while the others proceed to a collective; the caller passes it to
/// detail::collective_error. Call on all ranks, including nonwriters.
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
} // namespace madness::checkpoint
