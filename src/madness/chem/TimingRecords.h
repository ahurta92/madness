#pragma once

#include <madness/chem/IterationOutput.h>
#include <nlohmann/json.hpp>
#include <stdexcept>

namespace madness::iteration {

/// Driver-rank wall intervals at existing completion boundaries. Parent times
/// are measured independently, never reconstructed by adding child intervals.
/// This collector does not synchronize or retain work from previous attempts.
template<class Clock = std::chrono::steady_clock>
class TimingCollector {
    using json = nlohmann::json;
    struct Frame { json node; ElapsedTimer<Clock> timer; };
    std::vector<Frame> stack_;
    json finished_;

    static void require(bool ok, const char* message) {
        if (!ok) throw std::invalid_argument(message);
    }
    static void validate(const json& node) {
        require(node.is_object() && node.contains("name") && node.at("name").is_string() &&
                    !node.at("name").get<std::string>().empty(), "timing needs a step name");
        if (node.contains("seconds")) {
            require(node.at("seconds").is_number() &&
                        finite_nonnegative(node.at("seconds").get<double>()), "invalid wall interval");
            require(!node.contains("availability"), "measured timing cannot be unavailable");
        } else {
            require(node.contains("availability") &&
                        node.at("availability").at("status") == "unavailable" &&
                        !node.at("availability").at("reason").get<std::string>().empty(),
                    "missing wall interval needs a reason");
        }
        if (node.contains("children")) {
            require(node.at("children").is_array(), "timing children must be ordered");
            for (const auto& child : node.at("children")) validate(child);
        }
    }
    static json measured(const Frame& frame) {
        auto result = frame.node;
        result["seconds"] = frame.timer.seconds();
        validate(result);
        return result;
    }
public:
    void begin(const std::string& name, const json& context = json::object()) {
        require(context.is_object(), "timing context must be an object");
        for (const auto* key : {"name", "seconds", "availability", "children"})
            require(!context.contains(key), "timing context overrides a measured field");
        require(!name.empty() && finished_.is_null(), "timing attempt already ended or step name empty");
        auto node = context;
        node["name"] = name;
        stack_.push_back({std::move(node), {}});
    }
    void append(json node) {
        require(!stack_.empty(), "timing child needs an active parent");
        validate(node);
        auto& parent = stack_.back().node;
        if (!parent.contains("children")) parent["children"] = json::array();
        parent["children"].push_back(std::move(node));
    }
    void unavailable(const std::string& name, const std::string& reason) {
        append({{"name", name}, {"availability", {{"status", "unavailable"}, {"reason", reason}}}});
    }
    void end() {
        require(!stack_.empty(), "no timing scope to end");
        auto node = measured(stack_.back());
        stack_.pop_back();
        if (stack_.empty()) finished_ = std::move(node);
        else append(std::move(node));
    }
    /// A snapshot may include scopes still running; reading does not stop them.
    [[nodiscard]] json record() const {
        if (stack_.empty()) {
            require(!finished_.is_null(), "timing attempt has not started");
            return finished_;
        }
        auto node = measured(stack_.back());
        for (std::size_t i = stack_.size() - 1; i > 0; --i) {
            auto parent = measured(stack_[i - 1]);
            if (!parent.contains("children")) parent["children"] = json::array();
            parent["children"].push_back(std::move(node));
            node = std::move(parent);
        }
        return node;
    }
};

/// The process layout of one launch, supplied by the caller.
struct Execution {
    std::size_t mpi_ranks = 1;
    std::size_t threads_per_rank = 1;
};

/// The record for one launch, as write_timing_entry expects it; `started_at`
/// identifies the launch. `failure` is empty for a launch that ended normally.
/// Reads nothing from the runtime.
inline nlohmann::json timing_entry(const std::string& started_at, const Execution& execution,
                                   const nlohmann::json& timing_tree, const std::string& failure = {}) {
    nlohmann::json record{{"started_at", started_at},
        {"execution", {{"mpi_ranks", execution.mpi_ranks}, {"threads_per_rank", execution.threads_per_rank}}},
        {"timing", timing_tree}};
    if (!failure.empty()) record["failure"] = failure;
    return record;
}
} // namespace madness::iteration
