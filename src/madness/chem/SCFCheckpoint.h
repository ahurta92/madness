#pragma once

#include <madness/chem/NodeArchive.h>
#include <madness/chem/SCF.h>

namespace madness {

/// The node folders of a run directory.
std::filesystem::path ground_state_node(const std::filesystem::path& run);
std::filesystem::path virtuals_node(const std::filesystem::path& run);

/// Collective read of a node's description; rank 0 alone touches the folder.
std::optional<iteration::NodeDescription> read_node(World& world, const std::filesystem::path& folder);

/// What later calculations record to know whether a saved node has changed:
/// its settings, its level and its save time.
nlohmann::json node_identity(const iteration::NodeDescription& description);

/// Append one launch row (time, code version, host, processes, threads) to a
/// node's history.jsonl. Collective; rank 0 writes.
void write_launch_header(World& world, const std::filesystem::path& folder);

/// Application binding of SCF to one node folder (NodeArchive.h): the ground
/// state, or the virtuals of a frozen-occupied solve, whose lineage is the
/// ground-state identity. All operations are collective. The binding captures
/// SCF's scientific settings; the workflow owns protocol scheduling.
/// Keep this object alive while using its controls; it owns no orbital history.
class SCFCheckpoint {
public:
    struct Identity {
        std::string task_id, context_id;
    };
    struct Restart {
        const iteration::RestartPlan plan;
        const std::vector<iteration::Protocol> protocols;
        const std::optional<iteration::NodeDescription> saved;

    private:
        friend class SCFCheckpoint;
        const SCFCheckpoint* binding_;
        const std::optional<iteration::Request> display_request_;
        Restart(iteration::RestartPlan plan, std::vector<iteration::Protocol> protocols,
                std::optional<iteration::NodeDescription> saved, std::optional<iteration::Request> display_request,
                const SCFCheckpoint* binding)
            : plan(std::move(plan)), protocols(std::move(protocols)),
              saved(std::move(saved)), binding_(binding), display_request_(std::move(display_request)) {}
    };

    /// Capture effective geometry, model, operator and occupation settings,
    /// excluding tolerances, budgets and output controls. A frozen solve names
    /// its qualified ground state in `lineage`.
    SCFCheckpoint(World& world, SCF& scf, std::filesystem::path folder, Identity identity,
                  nlohmann::json lineage = nlohmann::json());
    SCFCheckpoint(const SCFCheckpoint&) = delete;
    SCFCheckpoint& operator=(const SCFCheckpoint&) = delete;

    const std::filesystem::path& folder() const { return folder_; }
    std::optional<iteration::NodeDescription> saved() const { return read_node(world_, folder_); }

    Restart plan(const std::vector<iteration::Protocol>& protocols,
                 const iteration::Request& request, bool require_existing = false);

    /// The caller first sets the chosen protocol and rebuilds its operators.
    /// Reuse loads the exact saved functions without evaluating or saving.
    /// Resume/guess discards convergence claims and starts fresh comparisons.
    /// Step-down reporting uses intermediate limits and leaves scheduling
    /// messages to the workflow. This flag changes no loading/restart decisions.
    void load(const Restart& restart, bool report_preparation = false);

    /// Load another node's orbitals (the ground state, for a frozen solve) as
    /// a guess. No saved metric or completion claim transfers.
    void load_guess(const std::filesystem::path& source);

    SCF::IterationControl control(const iteration::Request& request,
                                  iteration::Protocol protocol, bool final_protocol,
                                  std::size_t previous_updates = 0);

private:
    World& world_;
    SCF& scf_;
    std::filesystem::path folder_;
    Identity identity_;
    nlohmann::json settings_, lineage_, layout_;

    [[nodiscard]] std::string kind() const;
    [[nodiscard]] iteration::ExpectedTask expected() const;
    [[nodiscard]] nlohmann::json layout() const;
    [[nodiscard]] nlohmann::json capture_settings() const;
    void check_current() const;
    [[nodiscard]] bool loadable(const nlohmann::json& stored) const;
    void save(const iteration::IterationMetrics& metrics,
              const iteration::Decision& decision, const iteration::Progress& progress);
};

} // namespace madness
