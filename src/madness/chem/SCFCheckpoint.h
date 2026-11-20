#pragma once

#include <madness/chem/CheckpointStore.h>
#include <madness/chem/SCF.h>

namespace madness {

/// Application binding for shared SCF records and native orbital archives.
/// All operations are collective. The binding captures SCF's scientific settings;
/// the workflow owns attempts, protocol scheduling and qualification of occupied sources.
/// Keep this object alive while using its controls; it owns no orbital history.
class SCFCheckpoint {
public:
    struct Identity {
        std::string calculation_id, task_id, attempt_id, context_id;
    };
    struct Restart {
        const iteration::RestartPlan plan;
        const std::vector<iteration::Protocol> protocols;
        const std::optional<checkpoint::Manifest> checkpoint;

    private:
        friend class SCFCheckpoint;
        const SCFCheckpoint* binding_;
        const std::optional<iteration::Request> display_request_;
        Restart(iteration::RestartPlan plan, std::vector<iteration::Protocol> protocols,
                std::optional<checkpoint::Manifest> checkpoint, std::optional<iteration::Request> display_request,
                const SCFCheckpoint* binding)
            : plan(std::move(plan)), protocols(std::move(protocols)),
              checkpoint(std::move(checkpoint)), binding_(binding), display_request_(std::move(display_request)) {}
    };

    /// Capture effective geometry, model, operator and occupation settings,
    /// excluding tolerances, budgets and output controls.
    /// A frozen task requires its already-qualified occupied source in inputs.
    SCFCheckpoint(World& world, SCF& scf, checkpoint::CheckpointStore& store,
                  Identity identity,
                  checkpoint::json inputs = checkpoint::json::array());
    SCFCheckpoint(const SCFCheckpoint&) = delete;
    SCFCheckpoint& operator=(const SCFCheckpoint&) = delete;

    Restart plan(const std::vector<iteration::Protocol>& protocols,
                 const iteration::Request& request, bool require_existing = false);

    /// The result is immutable and valid only during this binding's lifetime.
    /// The caller first sets the chosen protocol and rebuilds its operators.
    /// Reuse loads the exact saved functions without evaluating or publishing.
    /// Resume/guess discards convergence claims and starts fresh comparisons.
    /// Preparation reporting uses intermediate limits and leaves scheduling
    /// messages to the workflow. This flag changes no loading/restart decisions.
    void load(const Restart& restart, bool report_preparation = false);

    /// Load verified physical functions as a guess for this task. No saved
    /// metric or completion claim transfers across task identities.
    void load_guess(const checkpoint::Manifest& source);

    SCF::IterationControl control(const iteration::Request& request,
                                  iteration::Protocol protocol, bool final_protocol,
                                  std::size_t previous_updates = 0);

private:
    World& world_;
    SCF& scf_;
    checkpoint::CheckpointStore& store_;
    Identity identity_;
    checkpoint::json settings_, inputs_, layout_;
    std::size_t publications_ = 0;

    [[nodiscard]] iteration::ExpectedTask expected() const;
    [[nodiscard]] checkpoint::json layout() const;
    [[nodiscard]] checkpoint::json capture_settings() const;
    void check_current() const;
    [[nodiscard]] bool loadable(const checkpoint::json& stored) const;
    void publish(const iteration::IterationMetrics& metrics,
                 const iteration::Decision& decision, const iteration::Progress& progress);
};

} // namespace madness
