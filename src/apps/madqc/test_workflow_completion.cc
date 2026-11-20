#include <madness/chem/Drivers.hpp>
#include <madness/chem/CalculationOutcome.h>
#include <madness/world/MADworld.h>

#include <cstdlib>
#include <fstream>
#include <iostream>

namespace fs = std::filesystem;
using namespace madness;
namespace {
int failures = 0;
void check(bool condition, const char* message) {
    if (!condition) { ++failures; std::cerr << "FAILED: " << message << '\n'; }
}

// Only the expensive numerical application is replaced; the real driver,
// context handoff, report writer and MPI failure propagation run unchanged.
class TestApplication : public Application {
public:
    explicit TestApplication(CalculationOutcome outcome, bool fail = false, bool invalid_text = false)
        : Application(Params{}), outcome_(std::move(outcome)), fail_(fail), invalid_text_(invalid_text) {}
    void run(const fs::path&) override {
        ++runs;
        if (fail_) throw std::runtime_error("numerical failure");
    }
    void publish_to_context(StepContext& ctx) override {
        ++publications;
        ctx.blob["qualified"] = true;
    }
    [[nodiscard]] nlohmann::json results() const override {
        return {{"type", "test"}, {"properties", {{"available", 2.0}}},
                {"description", invalid_text_ ? std::string(1, static_cast<char>(0xff)) : "valid"}};
    }
    [[nodiscard]] std::optional<CalculationOutcome> outcome() const override { return outcome_; }
    void print_parameters(World&) const override {}
    int runs = 0, publications = 0;
private:
    CalculationOutcome outcome_;
    bool fail_;
    bool invalid_text_;
};

nlohmann::json read(const fs::path& path) {
    std::ifstream input(path);
    return nlohmann::json::parse(input);
}

void workflow_case(World& world, const fs::path& root, bool complete, bool failure) {
    auto first = std::make_shared<TestApplication>(
        CalculationOutcome{complete, complete ? "final_criteria_met" : "iteration_limit"}, failure);
    auto second = std::make_shared<TestApplication>(CalculationOutcome{true, "computed"});
    qcapp::Workflow workflow;
    workflow.addDriver(std::make_unique<qcapp::SinglePointDriver>(first));
    workflow.addDriver(std::make_unique<qcapp::SinglePointDriver>(second));
    bool stopped = false, failed = false;
    try { workflow.run(root.string()); }
    catch (const qcapp::WorkflowStopped&) { stopped = true; }
    catch (const std::runtime_error&) { failed = true; }
    check(stopped == (!complete && !failure), "controlled stop is distinct from failure");
    check(failed == failure, "computation failure remains an error");
    check(first->publications == (complete && !failure ? 1 : 0), "unfinished reference not published");
    check(second->runs == (complete && !failure ? 1 : 0), "dependent execution gated");
    world.gop.fence();
    const auto records = read(root.string() + ".calc_info.json").at("tasks");
    check(records.size() == (complete && !failure ? 2 : 1), "only executed tasks reported");
    if (failure) check(records[0].at("type") == "task_failed", "failure record preserved");
    else {
        check(records[0].at("outcome").at("complete") == complete, "completion recorded");
        check(records[0].at("properties").at("available") == 2.0, "available results survive stop");
    }
}

void decision_cases() {
    iteration::Decision decision;
    decision.action = iteration::Action::complete;
    decision.reason = "saved_metrics_satisfy_request";
    check(CalculationOutcome::from_iteration(decision).complete, "saved completion releases work");
    decision.action = iteration::Action::stop;
    decision.reason = "iteration_limit";
    check(!CalculationOutcome::from_iteration(decision).complete, "cap does not release work");
    for (auto action : {iteration::Action::update, iteration::Action::finish_protocol}) {
        decision.action = action;
        bool rejected = false;
        try { (void)CalculationOutcome::from_iteration(decision); }
        catch (const std::invalid_argument&) { rejected = true; }
        check(rejected, "intermediate actions cannot complete an application");
    }
    RestartMetadata saved;
    saved.converged_for_thresh = 1e-4;
    saved.converged_for_dconv = 1e-3;
    check(saved.is_converged_to(1e-4, 1e-3), "valid precision stamps qualify");
    for (double invalid : {0.0, -1.0, -std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::quiet_NaN()}) {
        saved.converged_for_thresh = invalid;
        check(!saved.is_converged_to(1e-4, 1e-3), "invalid threshold cannot qualify");
        saved.converged_for_thresh = 1e-4;
        saved.converged_for_dconv = invalid;
        check(!saved.is_converged_to(1e-4, 1e-3), "invalid density precision cannot qualify");
        saved.converged_for_dconv = 1e-3;
    }
}

void nemo_repeat_case(World& world, const fs::path& root) {
    Params params(world);
    Molecule molecule;
    molecule.add_atom(0.0, 0.0, 0.0, 2.0, 2);
    params.set(molecule);
    auto& cp = params.get<CalculationParameters>();
    cp.set_user_defined_value("xc", std::string("hf"));
    cp.set_user_defined_value("protocol", std::vector<double>{1e-4});
    cp.set_user_defined_value("k", 6);
    cp.set_user_defined_value("econv", 1e-4);
    cp.set_user_defined_value("dconv", 1e-3);
    cp.set_user_defined_value("maxiter", 1);
    cp.set_user_defined_value("save", true);
    cp.set_user_defined_value("prefix", std::string("repeat"));
    class RepeatedNemo : public SCFApplication<nemo_lib> {
    public:
        using SCFApplication<nemo_lib>::SCFApplication;
        void continue_with_more_iterations() {
            params_.get<CalculationParameters>().set_user_defined_value("maxiter", 30);
            params_.get<CalculationParameters>().set_user_defined_value("restart", std::string("iterate"));
            calc()->get_calc()->param.set_user_defined_value("maxiter", 30);
            calc()->get_calc()->param.set_user_defined_value("restart", std::string("iterate"));
        }
    } app(world, params);
    app.run(root);
    check(!app.outcome()->complete, "first capped Nemo invocation is incomplete");
    app.continue_with_more_iterations();
    app.run(root);
    check(app.outcome()->complete, "same Nemo instance retries at unchanged geometry");
}

void report_failure_case(World& world, const fs::path& root, bool serialization_error = false) {
    const fs::path report = root.string() + ".calc_info.json";
    if (world.rank() == 0) {
        std::ofstream(report) << "{\"previous\":true}";
        if (!serialization_error) fs::create_directory(root.string() + ".calc_info.json.tmp");
    }
    world.gop.fence();
    auto first = std::make_shared<TestApplication>(CalculationOutcome{false, "iteration_limit"}, false, serialization_error);
    auto second = std::make_shared<TestApplication>(CalculationOutcome{true, "computed"});
    qcapp::Workflow workflow;
    workflow.addDriver(std::make_unique<qcapp::SinglePointDriver>(first));
    workflow.addDriver(std::make_unique<qcapp::SinglePointDriver>(second));
    bool failed = false;
    try { workflow.run(root.string()); }
    catch (const qcapp::WorkflowStopped&) { check(false, "write failure hidden by controlled stop"); }
    catch (const std::runtime_error&) { failed = true; }
    check(failed, "report failure reaches every rank");
    check(first->publications == 0 && second->runs == 0, "failure releases no dependent work");
    world.gop.fence();
    check(read(report).at("previous") == true, "failed report replacement preserves old snapshot");
}
} // namespace

int main(int argc, char** argv) {
    const bool repeat_nemo = argc > 1 && std::string(argv[1]) == "--nemo-repeat";
    initialize(argc, argv);
    {
        World world(SafeMPI::COMM_WORLD);
        std::string path;
        if (world.rank() == 0) {
            char name[] = "/tmp/madqc-workflow-XXXXXX";
            const char* created = mkdtemp(name);
            if (created) path = created;
        }
        world.gop.broadcast_serializable(path, 0);
        MADNESS_CHECK_THROW(!path.empty(), "test temporary directory creation failed");
        if (repeat_nemo) {
            startup(world, argc, argv, true);
            nemo_repeat_case(world, fs::path(path) / "nemo");
        } else {
            decision_cases();
            workflow_case(world, fs::path(path) / "complete", true, false);
            workflow_case(world, fs::path(path) / "stopped", false, false);
            workflow_case(world, fs::path(path) / "failed", true, true);
            report_failure_case(world, fs::path(path) / "write-failed");
            report_failure_case(world, fs::path(path) / "serialization-failed", true);
        }
        world.gop.fence();
        if (world.rank() == 0) fs::remove_all(path);
        world.gop.sum(failures);
    }
    finalize();
    return failures == 0 ? 0 : 1;
}
