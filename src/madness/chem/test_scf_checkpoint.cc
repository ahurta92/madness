#include <madness/chem/SCFCheckpoint.h>
#include <cstdlib>
#include <iostream>

using namespace madness;
namespace it = madness::iteration;
namespace cp = madness::checkpoint;
namespace {
int failures = 0;
void check(bool ok, const char* message) {
    if (!ok) { ++failures; std::cerr << "FAILED: " << message << '\n'; }
}
template<class F> void rejects(F&& action, const char* message) {
    bool rejected = false;
    try { action(); } catch (const std::exception&) { rejected = true; }
    check(rejected, message);
}
it::Request request() {
    it::Request result;
    // Floor factor zero keeps the written targets exact; this test is about
    // checkpoints and restart, not the protocol floor.
    result.targets = {{it::Metric::bsh_residual, 1e-3},
        {it::Metric::density_change, 1e-3}, {it::Metric::energy_change, 1e-5}};
    result.floor_factor = 0.0;
    result.max_updates = 25;
    return result;
}
void exercise(World& world) {
    Molecule molecule;
    molecule.add_atom(0.1, 0, 0, 2, 2);
    CalculationParameters parameters;
    parameters.set_user_defined_value("prefix", std::string("he"));
    parameters.set_user_defined_value("xc", std::string("hf"));
    parameters.set_user_defined_value("protocol", std::vector<double>{1e-4});
    parameters.set_user_defined_value("k", 6);
    parameters.set_user_defined_value("localize", std::string("canon"));
    parameters.set_user_defined_value("print_level", 0);
    parameters.set_user_defined_value("nio", world.size());
    SCF scf(world, parameters, molecule);
    const it::Protocol protocol{1e-4, 6};
    scf.set_protocol<3>(world, protocol.thresh);
    scf.make_nuclear_potential(world);
    scf.ao = scf.project_ao_basis(world, scf.aobasis);
    scf.initial_guess(world);

    cp::WorkflowRoots roots;
    cp::CheckpointStore store(world, "run", [&](const auto&, const auto&) { return roots; });
    store.initialize("a1");
    SCFCheckpoint saved(world, scf, store, {"ground_state", "scf", "a1", "he-context"});
    auto req = request();
    const auto fresh = saved.plan({protocol}, req);
    check(fresh.plan.action == it::RestartAction::fresh, "unpublished task starts fresh");
    rejects([&] { saved.plan({protocol}, req, true); }, "explicit restart requires a selected checkpoint");
    auto control = saved.control(req, protocol, true);
    int observations = 0;
    const auto publish = control.observer;
    control.observer = [&](World& w, const SCF& owner, const it::IterationMetrics& metrics,
                           const it::Decision& decision, const it::Progress& progress) {
        ++observations;
        publish(w, owner, metrics, decision, progress);
    };
    check(scf.solve(world, control).action == it::Action::complete, "occupied SCF source converges");
    const auto first = store.read_latest("ground_state", "scf");
    check(bool(first), "terminal decision publishes native functions and metrics together");
    const cp::GenerationRef source{"ground_state", "scf", first->generation};
    roots.generations.insert(source);
    const auto original = copy(world, scf.amo);
    const auto original_eps = copy(scf.aeps);
    const auto original_occ = copy(scf.aocc);
    const double energy = scf.current_energy;
    const auto metrics = cp::iteration_metrics_from_manifest(*first);
    auto relaxed = req;
    relaxed.max_updates = 0;
    for (auto& [metric, target] : relaxed.targets) target = 1.0;
    store.initialize("a2");
    SCFCheckpoint resumed(world, scf, store, {"ground_state", "scf", "a2", "he-context"});
    auto satisfied_tighter = req;
    for (auto& [metric, target] : satisfied_tighter.targets)
        target = (target + *metrics.states[0].metrics.at(metric).value) / 2;
    check(resumed.plan({protocol}, satisfied_tighter).plan.action == it::RestartAction::reuse,
          "tighter requests already met by saved metrics reuse without evaluation");
    const auto reuse = resumed.plan({protocol}, relaxed, true);
    check(reuse.plan.action == it::RestartAction::reuse, "saved evidence qualifies without another evaluation");
    const auto saved_parameters = scf.param;
    scf.param.set_user_defined_value("dconv", 1e-10);
    scf.param.set_user_defined_value("maxiter", 100);
    scf.param.set_user_defined_value("print_level", 1);
    check(resumed.plan({protocol}, relaxed).plan.action == it::RestartAction::reuse,
          "changing only tolerances, budget or reporting does not change scientific identity");
    scf.param = saved_parameters;
    auto changed_coords = molecule.get_all_coords();
    changed_coords(0, 0) = 0; changed_coords(0, 1) = 0.1;
    scf.molecule.set_all_coords(changed_coords);
    rejects([&] { resumed.plan({protocol}, relaxed); }, "live geometry changes invalidate an existing binding");
    SCFCheckpoint moved(world, scf, store, {"ground_state", "scf", "a2", "he-context"});
    rejects([&] { moved.plan({protocol}, relaxed); },
            "ordinary restart rejects a changed geometry even with equal coordinate norms");
    scf.molecule = molecule;
    if (world.rank() == 0) scf.molecule.update_rcut_with_eprec(1e-6);
    rejects([&] { resumed.load(reuse); }, "rank-local scientific changes fail collectively before loading");
    scf.molecule = molecule;
    const int before = observations;
    if (world.size() > 1) {
        if (world.rank() == 0) scf.param.set_user_defined_value("nmo_alpha", scf.param.nmo_alpha() + 1);
        rejects([&] { SCFCheckpoint mismatched(world, scf, store,
            {"ground_state", "rank_layout", "a2", "he-context"}); },
            "rank-local initial orbital layout is rejected collectively");
        scf.param = saved_parameters;
        if (world.rank() == 0) scf.param.set_user_defined_value("freeze_occupied", true);
        rejects([&] { SCFCheckpoint mismatched(world, scf, store,
            {"ground_state", "rank_layout", "a2", "he-context"}); },
            "rank-local frozen mode fails before divergent collective branches");
        scf.param = saved_parameters;
    }
    const auto saved_pcm = scf.pcm_param;
    scf.param.set_user_defined_value("pcm_data", std::string("explicit"));
    scf.pcm_param.set_derived_values(scf.param);
    rejects([&] { SCFCheckpoint unprepared(world, scf, store,
        {"ground_state", "pcm", "a2", "he-context"}); },
        "active PCM settings require a matching prepared operator");
    scf.param = saved_parameters;
    scf.pcm_param = saved_pcm;
    SCFCheckpoint other_context(world, scf, store, {"ground_state", "scf", "a2", "another-context"});
    rejects([&] { other_context.load(reuse); }, "reuse plan from another context is rejected before loading");
    resumed.load(reuse);
    check(observations == before, "loading a reusable checkpoint never runs the numerical evaluator");
    check(scf.last_calculation_outcome->complete, "reuse supplies a qualified workflow outcome");
    check(scf.param.prefix() == "he", "native loading restores the application prefix");
    check(std::abs(scf.current_energy - energy) < 1e-14, "loaded scalar energy is unchanged");
    for (std::size_t i = 0; i < original.size(); ++i)
        check((original[i] - scf.amo[i]).norm2() < 1e-12, "raw loaded functions match before reevaluation");
    check((original_eps - scf.aeps).normf() < 1e-14 && (original_occ - scf.aocc).normf() < 1e-14,
          "raw loaded eigenvalues and occupations match before reevaluation");
    check(store.read_latest("ground_state", "scf")->generation == first->generation,
          "reuse preserves the original source generation and attempt");

    auto tight = req;
    for (auto& [metric, target] : tight.targets) target = 0.0;
    tight.max_updates = 1;
    const auto more = resumed.plan({protocol}, tight);
    check(more.plan.action == it::RestartAction::resume && more.plan.main_updates == 0,
          "unmet tolerances receive a fresh update allowance");
    resumed.load(more);
    auto next = resumed.control(tight, protocol, true, more.plan.updates_at_protocol);
    const auto publish_next = next.observer;
    int fresh_checks = 0;
    next.observer = [&](World& w, const SCF& owner, const it::IterationMetrics& frame,
                        const it::Decision& decision, const it::Progress& progress) {
        if (fresh_checks++ == 0) {
            check(frame.info.iteration == 0, "restart iteration begins at zero");
            check(frame.states[0].metrics.at(it::Metric::density_change).status == it::Availability::unavailable,
                  "restart cannot invent zero density change from scalar history");
        }
        publish_next(w, owner, frame, decision, progress);
    };
    check(scf.solve(world, next).action == it::Action::stop && fresh_checks == 2,
          "one fresh update followed by a definitive endpoint check");
    check(store.read(source).generation == source.generation, "active exact source survives later publication");
    const auto requalified = resumed.plan({protocol}, relaxed);
    check(requalified.plan.action == it::RestartAction::reuse &&
          requalified.checkpoint->record.at("task").at("outcome").at("complete") == false,
          "new request qualifies an immutable checkpoint whose original attempt stopped incomplete");
    const cp::GenerationRef requalified_source{"ground_state", "scf", requalified.checkpoint->generation};
    roots.generations.insert(requalified_source);

    scf.molecule.update_rcut_with_eprec(1e-6);
    SCFCheckpoint different(world, scf, store, {"ground_state", "scf", "a3", "he-context"});
    const auto guess = different.plan({protocol}, relaxed);
    check(guess.plan.action == it::RestartAction::resume && !guess.plan.completion,
          "changed operator settings keep functions but discard completion");
    scf.molecule = molecule;
    const auto appended = resumed.plan({protocol, {1e-6, 8}}, relaxed);
    check(appended.plan.action != it::RestartAction::reuse, "a new final representation cannot reuse coarse evidence");
    const auto converted = resumed.plan({{1e-6, 8}}, relaxed);
    check(converted.plan.action == it::RestartAction::guess && converted.plan.start_protocol == 0,
          "removed saved protocol is a guess at the first requested protocol");
    rejects([&] { resumed.load(converted); }, "loading requires the caller to prepare the chosen representation");
    scf.param.set_user_defined_value("l", 40.0);
    SCFCheckpoint another_cell(world, scf, store, {"ground_state", "scf", "a3", "he-context"});
    rejects([&] { another_cell.plan({protocol}, relaxed); }, "unsupported cell conversion fails before native loading");
    scf.param.set_user_defined_value("l", 20.0);

    // Frozen virtual work has a separate task and pins its exact occupied source.
    auto vp = parameters;
    vp.set_user_defined_value("freeze_occupied", true);
    vp.set_user_defined_value("nvalpha", 1);
    SCF virtuals(world, vp, molecule);
    virtuals.set_protocol<3>(world, protocol.thresh);
    virtuals.make_nuclear_potential(world);
    virtuals.param.set_user_defined_value("prefix", (store.task_directory("ground_state", "scf") / source.generation / "scf").string());
    virtuals.load_mos(world);
    virtuals.param.set_user_defined_value("prefix", std::string("virtuals"));
    const cp::json inputs = cp::json::array({{{"role", "occupied"}, {"calculation_id", source.calculation_id},
        {"task_id", source.task_id}, {"generation", source.generation}}});
    auto requalified_inputs = inputs;
    requalified_inputs[0]["generation"] = requalified_source.generation;
    bool accepted_requalified = true;
    try {
        SCFCheckpoint from_requalified(world, virtuals, store,
            {"ground_state", "requalified_virtuals", "a3", "virtual-context"}, requalified_inputs);
    } catch (const std::exception&) { accepted_requalified = false; }
    check(accepted_requalified, "current occupied qualification overrides its historical incomplete outcome");
    rejects([&] { SCFCheckpoint missing(world, virtuals, store,
        {"ground_state", "virtuals", "a3", "virtual-context"}); },
        "frozen task cannot omit its occupied source");
    SCFCheckpoint frozen(world, virtuals, store, {"ground_state", "virtuals", "a3", "virtual-context"}, inputs);
    it::Request vr;
    vr.targets[it::Metric::bsh_residual] = 10.0;
    vr.floor_factor = 0.0;
    vr.max_updates = 0;
    check(virtuals.solve_virtuals(world, frozen.control(vr, protocol, true)).action == it::Action::complete,
          "frozen virtual task publishes through the same binding");
    auto vf = store.read_latest("ground_state", "virtuals");
    check(vf->record.at("task").at("inputs") == inputs, "virtual record retains exact occupied source");
    const auto vreuse = frozen.plan({protocol}, vr);
    frozen.load(vreuse);
    check(virtuals.current_energy == energy, "virtual reuse retains sourced occupied energy");
    virtuals.param.set_user_defined_value("aocc", std::vector<double>{0.0});
    SCFCheckpoint changed_frozen(world, virtuals, store,
        {"ground_state", "virtuals", "a4", "virtual-context"}, inputs);
    rejects([&] { changed_frozen.plan({protocol}, vr); },
            "changed frozen occupations require rebuilding from a qualified occupied source");
    virtuals.param.set_user_defined_value("aocc", std::vector<double>{});
    auto other_inputs = inputs;
    other_inputs[0]["role"] = "replacement_occupied";
    SCFCheckpoint rebound(world, virtuals, store, {"ground_state", "virtuals", "a4", "virtual-context"}, other_inputs);
    auto strict_virtual = vr;
    strict_virtual.targets[it::Metric::bsh_residual] = 0.0;
    const auto virtual_resume = frozen.plan({protocol}, strict_virtual);
    check(virtual_resume.plan.action == it::RestartAction::resume, "virtual fixture needs more iterations");
    rejects([&] { rebound.load(virtual_resume); }, "resume plan cannot cross a frozen-source binding");
    rejects([&] { rebound.plan({protocol}, vr); },
            "changed frozen dependencies cannot load the previous occupied bundle");
    SCFCheckpoint republished(world, virtuals, store, {"ground_state", "virtuals", "a3", "virtual-context"}, inputs);
    check(virtuals.solve_virtuals(world, republished.control(vr, protocol, true)).action == it::Action::complete,
          "reconstructed binding in the same attempt can publish a new generation");
    const auto next_vf = store.read_latest("ground_state", "virtuals");
    check(next_vf->generation != vf->generation, "same-attempt publication never overwrites an immutable generation");
    vf = next_vf;

    // Publish a genuine two-function active bundle, then request a one-function
    // guess for a model with a core potential. No Be numerical solve is needed:
    // this exercises the preflight boundary before the legacy reader can drop
    // an active occupied function at its nonzero core offset.
    const auto old_parameters = virtuals.param;
    const auto old_molecule = virtuals.molecule;
    virtuals.param.set_user_defined_value("freeze_occupied", false);
    // The full solve judges every default criterion: one update, loose targets.
    auto active_request = vr;
    active_request.targets[it::Metric::density_change] = 1000.0;
    active_request.targets[it::Metric::energy_change] = 1000.0;
    active_request.max_updates = 1;
    SCFCheckpoint active(world, virtuals, store, {"ground_state", "active", "a5", "active-context"});
    check(virtuals.solve(world, active.control(active_request, protocol, true)).action == it::Action::complete,
          "active-space source fixture publishes both occupied and virtual functions");
    const auto active_parameters = virtuals.param;
    const auto active_occupations = copy(virtuals.aocc);
    virtuals.param.set_user_defined_value("aocc", std::vector<double>{0.0, 1.0});
    SCFCheckpoint holes(world, virtuals, store, {"ground_state", "occupations", "a7", "occupation-context"});
    virtuals.solve(world, holes.control(active_request, protocol, true));
    check(bool(store.read_latest("ground_state", "occupations")), "occupation-hole fixture publishes real native functions");
    for (const auto& occupations : {std::vector<double>{}, std::vector<double>{1.0}}) {
        virtuals.param.set_user_defined_value("aocc", occupations);
        SCFCheckpoint changed(world, virtuals, store, {"ground_state", "occupations", "a8", "occupation-context"});
        rejects([&] { changed.plan({protocol}, active_request); },
                "ordinary restart rejects changed occupation requests");
    }
    virtuals.param = active_parameters;
    virtuals.aocc = active_occupations;
    Molecule core_molecule;
    core_molecule.add_atom(0, 0, 0, 4, 4);
    const char* previous_data = std::getenv("MRA_DATA_DIR");
    const std::optional<std::string> previous_data_directory = previous_data
        ? std::optional<std::string>(previous_data) : std::nullopt;
    setenv("MRA_DATA_DIR", MRA_CHEMDATA_DIR, 1);
    core_molecule.read_core_file("mcp");
    if (previous_data_directory) setenv("MRA_DATA_DIR", previous_data_directory->c_str(), 1);
    else unsetenv("MRA_DATA_DIR");
    check(core_molecule.n_core_orb_all() > 0, "core trimming fixture has a nonzero native reader offset");
    virtuals.molecule = core_molecule;
    virtuals.param.set_user_defined_value("nmo_alpha", 1);
    virtuals.param.set_user_defined_value("nvalpha", 0);
    SCFCheckpoint trimmed(world, virtuals, store, {"ground_state", "active", "a6", "core-context"});
    rejects([&] { trimmed.plan({protocol}, active_request); }, "core-potential virtual trimming is rejected before native reading");
    virtuals.param = old_parameters;
    virtuals.molecule = old_molecule;

    // Corruption is a strict failure, not a fallback to a fresh guess.
    if (world.rank() == 0) {
        const auto path = store.task_directory("ground_state", "virtuals") / vf->generation /
            vf->archives.front().shards.front().path;
        std::ofstream file(path, std::ios::binary | std::ios::app);
        file.put('x');
    }
    world.gop.fence();
    rejects([&] { frozen.plan({protocol}, vr); }, "selected corrupt native checkpoint is rejected collectively");
}
}
int main(int argc, char** argv) {
    initialize(argc, argv);
    {
        World world(SafeMPI::COMM_WORLD);
        startup(world, argc, argv, true);
        const auto original = std::filesystem::current_path();
        std::string directory;
        if (world.rank() == 0) {
            char pattern[] = "/tmp/madness-scf-checkpoint-XXXXXX";
            if (const auto* made = mkdtemp(pattern)) directory = made;
        }
        world.gop.broadcast_serializable(directory, 0);
        MADNESS_CHECK_THROW(!directory.empty(), "cannot create SCF checkpoint test directory");
        std::filesystem::current_path(directory);
        try { exercise(world); }
        catch (const std::exception& error) { ++failures; std::cerr << error.what() << '\n'; }
        world.gop.max(failures);
        std::filesystem::current_path(original);
        if (world.rank() == 0 && failures == 0) std::filesystem::remove_all(directory);
    }
    finalize();
    return failures ? 1 : 0;
}
