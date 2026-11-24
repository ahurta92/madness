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

    const auto ground = ground_state_node("run");
    const auto node = [&](const std::string& name) { return std::filesystem::path("run") / "states" / name; };
    SCFCheckpoint saved(world, scf, ground, {"scf", "he-context"});
    auto req = request();
    const auto fresh = saved.plan({protocol}, req);
    check(fresh.plan.action == it::RestartAction::fresh, "unpublished task starts fresh");
    rejects([&] { saved.plan({protocol}, req, true); }, "explicit restart requires a saved state");
    auto control = saved.control(req, protocol, true);
    int observations = 0;
    const auto publish = control.observer;
    control.observer = [&](World& w, const SCF& owner, const it::IterationMetrics& metrics,
                           const it::Decision& decision, const it::Progress& progress) {
        ++observations;
        publish(w, owner, metrics, decision, progress);
    };
    check(scf.solve(world, control).action == it::Action::complete, "occupied SCF source converges");
    const auto first = read_node(world, ground);
    check(bool(first), "terminal decision saves native functions and metrics together");
    const auto original = copy(world, scf.amo);
    const auto original_eps = copy(scf.aeps);
    const auto original_occ = copy(scf.aocc);
    const double energy = scf.current_energy;
    const auto metrics = it::stored_metrics(*first, {"scf", protocol, 0, "he-context"});
    auto relaxed = req;
    relaxed.max_updates = 0;
    for (auto& [metric, target] : relaxed.targets) target = 1.0;
    SCFCheckpoint resumed(world, scf, ground, {"scf", "he-context"});
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
    SCFCheckpoint moved(world, scf, ground, {"scf", "he-context"});
    rejects([&] { moved.plan({protocol}, relaxed); },
            "ordinary restart rejects a changed geometry even with equal coordinate norms");
    scf.molecule = molecule;
    if (world.rank() == 0) scf.molecule.update_rcut_with_eprec(1e-6);
    rejects([&] { resumed.load(reuse); }, "rank-local scientific changes fail collectively before loading");
    scf.molecule = molecule;
    const int before = observations;
    if (world.size() > 1) {
        if (world.rank() == 0) scf.param.set_user_defined_value("nmo_alpha", scf.param.nmo_alpha() + 1);
        rejects([&] { SCFCheckpoint mismatched(world, scf, node("rank_layout"),
            {"rank_layout", "he-context"}); },
            "rank-local initial orbital layout is rejected collectively");
        scf.param = saved_parameters;
        if (world.rank() == 0) scf.param.set_user_defined_value("freeze_occupied", true);
        rejects([&] { SCFCheckpoint mismatched(world, scf, node("rank_layout"),
            {"rank_layout", "he-context"}); },
            "rank-local frozen mode fails before divergent collective branches");
        scf.param = saved_parameters;
    }
    const auto saved_pcm = scf.pcm_param;
    scf.param.set_user_defined_value("pcm_data", std::string("explicit"));
    scf.pcm_param.set_derived_values(scf.param);
    rejects([&] { SCFCheckpoint unprepared(world, scf, node("pcm"),
        {"pcm", "he-context"}); },
        "active PCM settings require a matching prepared operator");
    scf.param = saved_parameters;
    scf.pcm_param = saved_pcm;
    SCFCheckpoint other_context(world, scf, ground, {"scf", "another-context"});
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
    check(read_node(world, ground) == first, "reuse saves nothing");

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
    const auto requalified = resumed.plan({protocol}, relaxed);
    check(requalified.plan.action == it::RestartAction::reuse && !requalified.saved->level_finished,
          "a new request qualifies a saved state whose solve stopped incomplete");

    scf.molecule.update_rcut_with_eprec(1e-6);
    SCFCheckpoint different(world, scf, ground, {"scf", "he-context"});
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
    SCFCheckpoint another_cell(world, scf, ground, {"scf", "he-context"});
    rejects([&] { another_cell.plan({protocol}, relaxed); }, "unsupported cell conversion fails before native loading");
    scf.param.set_user_defined_value("l", 20.0);

    // Frozen virtual work is a second node whose lineage is the ground-state identity.
    auto vp = parameters;
    vp.set_user_defined_value("freeze_occupied", true);
    vp.set_user_defined_value("nvalpha", 1);
    SCF virtuals(world, vp, molecule);
    virtuals.set_protocol<3>(world, protocol.thresh);
    virtuals.make_nuclear_potential(world);
    virtuals.param.set_user_defined_value("prefix", (it::Node(ground).archive() / "scf").string());
    virtuals.load_mos(world);
    virtuals.param.set_user_defined_value("prefix", std::string("virtuals"));
    const auto lineage = node_identity(*read_node(world, ground));
    const auto virtual_node = virtuals_node("run");
    rejects([&] { SCFCheckpoint missing(world, virtuals, virtual_node, {"virtuals", "virtual-context"}); },
        "frozen task cannot omit its occupied source");
    SCFCheckpoint frozen(world, virtuals, virtual_node, {"virtuals", "virtual-context"}, lineage);
    it::Request vr;
    vr.targets[it::Metric::bsh_residual] = 10.0;
    vr.floor_factor = 0.0;
    vr.max_updates = 0;
    check(virtuals.solve_virtuals(world, frozen.control(vr, protocol, true)).action == it::Action::complete,
          "frozen virtual task saves through the same binding");
    const auto vf = read_node(world, virtual_node);
    check(vf->kind == "virtuals" && vf->lineage == lineage, "virtuals node records its ground-state lineage");
    const auto vreuse = frozen.plan({protocol}, vr);
    frozen.load(vreuse);
    check(virtuals.current_energy == read_node(world, ground)->details.at("energy").get<double>(),
          "virtual reuse retains sourced occupied energy");
    virtuals.param.set_user_defined_value("aocc", std::vector<double>{0.0});
    SCFCheckpoint changed_frozen(world, virtuals, virtual_node, {"virtuals", "virtual-context"}, lineage);
    rejects([&] { changed_frozen.plan({protocol}, vr); },
            "changed frozen occupations require rebuilding from a qualified occupied source");
    virtuals.param.set_user_defined_value("aocc", std::vector<double>{});
    auto other_lineage = lineage;
    other_lineage["saved_at"] = "a later save";
    SCFCheckpoint rebound(world, virtuals, virtual_node, {"virtuals", "virtual-context"}, other_lineage);
    auto strict_virtual = vr;
    strict_virtual.targets[it::Metric::bsh_residual] = 0.0;
    const auto virtual_resume = frozen.plan({protocol}, strict_virtual);
    check(virtual_resume.plan.action == it::RestartAction::resume, "virtual fixture needs more iterations");
    rejects([&] { rebound.load(virtual_resume); }, "resume plan cannot cross a frozen-source binding");
    const auto changed_ground = rebound.plan({protocol}, vr);
    check(changed_ground.plan.action == it::RestartAction::resume && !changed_ground.plan.completion,
          "a ground state saved again turns the saved virtuals into a starting guess");
    SCFCheckpoint republished(world, virtuals, virtual_node, {"virtuals", "virtual-context"}, lineage);
    check(virtuals.solve_virtuals(world, republished.control(vr, protocol, true)).action == it::Action::complete,
          "reconstructed binding in the same launch can save again");

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
    SCFCheckpoint active(world, virtuals, node("active"), {"active", "active-context"});
    check(virtuals.solve(world, active.control(active_request, protocol, true)).action == it::Action::complete,
          "active-space source fixture publishes both occupied and virtual functions");
    const auto active_parameters = virtuals.param;
    const auto active_occupations = copy(virtuals.aocc);
    virtuals.param.set_user_defined_value("aocc", std::vector<double>{0.0, 1.0});
    SCFCheckpoint holes(world, virtuals, node("occupations"), {"occupations", "occupation-context"});
    virtuals.solve(world, holes.control(active_request, protocol, true));
    check(bool(read_node(world, node("occupations"))), "occupation-hole fixture saves real native functions");
    for (const auto& occupations : {std::vector<double>{}, std::vector<double>{1.0}}) {
        virtuals.param.set_user_defined_value("aocc", occupations);
        SCFCheckpoint changed(world, virtuals, node("occupations"), {"occupations", "occupation-context"});
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
    SCFCheckpoint trimmed(world, virtuals, node("active"), {"active", "core-context"});
    rejects([&] { trimmed.plan({protocol}, active_request); }, "core-potential virtual trimming is rejected before native reading");
    virtuals.param = old_parameters;
    virtuals.molecule = old_molecule;

    // A half-present node is an error, never a silent fresh start.
    const auto description = virtual_node / "description.json";
    std::string kept;
    if (world.rank() == 0) {
        std::ifstream in(description); kept.assign(std::istreambuf_iterator<char>(in), {});
        std::ofstream(description) << "{ not json";
    }
    world.gop.fence();
    rejects([&] { frozen.plan({protocol}, vr); }, "an unreadable description is rejected collectively");
    if (world.rank() == 0) std::filesystem::remove(description);
    world.gop.fence();
    rejects([&] { frozen.plan({protocol}, vr); }, "an archive without a description is rejected collectively");
    if (world.rank() == 0) std::ofstream(description) << kept;
    world.gop.fence();
    check(frozen.plan({protocol}, vr).plan.action == it::RestartAction::reuse, "the restored description reads again");
    if (world.rank() == 0) std::filesystem::rename(it::Node(virtual_node).archive(), virtual_node / "moved");
    world.gop.fence();
    rejects([&] { frozen.plan({protocol}, vr); }, "a description without its archive is rejected collectively");
    if (world.rank() == 0) std::filesystem::rename(virtual_node / "moved", it::Node(virtual_node).archive());
    world.gop.fence();
    check(!read_node(world, node("never_saved")), "a folder with nothing saved is no node");

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
