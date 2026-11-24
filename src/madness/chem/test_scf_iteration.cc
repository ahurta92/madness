#include <madness/chem/SCF.h>
#include <madness/chem/NodeArchive.h>
#include <madness/chem/RootCoordination.h>
#include <filesystem>
#include <cstdlib>
#include <iostream>

using namespace madness;
namespace it = madness::iteration;
namespace cp = madness::checkpoint;
namespace {
int failures = 0;
void check(bool ok, const char* why) {
    if (!ok) { ++failures; std::cerr << "FAILED: " << why << '\n'; }
}
void near(double a, double b, const char* why) { check(std::abs(a - b) < 1e-9, why); }
template<class F> void rejects(F&& action, const char* why) {
    bool rejected = false;
    try { action(); } catch (const std::exception&) { rejected = true; }
    check(rejected, why);
}
// Zero targets for all four metrics; the protocol floor sets the effective limit.
it::Request request(std::size_t cap) {
    it::Request value;
    for (const auto metric : {it::Metric::bsh_residual, it::Metric::density_change,
                              it::Metric::energy_change, it::Metric::kinetic_change})
        value.targets[metric] = 0.0;
    value.max_updates = cap;
    return value;
}
void exercise(World& world, int atom, const std::string& xc, bool preflight_only = false) {
    const std::string name = atom == 3 ? "li" : xc == "hf" ? "he" : "he_lda";
    Molecule molecule;
    molecule.add_atom(0.0, 0.0, 0.0, double(atom), atom);
    CalculationParameters p;
    p.set_user_defined_value("prefix", name);
    p.set_user_defined_value("xc", xc);
    p.set_user_defined_value("l", 20.0);
    p.set_user_defined_value("protocol", std::vector<double>{1e-4, 1e-6});
    p.set_user_defined_value("k", 6); // also test a threshold change at fixed k
    p.set_user_defined_value("localize", atom == 3 ? std::string("new") : std::string("canon"));
    p.set_user_defined_value("print_level", 0);
    p.set_user_defined_value("maxiter", 0); // shared settings, not this legacy allowance, control the overload
    if (atom == 3) {
        p.set_user_defined_value("nopen", 1);
        p.set_user_defined_value("nvalpha", 1);
        p.set_user_defined_value("nvbeta", 1);
    }
    SCF scf(world, p, molecule);
    scf.set_protocol<3>(world, 1e-4);
    scf.ao = scf.project_ao_basis(world, scf.aobasis);
    scf.make_nuclear_potential(world);
    scf.initial_guess(world);
    auto model = scf.iteration_model();
    check(!model.at(it::Metric::energy_change).floor_factor && !model.at(it::Metric::kinetic_change).floor_factor &&
          model.at(it::Metric::density_change).floor_factor == 5.0, // one atom: max(5, 1)
          "density keeps max(5, natom); energy and kinetic use the shared floor factor");
    check(!model.at(it::Metric::kinetic_change).default_criterion &&
          model.at(it::Metric::kinetic_change).default_target == 1e-5, "kinetic is a criterion only when written");
    const auto final_defaults = it::resolve({}, model, {1e-6, 6}, true);
    check(final_defaults.criteria.size() == 3, "three default criteria");
    near(final_defaults.criteria[0].effective, 1e-4, "residual default above the fine floor");
    near(final_defaults.criteria[1].effective, 1e-4, "physical spin density default");
    near(final_defaults.criteria[2].effective, 1e-5, "physical total energy default");
    const auto coarse_final = it::resolve({}, model, {1e-4, 6}, true);
    near(coarse_final.criteria[0].effective, 1e-4, "a coarse final level judges the BSH default as written");
    near(coarse_final.criteria[2].effective, 1e-5, "a coarse final level judges the energy default as written");
    const auto coarse_intermediate = it::resolve({}, model, {1e-4, 6}, false);
    near(coarse_intermediate.criteria[0].effective, 5e-4, "an intermediate level raises BSH to its floor");
    near(coarse_intermediate.criteria[2].effective, 1e-5, "energy change has no floor");
    auto r = request(1);
    auto strict = r; strict.floor_factor = 0.0; // never relaxed: zero targets stay zero
    const it::ExpectedTask task{"scf", name + "-scientific-context", {"ground"}};
    SCF::IterationControl control{it::resolve(r, model, {1e-4, 6}, false), false, task};
    near(control.settings.criteria[0].effective, 5e-4, "intermediate BSH floor is explicit");
    near(control.settings.criteria[1].effective, 5e-4, "intermediate spin density floor is explicit");
    std::vector<it::IterationMetrics> records;
    std::optional<SCF::OrbitalComparison> previous;
    // The endpoint is saved as a node, as SCFCheckpoint saves it.
    const auto folder = std::filesystem::absolute(name + "-run") / "node";
    cp::detail::on_root(world, [&] { std::filesystem::remove_all(folder.parent_path()); return std::string{}; });
    int writes = 0;
    control.observer = [&](World& w, const SCF& owner, const it::IterationMetrics& metrics,
                           const it::Decision& decision, const it::Progress& progress) {
        auto physical = scf.evaluate_orbitals(w);
        auto raw = scf.evaluate_residual(w, physical);
        const auto& values = metrics.states.at(0).metrics;
        check(metrics.states.size() == 1 && metrics.states[0].id == "ground", "SCF is one coupled subject");
        near(*values.at(it::Metric::bsh_residual).value, raw.maximum, "typed residual describes prepared orbitals");
        near(*values.at(it::Metric::energy_change).evidence.current, physical.energy, "typed energy describes prepared orbitals");
        near(*values.at(it::Metric::kinetic_change).evidence.current, physical.kinetic, "kinetic uses physical occupied energy");
        const auto& components = values.at(it::Metric::bsh_residual).evidence.components;
        check(components.size() == scf.amo.size() + (scf.param.have_beta() ? scf.bmo.size() : 0),
              "residual includes every iterated spin function");
        const auto& density = values.at(it::Metric::density_change);
        if (!previous) {
            for (const auto metric : {it::Metric::energy_change, it::Metric::density_change, it::Metric::kinetic_change})
                check(values.at(metric).status == it::Availability::unavailable && !values.at(metric).value &&
                      !values.at(metric).comparison, "first evaluation has no invented comparison");
        } else {
            near(*density.value, std::max((physical.alpha_density - previous->alpha_density).norm2(),
                                         (physical.beta_density - previous->beta_density).norm2()),
                 "density is the maximum physical spin change, not charge sum");
            near(*density.scale.value, std::max({physical.alpha_density.norm2(), physical.beta_density.norm2(),
                 previous->alpha_density.norm2(), previous->beta_density.norm2()}), "density relative scale uses both samples");
            near(*values.at(it::Metric::energy_change).value, std::abs(physical.energy - previous->energy),
                 "final rotations retain the preceding iteration energy");
            near(*values.at(it::Metric::kinetic_change).value, std::abs(physical.kinetic - previous->kinetic),
                 "kinetic change retains the preceding physical kinetic energy");
        }
        previous = SCF::OrbitalComparison{metrics.info, physical.alpha_density, physical.beta_density,
                                         physical.energy, physical.kinetic};
        records.push_back(metrics);
        if (decision.checkpoint) {
            ++writes;
            const auto staging = folder.parent_path() / "staging";
            cp::detail::on_root(w, [&] { std::filesystem::create_directories(staging); return std::string{}; });
            cp::detail::collective_error(w, owner.write_mos(w, (staging / "scf.restartdata").string(), owner.orbital_metadata()));
            w.gop.fence();
            cp::detail::on_root(w, [&] {
                it::NodeDescription d{"ground_state", metrics.info.protocol, progress.updates_at_protocol,
                    {{"xc", xc}}, it::node_metrics(metrics, SCF::iteration_metric_set()), nullptr, {}, "test", false, {}};
                it::Node(folder).save(d, [&](const std::filesystem::path& archive) {
                    for (const auto& file : std::filesystem::directory_iterator(staging))
                        std::filesystem::rename(file.path(), archive / file.path().filename());
                });
                std::filesystem::remove_all(staging);
                return std::string{};
            });
        }
    };
    auto wrong = control;
    wrong.settings.protocol.k = 8;
    rejects([&] { scf.solve(world, wrong); }, "different numerical representation is rejected before solving");
    wrong = control; wrong.task.states = {"wrong"};
    rejects([&] { scf.solve(world, wrong); }, "SCF cannot certify an arbitrary bundle");
    wrong = control;
    wrong.settings.max_updates = 0;
    int unexpected_observations = 0;
    wrong.observer = [&](World&, const SCF&, const it::IterationMetrics&, const it::Decision&, const it::Progress&) {
        ++unexpected_observations;
    };
    scf.param.set_user_defined_value("protocol", std::vector<double>{1e-5, 1e-6});
    rejects([&] { scf.solve(world, wrong); }, "a threshold absent from the requested sequence is rejected");
    check(unexpected_observations == 0, "unrequested protocol is rejected before numerical evaluation");
    // If the membership regression has already failed, do not follow it into
    // the known pre-fix empty-vector .back() failure while collecting RED evidence.
    if (failures == 0) {
        scf.param.set_user_defined_value("protocol", std::vector<double>{});
        rejects([&] { scf.solve(world, wrong); }, "empty requested protocol is rejected safely");
        check(unexpected_observations == 0, "empty protocol never reaches the observer");
    }
    scf.param.set_user_defined_value("protocol", std::vector<double>{1e-4, 1e-6});
    if (preflight_only) return;
    auto end = scf.solve(world, control);
    check(end.action == it::Action::finish_protocol && !scf.last_calculation_outcome->complete,
          "intermediate work does not establish final convergence");
    check(records.size() == 2 && records.back().info.iteration == 1 && writes == 1,
          "one allowed update, free terminal evaluation, endpoint-only publication");
    const auto loaded = cp::json::parse(cp::detail::on_root(world, [&] {
        const auto d = it::Node(folder).inspect();
        return d ? it::description_record(*d).dump() : std::string("null");
    }));
    check(!loaded.is_null(), "native payload and typed metrics published together");
    const auto restored = it::stored_metrics(it::read_description(loaded), {task.task_id, {}, 0, task.context_id});
    check(restored.info == records.back().info, "stored evidence preserves the exact iteration");
    auto loose = request(0);
    for (auto& [metric, target] : loose.targets) target = 10.0;
    auto final = it::resolve(loose, model, {1e-4, 6}, true);
    check(bool(it::saved_completion(final, task, restored, true)), "saved evidence can satisfy a new request without evaluation");
    near(it::resolve(r, model, {1e-4, 6}, true).criteria[0].effective, 0.0, "no floor at the final level");
    final = it::resolve(strict, model, {1e-4, 6}, true);
    near(final.criteria[0].effective, 0.0, "a zero floor factor leaves the zero target");
    check(!it::saved_completion(final, task, restored, true), "stricter unmet request must iterate");
    // Load the exact native payload with the existing SCF reader. It stays a
    // numerical guess for legacy readers: typed evidence is the new authority.
    scf.param.set_user_defined_value("prefix", (it::Node(folder).archive() / "scf").string());
    scf.load_mos(world);
    auto saved = scf.evaluate_orbitals(world);
    near(saved.energy, *restored.states[0].metrics.at(it::Metric::energy_change).evidence.current,
         "native loaded energy matches typed checkpoint evidence");
    scf.param.set_user_defined_value("prefix", name);
    previous.reset(); records.clear();
    control.previous_updates = 1;
    control.settings = final;
    control.final_protocol = true;
    // Use a new selected protocol below for helium; Li/DFT keep the same one to
    // verify fresh comparison history and allowance on restart.
    if (name == "he") {
        scf.set_protocol<3>(world, 1e-6);
        scf.project(world);
        scf.make_nuclear_potential(world);
        control.settings = it::resolve(strict, scf.iteration_model(), {1e-6, 6}, true);
        control.previous_updates = 0;
    } else scf.param.set_user_defined_value("protocol", std::vector<double>{1e-4});
    end = scf.solve(world, control);
    check(end.action == it::Action::stop && end.reason == "iteration_limit" && records.size() == 2,
          "strict final solve uses a fresh allowance and stops unconverged");
    check(records.front().info.iteration == 0 && records.back().info.iteration == 1,
          "restart counts main updates from zero");
    // Residual-only completion at iteration zero must not invoke KAIN. Kinetic
    // and density remain reported as unavailable diagnostics.
    control.settings = it::resolve(loose, scf.iteration_model(), control.settings.protocol, true);
    control.settings.criteria.resize(1);
    previous.reset(); records.clear();
    end = scf.solve(world, control);
    check(end.action == it::Action::complete && records.size() == 1 && end.info.iteration == 0,
          "satisfied residual-only request completes without another update");
    check(scf.last_calculation_outcome->complete, "shared final decision supplies workflow outcome");
}
}
int main(int argc, char** argv) {
    initialize(argc, argv);
    {
        World world(SafeMPI::COMM_WORLD);
        startup(world, argc, argv, true);
        const auto original = std::filesystem::current_path();
        std::string directory;
        if (world.rank() == 0) { char pattern[] = "/tmp/madness-scf-iteration-XXXXXX";
            if (const char* made = mkdtemp(pattern)) directory = made; }
        world.gop.broadcast_serializable(directory, 0);
        MADNESS_CHECK_THROW(!directory.empty(), "cannot create SCF iteration test directory");
        std::filesystem::current_path(directory);
        try {
            const bool preflight_only = argc > 1 && std::string(argv[1]) == "--preflight";
            exercise(world, 2, "hf", preflight_only);
            if (!preflight_only) { exercise(world, 3, "hf"); exercise(world, 2, "lda"); }
        }
        catch (const std::exception& error) { ++failures; std::cerr << error.what() << '\n'; }
        world.gop.fence();
        std::filesystem::current_path(original);
        if (world.rank() == 0 && failures == 0) std::filesystem::remove_all(directory);
        world.gop.max(failures);
    }
    finalize();
    return failures ? 1 : 0;
}
