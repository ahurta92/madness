#include <madness/chem/SCF.h>
#include <madness/chem/NodeArchive.h>
#include <madness/chem/RootCoordination.h>
#include <cstdlib>
#include <filesystem>
#include <iostream>

using namespace madness;
namespace it = madness::iteration;
namespace cp = madness::checkpoint;
namespace {
int failures = 0;
void check(bool ok, const char* why) {
    if (!ok) { ++failures; std::cerr << "FAILED: " << why << '\n'; }
}
void near(double a, double b, const char* why) { check(std::abs(a - b) < 1e-8, why); }
template<class F> void rejects(F&& f, const char* why) {
    bool rejected = false;
    try { f(); } catch (const std::exception&) { rejected = true; }
    check(rejected, why);
}

// Observe a live numerical object's lifetime at an existing printed boundary;
// no solver hook, allocator limit or assumed memory size is needed.
class LineProbe : public std::streambuf {
    std::streambuf* destination_;
    std::function<void(const std::string&)> observe_;
    std::string line_;
    void consume(char c) {
        if (c == '\n') { observe_(line_); line_.clear(); }
        else line_ += c;
    }
    std::streamsize xsputn(const char* text, std::streamsize n) override {
        for (std::streamsize i = 0; i < n; ++i) consume(text[i]);
        return destination_->sputn(text, n);
    }
    int overflow(int c) override {
        if (traits_type::eq_int_type(c, traits_type::eof())) return traits_type::not_eof(c);
        consume(traits_type::to_char_type(c));
        return destination_->sputc(traits_type::to_char_type(c));
    }
    int sync() override { return destination_->pubsync(); }
public:
    explicit LineProbe(std::function<void(const std::string&)> observer)
        : destination_(std::cout.rdbuf()), observe_(std::move(observer)) { std::cout.rdbuf(this); }
    ~LineProbe() override { std::cout.rdbuf(destination_); }
};

void review_regression(World& world, const std::string& mode) {
    const int atom = mode == "empty" ? 2 : 3;
    Molecule molecule;
    molecule.add_atom(0, 0, 0, double(atom), atom);
    CalculationParameters p;
    p.set_user_defined_value("prefix", std::string("review_source"));
    p.set_user_defined_value("xc", std::string("hf"));
    p.set_user_defined_value("protocol", std::vector<double>{1e-4});
    p.set_user_defined_value("k", 6);
    p.set_user_defined_value("localize", std::string("canon"));
    p.set_user_defined_value("print_level", 0);
    if (atom == 3) p.set_user_defined_value("nopen", 1);
    SCF source(world, p, molecule);
    source.set_protocol<3>(world, 1e-4);
    source.ao = source.project_ao_basis(world, source.aobasis);
    source.make_nuclear_potential(world);
    source.initial_guess(world);
    // This deliberately off-unit fixture tests numerical preservation, not a
    // claim that an unqualified native source is reusable as a converged state.
    if (mode == "frozen") {
        source.amo[0].scale(1.0 + 1e-6);
        source.bmo[0].scale(1.0 + 2e-6);
    }
    source.save_mos(world);
    p.set_user_defined_value("freeze_occupied", true);
    p.set_user_defined_value("nvalpha", mode == "empty" ? 0 : 1);
    p.set_user_defined_value("nvbeta", 1);
    SCF scf(world, p, molecule);
    scf.set_protocol<3>(world, 1e-4);
    scf.make_nuclear_potential(world);
    scf.load_mos(world);
    it::Request request;
    request.targets[it::Metric::bsh_residual] = 0.0;
    request.floor_factor = 0.0; // a strict, never-met request
    request.max_updates = 1;
    SCF::IterationControl control{it::resolve(request, scf.iteration_model(), {1e-4, 6}, true), true,
        {"virtuals", "review-fixed-context", {"virtuals"}}};
    if (mode == "empty" || mode == "occupations") {
        scf.last_calculation_outcome = CalculationOutcome{true, "preserve_previous_outcome"};
        scf.converged_for_dconv = 1e-5;
        scf.converged_for_thresh = 1e-6;
        if (mode == "occupations") scf.aocc = tensorT(scf.param.nalpha());
        rejects([&] { scf.solve_virtuals(world, control); }, "unusable virtual bundle is rejected");
        check(scf.last_calculation_outcome->complete &&
              scf.last_calculation_outcome->reason == "preserve_previous_outcome",
              "preflight rejection preserves the previous outcome");
        near(scf.converged_for_dconv, 1e-5, "preflight rejection preserves convergence stamp");
        near(scf.converged_for_thresh, 1e-6, "preflight rejection preserves source threshold");
    } else if (mode == "frozen") {
        near((scf.amo[0] - source.amo[0]).norm2(), 0, "frozen padding preserves exact alpha source");
        near((scf.bmo[0] - source.bmo[0]).norm2(), 0, "frozen padding preserves exact beta source");
        scf.amo[0].scale(1.0 + 1e-6);
        scf.bmo[0].scale(1.0 + 2e-6);
        const auto alpha = copy(world, scf.amo[0]), beta = copy(world, scf.bmo[0]);
        scf.solve_virtuals(world, control);
        near((scf.amo[0] - alpha).norm2(), 0, "KAIN update preserves off-unit alpha reference");
        near((scf.bmo[0] - beta).norm2(), 0, "KAIN update preserves off-unit beta reference");
    } else {
        MADNESS_CHECK_THROW(mode == "lifetime", "unknown review regression");
        if (world.size() != 1) return; // stdout probe is only a serial lifetime check
        const auto seed = std::weak_ptr(scf.amo.back().get_impl());
        scf.param.set_user_defined_value("print_level", 2);
        scf.param.set_user_defined_value("maxiter", 1);
        scf.param.set_user_defined_value("maxsub", 4);
        scf.param.set_user_defined_value("dconv", 0.0);
        int residual_lines = 0;
        bool released_before_beta = false;
        {
            LineProbe probe([&](const std::string& line) {
                if (line.starts_with("residuals") && ++residual_lines == 3)
                    released_before_beta = seed.expired();
            });
            scf.solve_virtuals(world);
        }
        check(residual_lines >= 3, "legacy lifetime fixture reaches the beta solve");
        check(released_before_beta, "legacy alpha KAIN history is released before beta evaluation");
    }
}

void exercise(World& world, int atom, const std::string& xc) {
    const std::string name = atom == 3 ? "li" : xc == "hf" ? "he" : "he_lda";
    Molecule molecule;
    molecule.add_atom(0.0, 0.0, 0.0, double(atom), atom);
    CalculationParameters p;
    p.set_user_defined_value("prefix", name);
    p.set_user_defined_value("xc", xc);
    p.set_user_defined_value("l", 20.0);
    p.set_user_defined_value("protocol", std::vector<double>{1e-4});
    p.set_user_defined_value("k", 6);
    p.set_user_defined_value("localize", std::string("canon"));
    p.set_user_defined_value("print_level", 0);
    if (atom == 3) p.set_user_defined_value("nopen", 1);
    SCF occupied(world, p, molecule);
    occupied.set_protocol<3>(world, 1e-4);
    occupied.ao = occupied.project_ao_basis(world, occupied.aobasis);
    occupied.make_nuclear_potential(world);
    occupied.initial_guess(world);
    it::Request ground_request;
    ground_request.targets = {{it::Metric::bsh_residual, 1e-3},
        {it::Metric::density_change, 1e-3}, {it::Metric::energy_change, 1e-5}};
    ground_request.max_updates = 30;
    const it::Protocol protocol{1e-4, 6};
    SCF::IterationControl ground_control{it::resolve(ground_request, occupied.iteration_model(), protocol, true), true,
        {"occupied", name + "-occupied-context", {"ground"}}};
    // Each endpoint is saved as a node, as SCFCheckpoint saves it.
    const auto run = std::filesystem::absolute(name + "-run");
    cp::detail::on_root(world, [&] { std::filesystem::remove_all(run); return std::string{}; });
    const auto publish = [&](const SCF& owner, const it::IterationMetrics& metrics,
                             const it::Decision&, const it::Progress& progress,
                             const std::filesystem::path& folder, const cp::json& lineage) {
        const auto staging = run / "staging";
        cp::detail::on_root(world, [&] { std::filesystem::create_directories(staging); return std::string{}; });
        cp::detail::collective_error(world, owner.write_mos(world, (staging / "scf.restartdata").string(),
                                                            owner.orbital_metadata()));
        world.gop.fence();
        cp::detail::on_root(world, [&] {
            const it::NodeDescription d{metrics.states.at(0).id == "virtuals" ? "virtuals" : "ground_state",
                metrics.info.protocol, progress.updates_at_protocol, {{"xc", xc}},
                it::node_metrics(metrics, SCF::iteration_metric_set()), lineage, {}, "test", false, {}};
            it::Node(folder).save(d, [&](const std::filesystem::path& archive) {
                for (const auto& file : std::filesystem::directory_iterator(staging))
                    std::filesystem::rename(file.path(), archive / file.path().filename());
            });
            std::filesystem::remove_all(staging);
            return std::string{};
        });
    };
    const auto saved_node = [&](const std::filesystem::path& folder) {
        const auto record = cp::json::parse(cp::detail::on_root(world, [&] {
            const auto d = it::Node(folder).inspect();
            return d ? it::description_record(*d).dump() : std::string("null");
        }));
        return record.is_null() ? std::nullopt : std::optional(it::read_description(record));
    };
    ground_control.observer = [&](World&, const SCF& owner, const it::IterationMetrics& metrics,
                                  const it::Decision& decision, const it::Progress& progress) {
        if (decision.checkpoint) publish(owner, metrics, decision, progress, run / "occupied", nullptr);
    };
    check(occupied.solve(world, ground_control).action == it::Action::complete,
          "fixed occupied source satisfies its explicitly selected request");

    // The existing native reader supplies padding; the source itself stays intact.
    p.set_user_defined_value("nvalpha", atom == 3 ? 2 : 1);
    if (atom == 3) p.set_user_defined_value("nvbeta", 1);
    p.set_user_defined_value("freeze_occupied", true);
    p.set_user_defined_value("maxiter", 0); // the typed request owns its allowance
    p.set_user_defined_value("prefix", (it::Node(run / "occupied").archive() / "scf").string());
    SCF scf(world, p, molecule);
    scf.set_protocol<3>(world, 1e-4);
    scf.make_nuclear_potential(world);
    scf.load_mos(world);
    for (int i = 0; i < scf.param.nalpha(); ++i)
        near((scf.amo[i] - occupied.amo[i]).norm2(), 0, "padded occupied alpha matches its published source");
    if (scf.param.have_beta())
        for (int i = 0; i < scf.param.nbeta(); ++i)
            near((scf.bmo[i] - occupied.bmo[i]).norm2(), 0, "padded occupied beta matches its published source");
    const auto alpha = copy(world, scf.amo);
    const auto beta = scf.param.have_beta() ? copy(world, scf.bmo) : vecfuncT{};
    const auto alpha_eps = copy(scf.aeps), beta_eps = copy(scf.beps);
    const double occupied_energy = scf.current_energy;
    const auto model = scf.iteration_model();
    // Targets for the fixed changes are written but cannot be judged here:
    // they are reported as not applicable, and only the residual decides.
    it::Request request;
    for (const auto metric : {it::Metric::bsh_residual, it::Metric::density_change,
                              it::Metric::energy_change, it::Metric::kinetic_change})
        request.targets[metric] = 0.0;
    request.floor_factor = 0.0; // a strict, never-met request
    request.max_updates = 1;
    request.checkpoint_interval = 1;
    for (const auto metric : {it::Metric::density_change, it::Metric::energy_change, it::Metric::kinetic_change})
        check(model.at(metric).support == it::Support::not_applicable, "fixed occupied changes are inapplicable");
    auto inapplicable = request;
    inapplicable.targets.erase(it::Metric::bsh_residual);
    const auto fixed_only = it::resolve(inapplicable, model, protocol, true);
    check(fixed_only.criteria.size() == 1 && fixed_only.criteria[0].metric == it::Metric::bsh_residual &&
          fixed_only.criteria[0].requested == 1e-4 && fixed_only.not_applicable.size() == 3,
          "fixed changes alone never converge virtuals: the residual default still decides");
    SCF::IterationControl control{it::resolve(request, model, protocol, true), true,
        {"virtuals", name + "-fixed-source-g1", {"virtuals"}}};
    int writes = 0;
    std::vector<it::IterationMetrics> records;
    std::vector<it::IterationTiming> intervals;
    control.timing_observer = [&](const it::IterationInfo& info, const it::IterationTiming& interval) {
        check(!records.empty() && info == records.back().info, "timing describes the evaluated virtual bundle");
        intervals.push_back(interval);
    };
    const cp::json inputs{{"fixed_occupied_reference", "occupied"}};
    control.observer = [&](World& w, const SCF& owner, const it::IterationMetrics& metrics,
                           const it::Decision& decision, const it::Progress& progress) {
        check(metrics.states.size() == 1 && metrics.states[0].id == "virtuals", "one task covers both virtual spin blocks");
        const auto& values = metrics.states[0].metrics;
        for (const auto metric : {it::Metric::density_change, it::Metric::energy_change, it::Metric::kinetic_change})
            check(values.at(metric).status == it::Availability::not_applicable && !values.at(metric).value,
                  "fixed occupied change is not reported as a measured zero");
        auto physical = scf.evaluate_orbitals(w);
        std::vector<double> independently_measured;
        const auto inspect = [&](vecfuncT& mo, tensorT& eps, const vecfuncT& original,
                                 const tensorT& original_eps, tensorT& full_fock,
                                 const vecfuncT& potential, int nocc) {
            for (int i = 0; i < nocc; ++i) {
                near((mo[i] - original[i]).norm2(), 0, "occupied functions are unchanged");
                near(eps[i], original_eps[i], "occupied eigenvalues are unchanged");
            }
            const auto nv = static_cast<int>(mo.size()) - nocc;
            if (nv == 0) return;
            vecfuncT functions(mo.begin() + nocc, mo.end());
            vecfuncT applied = copy(w, vecfuncT(potential.begin() + nocc, potential.end()));
            tensorT fock = copy(full_fock(Slice(nocc, -1), Slice(nocc, -1)));
            for (int i = 0; i < nv; ++i)
                near(eps[nocc + i], fock(i, i), "saved eigenvalues describe the prepared virtual functions");
            tensorT occupations(nv);
            double error = 0;
            std::vector<double> norms;
            scf.compute_residual(w, occupations, fock, functions, applied, error, &norms);
            independently_measured.insert(independently_measured.end(), norms.begin(), norms.end());
        };
        inspect(scf.amo, scf.aeps, alpha, alpha_eps, physical.fock_alpha, physical.potential_alpha, scf.param.nalpha());
        if (scf.param.have_beta())
            inspect(scf.bmo, scf.beps, beta, beta_eps, physical.fock_beta, physical.potential_beta, scf.param.nbeta());
        const auto& residual = values.at(it::Metric::bsh_residual);
        check(residual.evidence.components.size() == independently_measured.size(), "all active virtual functions have evidence");
        for (std::size_t i = 0; i < independently_measured.size(); ++i)
            near(residual.evidence.components[i].residual_norm, independently_measured[i], "residual evidence matches independent evaluation");
        near(*residual.value, *std::max_element(independently_measured.begin(), independently_measured.end()),
             "task residual is the maximum over both spins");
        near(owner.current_energy, occupied_energy, "virtual iteration preserves source energy");
        check(progress.main_updates == metrics.info.iteration &&
              progress.updates_at_protocol == control.previous_updates + metrics.info.iteration, "update counts preserve restart cadence");
        records.push_back(metrics);
        if (decision.checkpoint) { ++writes; publish(owner, metrics, decision, progress, run / "virtuals", inputs); }
    };
    auto bad = control;
    bad.task.states = {"ground"};
    rejects([&] { scf.solve_virtuals(world, bad); }, "frozen task cannot claim ordinary ground evidence");
    bad = control;
    bad.settings.criteria.push_back({it::Metric::density_change, 0.0, 0.0, 0.0, it::LimitRule::plain});
    rejects([&] { scf.solve_virtuals(world, bad); }, "caller cannot make a fixed density a criterion");
    bad = control; bad.settings.protocol.k = 8;
    rejects([&] { scf.solve_virtuals(world, bad); }, "wrong representation rejected");
    if (scf.param.have_beta()) {
        const auto complete_beta = scf.bmo;
        scf.bmo.resize(scf.param.nbeta());
        rejects([&] { scf.solve_virtuals(world, control); }, "missing requested beta virtuals cannot qualify a partial bundle");
        scf.bmo = complete_beta;
    }
    rejects([&] { scf.solve(world, ground_control); }, "ordinary entry cannot silently solve frozen task");
    check(records.empty(), "invalid controls never reach numerical observation");
    auto end = scf.solve_virtuals(world, control);
    check(end.action == it::Action::stop && end.reason == "iteration_limit", "strict virtual request respects update cap");
    check(records.size() == 2 && records.back().info.iteration == 1, "one KAIN update followed by terminal evaluation");
    check(intervals.size() == 2, "quiet frozen SCF forwards both measured intervals");
    if (intervals.size() == 2) {
        const auto has_update = [](const it::IterationTiming& interval) {
            return std::any_of(interval.steps.begin(), interval.steps.end(),
                               [](const auto& step) { return step.first == "advance"; });
        };
        check(has_update(intervals.front()) && !has_update(intervals.back()),
              "frozen terminal evaluation has no update timing");
    }
    check(!scf.last_calculation_outcome->complete && scf.converged_for_dconv > 1.0,
          "typed solve cannot leave a legacy success stamp");
    const auto loaded = saved_node(run / "virtuals");
    check(bool(loaded), "native virtual bundle and metrics are published");
    auto saved = it::stored_metrics(*loaded, {control.task.task_id, {}, 0, control.task.context_id});
    check(saved.info == records.back().info, "native evidence identifies exact virtual iteration");
    check(!it::saved_completion(control.settings, control.task, saved, true), "unmet saved request must iterate");
    scf.param.set_user_defined_value("prefix", (it::Node(run / "virtuals").archive() / "scf").string());
    scf.load_mos(world);
    control.previous_updates = 1;
    records.clear();
    end = scf.solve_virtuals(world, control);
    check(end.action == it::Action::stop && records.size() == 2 && records.front().info.iteration == 0,
          "native restart has a fresh allowance and fresh KAIN history");
    // A zero allowance still evaluates once; at an intermediate protocol it
    // permits advancement, but does not claim final completion.
    auto capped = request;
    capped.max_updates = 0;
    control.settings = it::resolve(capped, model, protocol, true);
    records.clear();
    end = scf.solve_virtuals(world, control);
    check(end.action == it::Action::stop && records.size() == 1 && end.info.iteration == 0,
          "strict zero allowance evaluates without applying KAIN");
    scf.param.set_user_defined_value("protocol", std::vector<double>{1e-4, 1e-6});
    control.settings = it::resolve(capped, model, protocol, true);
    control.final_protocol = false;
    records.clear();
    end = scf.solve_virtuals(world, control);
    check(end.action == it::Action::finish_protocol && records.size() == 1 && !scf.last_calculation_outcome->complete,
          "intermediate virtual endpoint can advance without final completion");
    scf.param.set_user_defined_value("protocol", std::vector<double>{1e-4});
    auto relaxed = request;
    relaxed.max_updates = 0;
    relaxed.targets[it::Metric::bsh_residual] = 2.0 * *records.back().states[0].metrics.at(it::Metric::bsh_residual).value;
    control.settings = it::resolve(relaxed, model, protocol, true);
    control.final_protocol = true;
    check(bool(it::saved_completion(control.settings, control.task, records.back(), true)), "sufficient saved virtual evidence reuses without evaluation");
    records.clear();
    end = scf.solve_virtuals(world, control);
    check(end.action == it::Action::complete && records.size() == 1 && end.info.iteration == 0,
          "satisfied virtual request completes without KAIN even with zero allowance");
    check(scf.last_calculation_outcome->complete && scf.converged_for_dconv > 1.0,
          "shared outcome succeeds without promoting a legacy stamp");
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
            char pattern[] = "/tmp/madness-scf-frozen-XXXXXX";
            if (const char* made = mkdtemp(pattern)) directory = made;
        }
        world.gop.broadcast_serializable(directory, 0);
        MADNESS_CHECK_THROW(!directory.empty(), "cannot create frozen SCF test directory");
        std::filesystem::current_path(directory);
        try {
            if (argc == 3 && std::string(argv[1]) == "--review") review_regression(world, argv[2]);
            else {
                exercise(world, 2, "hf");
                exercise(world, 3, "hf");
                exercise(world, 2, "lda");
            }
        } catch (const std::exception& e) { ++failures; std::cerr << e.what() << '\n'; }
        world.gop.fence();
        std::filesystem::current_path(original);
        world.gop.max(failures);
        if (world.rank() == 0 && failures == 0) std::filesystem::remove_all(directory);
    }
    finalize();
    return failures ? 1 : 0;
}
