#include <madness/chem/SCFRun.h>
#include <madness/chem/SCFCheckpoint.h>
#include <madness/chem/IterationOutput.h>
#include <madness/chem/SCFScientificSettings.h>
#include <madness/chem/TimingRecords.h>
#include <chrono>
#include <ctime>
#include <exception>
#include <type_traits>

namespace madness {
namespace {
namespace cp = checkpoint;
namespace it = iteration;
using json = cp::json;

void prepare(World& world, SCF& scf, it::Protocol protocol) {
    scf.set_protocol<3>(world, protocol.thresh);
    cp::records::require(FunctionDefaults<3>::get_k() == protocol.k, "inconsistent SCF protocol order");
    scf.make_nuclear_potential(world);
    scf.pcm_param.set_derived_values(scf.param);
    if (scf.param.pcm_data() != "none") scf.pcm = PCM(world, scf.molecule, scf.pcm_param, false);
}

std::vector<it::Protocol> protocols(const SCF& scf) {
    std::vector<it::Protocol> result;
    for (const double t : scf.param.protocol()) {
        const int k = scf.param.k() == -1 ? (t >= 0.9e-2 ? 4 : t >= 0.9e-4 ? 6 :
            t >= 0.9e-6 ? 8 : t >= 0.9e-8 ? 10 : 12) : scf.param.k();
        result.push_back({t, k});
    }
    cp::records::require(!result.empty(), "SCF protocol cannot be empty");
    return result;
}

std::string attempt(World& world) {
    return cp::detail::on_root(world, [] {
        static std::size_t sequence = 0;
        return "a" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) +
            "-" + std::to_string(++sequence);
    });
}

json read_preparation(World& world, const std::filesystem::path& path) {
    return json::parse(cp::detail::on_root(world, [&] {
        if (!std::filesystem::exists(path)) return json().dump();
        const auto value = cp::detail::read_json(path);
        cp::records::require(value.at("version") == 1, "unsupported SCF preparation record");
        return value.dump();
    }));
}
void write_preparation(World& world, const std::filesystem::path& path, const json& record) {
    cp::detail::on_root(world, [&] {
        const auto temporary = path.string() + ".tmp";
        std::ofstream file(temporary);
        file.exceptions(std::ios::badbit | std::ios::failbit);
        file << record.dump(2); file.close();
        std::filesystem::rename(temporary, path);
        return std::string{};
    });
}

void retire_preparation(World& world, cp::CheckpointStore& store, const std::filesystem::path& path,
                        const std::string& id) {
    if (id.empty()) return;
    // The marker makes a crash between inventory retirement and pointer removal
    // recoverable. Retirement and cleanup are idempotent.
    write_preparation(world, path, {{"version", 1}, {"calculation_id", id}, {"retiring", true}});
    store.retire_calculation(id);
    const auto warning = store.cleanup();
    if (!warning.empty() && world.rank() == 0) print(warning);
    write_preparation(world, path, {{"version", 1}});
}

void resize_bundle(World& world, SCF& scf, int alpha, int beta) {
    scf.param.set_user_defined_value("nmo_alpha", alpha);
    if (scf.param.have_beta()) scf.param.set_user_defined_value("nmo_beta", beta);
    const auto resize = [&](vecfuncT& functions, tensorT& eps, tensorT& occ, std::vector<int>& sets, int n) {
        if (static_cast<int>(functions.size()) > n) {
            functions.resize(n); sets.resize(n);
            eps = copy(eps(Slice(0, n - 1))); occ = copy(occ(Slice(0, n - 1)));
        }
    };
    resize(scf.amo, scf.aeps, scf.aocc, scf.aset, alpha);
    if (scf.param.have_beta()) resize(scf.bmo, scf.beps, scf.bocc, scf.bset, beta);
    scf.pad_virtuals_from_guess(world);
}

void report(World& world, cp::CheckpointStore& store, const std::string& owner, const std::string& calculation,
            const cp::Manifest& manifest, const SCF& scf, const std::optional<cp::Manifest>& occupied) {
    json states = json::object(), metrics = json::object();
    const auto add = [&](const cp::Manifest& saved) {
        const auto frame = cp::iteration_metrics_from_manifest(saved);
        for (const auto& state : frame.states) {
            states[state.id] = {{"task_id", frame.info.task_id}, {"subject", state.id},
                {"attempt_id", frame.info.attempt_id}, {"generation", saved.generation},
                {"protocol", {{"threshold", frame.info.protocol.thresh}, {"k", frame.info.protocol.k}}},
                {"iteration", frame.info.iteration}};
            json values = json::object(), unavailable = json::object();
            for (const auto& [metric, value] : state.metrics) {
                if (value.status == it::Availability::measured) values[it::metric_name(metric)] = *value.value;
                else unavailable[it::metric_name(metric)] = {{"status", value.status == it::Availability::not_applicable ?
                    "not_applicable" : "unavailable"}, {"reason", value.reason}};
            }
            metrics[state.id] = {{"values", values}, {"availability", unavailable}};
        }
    };
    if (occupied) add(*occupied);
    add(manifest);
    store.publish_source_map({{"version", 2}, {"calculation_id", calculation}, {"attempt_id", owner},
        {"states", states}, {"tasks", json::array()}});
    const bool complete = scf.last_calculation_outcome && scf.last_calculation_outcome->complete;
    json energy{{"units", "hartree"}, {"value_type", "real"}, {"status", complete ? "complete" : "unavailable"}};
    if (complete) energy["value"] = scf.current_energy;
    else energy["reason"] = "criteria_unmet";
    const auto warning = store.publish_results({{"version", 2}, {"attempt_id", owner}, {"calculations", json::array({
        {{"calculation_id", calculation}, {"application", "moldft"}, {"model", scf.param.xc()},
         {"protocol", manifest.record.at("task").at("definition").at("protocol")},
         {"properties", {{"total_energy", energy}}},
         {"state_metrics", json::array({{{"states", metrics}}})}}})}});
    if (!warning.empty() && world.rank() == 0) print(warning);
}
} // namespace

std::filesystem::path scf_checkpoint_directory(const SCF& scf) {
    if (!scf.checkpoint_root.empty()) return scf.checkpoint_root;
    return scf.work_dir.empty() ? std::filesystem::current_path() : scf.work_dir;
}

std::filesystem::path scf_checkpoint_archive(World& world, const std::filesystem::path& root, bool frozen) {
    cp::CheckpointStore store(world, root);
    const std::string task = frozen ? "virtuals" : "scf";
    const auto saved = store.read_latest("ground_state", task);
    cp::records::require(saved.has_value(), "no selected SCF checkpoint");
    return store.task_directory("ground_state", task) / saved->generation / "scf.restartdata";
}

std::filesystem::path scf_checkpoint_archive(World& world, const SCF& scf) {
    cp::records::require(scf.last_checkpoint.has_value(), "SCF has no qualified checkpoint identity");
    cp::records::require(scf.last_calculation_outcome && scf.last_calculation_outcome->complete &&
        scf.last_checkpoint->calculation_id == "ground_state", "SCF reference is not a completed requested bundle");
    cp::CheckpointStore store(world, scf_checkpoint_directory(scf));
    const auto& ref = *scf.last_checkpoint;
    const auto saved = store.read(ref); // exact generation; never follow a newer latest pointer
    return store.task_directory(ref.calculation_id, ref.task_id) / saved.generation / "scf.restartdata";
}

GroundReference scf_ground_reference(World& world, const SCF& scf) {
    const auto archive = scf_checkpoint_archive(world, scf);
    cp::CheckpointStore store(world, scf_checkpoint_directory(scf));
    const auto saved = store.read(*scf.last_checkpoint);
    const auto& definition = saved.record.at("task").at("definition");
    const auto protocol = cp::records::read_protocol(definition.at("protocol"));
    // The store keeps no save time; the manifest is written once, when the
    // generation is published, so its modification time is that time.
    const auto saved_at = cp::detail::on_root(world, [&] {
        const auto written = std::filesystem::last_write_time(archive.parent_path() / "manifest.json");
        const auto seconds = std::chrono::system_clock::to_time_t(
            std::chrono::time_point_cast<std::chrono::system_clock::duration>(
                std::chrono::file_clock::to_sys(written)));
        std::tm utc{};
        gmtime_r(&seconds, &utc);
        char text[32];
        std::strftime(text, sizeof text, "%Y-%m-%dT%H:%M:%SZ", &utc);
        return std::string(text);
    });
    return {archive, protocol, {{"settings", definition.at("scientific_settings")},
                                {"protocol", cp::records::protocol_record(protocol)}, {"saved_at", saved_at}}};
}

void require_scf_restart(World& world, const std::filesystem::path& root, bool frozen) {
    cp::CheckpointStore store(world, root);
    if (frozen && store.read_latest("ground_state", "scf")) return;
    const std::string task = frozen ? "virtuals" : "scf";
    if (store.read_latest("ground_state", task)) return;
    const auto progress = read_preparation(world, root / "checkpoints" / (task + "_preparation.json"));
    cp::records::require(!progress.is_null() && progress.contains("calculation_id") &&
        !progress.value("retiring", false), "no selected SCF or preparation checkpoint");
    cp::records::require(store.read_latest(progress.at("calculation_id"), task).has_value(),
        "no selected SCF preparation checkpoint");
}

void inspect_scf(World& world, const std::filesystem::path& root, const CalculationParameters& plotting) {
    cp::CheckpointStore store(world, root);
    const std::string task = plotting.freeze_occupied() ? "virtuals" : "scf";
    const auto saved = store.read_latest("ground_state", task);
    cp::records::require(saved.has_value(), "no selected SCF checkpoint");
    cp::records::require(saved->archives.size() == 1 && saved->archives.front().role == "orbitals" &&
        saved->archives.front().basename == "scf.restartdata", "unsupported SCF native payload");
    if (world.rank() == 0) print(saved->record.dump(2));
    if (!plotting.get<bool>("plotdens") && !plotting.get<bool>("plotcoul") &&
        plotting.get<int>("plothi") < plotting.get<int>("plotlo")) return;
    cp::records::require(saved->archives.front().writers <= world.size(), "insufficient native archive reader ranks");
    const auto directory = store.task_directory("ground_state", task) / saved->generation;
    const auto metadata = peek_restartdata(world, (directory / "scf.restartdata").string());
    cp::records::require(metadata.has_value(), "missing native SCF metadata");
    const auto& definition = saved->record.at("task").at("definition");
    const auto protocol = cp::records::read_protocol(definition.at("protocol"));
    auto parameters = plotting;
    parameters.from_json(definition.at("scientific_settings").at("dft"));
    parameters.set_user_defined_value("k", protocol.k);
    parameters.set_user_defined_value("protocol", std::vector<double>{protocol.thresh});
    const auto& layout = definition.at("orbitals");
    for (const auto* key : {"nmo_alpha", "nmo_beta", "nalpha", "nbeta"})
        parameters.set_user_defined_value(key, layout.at(key).get<int>());
    parameters.set_user_defined_value("nvalpha", layout.at("nmo_alpha").get<int>() - layout.at("nalpha").get<int>());
    parameters.set_user_defined_value("nvbeta", std::max(0, layout.at("nmo_beta").get<int>() - layout.at("nbeta").get<int>()));
    parameters.set_user_defined_value("freeze_occupied", layout.at("freeze_occupied").get<bool>());
    // Inspection loads the saved representation; current convergence requests
    // are not assessed and no operator or metric evaluation is needed.
    SCF scf(world, parameters, metadata->molecule);
    scf.set_protocol<3>(world, protocol.thresh);
    scf.make_nuclear_potential(world);
    scf.param.set_user_defined_value("prefix", (directory / "scf").string());
    scf.load_mos(world, false);
    scf.param.set_user_defined_value("prefix", plotting.prefix());
    scf.do_plots(world);
}

// The single recursive call qualifies an ordinary occupied source; that branch
// explicitly clears freeze_occupied, so recursion is bounded to one level.
// NOLINTNEXTLINE(misc-no-recursion)
double run_scf(World& world, SCF& scf) {
    scf.last_checkpoint.reset();
    scf.last_calculation_outcome.reset();
    cp::records::require(scf.param.restart() != "read_only", "read_only is inspection, not a solver request");
    const auto request = scf_iteration_request(scf.param);
    const auto ladder = protocols(scf);
    // Validate the request at every level before any work; warn once about
    // final-level targets below the protocol floor.
    for (std::size_t i = 0; i < ladder.size(); ++i) {
        const bool final_level = i + 1 == ladder.size();
        const auto settings = it::resolve(request, scf.iteration_model(), ladder[i], final_level, "ground state");
        if (final_level && world.rank() == 0 &&
            it::prints(it::verbosity(scf.param.print_level()), it::needed::targets_below_floor))
            it::write_targets_below_floor(std::cout, "SCF", scf.param.freeze_occupied() ? "virtuals" : "ground_state",
                                     settings);
    }
    const auto root = scf_checkpoint_directory(scf);
    cp::CheckpointStore store(world, root);
    const auto owner = attempt(world);
    store.initialize(owner);
    const std::string calculation = "ground_state";
    it::TimingCollector<> timing;
    timing.begin("run", {{"clock_scope", "driver_rank"}});
    timing.begin("scf", {{"calculation_id", calculation}});
    // Occupied qualification recurses once with freeze_occupied cleared.
    // NOLINTNEXTLINE(misc-no-recursion)
    const auto measured = [&](const char* name, auto&& operation, const json& context = json::object()) {
        timing.begin(name, context);
        // A failed operation leaves its scope open for the failure snapshot.
        if constexpr (std::is_void_v<std::invoke_result_t<decltype(operation)>>) {
            operation();
            timing.end();
        } else {
            auto result = operation();
            timing.end();
            return result;
        }
    };
    const auto observe_iteration = [&](const it::IterationInfo& info, const it::IterationTiming& interval) {
        json children = json::array();
        for (const auto& [name, seconds] : interval.steps)
            children.push_back({{"name", name}, {"seconds", seconds}});
        timing.append({{"name", "iteration"}, {"iteration", info.iteration},
            {"task_id", info.task_id}, {"attempt_id", info.attempt_id},
            {"seconds", interval.seconds}, {"children", children}});
    };
    const auto publish_timing = [&](const std::string& failure = {}) {
        const it::Execution execution{static_cast<std::size_t>(world.size()),
                                      static_cast<std::size_t>(ThreadPool::size() + 1)};
        store.publish_timings(it::timing_entry(owner, execution, timing.record(), failure));
    };
    // NOLINTNEXTLINE(misc-no-recursion)
    const auto execute = [&]() -> double {
        const std::string task = scf.param.freeze_occupied() ? "virtuals" : "scf";
        const std::string context = "scf-physical-v1";
        const auto mode = scf.param.restart();
        const auto preparation_file = root / "checkpoints" / (task + "_preparation.json");
        auto progress = read_preparation(world, preparation_file);
        std::string preparation_id;
        if (!progress.is_null() && progress.contains("calculation_id")) preparation_id = progress.at("calculation_id");
        if (!preparation_id.empty() && progress.value("retiring", false)) {
            retire_preparation(world, store, preparation_file, preparation_id);
            preparation_id.clear(); progress = nullptr;
        }
        timing.begin("prepare");
        scf.xc.initialize(scf.param.xc(), !scf.param.spin_restricted(), world, false);
        scf.dispersion = DispersionCorrection(scf.param.dispersion(), scf.param.xc(),
                                            scf.param.dispersion_functional(), scf.param.dispersion_atm());

        json inputs = json::array();
        std::optional<cp::Manifest> occupied_manifest;
        std::shared_ptr<SCF> occupied;
        if (scf.param.freeze_occupied()) {
            cp::records::require(mode != "none" && mode != "ao" && mode != "nwchem", "frozen SCF requires an occupied checkpoint");
            auto p = scf.param;
            p.set_user_defined_value("freeze_occupied", false);
            p.set_user_defined_value("nvalpha", 0); p.set_user_defined_value("nvbeta", 0);
            p.set_user_defined_value("nmo_alpha", p.nalpha()); p.set_user_defined_value("nmo_beta", p.nbeta());
            p.set_user_defined_value("nv_extra", 0);
            p.set_user_defined_value("restart", std::string("iterate"));
            occupied = std::make_shared<SCF>(world, p, scf.molecule, scf.pcm_param);
            occupied->work_dir = scf.work_dir;
            occupied->checkpoint_root = root;
            // NOLINTNEXTLINE(misc-no-recursion)
            measured("occupied_reference", [&] { run_scf(world, *occupied); });
            cp::records::require(occupied->last_calculation_outcome->complete, "occupied SCF does not satisfy the current request");
            cp::records::require(occupied->last_checkpoint.has_value(), "occupied SCF did not select a checkpoint");
            occupied_manifest = store.read(*occupied->last_checkpoint);
            inputs.push_back({{"role", "occupied"}, {"calculation_id", calculation}, {"task_id", "scf"},
                              {"generation", occupied_manifest->generation}});
            store.initialize(owner);
        }
        prepare(world, scf, ladder.front());
        timing.end();
        const int target_alpha = scf.param.nmo_alpha(), target_beta = scf.param.nmo_beta();
        std::optional<cp::Manifest> preparation_saved;
        if (!preparation_id.empty()) preparation_saved = store.read_latest(preparation_id, task);
        const bool use_preparation = preparation_saved && (mode == "auto" || mode == "iterate") &&
            preparation_saved->record.at("task").at("inputs") == inputs;
        std::size_t start = 0;
        std::size_t previous_updates = 0;
        bool padded_or_new = false;
        {
            SCFCheckpoint binding(world, scf, store, {calculation, task, owner, context}, inputs);
            const auto saved = store.read_latest(calculation, task);
            const bool source_changed = saved && saved->record.at("task").at("inputs") != inputs;
            const bool seed = mode == "none" || mode == "ao" || mode == "nwchem" || source_changed;
            if (!seed) {
                const auto restart = measured("checkpoint_planning", [&] { return binding.plan(ladder, request); });
                if (restart.plan.action != it::RestartAction::fresh) {
                    start = restart.plan.start_protocol;
                    measured("prepare", [&] { prepare(world, scf, ladder[start]); });
                    measured("load", [&] { binding.load(restart); }, {{"task_id", task}});
                    previous_updates = restart.plan.updates_at_protocol;
                    padded_or_new = restart.checkpoint->record.at("task").at("definition").at("orbitals").at("nmo_alpha").get<int>() < target_alpha ||
                        (scf.param.have_beta() && restart.checkpoint->record.at("task").at("definition").at("orbitals").at("nmo_beta").get<int>() < target_beta);
                    if (restart.plan.action == it::RestartAction::reuse) {
                        measured("publish_results", [&] {
                            report(world, store, owner, calculation, *restart.checkpoint, scf, occupied_manifest);
                            retire_preparation(world, store, preparation_file, preparation_id);
                        });
                        return scf.current_energy;
                    }
                } else padded_or_new = true;
            } else padded_or_new = true;
            if (padded_or_new && scf.param.freeze_occupied())
                measured("load_guess", [&] { binding.load_guess(*occupied_manifest); });
            else if (padded_or_new && (seed || !saved) && !use_preparation) {
                timing.begin("initial_guess");
                auto legacy_mode = restart_mode_from_string(mode);
                if (mode == "auto") {
                    const auto exists = cp::detail::on_root(world, [&] {
                        return std::filesystem::exists(scf.param.prefix() + ".restartdata.00000") ? "yes" : "no";
                    });
                    legacy_mode = restart_mode_from_string(exists == "yes" ? "iterate" : "none");
                }
                auto guess = make_restart_plan(world, legacy_mode, scf.param, scf.molecule,
                    scf.restart_representation, RestartCapabilities::all(), scf.restart_ncf);
                guess.protocol_start = 0; guess.iterate = true;
                scf.ao = scf.project_ao_basis(world, scf.aobasis);
                scf.get_initial_orbitals(world, guess);
                if (guess.source == RestartSource::restartdata) {
                    const auto defaults = [](int n, int occupied_count) {
                        tensorT values(n);
                        for (int i = 0; i < occupied_count; ++i) values[i] = 1.0;
                        return values;
                    };
                    scf.aocc = defaults(scf.param.nmo_alpha(), scf.param.nalpha());
                    if (scf.param.have_beta()) scf.bocc = defaults(scf.param.nmo_beta(), scf.param.nbeta());
                    scf.apply_explicit_occupations(world);
                    scf.aset = scf.group_orbital_sets(world, scf.aeps, scf.aocc, scf.param.nmo_alpha());
                    if (scf.param.have_beta()) scf.bset = scf.group_orbital_sets(world, scf.beps, scf.bocc, scf.param.nmo_beta());
                }
                padded_or_new = guess.source == RestartSource::initial_guess || scf.last_load_padded;
                scf.converged_for_thresh = scf.converged_for_dconv = 1e10;
                scf.last_calculation_outcome = CalculationOutcome{false, "fresh_evaluation_required"};
                timing.end();
            }
        }
        const int extra = (padded_or_new || !preparation_id.empty()) ? scf.param.nv_extra() : 0;
        if ((extra > 0 || !preparation_id.empty()) && start == 0) {
            cp::records::require(scf.param.nv_its() >= 0, "negative virtual preparation budget");
            auto prep_request = request; prep_request.max_updates = scf.param.nv_its();
            const int step = scf.param.nv_step() > 0 ? scf.param.nv_step() : std::max(extra, 1);
            int count = target_alpha + extra;
            auto saved = use_preparation ? preparation_saved : std::nullopt;
            const json schedule{{"nv_extra", scf.param.nv_extra()}, {"nv_step", scf.param.nv_step()},
                {"nv_its", scf.param.nv_its()}, {"target_alpha", target_alpha}, {"target_beta", target_beta},
                {"protocol", cp::records::protocol_record(ladder.front())}};
            // Embed the schedule in the saved evidence as well as the pointer. A
            // crash after writing a new pointer cannot relabel an old completion.
            const auto preparation_context = context + ":preparation:" + schedule.dump();
            const bool resume = saved && progress.value("schedule", json()) == schedule &&
                saved->record.at("task").at("definition").at("context_id") == preparation_context;
            if (resume)
                count = saved->record.at("task").at("definition").at("orbitals").at("nmo_alpha");
            else {
                if (preparation_id.empty()) preparation_id = calculation + "-preparation-" + owner;
            }
            scf.param.set_user_defined_value("nmo_alpha", count);
            if (scf.param.have_beta()) scf.param.set_user_defined_value("nmo_beta", target_beta + count - target_alpha);
            if (!saved) resize_bundle(world, scf, count, target_beta + count - target_alpha);
            // Each membership is a separate solve attempt: iteration zero at a
            // smaller bundle must not share an identity with the previous bundle.
            const auto preparation_attempt = [&](int n) { return owner + "-n" + std::to_string(n); };
            SCFCheckpoint first(world, scf, store, {preparation_id, task, preparation_attempt(count), preparation_context}, inputs);
            write_preparation(world, preparation_file, {{"version", 1}, {"calculation_id", preparation_id},
                {"schedule", schedule}});
            std::size_t preparation_updates = 0;
            if (saved) {
                timing.begin("orbital_preparation", {{"calculation_id", preparation_id}, {"task_id", task},
                    {"protocol", cp::records::protocol_record(ladder.front())}});
                const auto restart = measured("checkpoint_planning", [&] { return first.plan({ladder.front()}, prep_request); });
                measured("load", [&] { first.load(restart, true); }); // intermediate preparation reporting only
                preparation_updates = restart.plan.updates_at_protocol;
                const auto settings = scf_scientific_settings(scf);
                if (resume && saved->record.at("task").at("iteration").at("protocol_complete") == true &&
                    saved->record.at("task").at("definition").at("scientific_settings") == settings &&
                    scf_settings_allow_reuse(settings)) {
                    if (world.rank() == 0 && scf.param.print_level() > 0)
                        print("Completed preparation bundle; continue with a smaller orbital bundle.");
                    count = std::max(target_alpha, count - step);
                    preparation_updates = 0;
                } else if (world.rank() == 0 && scf.param.print_level() > 0) {
                    print("Continue preparation on the loaded orbital bundle.");
                }
                timing.end();
            }
            for (; count > target_alpha; count = std::max(target_alpha, count - step)) {
                timing.begin("orbital_preparation", {{"calculation_id", preparation_id}, {"orbitals_alpha", count},
                    {"protocol", cp::records::protocol_record(ladder.front())}});
                timing.begin("prepare");
                resize_bundle(world, scf, count, target_beta + count - target_alpha);
                scf.reset_aobasis("sto-3g"); scf.ao = scf.project_ao_basis(world, scf.aobasis);
                SCFCheckpoint binding(world, scf, store, {preparation_id, task, preparation_attempt(count), preparation_context}, inputs);
                auto control = binding.control(prep_request, ladder.front(), false, preparation_updates);
                control.preparing_virtuals = true;
                control.timing_observer = observe_iteration;
                timing.end();
                const auto decision = measured("task", [&] {
                    return scf.param.freeze_occupied() ? scf.solve_virtuals(world, control) : scf.solve(world, control);
                }, {{"task_id", task}});
                timing.end();
                cp::records::require(decision.action == it::Action::finish_protocol, "SCF virtual preparation stopped; endpoint saved");
                preparation_updates = 0;
            }
            resize_bundle(world, scf, target_alpha, target_beta);
        }
        for (std::size_t stage = start; stage < ladder.size(); ++stage) {
            timing.begin("protocol", {{"protocol", cp::records::protocol_record(ladder[stage])}});
            timing.begin("prepare");
            prepare(world, scf, ladder[stage]);
            if (stage > start) scf.project(world);
            if (occupied) {
                const auto restore = [&](const vecfuncT& source, vecfuncT& target,
                                         const tensorT& source_eps, tensorT& target_eps) {
                    auto functions = copy(world, source);
                    if (ladder[stage] != ladder.back()) {
                        for (auto& function : functions)
                            function = madness::project(function, ladder[stage].k, ladder[stage].thresh, false);
                        world.gop.fence();
                    }
                    std::copy(functions.begin(), functions.end(), target.begin());
                    for (std::size_t i = 0; i < source_eps.size(); ++i)
                        target_eps[static_cast<long>(i)] = source_eps[static_cast<long>(i)];
                };
                restore(occupied->amo, scf.amo, occupied->aeps, scf.aeps);
                if (scf.param.have_beta()) {
                    restore(occupied->bmo, scf.bmo, occupied->beps, scf.beps);
                }
                scf.current_energy = occupied->current_energy;
            }
            scf.reset_aobasis("sto-3g"); scf.ao.clear(); world.gop.fence();
            scf.ao = scf.project_ao_basis(world, scf.aobasis);
            SCFCheckpoint binding(world, scf, store, {calculation, task, owner, context}, inputs);
            auto control = binding.control(request, ladder[stage], stage + 1 == ladder.size(),
                                                 stage == start ? previous_updates : 0);
            control.timing_observer = observe_iteration;
            timing.end();
            const auto decision = measured("task", [&] {
                return scf.param.freeze_occupied() ? scf.solve_virtuals(world, control) : scf.solve(world, control);
            }, {{"task_id", task}});
            timing.end();
            if (decision.action != it::Action::finish_protocol) break;
        }
        cp::records::require(scf.last_checkpoint.has_value(), "SCF ended without a shared endpoint checkpoint");
        const auto last = store.read(*scf.last_checkpoint);
        measured("publish_results", [&] {
            report(world, store, owner, calculation, last, scf, occupied_manifest);
            retire_preparation(world, store, preparation_file, preparation_id);
        });
        return scf.current_energy;
    };
    try {
        const auto energy = execute();
        timing.end(); // scf
        timing.end(); // run
        publish_timing();
        return energy;
    } catch (...) {
        const auto original = std::current_exception();
        std::string reason = "unknown SCF failure";
        try { std::rethrow_exception(original); }
        catch (const std::exception& error) { reason = error.what(); }
        catch (...) {}
        try { publish_timing(reason); } catch (...) {} // preserve the solver/publication failure
        std::rethrow_exception(original);
    }
}

double MolecularEnergy::value(const Tensor<double>& x) {
    if (!calc.param.scf_input_owner || calc.checkpoint_observer) return legacy_value(x);
    calc.molecule.set_all_coords(x.reshape(static_cast<long>(calc.molecule.natom()), 3));
    return run_scf(world, calc);
}
} // namespace madness
