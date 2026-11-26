#include <madness/chem/SCFRun.h>
#include <madness/chem/SCFCheckpoint.h>
#include <madness/chem/RunFiles.h>
#include <madness/chem/IterationOutput.h>
#include <madness/chem/SCFScientificSettings.h>
#include <madness/chem/TimingRecords.h>
#include <madness/chem/IterationHistory.h>
#include <madness/misc/info.h>
#include <unistd.h>
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

// One node's measured values, keyed in results.json by its folder name.
json state_metrics(const it::NodeDescription& saved) {
    json subjects = json::object();
    for (const auto& state : it::stored_metrics(saved, {"scf", saved.protocol, 0, "saved"}).states) {
        json values = json::object(), unavailable = json::object();
        for (const auto& [metric, value] : state.metrics) {
            if (value.status == it::Availability::measured) values[it::metric_name(metric)] = *value.value;
            else unavailable[it::metric_name(metric)] = {{"status", value.status == it::Availability::not_applicable ?
                "not_applicable" : "unavailable"}, {"reason", value.reason}};
        }
        subjects[state.id] = {{"values", values}, {"availability", unavailable}};
    }
    return {{"protocol", cp::records::protocol_record(saved.protocol)}, {"subjects", subjects}};
}

// A tensor's elements in row-major order.
json tensor_values(const Tensor<double>& t) {
    const Tensor<double> flat = copy(t);
    json values = json::array();
    for (long i = 0; i < flat.size(); ++i) values.push_back(flat.ptr()[i]);
    return values;
}

void report(World& world, const std::filesystem::path& root, const std::filesystem::path& folder, const SCF& scf,
            const std::optional<it::NodeDescription>& occupied) {
    const auto saved = read_node(world, folder);
    cp::records::require(saved.has_value(), "SCF ended without a saved state");
    const auto name = folder.filename().string();
    json metrics = json::object();
    if (occupied) metrics[ground_state_node(root).filename().string()] = state_metrics(*occupied);
    metrics[name] = state_metrics(*saved);
    const bool complete = scf.last_calculation_outcome && scf.last_calculation_outcome->complete;
    const auto property = [&](const char* units, const json& value) {
        json p{{"units", units}, {"value_type", "real"}, {"status", complete ? "complete" : "unavailable"}};
        if (complete) p["value"] = value;
        else p["availability"] = {{"value", name + " has not met the requested criteria"}};
        return p;
    };
    json orbitals{{"alpha", tensor_values(scf.aeps)}};
    if (!scf.param.spin_restricted() && scf.param.have_beta()) orbitals["beta"] = tensor_values(scf.beps);
    json properties{{"total_energy", property("hartree", scf.current_energy)},
                    {"orbital_energies", property("hartree", orbitals)}};
    // The final energy terms, when this run evaluated the energy.
    const auto terms = scf.e_data.final_values();
    if (!terms.empty()) {
        json components{{"kinetic", terms.at("e_kinetic")}, {"nuclear_attraction", terms.at("e_nuclear")},
                        {"coulomb", terms.at("e_coulomb")}, {"exchange_correlation", terms.at("e_xc")},
                        {"nuclear_repulsion", terms.at("e_nrep")}};
        if (scf.param.pcm_data() != "none") components["pcm"] = terms.at("e_pcm");
        properties["energy_components"] = property("hartree", components);
    }
    if (scf.dispersion.active())
        properties["dispersion_energy"] = property("hartree", scf.dispersion.energy(world, scf.molecule));
    const auto request = scf_iteration_request(scf.param);
    const auto final = protocols(scf).back();
    // Every ground-state value comes from the one saved state: the entry
    // carries its convergence once, and the properties carry none.
    json entry{{"calculation_id", "ground_state"}, {"application", "moldft"},
        {"model", scf.param.xc()},
        {"protocol", cp::records::protocol_record(saved->protocol)},
        {"targets", it::results_targets(it::resolve(request, scf.iteration_model(), final, true, "ground state"))},
        {"properties", properties},
        {"state_metrics", metrics}};
    it::add_results_source(entry, name);
    it::write_results_entry(world, root, std::move(entry));
}
} // namespace

std::filesystem::path scf_checkpoint_directory(const SCF& scf) {
    if (!scf.checkpoint_root.empty()) return scf.checkpoint_root;
    return scf.work_dir.empty() ? std::filesystem::current_path() : scf.work_dir;
}

std::filesystem::path scf_checkpoint_archive(World& world, const std::filesystem::path& root, bool frozen) {
    const auto folder = frozen ? virtuals_node(root) : ground_state_node(root);
    cp::records::require(read_node(world, folder).has_value(), "no saved SCF state");
    return it::Node(folder).archive() / "scf.restartdata";
}

namespace {
// The qualified state of a completed SCF, never a step-down state: its bundle
// is exactly the requested one.
it::NodeDescription qualified(World& world, const SCF& scf) {
    cp::records::require(scf.last_node.has_value(), "SCF has no qualified saved state");
    cp::records::require(scf.last_calculation_outcome && scf.last_calculation_outcome->complete,
                         "SCF reference is not a completed requested bundle");
    const auto saved = read_node(world, *scf.last_node);
    cp::records::require(saved.has_value(), "the qualified SCF state is missing");
    const auto& orbitals = saved->details.at("orbitals");
    cp::records::require(orbitals.at("nmo_alpha") == scf.param.nmo_alpha() &&
        orbitals.at("nmo_beta") == (scf.param.have_beta() ? scf.param.nmo_beta() : 0),
        "the saved SCF state is not the requested bundle");
    return *saved;
}
} // namespace

std::filesystem::path scf_checkpoint_archive(World& world, const SCF& scf) {
    (void)qualified(world, scf);
    return it::Node(*scf.last_node).archive() / "scf.restartdata";
}

GroundReference scf_ground_reference(World& world, const SCF& scf) {
    const auto saved = qualified(world, scf);
    return {it::Node(*scf.last_node).archive() / "scf.restartdata", saved.protocol, node_identity(saved)};
}

void report_scf_properties(World& world, const SCF& scf, const Tensor<double>& dipole,
                           const Tensor<double>& gradient, const json& populations) {
    // A run that saved no node (the legacy solve) writes no results.json.
    if ((dipole.size() == 0 && gradient.size() == 0 && populations.empty()) || !scf.last_node) return;
    // Their convergence is the entry's: they come from its one saved state.
    const auto property = [](const char* units, const json& value) {
        return json{{"units", units}, {"value_type", "real"}, {"status", "complete"}, {"value", value}};
    };
    json properties = json::object();
    if (dipole.size() > 0) properties["dipole"] = property("atomic_units", tensor_values(dipole));
    if (gradient.size() > 0) properties["gradient"] = property("hartree/bohr", tensor_values(gradient));
    if (!populations.empty()) properties["populations"] = property("e", populations);
    it::add_results_properties(world, scf_checkpoint_directory(scf), "ground_state", properties);
}

bool has_saved_ground_state(World& world, const std::filesystem::path& root) {
    return read_node(world, ground_state_node(root)).has_value();
}

void require_scf_restart(World& world, const std::filesystem::path& root) {
    cp::records::require(has_saved_ground_state(world, root), "--continue needs a saved ground state");
}

void inspect_scf(World&, const std::filesystem::path&, const CalculationParameters&) {
    // It returns as a loop over the node descriptions of a run.
    throw std::runtime_error("inspection of saved states is not available in this version");
}

// The single recursive call qualifies an ordinary occupied source; that branch
// explicitly clears freeze_occupied, so recursion is bounded to one level.
// NOLINTNEXTLINE(misc-no-recursion)
double run_scf(World& world, SCF& scf) {
    scf.last_node.reset();
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
    it::require_current_run_folder(world, root);
    it::TimingCollector<> timing;
    timing.begin("run", {{"clock_scope", "driver_rank"}});
    timing.begin("scf", {{"calculation_id", "ground_state"}});
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
            {"task_id", info.task_id},
            {"seconds", interval.seconds}, {"children", children}});
    };
    const auto publish_timing = [&](const std::string& failure = {}) {
        const it::Execution execution{static_cast<std::size_t>(world.size()),
                                      static_cast<std::size_t>(ThreadPool::size() + 1)};
        it::write_timing_entry(world, root, it::timing_entry(it::launch_started_at(), execution, timing.record(), failure));
    };
    // NOLINTNEXTLINE(misc-no-recursion)
    const auto execute = [&]() -> double {
        const bool frozen = scf.param.freeze_occupied();
        const std::string task = frozen ? "virtuals" : "scf";
        const std::string context = "scf-physical-v1";
        const auto mode = scf.param.restart();
        const auto folder = frozen ? virtuals_node(root) : ground_state_node(root);
        write_launch_header(world, folder);
        timing.begin("prepare");
        scf.xc.initialize(scf.param.xc(), !scf.param.spin_restricted(), world, false);
        scf.dispersion = DispersionCorrection(scf.param.dispersion(), scf.param.xc(),
                                            scf.param.dispersion_functional(), scf.param.dispersion_atm());

        // The virtuals node's lineage is the identity of the ground state it was solved in.
        json lineage;
        std::optional<it::NodeDescription> occupied_saved;
        std::shared_ptr<SCF> occupied;
        if (frozen) {
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
            occupied_saved = qualified(world, *occupied);
            lineage = node_identity(*occupied_saved);
        }
        prepare(world, scf, ladder.front());
        timing.end();
        // A ground state started from a seed names it as its lineage; a restart
        // keeps the lineage of the state it continues.
        const bool seed = mode == "none" || mode == "ao" || mode == "nwchem";
        if (!frozen && !seed)
            if (const auto previous = read_node(world, folder)) lineage = previous->lineage;
        const int target_alpha = scf.param.nmo_alpha(), target_beta = scf.param.nmo_beta();
        std::size_t start = 0;
        std::size_t previous_updates = 0;
        bool padded_or_new = false;
        // A saved state with more orbitals than requested is an unfinished
        // step-down at the first level: continue the step-down from it.
        std::optional<it::NodeDescription> step_down;
        {
            SCFCheckpoint binding(world, scf, folder, {task, context}, lineage);
            const auto saved = binding.saved();
            if (!seed && saved && saved->details.at("orbitals").at("nmo_alpha").get<int>() > target_alpha &&
                saved->protocol == ladder.front() && saved->lineage == lineage)
                step_down = saved;
            else if (!seed) {
                const auto restart = measured("checkpoint_planning", [&] { return binding.plan(ladder, request); });
                if (restart.plan.action != it::RestartAction::fresh) {
                    start = restart.plan.start_protocol;
                    measured("prepare", [&] { prepare(world, scf, ladder[start]); });
                    measured("load", [&] { binding.load(restart); }, {{"task_id", task}});
                    previous_updates = restart.plan.updates_at_protocol;
                    const auto& orbitals = restart.saved->details.at("orbitals");
                    padded_or_new = orbitals.at("nmo_alpha").get<int>() < target_alpha ||
                        (scf.param.have_beta() && orbitals.at("nmo_beta").get<int>() < target_beta);
                    if (restart.plan.action == it::RestartAction::reuse) {
                        measured("publish_results", [&] { report(world, root, folder, scf, occupied_saved); });
                        return scf.current_energy;
                    }
                } else padded_or_new = true;
            } else padded_or_new = true;
            if (padded_or_new && frozen)
                measured("load_guess", [&] { binding.load_guess(ground_state_node(root)); });
            else if (padded_or_new && (seed || !saved)) {
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
                if (guess.source == RestartSource::restartdata)
                    lineage = {{"seed", scf.param.prefix() + ".restartdata"}};
                else if (guess.source == RestartSource::restartao || guess.source == RestartSource::nwchem)
                    lineage = {{"seed", guess.source == RestartSource::nwchem ? "nwchem" : scf.param.prefix() + ".restartaodata"}};
                else lineage = nullptr;
                padded_or_new = guess.source == RestartSource::initial_guess || scf.last_load_padded;
                scf.converged_for_thresh = scf.converged_for_dconv = 1e10;
                scf.last_calculation_outcome = CalculationOutcome{false, "fresh_evaluation_required"};
                timing.end();
            }
        }
        // The nv_extra step-down: converge extra virtuals first, then drop them
        // stepwise. Each stage saves into the same node as an unfinished state
        // with its larger bundle; the final save replaces it.
        const int extra = padded_or_new ? scf.param.nv_extra() : 0;
        if ((extra > 0 || step_down) && start == 0) {
            cp::records::require(scf.param.nv_its() >= 0, "negative virtual preparation budget");
            auto prep_request = request; prep_request.max_updates = scf.param.nv_its();
            const int step = scf.param.nv_step() > 0 ? scf.param.nv_step() : std::max(scf.param.nv_extra(), 1);
            int count = step_down ? step_down->details.at("orbitals").at("nmo_alpha").get<int>() : target_alpha + extra;
            const json schedule{{"nv_extra", scf.param.nv_extra()}, {"nv_step", scf.param.nv_step()},
                {"nv_its", scf.param.nv_its()}, {"target_alpha", target_alpha}, {"target_beta", target_beta},
                {"protocol", cp::records::protocol_record(ladder.front())}};
            // In memory only: a step-down stage never compares with the final bundle.
            const auto preparation_context = context + ":preparation:" + schedule.dump();
            // Each membership is a separate solve, so a change is never measured
            // against the state of another bundle.
            scf.param.set_user_defined_value("nmo_alpha", count);
            if (scf.param.have_beta()) scf.param.set_user_defined_value("nmo_beta", target_beta + count - target_alpha);
            if (!step_down) resize_bundle(world, scf, count, target_beta + count - target_alpha);
            std::size_t preparation_updates = 0;
            if (step_down) {
                timing.begin("orbital_preparation", {{"orbitals_alpha", count}, {"task_id", task},
                    {"protocol", cp::records::protocol_record(ladder.front())}});
                SCFCheckpoint first(world, scf, folder, {task, preparation_context}, lineage);
                const auto restart = measured("checkpoint_planning", [&] { return first.plan({ladder.front()}, prep_request); });
                measured("load", [&] { first.load(restart, true); }); // intermediate step-down reporting only
                // A step-down state is never a qualified result.
                scf.last_calculation_outcome = CalculationOutcome{false, "fresh_evaluation_required"};
                scf.last_node.reset();
                preparation_updates = restart.plan.updates_at_protocol;
                const auto settings = scf_scientific_settings(scf);
                // level_finished only chooses the next, smaller bundle; a step-down state is never reused.
                if (step_down->level_finished && step_down->settings == settings && scf_settings_allow_reuse(settings)) {
                    if (world.rank() == 0 && scf.param.print_level() > 0)
                        print("Completed step-down bundle; continue with a smaller orbital bundle.");
                    count = std::max(target_alpha, count - step);
                    preparation_updates = 0;
                } else if (world.rank() == 0 && scf.param.print_level() > 0) {
                    print("Continue the step-down on the loaded orbital bundle.");
                }
                timing.end();
            }
            for (; count > target_alpha; count = std::max(target_alpha, count - step)) {
                timing.begin("orbital_preparation", {{"orbitals_alpha", count},
                    {"protocol", cp::records::protocol_record(ladder.front())}});
                timing.begin("prepare");
                resize_bundle(world, scf, count, target_beta + count - target_alpha);
                scf.reset_aobasis("sto-3g"); scf.ao = scf.project_ao_basis(world, scf.aobasis);
                SCFCheckpoint binding(world, scf, folder, {task, preparation_context}, lineage);
                auto control = binding.control(prep_request, ladder.front(), false, preparation_updates);
                control.preparing_virtuals = true;
                control.timing_observer = observe_iteration;
                timing.end();
                const auto decision = measured("task", [&] {
                    return frozen ? scf.solve_virtuals(world, control) : scf.solve(world, control);
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
            SCFCheckpoint binding(world, scf, folder, {task, context}, lineage);
            auto control = binding.control(request, ladder[stage], stage + 1 == ladder.size(),
                                                 stage == start ? previous_updates : 0);
            control.timing_observer = observe_iteration;
            timing.end();
            const auto decision = measured("task", [&] {
                return frozen ? scf.solve_virtuals(world, control) : scf.solve(world, control);
            }, {{"task_id", task}});
            timing.end();
            if (decision.action != it::Action::finish_protocol) break;
        }
        cp::records::require(scf.last_node.has_value(), "SCF ended without a saved endpoint");
        measured("publish_results", [&] { report(world, root, folder, scf, occupied_saved); });
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
