#include <madness/chem/SCFCheckpoint.h>
#include <madness/chem/IterationOutput.h>
#include <madness/chem/SCFScientificSettings.h>

namespace madness {
namespace cp = checkpoint;
namespace it = iteration;
namespace {
it::MetricLabels metric_labels() {
    return {{it::Metric::bsh_residual, "BSH residual, maximum"},
            {it::Metric::density_change, "Maximum spin-density change"},
            {it::Metric::energy_change, "Energy change (Ha)"},
            {it::Metric::kinetic_change, "Kinetic change (Ha)"}};
}
bool same_occupation_request(const cp::json& stored, const cp::json& current) {
    if (!stored.contains("dft")) return false;
    const auto keys = {"aocc", "bocc", "smear"};
    return std::all_of(keys.begin(), keys.end(), [&](const char* key) {
        return stored.at("dft").contains(key) && stored.at("dft").at(key) == current.at("dft").at(key);
    });
}
} // namespace

SCFCheckpoint::SCFCheckpoint(World& world, SCF& scf, cp::CheckpointStore& store,
                             Identity identity, cp::json inputs)
    : world_(world), scf_(scf), store_(store), identity_(std::move(identity)),
      settings_(capture_settings()), inputs_(std::move(inputs)), layout_(layout()) {
    cp::records::identifier(identity_.calculation_id);
    cp::records::identifier(identity_.task_id);
    cp::records::identifier(identity_.attempt_id);
    cp::records::require(!identity_.context_id.empty() && settings_.is_object() && !settings_.empty(),
                         "SCF needs a complete scientific settings snapshot and context identity");
    cp::records::finite_numbers(settings_);
    cp::records::require(scf_.restart_representation == Representation::mo,
                         "SCF checkpoint binding requires physical orbitals");
    cp::records::require(!scf_.checkpoint_observer, "cannot combine shared and legacy checkpoint observers");
    const auto sources = cp::records::inputs(inputs_);
    if (scf_.param.freeze_occupied()) {
        cp::records::require(sources.size() == 1, "frozen SCF needs one exact occupied source");
        const auto source = store_.read(sources.front());
        // Current qualification belongs to the workflow. An immutable source
        // may satisfy today's request even when its original attempt stopped.
        cp::records::require(source.record.at("stage").at("application") == "scf" &&
                            source.record.at("task").at("kind") == "iterative" &&
                            source.record.at("task").at("definition").value("subject", "") == "ground",
                            "frozen SCF source must be an occupied SCF task");
    }
    store_.register_task(identity_.calculation_id, "scf", identity_.task_id, "iterative", true, inputs_);
}

cp::json SCFCheckpoint::capture_settings() const {
    cp::json result;
    std::string error, agreement;
    try {
        result = scf_scientific_settings(scf_);
        cp::records::require(scf_.param.pcm_data() == "none" ||
            scf_.pcm.model_settings() == PCM::model_settings(scf_.molecule, scf_.pcm_param),
            "prepare the PCM operator for the current settings and geometry before binding SCF");
        agreement = cp::json{{"scientific_settings", result}, {"orbitals", layout()}}.dump();
    } catch (const std::exception& exception) {
        error = exception.what();
    }
    cp::detail::collective_error(world_, error);
    auto root_settings = agreement;
    world_.gop.broadcast_serializable(root_settings, 0);
    cp::detail::collective_error(world_, root_settings == agreement ? std::string{} :
        std::string("SCF scientific settings or orbital layout differ between ranks"));
    return result;
}

void SCFCheckpoint::check_current() const {
    const auto current = capture_settings();
    cp::detail::collective_error(world_, current == settings_ && layout() == layout_ ? std::string{} :
        std::string("SCF settings or bundle changed; construct a new checkpoint binding"));
}

it::ExpectedTask SCFCheckpoint::expected() const {
    return {identity_.task_id, identity_.context_id,
            {scf_.param.freeze_occupied() ? "virtuals" : "ground"}};
}

cp::json SCFCheckpoint::layout() const {
    return {{"representation", to_string(scf_.restart_representation)}, {"L", scf_.param.L()},
            {"spin_restricted", scf_.param.spin_restricted()}, {"nalpha", scf_.param.nalpha()},
            {"nbeta", scf_.param.nbeta()}, {"nmo_alpha", scf_.param.nmo_alpha()},
            {"nmo_beta", scf_.param.have_beta() ? scf_.param.nmo_beta() : 0},
            {"freeze_occupied", scf_.param.freeze_occupied()}};
}

bool SCFCheckpoint::loadable(const cp::json& stored) const {
    // The existing native reader supports k/threshold conversion and virtual
    // padding/trimming. It does not convert spin layouts, cells or representations.
    for (const auto* key : {"representation", "L", "spin_restricted", "nalpha", "nbeta", "freeze_occupied"})
        if (stored.at(key) != layout_.at(key)) return false;
    // Shared SCF archives already contain the active orbitals. The legacy
    // reader instead removes a core prefix whenever it trims an oversized
    // bundle. Do not let it drop active occupied functions (or overrun a slice).
    if (scf_.molecule.n_core_orb_all() != 0 &&
        (stored.at("nmo_alpha").get<int>() > scf_.param.nmo_alpha() ||
         (scf_.param.have_beta() && stored.at("nmo_beta").get<int>() > scf_.param.nmo_beta())))
        return false;
    return stored.at("nmo_alpha").get<int>() >= scf_.param.nalpha() &&
        (!scf_.param.have_beta() || stored.at("nmo_beta").get<int>() >= scf_.param.nbeta());
}

SCFCheckpoint::Restart SCFCheckpoint::plan(const std::vector<it::Protocol>& protocols,
                                           const it::Request& request, bool require_existing) {
    check_current();
    auto checkpoint = store_.read_latest(identity_.calculation_id, identity_.task_id);
    cp::records::require(!require_existing || checkpoint.has_value(), "no selected SCF checkpoint");
    std::optional<it::SavedTask> saved;
    if (checkpoint) {
        const auto& record = checkpoint->record;
        const auto& task = record.at("task");
        const auto& definition = task.at("definition");
        cp::records::require(definition.at("scf_checkpoint_version") == 1,
                             "unsupported SCF checkpoint binding version");
        cp::records::require(scf_model_settings(definition.at("scientific_settings")) == scf_model_settings(settings_),
                             "SCF model differs from the saved checkpoint; use a new output prefix or an explicit seed");
        cp::records::require(loadable(definition.at("orbitals")), "unsupported SCF archive conversion");
        cp::records::require(!scf_.param.freeze_occupied() || task.at("inputs") == inputs_,
                             "frozen source changed; rebuild the bundle from the qualified occupied source");
        cp::records::require(!scf_.param.freeze_occupied() ||
            same_occupation_request(definition.at("scientific_settings"), settings_),
            "frozen occupation request changed; rebuild from the qualified occupied source");
        cp::records::require(record.at("metric_definitions") == SCF::iteration_metric_definitions(),
                             "unsupported SCF metric definitions");
        saved = it::SavedTask{};
        saved->metrics = cp::iteration_metrics_from_manifest(*checkpoint);
        saved->protocol_complete = task.at("iteration").at("protocol_complete");
        saved->updates_at_protocol = cp::records::count(task.at("iteration").at("updates_at_protocol"));
        saved->conversion_supported = true;
        if (!scf_settings_allow_reuse(settings_) ||
            definition.at("scientific_settings") != settings_ || definition.at("orbitals") != layout_ ||
            task.at("inputs") != inputs_ || saved->metrics.info.context_id != identity_.context_id)
            saved->use = it::SavedEvidenceUse::guess_only;
    }
    const auto plan = it::plan_restart(protocols, request, scf_.iteration_model(), expected(), saved);
    // Resolve display limits after the selected protocol has been prepared;
    // planning may still have the first protocol's FunctionDefaults in scope.
    const auto display_request = saved && saved->use == it::SavedEvidenceUse::compatible
        ? std::optional<it::Request>(request) : std::nullopt;
    return {plan, protocols, std::move(checkpoint), display_request, this};
}

void SCFCheckpoint::load(const Restart& restart, bool report_preparation) {
    it::ElapsedTimer<> load_timer;
    cp::records::require(restart.binding_ == this, "restart plan belongs to another SCF checkpoint binding");
    cp::records::require(restart.checkpoint.has_value() && restart.plan.action != it::RestartAction::fresh,
                         "no SCF checkpoint selected for loading");
    check_current();
    const auto target = restart.protocols.at(restart.plan.start_protocol);
    cp::records::require(target == it::Protocol{FunctionDefaults<3>::get_thresh(), FunctionDefaults<3>::get_k()},
                         "prepare the chosen SCF protocol before loading");
    const auto& original = *restart.checkpoint;
    cp::records::require(original.record.at("calculation_id") == identity_.calculation_id &&
                        original.record.at("task").at("task_id") == identity_.task_id,
                         "checkpoint belongs to another SCF task");
    // Validate exact files again: planning does not grant a stale integrity check.
    const auto manifest = store_.read({identity_.calculation_id, identity_.task_id, original.generation});
    cp::records::require(manifest.archives.size() == 1 && manifest.archives.front().role == "orbitals" &&
                        manifest.archives.front().basename == "scf.restartdata",
                         "unsupported SCF native payload");
    cp::records::require(manifest.archives.front().writers <= world_.size(),
                         "native SCF archive needs at least as many ranks as its I/O writers");
    const auto& definition = manifest.record.at("task").at("definition");
    cp::records::require(loadable(definition.at("orbitals")), "unsupported SCF archive conversion");
    const bool reset_occupations = !same_occupation_request(definition.at("scientific_settings"), settings_);
    cp::records::require(!scf_.param.freeze_occupied() || !reset_occupations,
                         "frozen occupation request changed; rebuild from the qualified occupied source");
    const bool reuse = restart.plan.action == it::RestartAction::reuse;
    if (reuse)
        cp::records::require(definition.at("scientific_settings") == settings_ &&
                             definition.at("orbitals") == layout_ &&
                             manifest.record.at("task").at("inputs") == inputs_ &&
                             cp::records::read_protocol(definition.at("protocol")) == target,
                             "saved completion does not describe the requested SCF bundle");
    // The legacy reader takes its basename from prefix. Restore the entire
    // parameter object so this temporary path never becomes user input/output.
    struct RestoreParameters {
        CalculationParameters& current;
        CalculationParameters original;
        ~RestoreParameters() { current = original; }
    } restore{scf_.param, scf_.param};
    scf_.param.set_user_defined_value("prefix", (store_.task_directory(identity_.calculation_id,
        identity_.task_id) / manifest.generation / "scf").string());
    scf_.load_mos(world_);
    if (reset_occupations) {
        // Archived holes are properties of the source state. Start the target
        // request from the same complete default occupations as a fresh guess,
        // then apply its explicit prefix; no source suffix may survive.
        const auto defaults = [](int count, int occupied) {
            tensorT result(count);
            for (int i = 0; i < occupied; ++i) result[i] = 1.0;
            return result;
        };
        scf_.aocc = defaults(scf_.param.nmo_alpha(), scf_.param.nalpha());
        if (scf_.param.have_beta()) scf_.bocc = defaults(scf_.param.nmo_beta(), scf_.param.nbeta());
        scf_.apply_explicit_occupations(world_);
        scf_.aset = scf_.group_orbital_sets(world_, scf_.aeps, scf_.aocc, scf_.param.nmo_alpha());
        if (scf_.param.have_beta())
            scf_.bset = scf_.group_orbital_sets(world_, scf_.beps, scf_.bocc, scf_.param.nmo_beta());
    }
    if (reuse) {
        const auto exact = [&](const vecfuncT& functions, int count) {
            return functions.size() == static_cast<std::size_t>(count) &&
                std::all_of(functions.begin(), functions.end(), [&](const auto& function) {
                    return function.k() == target.k && function.thresh() == target.thresh;
                });
        };
        cp::records::require(exact(scf_.amo, scf_.param.nmo_alpha()) &&
                             (!scf_.param.have_beta() || exact(scf_.bmo, scf_.param.nmo_beta())) &&
                             scf_.current_energy == definition.at("energy").get<double>(),
                             "native SCF payload changed while loading reusable evidence");
        scf_.last_calculation_outcome = CalculationOutcome{true, "saved_criteria_met"};
        scf_.last_checkpoint = cp::GenerationRef{identity_.calculation_id, identity_.task_id, manifest.generation};
    } else {
        set_thresh(world_, scf_.amo, target.thresh);
        if (scf_.param.have_beta()) set_thresh(world_, scf_.bmo, target.thresh);
        scf_.last_calculation_outcome = CalculationOutcome{false, "fresh_evaluation_required"};
        // Frozen energy belongs to its qualified occupied source, even when k
        // conversion invalidates the combined occupied+virtual legacy stamp.
        if (scf_.param.freeze_occupied()) scf_.current_energy = definition.at("energy").get<double>();
    }
    scf_.converged_for_thresh = 1e10;
    scf_.converged_for_dconv = 1e10;
    if (world_.rank() == 0 && scf_.param.print_level() > 0) {
        const double load_seconds = load_timer.seconds();
        const auto metrics = cp::iteration_metrics_from_manifest(manifest);
        it::Assessment assessment;
        if (restart.display_request_)
            assessment = it::assess(it::resolve(*restart.display_request_, scf_.iteration_model(), target,
                target == restart.protocols.back(), "ground state"), expected(), metrics);
        it::write_restart_summary(std::cout, report_preparation ? "SCF preparation" : "SCF", metrics,
            target, assessment, reuse, metric_labels(),
            report_preparation ? "Loaded saved preparation functions." : "");
        it::write_checkpoint_timing(std::cout, "SCF", identity_.task_id, load_seconds);
    }

}

void SCFCheckpoint::load_guess(const cp::Manifest& source) {
    check_current();
    const cp::GenerationRef ref{source.record.at("calculation_id"),
        source.record.at("task").at("task_id"), source.generation};
    const auto saved = store_.read(ref);
    const auto& definition = saved.record.at("task").at("definition");
    auto orbitals = definition.at("orbitals");
    orbitals["freeze_occupied"] = layout_.at("freeze_occupied");
    cp::records::require(loadable(orbitals), "unsupported SCF guess conversion");
    cp::records::require(saved.archives.size() == 1 && saved.archives.front().role == "orbitals" &&
        saved.archives.front().basename == "scf.restartdata" && saved.archives.front().writers <= world_.size(),
        "unsupported SCF guess archive");
    if (scf_.param.freeze_occupied()) {
        const auto sources = cp::records::inputs(inputs_);
        cp::records::require(std::find(sources.begin(), sources.end(), ref) != sources.end(),
                             "frozen guess must be the qualified occupied source");
    }
    const auto parameters = scf_.param;
    try {
        scf_.param.set_user_defined_value("prefix", (store_.task_directory(ref.calculation_id, ref.task_id) /
            ref.generation / "scf").string());
        scf_.load_mos(world_);
    } catch (...) { scf_.param = parameters; throw; }
    scf_.param = parameters;
    if (!scf_.param.freeze_occupied()) {
        scf_.aocc = tensorT(scf_.param.nmo_alpha());
        for (int i = 0; i < scf_.param.nalpha(); ++i) scf_.aocc[i] = 1.0;
        if (scf_.param.have_beta()) {
            scf_.bocc = tensorT(scf_.param.nmo_beta());
            for (int i = 0; i < scf_.param.nbeta(); ++i) scf_.bocc[i] = 1.0;
        }
        scf_.apply_explicit_occupations(world_);
        scf_.aset = scf_.group_orbital_sets(world_, scf_.aeps, scf_.aocc, scf_.param.nmo_alpha());
        if (scf_.param.have_beta()) scf_.bset = scf_.group_orbital_sets(world_, scf_.beps, scf_.bocc, scf_.param.nmo_beta());
    } else scf_.current_energy = definition.at("energy").get<double>();
    scf_.converged_for_thresh = scf_.converged_for_dconv = 1e10;
    scf_.last_calculation_outcome = CalculationOutcome{false, "fresh_evaluation_required"};
}

SCF::IterationControl SCFCheckpoint::control(const it::Request& request, it::Protocol protocol,
                                             bool final_protocol, std::size_t previous_updates) {
    check_current();
    SCF::IterationControl result{it::resolve(request, scf_.iteration_model(), protocol, final_protocol, "ground state"),
                                final_protocol, expected(), identity_.attempt_id, previous_updates};
    result.observer = [this, final_protocol](World& world, const SCF& owner, const it::IterationMetrics& metrics,
                             const it::Decision& decision, const it::Progress& progress) {
        cp::records::require(&world == &world_ && &owner == &scf_, "SCF control used with another application");
        if (world.rank() == 0 && scf_.param.print_level() > 0)
            it::write_iteration_metrics(std::cout, "SCF", metrics, decision, metric_labels(), final_protocol);
        if (decision.checkpoint) {
            publish(metrics, decision, progress);
            if (world.rank() == 0 && scf_.param.print_level() > 0)
                it::write_checkpoint_saved(std::cout, decision);
        }
    };
    return result;
}

void SCFCheckpoint::publish(const it::IterationMetrics& metrics, const it::Decision& decision,
                             const it::Progress& progress) {
    check_current();
    cp::records::require(scf_.amo.size() == static_cast<std::size_t>(scf_.param.nmo_alpha()) &&
                         (!scf_.param.have_beta() || scf_.bmo.size() == static_cast<std::size_t>(scf_.param.nmo_beta())),
                         "SCF payload does not match its declared bundle");
    cp::json stage{{"application", "scf"}, {"settings", settings_},
        {"requested_properties", cp::json::array()}, {"inputs", inputs_},
        {"outcome", {{"complete", decision.action == it::Action::complete}, {"reason", decision.reason}}}};
    // A workflow may reconstruct a binding within one attempt. Published
    // identities remain reserved even after their numerical payload is reclaimed.
    const auto generation = cp::detail::on_root(world_, [&] {
        const auto directory = store_.task_directory(identity_.calculation_id, identity_.task_id);
        for (;;) {
            cp::records::require(publications_ < std::numeric_limits<std::size_t>::max(),
                                 "SCF generation counter overflow");
            const auto candidate = identity_.attempt_id + "-" + std::to_string(++publications_);
            const auto exists = [&](const std::string& name) {
                return std::filesystem::exists(std::filesystem::symlink_status(directory / name));
            };
            if (!exists(candidate) && !exists(".published-" + candidate + ".json") &&
                !exists(".deleting-" + candidate)) return candidate;
        }
    });
    auto manifest = cp::iterative_manifest(identity_.calculation_id, generation, stage,
        {{"subject", expected().states.front()}, {"scf_checkpoint_version", 1},
         {"scientific_settings", settings_}, {"orbitals", layout_}, {"energy", scf_.current_energy}},
        inputs_, SCF::iteration_metric_definitions(), SCF::iteration_metric_set(), metrics, decision, progress);
    const auto warning = store_.publish(manifest, [&](const std::filesystem::path& directory) {
        cp::Payload payload;
        payload.local_error = scf_.write_mos(world_, (directory / "scf.restartdata").string(), scf_.orbital_metadata());
        world_.gop.fence();
        if (world_.rank() == 0 && payload.local_error.empty()) {
            // ParallelOutputArchive clamps the requested writer count to ranks
            // and its native maximum. Read the actual count from shard zero.
            archive::BinaryFstreamInputArchive header((directory / "scf.restartdata.00000").string());
            int writers = 0;
            header & writers;
            payload.archives.push_back(cp::inspect_archive(directory, "orbitals", "scf.restartdata", {}, writers));
        }
        return payload;
    });
    scf_.last_checkpoint = cp::GenerationRef{identity_.calculation_id, identity_.task_id, generation};
    if (!warning.empty() && world_.rank() == 0) print(warning);
}

} // namespace madness
