#include <madness/chem/SCFCheckpoint.h>
#include <madness/chem/IterationHistory.h>
#include <madness/chem/IterationOutput.h>
#include <madness/chem/SCFScientificSettings.h>
#include <madness/chem/StageCheckpoint.h>
#include <madness/misc/info.h>

namespace madness {
namespace cp = checkpoint;
namespace it = iteration;
using json = nlohmann::json;
namespace {
it::MetricLabels metric_labels() {
    return {{it::Metric::bsh_residual, "BSH residual, maximum"},
            {it::Metric::density_change, "Maximum spin-density change"},
            {it::Metric::energy_change, "Energy change (Ha)"},
            {it::Metric::kinetic_change, "Kinetic change (Ha)"}};
}
bool same_occupation_request(const json& stored, const json& current) {
    if (!stored.contains("dft")) return false;
    const auto keys = {"aocc", "bocc", "smear"};
    return std::all_of(keys.begin(), keys.end(), [&](const char* key) {
        return stored.at("dft").contains(key) && stored.at("dft").at(key) == current.at("dft").at(key);
    });
}
// The archive's facts this binding checks before loading it.
const json& orbitals(const it::NodeDescription& d) { return d.details.at("orbitals"); }
} // namespace

std::filesystem::path ground_state_node(const std::filesystem::path& run) { return run / "states" / "ground_state"; }
std::filesystem::path virtuals_node(const std::filesystem::path& run) { return run / "states" / "virtuals"; }

std::optional<it::NodeDescription> read_node(World& world, const std::filesystem::path& folder) {
    // No saved state only when neither a description nor an archive is there
    // (a folder with just the history is a launch that has not saved yet).
    // Anything half-present is an error, never a silent fresh start.
    const auto text = cp::detail::on_root(world, [&] {
        namespace fs = std::filesystem;
        const it::Node node(folder);
        const auto saved = node.inspect(); // finishes or discards an interrupted save first
        const bool described = fs::exists(fs::symlink_status(folder / "description.json"));
        const bool archived = fs::is_directory(node.archive());
        if (!saved && described) throw std::runtime_error("unreadable node description: " + (folder / "description.json").string());
        if (!saved && archived) throw std::runtime_error("node archive without a description: " + folder.string());
        if (saved && !archived) throw std::runtime_error("node description without an archive: " + folder.string());
        return saved ? it::description_record(*saved).dump() : std::string{};
    });
    if (text.empty()) return std::nullopt;
    return it::read_description(json::parse(text));
}

json node_identity(const it::NodeDescription& d) {
    return {{"settings", d.settings}, {"protocol", cp::records::protocol_record(d.protocol)},
            {"saved_at", d.saved_at}};
}

SCFCheckpoint::SCFCheckpoint(World& world, SCF& scf, std::filesystem::path folder, Identity identity, json lineage)
    : world_(world), scf_(scf), folder_(std::move(folder)), identity_(std::move(identity)),
      settings_(capture_settings()), lineage_(std::move(lineage)), layout_(layout()) {
    cp::records::identifier(identity_.task_id);
    cp::records::identifier(identity_.attempt_id);
    cp::records::require(!identity_.context_id.empty() && settings_.is_object() && !settings_.empty(),
                         "SCF needs a complete scientific settings snapshot and context identity");
    cp::records::finite_numbers(settings_);
    cp::records::require(scf_.restart_representation == Representation::mo,
                         "SCF checkpoint binding requires physical orbitals");
    cp::records::require(!scf_.checkpoint_observer, "cannot combine shared and legacy checkpoint observers");
    cp::records::require(!scf_.param.freeze_occupied() || !lineage_.is_null(),
                         "frozen SCF needs the identity of its qualified ground state");
}

std::string SCFCheckpoint::kind() const { return scf_.param.freeze_occupied() ? "virtuals" : "ground_state"; }

json SCFCheckpoint::capture_settings() const {
    json result;
    std::string error, agreement;
    try {
        result = scf_scientific_settings(scf_);
        cp::records::require(scf_.param.pcm_data() == "none" ||
            scf_.pcm.model_settings() == PCM::model_settings(scf_.molecule, scf_.pcm_param),
            "prepare the PCM operator for the current settings and geometry before binding SCF");
        agreement = json{{"scientific_settings", result}, {"orbitals", layout()}}.dump();
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
    return {identity_.task_id, identity_.context_id, {scf_.param.freeze_occupied() ? "virtuals" : "ground"}};
}

json SCFCheckpoint::layout() const {
    return {{"representation", to_string(scf_.restart_representation)}, {"L", scf_.param.L()},
            {"spin_restricted", scf_.param.spin_restricted()}, {"nalpha", scf_.param.nalpha()},
            {"nbeta", scf_.param.nbeta()}, {"nmo_alpha", scf_.param.nmo_alpha()},
            {"nmo_beta", scf_.param.have_beta() ? scf_.param.nmo_beta() : 0},
            {"freeze_occupied", scf_.param.freeze_occupied()}};
}

bool SCFCheckpoint::loadable(const json& stored) const {
    // The existing native reader supports k/threshold conversion and virtual
    // padding/trimming. It does not convert spin layouts, cells or representations.
    for (const auto* key : {"representation", "L", "spin_restricted", "nalpha", "nbeta", "freeze_occupied"})
        if (stored.at(key) != layout_.at(key)) return false;
    // Saved archives already contain the active orbitals. The legacy reader
    // instead removes a core prefix whenever it trims an oversized bundle. Do
    // not let it drop active occupied functions (or overrun a slice).
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
    auto saved = read_node(world_, folder_);
    cp::records::require(!require_existing || saved.has_value(), "no saved SCF state");
    std::optional<it::SavedTask> task;
    if (saved) {
        const auto& d = *saved;
        cp::records::require(d.kind == kind(), "the SCF node holds another kind of state");
        cp::records::require(d.details.at("scf_node_version") == 1, "unsupported SCF node version");
        cp::records::require(scf_model_settings(d.settings) == scf_model_settings(settings_),
                             "SCF model differs from the saved state; use a new output prefix or an explicit seed");
        cp::records::require(loadable(orbitals(d)), "unsupported SCF archive conversion");
        cp::records::require(!scf_.param.freeze_occupied() || same_occupation_request(d.settings, settings_),
            "frozen occupation request changed; rebuild from the qualified occupied source");
        task = it::SavedTask{};
        task->metrics = it::stored_metrics(d, {identity_.task_id, identity_.attempt_id, d.protocol, 0,
                                               identity_.context_id});
        // level_finished never decides reuse. plan_restart reads it only to
        // advance past a finished intermediate level; reuse at the final level
        // is decided by saved_completion on the stored metrics alone.
        task->protocol_complete = d.level_finished;
        task->updates_at_protocol = d.updates_at_protocol;
        task->conversion_supported = true;
        // A different bundle (a step-down state with extra orbitals among them)
        // or a changed ground state for the virtuals: the functions are a guess.
        if (!scf_settings_allow_reuse(settings_) || d.settings != settings_ || orbitals(d) != layout_ ||
            d.lineage != lineage_)
            task->use = it::SavedEvidenceUse::guess_only;
    }
    const auto plan = it::plan_restart(protocols, request, scf_.iteration_model(), expected(), task);
    // Resolve display limits after the selected protocol has been prepared;
    // planning may still have the first protocol's FunctionDefaults in scope.
    const auto display_request = task && task->use == it::SavedEvidenceUse::compatible
        ? std::optional<it::Request>(request) : std::nullopt;
    return {plan, protocols, std::move(saved), display_request, this};
}

void SCFCheckpoint::load(const Restart& restart, bool report_preparation) {
    it::ElapsedTimer<> load_timer;
    cp::records::require(restart.binding_ == this, "restart plan belongs to another SCF checkpoint binding");
    cp::records::require(restart.saved.has_value() && restart.plan.action != it::RestartAction::fresh,
                         "no saved SCF state selected for loading");
    check_current();
    const auto target = restart.protocols.at(restart.plan.start_protocol);
    cp::records::require(target == it::Protocol{FunctionDefaults<3>::get_thresh(), FunctionDefaults<3>::get_k()},
                         "prepare the chosen SCF protocol before loading");
    // Read again: planning does not grant a stale description.
    const auto current = read_node(world_, folder_);
    cp::records::require(current == restart.saved, "the saved SCF state changed after planning");
    const auto& d = *current;
    cp::records::require(d.details.at("writers").get<int>() <= world_.size(),
                         "native SCF archive needs at least as many ranks as its I/O writers");
    cp::records::require(loadable(orbitals(d)), "unsupported SCF archive conversion");
    const bool reset_occupations = !same_occupation_request(d.settings, settings_);
    cp::records::require(!scf_.param.freeze_occupied() || !reset_occupations,
                         "frozen occupation request changed; rebuild from the qualified occupied source");
    const bool reuse = restart.plan.action == it::RestartAction::reuse;
    // A state with a different bundle (extra step-down orbitals) never qualifies.
    if (reuse)
        cp::records::require(d.settings == settings_ && orbitals(d) == layout_ && d.lineage == lineage_ &&
                             d.protocol == target, "saved completion does not describe the requested SCF bundle");
    // The legacy reader takes its basename from prefix. Restore the entire
    // parameter object so this temporary path never becomes user input/output.
    struct RestoreParameters {
        CalculationParameters& current;
        CalculationParameters original;
        ~RestoreParameters() { current = original; }
    } restore{scf_.param, scf_.param};
    scf_.param.set_user_defined_value("prefix", (it::Node(folder_).archive() / "scf").string());
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
    const double energy = d.details.at("energy").get<double>();
    if (reuse) {
        const auto exact = [&](const vecfuncT& functions, int count) {
            return functions.size() == static_cast<std::size_t>(count) &&
                std::all_of(functions.begin(), functions.end(), [&](const auto& function) {
                    return function.k() == target.k && function.thresh() == target.thresh;
                });
        };
        cp::records::require(exact(scf_.amo, scf_.param.nmo_alpha()) &&
                             (!scf_.param.have_beta() || exact(scf_.bmo, scf_.param.nmo_beta())) &&
                             scf_.current_energy == energy,
                             "native SCF payload changed while loading reusable evidence");
        scf_.last_calculation_outcome = CalculationOutcome{true, "saved_criteria_met"};
        scf_.last_node = folder_;
    } else {
        set_thresh(world_, scf_.amo, target.thresh);
        if (scf_.param.have_beta()) set_thresh(world_, scf_.bmo, target.thresh);
        scf_.last_calculation_outcome = CalculationOutcome{false, "fresh_evaluation_required"};
        // Frozen energy belongs to its qualified occupied source, even when k
        // conversion invalidates the combined occupied+virtual legacy stamp.
        if (scf_.param.freeze_occupied()) scf_.current_energy = energy;
    }
    scf_.converged_for_thresh = 1e10;
    scf_.converged_for_dconv = 1e10;
    if (world_.rank() == 0 && it::prints(it::verbosity(scf_.param.print_level()), it::needed::restart_summary)) {
        const double load_seconds = load_timer.seconds();
        const auto metrics = it::stored_metrics(d, {identity_.task_id, identity_.attempt_id, d.protocol, 0,
                                                    identity_.context_id});
        it::Assessment assessment;
        if (restart.display_request_)
            assessment = it::assess(it::resolve(*restart.display_request_, scf_.iteration_model(), target,
                target == restart.protocols.back(), "ground state"), expected(), metrics);
        it::write_restart_summary(std::cout, report_preparation ? "SCF preparation" : "SCF", metrics,
            target, assessment, reuse, metric_labels(),
            report_preparation ? "Loaded saved step-down functions." : "");
        if (it::prints(it::verbosity(scf_.param.print_level()), it::needed::checkpoint_timing))
            it::write_checkpoint_timing(std::cout, "SCF", identity_.task_id, load_seconds);
    }
}

void SCFCheckpoint::load_guess(const std::filesystem::path& source) {
    check_current();
    const auto saved = read_node(world_, source);
    cp::records::require(saved.has_value(), "no saved SCF state to load as a guess");
    auto stored = orbitals(*saved);
    stored["freeze_occupied"] = layout_.at("freeze_occupied");
    cp::records::require(loadable(stored), "unsupported SCF guess conversion");
    cp::records::require(saved->details.at("writers").get<int>() <= world_.size(), "unsupported SCF guess archive");
    if (scf_.param.freeze_occupied())
        cp::records::require(saved->kind == "ground_state" && node_identity(*saved) == lineage_,
                             "frozen guess must be the qualified occupied source");
    const auto parameters = scf_.param;
    try {
        scf_.param.set_user_defined_value("prefix", (it::Node(source).archive() / "scf").string());
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
    } else scf_.current_energy = saved->details.at("energy").get<double>();
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
        if (world.rank() == 0 && it::prints(it::verbosity(scf_.param.print_level()), it::needed::iteration_metrics))
            it::write_iteration_metrics(std::cout, "SCF", metrics, decision, metric_labels(), final_protocol);
        // One history row per decision, before any save can fail.
        cp::detail::on_root(world, [&] {
            it::append_history(it::Node(folder_).history(),
                               it::history_row(metrics, decision, SCF::iteration_metric_set()));
            return std::string{};
        });
        if (decision.checkpoint) {
            save(metrics, decision, progress);
            if (world.rank() == 0 && it::prints(it::verbosity(scf_.param.print_level()), it::needed::checkpoint_saved))
                it::write_checkpoint_saved(std::cout, decision);
        }
    };
    return result;
}

void SCFCheckpoint::save(const it::IterationMetrics& metrics, const it::Decision& decision,
                         const it::Progress& progress) {
    check_current();
    cp::records::require(progress.state_valid && decision.info == metrics.info &&
                         metrics.info.iteration == progress.main_updates, "sample is not eligible for saving");
    cp::records::require(scf_.amo.size() == static_cast<std::size_t>(scf_.param.nmo_alpha()) &&
                         (!scf_.param.have_beta() || scf_.bmo.size() == static_cast<std::size_t>(scf_.param.nmo_beta())),
                         "SCF payload does not match its declared bundle");
    // The native archive is written by every rank, so it is staged beside the
    // node first; the one save then moves it into place.
    const auto staging = folder_ / "archive.staging";
    cp::detail::on_root(world_, [&] {
        std::filesystem::remove_all(staging);
        std::filesystem::create_directories(staging);
        return std::string{};
    });
    cp::detail::collective_error(world_, scf_.write_mos(world_, (staging / "scf.restartdata").string(),
                                                        scf_.orbital_metadata()));
    world_.gop.fence();
    cp::detail::on_root(world_, [&] {
        // ParallelOutputArchive clamps the requested writer count to ranks
        // and its native maximum. Read the actual count from shard zero.
        archive::BinaryFstreamInputArchive header((staging / "scf.restartdata.00000").string());
        int writers = 0;
        header & writers;
        it::NodeDescription d;
        d.kind = kind();
        d.protocol = metrics.info.protocol;
        d.updates_at_protocol = progress.updates_at_protocol;
        d.settings = settings_;
        d.metrics = it::node_metrics(metrics, SCF::iteration_metric_set());
        d.lineage = lineage_;
        d.code_version = info::git_commit();
        d.level_finished = decision.action == it::Action::finish_protocol || decision.action == it::Action::complete;
        d.details = {{"scf_node_version", 1}, {"orbitals", layout_}, {"energy", scf_.current_energy},
                     {"writers", writers}};
        it::Node(folder_).save(d, [&](const std::filesystem::path& archive) {
            for (const auto& file : std::filesystem::directory_iterator(staging))
                std::filesystem::rename(file.path(), archive / file.path().filename());
        });
        std::filesystem::remove_all(staging);
        return std::string{};
    });
    scf_.last_node = folder_;
}

} // namespace madness
