#include <madness/chem/SCFScientificSettings.h>
#include <madness/chem/CheckpointRecords.h>
#include <madness/chem/SCF.h>
#include <sstream>
#include <cctype>

namespace madness {
namespace {
using json = nlohmann::json;

json tensor_settings(const Tensor<double>& tensor) {
    json dimensions = json::array(), values = json::array();
    for (int axis = 0; axis < tensor.ndim(); ++axis) dimensions.push_back(tensor.dim(axis));
    if (tensor.size() > 0) {
        const auto flat = copy(tensor).flat();
        for (long i = 0; i < flat.size(); ++i) values.push_back(flat[i]);
    }
    return {{"dimensions", dimensions}, {"values", values}};
}

json equation_parameters(const CalculationParameters& parameters) {
    auto result = parameters.to_json().at("parameters");
    // Only known non-scientific controls are removed. Adding a new parameter
    // therefore invalidates old evidence until its role has been considered.
    for (const auto* key : {
             "prefix", "hfexalg", "hfex_granularity", "hfex_cost_aware", "hfex_accumulation",
             "bsh_apply", "memory", "econv", "dconv", "convergence_criteria", "conv_only_dens",
             "convergence.use_default_criteria", "accuracy.bsh_residual", "accuracy.kinetic_trace_change",
             "convergence.criteria", "convergence.bsh_residual", "convergence.bsh_residual.relative",
             "convergence.density_change", "convergence.density_change.relative", "convergence.energy_change",
             "convergence.energy_change.relative", "convergence.kinetic_change", "convergence.kinetic_change.relative",
             "k", "protocol", "maxiter", "maxiter_per_rung", "restart.start_rung", "checkpoint.every_iterations",
             "diagnostic.iterations_per_rung", "progress.window", "progress.minimum_relative_improvement",
             "maxsub", "maxrotn", "nio", "save", "restart", "restartao", "no_compute",
             "nvalpha", "nvbeta", "nmo_alpha", "nmo_beta", "freeze_occupied", "nv_extra", "nv_step", "nv_its", "nv_factor",
             "plotlo", "plothi", "plotdens", "plotcoul", "plotcube", "npt_plot", "plot_cell",
             "aobasis", "nwfile", "molecular_structure", "derivatives", "dipole", "print_level",
             "population", "population_basis", "population_minbasis", "vnucextra", "loadbalparts",
             "gopt", "gtol", "gtest", "gval", "gprec", "gmaxiter", "ginitial_hessian", "algopt"})
        result.erase(key);
    return result;
}
} // namespace

nlohmann::json scf_scientific_settings(const SCF& scf) {
    json unavailable = json::array();
    json pcm = {{"enabled", scf.param.pcm_data() != "none"}};
    if (scf.param.pcm_data() != "none") {
        pcm["parameters"] = scf.pcm_param.to_json().at("parameters");
        pcm["reader_own"] = scf.pcm_param.reader_own();
        pcm["prepared_operator"] = scf.pcm.model_settings();
        if (scf.pcm.model_settings() != PCM::model_settings(scf.molecule, scf.pcm_param))
            unavailable.push_back("unprepared_pcm_operator");
        if (scf.pcm_param.reader_own() || scf.param.pcm_data() == "pcmsolver_reader_own")
            unavailable.push_back("external_pcm_reader_settings");
        else if (scf.pcm_param.cavity_type() == "restart")
            unavailable.push_back("external_pcm_cavity_data");
        else
            pcm["parameters"].erase("restart_name");
    }
    json gth = json::array();
    std::set<unsigned int> recorded;
    for (const auto& atom : scf.molecule.get_atoms()) {
        if (!atom.pseudo_atom || !recorded.insert(atom.atomic_number).second) continue;
        const auto number = atom.atomic_number;
        checkpoint::records::require(number > 0 && number <= 118, "invalid GTH atomic number");
        if (!scf.gthpseudopotential || scf.gthpseudopotential->localp[number - 1].size() != 6) {
            unavailable.push_back("unloaded_gth_data_" + std::to_string(number));
            continue;
        }
        const auto& operator_data = *scf.gthpseudopotential;
        gth.push_back({{"atomic_number", number}, {"local", tensor_settings(operator_data.localp[number - 1])},
                       {"radii", tensor_settings(operator_data.radii[number - 1])},
                       {"h", tensor_settings(operator_data.hlij[number - 1])},
                       {"k", tensor_settings(operator_data.klij[number - 1])}});
    }
    json boundaries = json::array();
    for (std::size_t axis = 0; axis < 3; ++axis)
        boundaries.push_back({static_cast<int>(FunctionDefaults<3>::get_bc()(axis, 0)),
                              static_cast<int>(FunctionDefaults<3>::get_bc()(axis, 1))});
    json result{{"version", 1}, {"dft", equation_parameters(scf.param)},
        {"molecule", scf.molecule.electronic_settings()},
        {"xc", {{"specification", scf.xc.specification()}, {"polarized", scf.xc.is_spin_polarized()},
                 {"hf_exchange", scf.xc.hf_exchange_coefficient()}, {"rhomin", scf.xc.get_rhomin()},
                 {"rhotol", scf.xc.get_rhotol()}, {"tautol", scf.xc.get_tautol()},
                 {"tauwmargin", scf.xc.get_tauwmargin()}}},
        {"dispersion", {{"model", scf.dispersion.description()}, {"atm", scf.dispersion.includes_three_body()},
                         {"library_version", DispersionCorrection::library_version()}}},
        {"pcm", pcm}, {"gth", gth},
        {"representation_settings", {{"cell", tensor_settings(FunctionDefaults<3>::get_cell())},
                                     {"boundary_conditions", boundaries},
                                     {"truncate_mode", FunctionDefaults<3>::get_truncate_mode()}}},
        {"unavailable", unavailable}};
#ifdef MADNESS_HAS_LIBXC
    result["xc"]["backend"] = "libxc";
#else
    result["xc"]["backend"] = "ldaonly";
#endif
    checkpoint::records::finite_numbers(result);
    return result;
}

nlohmann::json scf_model_settings(const nlohmann::json& settings) {
    checkpoint::records::require(settings.at("version") == 1,
                                 "unsupported SCF scientific settings snapshot");
    auto model = settings;
    auto& dft = model.at("dft");
    for (const auto* key : {"lo", "deriv", "dft_deriv", "orbitalshift"}) dft.erase(key);
    auto& molecule = model.at("molecule");
    molecule.at("parameters").erase("eprec");
    molecule.erase("smoothing_radii_inverse");
    for (auto& core : molecule.at("core_potentials"))
        for (const auto* key : {"eprec", "rcut0", "rcut"}) core.at("potential").erase(key);
    auto& xc = model.at("xc");
    for (const auto* key : {"rhomin", "rhotol", "tautol", "tauwmargin"}) xc.erase(key);
    // XC input embeds both the functional and numerical screening controls.
    // Preserve every unknown token and weight; only strip supported screens.
    const auto functional = [](const std::string& specification) {
        std::istringstream in(specification);
        std::string token, result;
        while (in >> token) {
            for (auto& c : token) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            if (token == "RHOMIN" || token == "RHOTOL" || token == "TAUTOL") {
                std::string value;
                checkpoint::records::require(bool(in >> value), "missing XC screening value");
                continue;
            }
            if (!result.empty()) result += ' ';
            result += token;
        }
        return result;
    };
    dft.at("xc") = functional(dft.at("xc").get<std::string>());
    xc.at("specification") = functional(xc.at("specification").get<std::string>());
    return model;
}

bool scf_settings_allow_reuse(const nlohmann::json& settings) {
    checkpoint::records::require(settings.at("version") == 1 && settings.at("unavailable").is_array(),
                                 "unsupported SCF scientific settings snapshot");
    return settings.at("unavailable").empty();
}
} // namespace madness
