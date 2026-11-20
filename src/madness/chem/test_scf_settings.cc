#include <madness/chem/SCFScientificSettings.h>
#include <madness/chem/SCF.h>
#include <cstdlib>
#include <iostream>

using namespace madness;
namespace {
int failures = 0;
void check(bool ok, const char* message) {
    if (!ok) { ++failures; std::cerr << "FAILED: " << message << '\n'; }
}
struct FutureParameters : CalculationParameters {
    FutureParameters() { initialize<double>("future_operator_threshold", 1e-7); }
};
void exercise(World& world) {
    Molecule molecule;
    molecule.add_atom(0.1, 0, 0, 2, 2);
    CalculationParameters p;
    p.set_user_defined_value("localize", std::string("canon"));
    p.set_user_defined_value("print_level", 0);
    p.set_user_defined_value("protocol", std::vector<double>{1e-4});
    p.set_user_defined_value("k", 6);
    SCF scf(world, p, molecule);
    scf.set_protocol<3>(world, 1e-4);
    const auto base = scf_scientific_settings(scf);
    check(scf_settings_allow_reuse(base), "ordinary all-electron settings are available for reuse");

    const auto model = scf_model_settings(base);
    auto numerical = base;
    numerical["dft"]["lo"] = 1e-8;
    numerical["dft"]["deriv"] = "ble";
    numerical["molecule"]["parameters"]["eprec"] = 1e-8;
    numerical["molecule"]["smoothing_radii_inverse"] = {100.0};
    numerical["xc"]["rhomin"] = 1e-15;
    numerical["xc"]["specification"] = "HF RHOMIN 1e-15 RHOTOL 1e-10";
    numerical["dft"]["xc"] = "HF RHOMIN 1e-15 RHOTOL 1e-10";
    check(scf_model_settings(numerical) == model, "numerical refinement preserves the physical model");
    check(numerical != base, "numerical refinement still invalidates saved metric evidence");
    for (const auto* key : {"localize", "xc", "aocc", "smear", "future_parameter"}) {
        auto different = base;
        different["dft"][key] = "changed";
        check(scf_model_settings(different) != model, "physical and unknown parameter changes remain model differences");
    }
    const auto old_parameters = scf.param;
    scf.param.set_user_defined_value("prefix", std::string("somewhere_else"));
    scf.param.set_user_defined_value("maxiter", 100);
    scf.param.set_user_defined_value("dconv", 1e-9);
    scf.param.set_user_defined_value("econv", 1e-10);
    scf.param.set_user_defined_value("maxsub", 0);
    scf.param.set_user_defined_value("maxrotn", 0.1);
    scf.param.set_user_defined_value("print_level", 10);
    scf.param.set_user_defined_value("nio", 20);
    scf.param.set_user_defined_value("save", false);
    scf.param.set_user_defined_value("aobasis", std::string("sto-3g"));
    scf.param.set_user_defined_value("nv_extra", 3);
    scf.param.set_user_defined_value("derivatives", true);
    scf.param.set_user_defined_value("population", std::vector<std::string>{"iao"});
    check(scf_scientific_settings(scf) == base, "budgets, tolerances, acceleration, seed and output options do not change the equations");
    scf.param = old_parameters;
    scf.param.set_user_defined_value("protocol", std::vector<double>{1e-4, 1e-6});
    scf.param.set_user_defined_value("k", 8);
    scf.set_protocol<3>(world, 1e-6);
    check(scf_scientific_settings(scf) == base, "protocol representation is checked by the shared planner separately");
    scf.param = old_parameters;
    scf.set_protocol<3>(world, 1e-4);

    auto coords = scf.molecule.get_all_coords();
    coords(0, 0) = 0; coords(0, 1) = 0.1;
    scf.molecule.set_all_coords(coords);
    check(scf.molecule.get_all_coords().sumsq() == molecule.get_all_coords().sumsq(), "geometry fixture has equal coordinate norms");
    check(scf_scientific_settings(scf) != base, "distinct geometries with equal coordinate norms must differ");
    scf.molecule = molecule;
    coords = scf.molecule.get_all_coords(); coords(0, 0) += 1e-12;
    scf.molecule.set_all_coords(coords);
    check(scf_scientific_settings(scf) != base, "scientific geometry comparison has no hidden displacement tolerance");
    scf.molecule = molecule;
    scf.molecule.parameters.set_user_defined_value("source_name", std::string("another.xyz"));
    scf.molecule.parameters.set_user_defined_value("units", std::string("angstrom"));
    scf.molecule.parameters.set_user_defined_value("no_orient", true);
    check(scf_scientific_settings(scf) == base, "input provenance does not replace effective coordinates");
    scf.molecule = molecule;
    Tensor<double> field(3L); field(2) = 0.01;
    scf.molecule = Molecule(molecule.get_atoms(), molecule.get_eprec(), CorePotentialManager{}, field);
    check(scf.molecule.to_json() == molecule.to_json(), "legacy molecule JSON omits the effective constructor field");
    check(scf_scientific_settings(scf) != base, "effective electric field participates independently of input provenance");
    scf.molecule = molecule;
    scf.molecule.set_rcut(100.0);
    check(scf_scientific_settings(scf) != base, "effective smoothing-radius override participates");
    scf.molecule = molecule;
    scf.molecule.update_rcut_with_eprec(1e-7);
    check(scf_scientific_settings(scf) != base, "nuclear smoothing precision participates");
    scf.molecule = molecule;

    for (const auto& key : {"lo", "smear", "orbitalshift"}) {
        scf.param.set_user_defined_value(key, 0.01);
        check(scf_scientific_settings(scf) != base, "equation and residual operator parameters participate");
        scf.param = old_parameters;
    }
    scf.param.set_user_defined_value("deriv", std::string("ble"));
    check(scf_scientific_settings(scf) != base, "derivative convention participates");
    scf.param = old_parameters;
    scf.param.set_user_defined_value("aocc", std::vector<double>{0.0});
    check(scf_scientific_settings(scf) != base, "requested occupation holes participate");
    scf.param = old_parameters;
    scf.param.set_user_defined_value("dispersion", std::string("d3bj"));
    const auto dispersion = scf_scientific_settings(scf);
    check(dispersion != base, "dispersion model participates even before operator preparation");
    scf.param.set_user_defined_value("dispersion_atm", true);
    check(scf_scientific_settings(scf) != dispersion, "three-body dispersion choice participates");
    scf.param = old_parameters;
    if (DispersionCorrection::available()) {
        scf.dispersion = DispersionCorrection("d3bj", "pbe", "pbe", false);
        const auto two_body = scf_scientific_settings(scf);
        scf.dispersion = DispersionCorrection("d3bj", "pbe", "pbe", true);
        check(scf_scientific_settings(scf) != two_body, "actual initialized dispersion includes the three-body choice");
        scf.dispersion = DispersionCorrection{};
    }
    scf.xc.initialize("HF RHOTOL 1e-9 RHOMIN 1e-14", false, world, false);
    check(scf_scientific_settings(scf) != base, "actual initialized XC screening participates even if param.xc is unchanged");
    scf.xc.initialize(scf.param.xc(), false, world, false);
    check(scf_scientific_settings(scf) == base, "restored runtime XC returns the same scientific settings");
    scf.pcm_param.set_user_defined_value("epsilon", 12.0);
    check(scf_scientific_settings(scf) == base, "inactive PCM settings do not affect SCF reuse");
    scf.param.set_user_defined_value("pcm_data", std::string("explicit"));
    scf.pcm_param.set_derived_values(scf.param);
    const auto pcm = scf_scientific_settings(scf);
    check(!scf_settings_allow_reuse(pcm), "PCM requires a prepared solver with matching settings");
#ifdef MADNESS_HAS_PCM
    scf.pcm = PCM(world, scf.molecule, scf.pcm_param, false);
    check(scf_settings_allow_reuse(scf_scientific_settings(scf)), "prepared host PCM can qualify matching settings");
#endif
    scf.pcm_param.set_user_defined_value("epsilon", 15.0);
    check(scf_scientific_settings(scf) != pcm, "active resolved PCM dielectric participates");
    check(!scf_settings_allow_reuse(scf_scientific_settings(scf)), "changed PCM inputs cannot relabel a prepared operator");
#ifdef MADNESS_HAS_PCM
    scf.pcm = PCM(world, scf.molecule, scf.pcm_param, false);
    check(scf_settings_allow_reuse(scf_scientific_settings(scf)), "rebuilt PCM adopts the changed dielectric");
#endif
    scf.pcm_param.set_user_defined_value("cavity_type", std::string("restart"));
    check(!scf_settings_allow_reuse(scf_scientific_settings(scf)), "external cavity filename is insufficient for reuse");
    scf.param.set_user_defined_value("pcm_data", std::string("pcmsolver_reader_own"));
    scf.pcm_param.set_derived_values(scf.param);
    check(!scf_settings_allow_reuse(scf_scientific_settings(scf)), "unexposed external PCM settings require fresh evidence");
    scf.param = old_parameters;
    scf.pcm = PCM{};

    FutureParameters future;
    future.set_user_defined_value("print_level", 0);
    SCF future_scf(world, future, molecule);
    const auto before_future = scf_scientific_settings(future_scf);
    future_scf.param.set_user_defined_value("future_operator_threshold", 1e-8);
    check(scf_scientific_settings(future_scf) != before_future, "unclassified new parameters participate conservatively");

    Molecule core;
    core.add_atom(0, 0, 0, 4, 4);
    const auto* previous = std::getenv("MRA_DATA_DIR");
    const std::optional<std::string> directory = previous ? std::optional<std::string>(previous) : std::nullopt;
    setenv("MRA_DATA_DIR", MRA_CHEMDATA_DIR, 1);
    core.read_core_file("mcp");
    if (directory) setenv("MRA_DATA_DIR", directory->c_str(), 1); else unsetenv("MRA_DATA_DIR");
    const auto core_settings = core.electronic_settings();
    core.set_core_rcut(120.0);
    check(core.electronic_settings() != core_settings, "loaded core-potential smoothing participates");
    check(!core_settings.at("core_potentials").empty(), "snapshot contains loaded core coefficients, not their filename");

    Molecule pseudo;
    pseudo.add_atom(0, 0, 0, 2, 2, true);
    scf.molecule = pseudo;
    check(!scf_settings_allow_reuse(scf_scientific_settings(scf)), "unloaded GTH coefficients cannot qualify reuse");
    scf.gthpseudopotential = std::make_shared<GTHPseudopotential<double>>(world, pseudo);
    scf.gthpseudopotential->load_pseudo_from_file(world, std::string(MAD_ROOT_DIR) + "/src/apps/moldft/gth.xml");
    const auto gth = scf_scientific_settings(scf);
    check(scf_settings_allow_reuse(gth), "loaded GTH coefficients support an exact snapshot");
    scf.gthpseudopotential->localp[1](0L) += 0.125;
    check(scf_scientific_settings(scf) != gth, "changed GTH coefficients invalidate scientific equality");
}
}
int main(int argc, char** argv) {
    initialize(argc, argv);
    {
        World world(SafeMPI::COMM_WORLD);
        startup(world, argc, argv, true);
        try { exercise(world); }
        catch (const std::exception& error) { ++failures; std::cerr << error.what() << '\n'; }
        world.gop.max(failures);
    }
    finalize();
    return failures ? 1 : 0;
}
