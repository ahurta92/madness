// The text summary (summary.txt) of hand-written results.json objects: the
// whole text must equal the expected text stored here.
#include <madness/chem/RunSummary.h>

#include <iostream>
#include <string>

namespace {

using json = nlohmann::json;

// Everything complete: a ground state and a response with alpha and one root.
const char* all_complete = R"json({
  "version": 5,
  "provenance": {"code_version": "900a79af8c1e2b3d4f5a6b7c8d9e0f1a2b3c4d5e-dirty",
                 "launched_at": "2026-11-26T01:51:07.123Z", "mpi_processes": 1,
                 "threads_per_process": 16, "archive_io": "native"},
  "calculations": [
    {"calculation_id": "ground_state", "application": "moldft", "model": "hf",
     "protocol": {"threshold": 0.0001, "k": 6},
     "targets": {"bsh_residual": 0.0005, "density_change": 0.0001, "energy_change": 1e-05},
     "convergence": {"bsh_residual": 2.16e-05, "density_change": 2.79e-05, "energy_change": 2.73e-08},
     "properties": {
       "total_energy": {"units": "hartree", "value_type": "real", "value": -1.1335935503913517, "status": "complete"},
       "orbital_energies": {"units": "hartree", "value_type": "real", "value": {"alpha": [-0.594894213856056]},
                            "status": "complete"},
       "energy_components": {"units": "hartree", "value_type": "real", "status": "complete",
                             "value": {"coulomb": 1.3178222663, "exchange_correlation": -0.6589122039,
                                       "kinetic": 1.1269592249, "nuclear_attraction": -3.6345671553,
                                       "nuclear_repulsion": 0.7151043390}},
       "dipole": {"units": "atomic_units", "value_type": "real", "value": [1.2e-14, -3.0e-15, 0.0],
                  "status": "complete"}
     },
     "state_metrics": {}},
    {"calculation_id": "response", "application": "molresponse",
     "protocol": {"threshold": 0.0001, "k": 6},
     "targets": {"driven_states": {"bsh_residual": 0.0005, "density_change": 0.0005},
                 "excited_states": {"bsh_residual": 0.0005, "density_change": 0.0005, "energy_change": 0.0005}},
     "properties": {
       "polarizability": {"units": "bohr^3", "value_type": "real", "values": [
         {"frequency": 0.0, "components": {"zz": 6.436397286910841}, "status": "complete",
          "convergence": {"bsh_residual": 2.25e-06, "density_change": 3.42e-05}}]},
       "excited_states": {"value_type": "real",
         "units": {"excitation_energy": "hartree", "oscillator_strength": "dimensionless",
                   "transition_dipole": "atomic_units", "transition_quadrupole": "atomic_units"},
         "values": [
           {"root": "es_root_0000", "excitation_energy": 0.4658334146727718,
            "oscillator_strength": 0.2854800950054045,
            "transition_dipole": {"x": 1.8e-15, "y": 1.9e-15, "z": -0.9587783682178349},
            "transition_quadrupole": {"xx": -9.9e-14, "zz": -1.1e-13},
            "status": "complete",
            "convergence": {"bsh_residual": 3.85e-05, "density_change": 0.000109, "energy_change": 1e-07}}]}
     },
     "run": {"stop_reason": "complete", "dropped_work": [], "stalled": []},
     "state_metrics": {}}
  ]
})json";

const char* all_complete_text = R"text(madqc summary
  code      900a79af8-dirty
  launched  2026-11-26 01:51 UTC
  run       1 process x 16 threads

Ground state (moldft, hf), level 1.0e-04 k 6
  reached   bsh_residual 2.16e-05   density_change 2.79e-05   energy_change 2.73e-08
  targets   bsh_residual 5.00e-04   density_change 1.00e-04   energy_change 1.00e-05

  property      value        units
  total_energy  -1.13359355  hartree

  dipole (atomic_units)
  x  y  z
  0  0  0

  energy_components (hartree)
  term                  value
  kinetic               1.12695922
  nuclear_attraction    -3.63456716
  coulomb               1.31782227
  exchange_correlation  -0.65891220
  nuclear_repulsion     0.71510434

  orbital_energies (hartree)
  spin   orbital  value
  alpha  0        -0.59489421

Response (molresponse), level 1.0e-04 k 6
  targets   driven:  bsh_residual 5.00e-04   density_change 5.00e-04
            excited: bsh_residual 5.00e-04   density_change 5.00e-04   energy_change 5.00e-04

  excited states (energy in hartree)
  root          energy      osc. strength  bsh_residual  density_change  energy_change
  es_root_0000  0.46583341  0.285480       3.85e-05      1.09e-04        1.00e-07

  transition dipole (atomic_units)
  root          x  y  z
  es_root_0000  0  0  -0.958778

  polarizability (bohr^3)
  frequency  component  value     bsh_residual  density_change
  0.00000    zz         6.436397  2.25e-06      3.42e-05

  not printed (see results.json): transition_quadrupole

  run: complete. Nothing dropped. Nothing stalled.
)text";

// An excited root that never converged, and beta work that was dropped.
const char* unavailable_and_dropped = R"json({
  "version": 5,
  "provenance": {"code_version": "unknown", "launched_at": "2026-11-26T02:00:00.000Z",
                 "mpi_processes": 2, "threads_per_process": 8, "archive_io": "native"},
  "calculations": [
    {"calculation_id": "response", "application": "molresponse",
     "protocol": {"threshold": 0.0001, "k": 6},
     "targets": {"bsh_residual": 0.0005, "density_change": 0.0005},
     "properties": {
       "polarizability": {"units": "bohr^3", "value_type": "real", "values": [
         {"frequency": 0.0, "components": {"xx": 4.5, "zz": 6.4}, "status": "complete",
          "convergence": {"bsh_residual": 2.0e-06, "density_change": 3.0e-05}},
         {"frequency": 0.05, "components": {"xx": 4.6}, "status": "partial",
          "availability": {"zz": "dipole_z-f0.05000 has not met the requested criteria"},
          "convergence": {"bsh_residual": 0.0031, "density_change": 0.012}}]},
       "excited_states": {"value_type": "real",
         "units": {"excitation_energy": "hartree", "oscillator_strength": "dimensionless"},
         "values": [
           {"root": "es_root_0000", "status": "unavailable",
            "availability": {"excitation_energy": "excited-full has not met the requested criteria",
                             "oscillator_strength": "excited-full has not met the requested criteria"},
            "convergence": {"bsh_residual": 0.021, "density_change": 0.0043}}]},
       "hyperpolarizability": {"units": "atomic_units", "value_type": "real", "values": [
         {"frequency_b": 0.0, "frequency_c": 0.05, "components": {}, "status": "unavailable",
          "availability": {"zzz": "vbc_dipole_z-dipole_z-f0.00000-f0.05000 has not met the requested criteria"},
          "convergence": {}}]}
     },
     "run": {"stop_reason": "complete_with_dropped_beta",
             "dropped_work": [{"id": "vbc_dipole_z-dipole_z-f0.00000-f0.05000", "top_protocol_key": "1e-04_k6",
                               "reason": "prerequisite dipole_z-f0.05000 did not converge at the top rung"}],
             "stalled": []},
     "state_metrics": {}}
  ]
})json";

const char* unavailable_and_dropped_text = R"text(NOT COMPLETE
  response: run stopped: complete_with_dropped_beta
  response: dropped: vbc_dipole_z-dipole_z-f0.00000-f0.05000: prerequisite dipole_z-f0.05000 did not
        converge at the top rung
  response excited_states root es_root_0000: excited-full has not met the requested criteria
  response hyperpolarizability frequency_b 0.00000 frequency_c 0.05000:
        vbc_dipole_z-dipole_z-f0.00000-f0.05000 has not met the requested criteria
  response polarizability frequency 0.05000: dipole_z-f0.05000 has not met the requested criteria

madqc summary
  code      unknown
  launched  2026-11-26 02:00 UTC
  run       2 processes x 8 threads

Response (molresponse), level 1.0e-04 k 6
  targets   bsh_residual 5.00e-04   density_change 5.00e-04

  excited states (energy in hartree)
  root          energy  osc. strength  bsh_residual  density_change
  es_root_0000  -       -              2.10e-02*     4.30e-03*

  hyperpolarizability (atomic_units)
  frequency_b  frequency_c  component  value
  0.00000      0.05000      zzz        -

  polarizability (bohr^3)
  frequency  component  value     bsh_residual  density_change
  0.00000    xx         4.500000  2.00e-06      3.00e-05
             zz         6.400000
  0.05000    xx         4.600000  3.10e-03*     1.20e-02*
             zz         -

  run: complete_with_dropped_beta. 1 dropped. Nothing stalled.

  * above target (driven states can pass under the scaled rule)
)text";

// Quadratic response, two-photon absorption and Raman, with one driven value
// above its target that the scaled rule still passed.
const char* quadratic_two_photon_raman = R"json({
  "version": 5,
  "provenance": {"code_version": "0123456789abcdef0123456789abcdef01234567", "launched_at": "2026-11-26T03:04:05.000Z",
                 "mpi_processes": 1, "threads_per_process": 4, "archive_io": "hdf5"},
  "calculations": [
    {"calculation_id": "response", "application": "molresponse",
     "protocol": {"threshold": 1e-06, "k": 8},
     "targets": {"driven_states": {"bsh_residual": 0.0001, "density_change": 0.0001},
                 "excited_states": {"bsh_residual": 0.0001, "density_change": 0.0001, "energy_change": 1e-05}},
     "properties": {
       "hyperpolarizability": {"units": "atomic_units", "value_type": "real", "values": [
         {"frequency_b": 0.02, "frequency_c": 0.02, "components": {"zzz": -2.3375398838787476e-16, "xxz": 12.3456789},
          "status": "complete", "convergence": {"bsh_residual": 0.000104, "density_change": 4.1e-05}}]},
       "two_photon_absorption": {"value_type": "real", "scientific_status": "preliminary",
         "units": {"energy": "hartree", "transition_tensor": "atomic_units"},
         "derived_units": {"D_linear": "atomic_units", "D_circular": "atomic_units",
                           "sigma_linear_gm": "GM", "sigma_circular_gm": "GM"},
         "values": [
           {"root": "es_root_0000", "energy": 0.32410279410481674, "frequency": 0.16205139705240837,
            "components": {"xx": -0.004145304270079069, "xz": -2.9349669604664177, "zz": -0.008598678300069702},
            "derived": {"D_linear": 2.297096666106714, "D_circular": 3.4456213046577853,
                        "sigma_linear_gm": 0.13091079971237918, "sigma_circular_gm": 0.19636484922650924},
            "status": "complete",
            "convergence": {"bsh_residual": 8.0e-05, "density_change": 4.29e-05, "energy_change": 2.27e-06}}]},
       "polarizability_gradient": {"units": "atomic_units", "value_type": "real", "values": [
         {"frequency": 0.0, "components": {"0x": {"xx": 1.0}, "0y": {"xx": 2.0}}, "status": "complete",
          "convergence": {"bsh_residual": 5.0e-05, "density_change": 6.0e-05}}]},
       "vibrational_analysis": {"value_type": "real", "status": "complete",
         "units": {"hessian": "Eh/bohr^2", "frequencies": "hartree", "reducedmass": "m_e"},
         "value": {"hessian": [[0.5, -0.5], [-0.5, 0.5]], "frequencies": [1.2e-12, 0.0193218],
                   "normal_modes": [[0.7, 0.7], [-0.7, 0.7]], "reduced_mass": [1.0, 1604.27]},
         "convergence": {"bsh_residual": 5.0e-05, "density_change": 6.0e-05}},
       "raman_spectrum": {"value_type": "real",
         "units": {"freq_cm1": "cm^-1", "pol_int": "Angstrom^4/amu (activity, linear polarization)",
                   "dep_ratio": "dimensionless (3 beta'^2 / activity)"},
         "values": [
           {"frequency": 0.0, "status": "complete",
            "value": [{"mode": 1, "freq_cm1": 1405.6271, "pol_int": 42.98765, "dep_ratio": 0.2499998, "converged": true}],
            "convergence": {"bsh_residual": 5.0e-05, "density_change": 6.0e-05}}]}
     },
     "run": {"stop_reason": "complete", "dropped_work": [], "stalled": []},
     "state_metrics": {}}
  ]
})json";

const char* quadratic_two_photon_raman_text = R"text(madqc summary
  code      012345678
  launched  2026-11-26 03:04 UTC
  run       1 process x 4 threads

Response (molresponse), level 1.0e-06 k 8
  targets   driven:  bsh_residual 1.00e-04   density_change 1.00e-04
            excited: bsh_residual 1.00e-04   density_change 1.00e-04   energy_change 1.00e-05

  hyperpolarizability (atomic_units)
  frequency_b  frequency_c  component  value      bsh_residual  density_change
  0.02000      0.02000      xxz        12.345679  1.04e-04*     4.10e-05
                            zzz        0

  polarizability_gradient (atomic_units)
  frequency  bsh_residual  density_change
  0.00000    5.00e-05      6.00e-05

  raman_spectrum (freq_cm1 cm^-1, pol_int Angstrom^4/amu, dep_ratio dimensionless)
  frequency  mode  freq_cm1     pol_int    dep_ratio  bsh_residual  density_change
  0.00000    1     1405.627100  42.987650  0.250000   5.00e-05      6.00e-05

  two_photon_absorption (energy hartree, transition_tensor atomic_units), preliminary
  root          frequency  component        value
  es_root_0000  0.16205    energy           0.32410279
                           xx               -0.004145
                           xz               -2.934967
                           zz               -0.008599
                           D_linear         2.297097 atomic_units
                           sigma_linear_gm  0.130911 GM

  root          frequency  bsh_residual  density_change  energy_change
  es_root_0000  0.16205    8.00e-05      4.29e-05        2.27e-06

  vibrational_analysis (frequency in hartree)
  mode  frequency   reduced_mass
  0     0           1.000000
  1     0.01932180  1604.270000

  not printed (see results.json): polarizability_gradient components, vibrational_analysis hessian,
    vibrational_analysis normal_modes

  run: complete. Nothing dropped. Nothing stalled.

  * above target (driven states can pass under the scaled rule)
)text";

// Derived quantities present: a full polarizability, a full
// hyperpolarizability (its 27 components are not printed, its derived
// quantities are), an excited root, the HOMO and the dipole magnitude.
const char* derived_present = R"json({
  "version": 5,
  "provenance": {"code_version": "unknown", "launched_at": "2026-11-26T04:00:00.000Z",
                 "mpi_processes": 1, "threads_per_process": 2, "archive_io": "native"},
  "calculations": [
    {"calculation_id": "ground_state", "application": "moldft", "model": "hf",
     "protocol": {"threshold": 0.0001, "k": 6},
     "convergence": {"bsh_residual": 2.0e-05},
     "properties": {
       "orbital_energies": {"units": "hartree", "value_type": "real", "status": "complete",
                            "value": {"alpha": [-20.5, -1.3, -0.49]}, "occupied": {"alpha": 2},
                            "orbitals": "canonical",
                            "derived": {"homo_energy": -1.3},
                            "definitions": {"homo_energy": "highest occupied"},
                            "derived_units": {"homo_energy": "hartree"}},
       "dipole": {"units": "atomic_units", "value_type": "real", "status": "complete", "value": [0.0, 0.0, 0.7],
                  "derived": {"dipole_magnitude": 0.7, "dipole_magnitude_debye": 1.7792225},
                  "derived_units": {"dipole_magnitude": "atomic_units", "dipole_magnitude_debye": "Debye"}}
     },
     "state_metrics": {}},
    {"calculation_id": "response", "application": "molresponse",
     "protocol": {"threshold": 0.0001, "k": 6},
     "targets": {"bsh_residual": 0.0005, "density_change": 0.0005},
     "properties": {
       "polarizability": {"units": "bohr^3", "value_type": "real",
         "derived_units": {"isotropic": "bohr^3", "anisotropy": "bohr^3"},
         "values": [
         {"frequency": 0.0, "components": {"xx": 1.0, "xy": 0.0, "xz": 0.0, "yx": 0.0, "yy": 2.0, "yz": 0.0,
                                            "zx": 0.0, "zy": 0.0, "zz": 3.0},
          "derived": {"isotropic": 2.0, "anisotropy": 1.7320508075688772},
          "status": "complete", "convergence": {"bsh_residual": 2.0e-06, "density_change": 3.0e-05}}]},
       "hyperpolarizability": {"units": "atomic_units", "value_type": "real",
         "derived_units": {"vector": "atomic_units", "vector_magnitude": "atomic_units", "parallel": "atomic_units",
                           "hrs": "atomic_units", "hrs_depolarization": "dimensionless"},
         "values": [
         {"frequency_b": 0.0, "frequency_c": 0.0,
          "components": {"xxx": 0, "xxy": 0, "xxz": 0, "xyx": 0, "xyy": 0, "xyz": 0, "xzx": 0, "xzy": 0, "xzz": 0,
                         "yxx": 0, "yxy": 0, "yxz": 0, "yyx": 0, "yyy": 0, "yyz": 0, "yzx": 0, "yzy": 0, "yzz": 0,
                         "zxx": 0, "zxy": 0, "zxz": 0, "zyx": 0, "zyy": 0, "zyz": 0, "zzx": 0, "zzy": 0, "zzz": 1.0},
          "derived": {"vector": {"x": 0.0, "y": 0.0, "z": 1.0}, "vector_magnitude": 1.0, "parallel": 0.6,
                      "hrs": 0.4140393356054125, "hrs_depolarization": 5.0},
          "status": "complete", "convergence": {"bsh_residual": 4.0e-06, "density_change": 5.0e-05}}]},
       "excited_states": {"value_type": "real",
         "units": {"excitation_energy": "hartree", "oscillator_strength": "dimensionless"},
         "derived_units": {"energy_ev": "eV", "wavelength_nm": "nm"},
         "values": [
           {"root": "es_root_0000", "excitation_energy": 0.4658334146727718, "oscillator_strength": 0.2854800950054045,
            "derived": {"energy_ev": 12.675858, "wavelength_nm": 97.810102},
            "status": "complete", "convergence": {"bsh_residual": 3.85e-05, "density_change": 0.000109}}]}
     },
     "run": {"stop_reason": "complete", "dropped_work": [], "stalled": [], "unconverged": []},
     "state_metrics": {}}
  ]
})json";

const char* derived_present_text = R"text(madqc summary
  code      unknown
  launched  2026-11-26 04:00 UTC
  run       1 process x 2 threads

Ground state (moldft, hf), level 1.0e-04 k 6
  reached   bsh_residual 2.00e-05

  dipole (atomic_units)
  x  y  z
  0  0  0.700000
  dipole_magnitude        0.700000
  dipole_magnitude_debye  1.77922 Debye

  orbital_energies (hartree)
  spin   orbital  value
  alpha  0        -20.50000000
  alpha  1        -1.30000000
  alpha  2        -0.49000000
  homo_energy  -1.30000000

Response (molresponse), level 1.0e-04 k 6
  targets   bsh_residual 5.00e-04   density_change 5.00e-04

  excited states (energy in hartree)
  root          energy      osc. strength  bsh_residual  density_change
  es_root_0000  0.46583341  0.285480       3.85e-05      1.09e-04
  es_root_0000  energy_ev      12.675858 eV
  es_root_0000  wavelength_nm  97.810102 nm

  hyperpolarizability (atomic_units)
  frequency_b  frequency_c  quantity            value             bsh_residual  density_change
  0.00000      0.00000      vector              (0, 0, 1.000000)  4.00e-06      5.00e-05
                            vector_magnitude    1.000000
                            parallel            0.600000
                            hrs                 0.414039
                            hrs_depolarization  5.000000

  polarizability (bohr^3)
  frequency  component   value     bsh_residual  density_change
  0.00000    xx          1.000000  2.00e-06      3.00e-05
             xy          0
             xz          0
             yx          0
             yy          2.000000
             yz          0
             zx          0
             zy          0
             zz          3.000000
             isotropic   2.000000
             anisotropy  1.732051

  not printed (see results.json): hyperpolarizability components

  run: complete. Nothing dropped. Nothing stalled.
)text";

// A run whose states did not converge: the run line counts them.
const char* unconverged = R"json({
  "version": 5,
  "provenance": {"code_version": "unknown", "launched_at": "2026-11-26T05:00:00.000Z",
                 "mpi_processes": 1, "threads_per_process": 2, "archive_io": "native"},
  "calculations": [
    {"calculation_id": "response", "application": "molresponse",
     "protocol": {"threshold": 0.0001, "k": 6},
     "targets": {"bsh_residual": 0.0005, "density_change": 0.0005},
     "properties": {
       "polarizability": {"units": "bohr^3", "value_type": "real", "values": [
         {"frequency": 0.0, "components": {}, "status": "unavailable",
          "availability": {"zz": "dipole_z-f0.00000 has not met the requested criteria"},
          "convergence": {"bsh_residual": 0.39, "density_change": 0.18}}]}
     },
     "run": {"stop_reason": "unconverged", "dropped_work": [], "stalled": [],
             "unconverged": ["dipole_z-f0.00000", "excited-full"]},
     "state_metrics": {}}
  ]
})json";

const char* unconverged_text = R"text(NOT COMPLETE
  response: run stopped: unconverged
  response polarizability frequency 0.00000: dipole_z-f0.00000 has not met the requested criteria

madqc summary
  code      unknown
  launched  2026-11-26 05:00 UTC
  run       1 process x 2 threads

Response (molresponse), level 1.0e-04 k 6
  targets   bsh_residual 5.00e-04   density_change 5.00e-04

  polarizability (bohr^3)
  frequency  component  value  bsh_residual  density_change
  0.00000    zz         -      3.90e-01*     1.80e-01*

  run: unconverged (2 states). Nothing dropped. Nothing stalled.

  * above target (driven states can pass under the scaled rule)
)text";

int check(const char* name, const char* input, const char* expected) {
    const auto text = madness::iteration::results_summary(json::parse(input));
    if (text == expected) {
        std::cout << "PASS  " << name << "\n";
        return 0;
    }
    std::cout << "FAIL  " << name << "\n--- expected\n" << expected << "--- got\n" << text << "---\n";
    return 1;
}

} // namespace

int main() {
    int failures = 0;
    failures += check("all complete", all_complete, all_complete_text);
    failures += check("unavailable item and dropped work", unavailable_and_dropped, unavailable_and_dropped_text);
    failures += check("hyperpolarizability, two-photon and Raman", quadratic_two_photon_raman,
                      quadratic_two_photon_raman_text);
    failures += check("derived quantities present", derived_present, derived_present_text);
    failures += check("unconverged states", unconverged, unconverged_text);
    // Every line fits the page.
    for (const auto* input : {all_complete, unavailable_and_dropped, quadratic_two_photon_raman, derived_present,
                              unconverged}) {
        const auto text = madness::iteration::results_summary(json::parse(input));
        std::size_t start = 0;
        while (start < text.size()) {
            const auto end = text.find('\n', start);
            if (end - start > 100) {
                std::cout << "FAIL  line longer than 100 characters: " << text.substr(start, end - start) << "\n";
                ++failures;
            }
            start = end + 1;
        }
    }
    return failures == 0 ? 0 : 1;
}
