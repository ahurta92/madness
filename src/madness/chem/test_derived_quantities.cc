// The formulas of DerivedQuantities.h against values known by hand, and the
// results.json layer of DerivedResults.h on hand-written entries.
#include <madness/chem/DerivedResults.h>

#include <cmath>
#include <iostream>
#include <random>
#include <string>

namespace {

using json = nlohmann::json;
namespace dq = madness::derived;

int failures = 0;

void expect(bool ok, const std::string& what) {
    std::cout << (ok ? "PASS  " : "FAIL  ") << what << "\n";
    if (!ok) ++failures;
}

bool near(double a, double b, double tolerance) { return std::abs(a - b) <= tolerance; }

dq::Tensor3 zero3() { return dq::Tensor3{}; }

/// A rotation matrix from three angles (ZYZ Euler).
dq::Tensor2 rotation(double alpha, double beta, double gamma) {
    const double ca = std::cos(alpha), sa = std::sin(alpha), cb = std::cos(beta), sb = std::sin(beta),
                 cg = std::cos(gamma), sg = std::sin(gamma);
    return {{{ca * cb * cg - sa * sg, -ca * cb * sg - sa * cg, ca * sb},
             {sa * cb * cg + ca * sg, -sa * cb * sg + ca * cg, sa * sb},
             {-sb * cg, sb * sg, cb}}};
}

dq::Tensor3 rotate(const dq::Tensor3& b, const dq::Tensor2& r) {
    dq::Tensor3 out{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) {
                double sum = 0.0;
                for (int l = 0; l < 3; ++l)
                    for (int m = 0; m < 3; ++m)
                        for (int n = 0; n < 3; ++n) sum += r[i][l] * r[j][m] * r[k][n] * b[l][m][n];
                out[i][j][k] = sum;
            }
    return out;
}

/// The same averages on a far finer grid (Gauss-Legendre with 16 points in
/// cos(beta), 32 points in alpha and gamma): the check that the 4 x 8 x 8 rule
/// is exact.
dq::HyperRayleighAverages fine_averages(const dq::Tensor3& b) {
    // Gauss-Legendre nodes and weights by Newton iteration on P_16.
    constexpr int ng = 16, nu = 32;
    std::array<double, ng> x{}, w{};
    for (int i = 0; i < ng; ++i) {
        double z = std::cos(M_PI * (i + 0.75) / (ng + 0.5)), dp = 0.0;
        for (int it = 0; it < 100; ++it) {
            double p0 = 1.0, p1 = z;
            for (int n = 2; n <= ng; ++n) {
                const double p2 = ((2.0 * n - 1.0) * z * p1 - (n - 1.0) * p0) / n;
                p0 = p1;
                p1 = p2;
            }
            dp = ng * (z * p1 - p0) / (z * z - 1.0);
            const double step = p1 / dp;
            z -= step;
            if (std::abs(step) < 1e-16) break;
        }
        x[i] = z;
        w[i] = 2.0 / ((1.0 - z * z) * dp * dp);
    }
    dq::HyperRayleighAverages sum;
    double total = 0.0;
    for (int a = 0; a < nu; ++a)
        for (int c = 0; c < nu; ++c)
            for (int g = 0; g < ng; ++g) {
                const auto r = rotation(2 * M_PI * a / nu, std::acos(x[g]), 2 * M_PI * c / nu);
                double zzz = 0.0, zxx = 0.0;
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                        for (int k = 0; k < 3; ++k) {
                            zzz += r[2][i] * r[2][j] * r[2][k] * b[i][j][k];
                            zxx += r[2][i] * r[0][j] * r[0][k] * b[i][j][k];
                        }
                sum.zzz += w[g] * zzz * zzz;
                sum.zxx += w[g] * zxx * zxx;
                total += w[g];
            }
    sum.zzz /= total;
    sum.zxx /= total;
    return sum;
}

void test_polarizability() {
    expect(near(dq::polarizability_isotropic(1, 2, 3), 2.0, 1e-15), "polarizability (1,2,3): isotropic 2");
    expect(near(dq::polarizability_anisotropy(1, 2, 3, 0, 0, 0), std::sqrt(3.0), 1e-15),
           "polarizability (1,2,3): anisotropy sqrt(3)");
}

void test_hyperpolarizability() {
    auto b = zero3();
    b[2][2][2] = 1.0;
    const auto v = dq::hyperpolarizability_vector(b);
    expect(v[0] == 0.0 && v[1] == 0.0 && near(v[2], 1.0, 1e-15), "beta_zzz = 1: vector (0, 0, 1)");
    const auto a = dq::hyper_rayleigh_averages(b);
    expect(near(a.zzz, 1.0 / 7.0, 1e-14), "beta_zzz = 1: <beta_ZZZ^2> = 1/7");
    expect(near(a.zxx, 1.0 / 35.0, 1e-14), "beta_zzz = 1: <beta_ZXX^2> = 1/35");
    expect(near(dq::hyper_rayleigh_depolarization(a), 5.0, 1e-12), "beta_zzz = 1: depolarization 5");
    expect(near(dq::hyper_rayleigh(a), std::sqrt(6.0 / 35.0), 1e-14), "beta_zzz = 1: hrs sqrt(1/7 + 1/35)");
    expect(near(dq::hyperpolarizability_parallel(v, {0.0, 0.0, -2.0}), -0.6, 1e-15),
           "beta_zzz = 1, dipole along -z: parallel -3/5");

    auto o = zero3();
    o[0][0][0] = 1.0;
    o[0][1][1] = o[1][0][1] = o[1][1][0] = -1.0;
    const auto ov = dq::hyperpolarizability_vector(o);
    expect(dq::norm(ov) < 1e-15, "octupolar: vector 0");
    expect(near(dq::hyper_rayleigh_depolarization(dq::hyper_rayleigh_averages(o)), 1.5, 1e-12),
           "octupolar: depolarization 1.5");

    // A random tensor with no index symmetry, and the same tensor rotated.
    std::mt19937 generator(20261126); // NOLINT(cert-msc51-cpp): a fixed seed keeps the test reproducible
    std::uniform_real_distribution<double> uniform(-1.0, 1.0);
    dq::Tensor3 r{};
    for (auto& plane : r)
        for (auto& row : plane)
            for (auto& x : row) x = uniform(generator);
    const auto R = rotation(2.0 * M_PI * (0.5 + 0.5 * uniform(generator)), M_PI * (0.5 + 0.5 * uniform(generator)),
                            2.0 * M_PI * (0.5 + 0.5 * uniform(generator)));
    const auto rr = rotate(r, R);
    const auto ar = dq::hyper_rayleigh_averages(r), arr = dq::hyper_rayleigh_averages(rr);
    expect(near(dq::norm(dq::hyperpolarizability_vector(r)), dq::norm(dq::hyperpolarizability_vector(rr)), 1e-12),
           "random tensor rotated: same vector_magnitude");
    expect(near(dq::hyper_rayleigh(ar), dq::hyper_rayleigh(arr), 1e-12), "random tensor rotated: same hrs");
    expect(near(dq::hyper_rayleigh_depolarization(ar), dq::hyper_rayleigh_depolarization(arr), 1e-12),
           "random tensor rotated: same hrs_depolarization");
    const auto fine = fine_averages(r);
    expect(near(ar.zzz, fine.zzz, 1e-13) && near(ar.zxx, fine.zxx, 1e-13),
           "random tensor: the 4 x 8 x 8 rule equals a 16 x 32 x 32 grid");
}

/// The strength written with all three terms:
/// delta = F dF + G dG + H dH, dF = sum S_aa S_bb / 30, dG = sum S_ab S_ab / 30,
/// dH = sum S_ab S_ba / 30.
double stated_strength(const dq::Tensor2& S, double F, double G, double H) {
    double f = 0, g = 0, h = 0;
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) {
            f += S[a][a] * S[b][b];
            g += S[a][b] * S[a][b];
            h += S[a][b] * S[b][a];
        }
    return (F * f + G * g + H * h) / 30.0;
}

void test_two_photon() {
    const dq::Tensor2 diagonal{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    const auto d = dq::two_photon_strengths(diagonal);
    expect(near(d.Df, 0.3, 1e-15) && near(d.Dg, 0.1, 1e-15), "S = (1,1,1): Df 3/10, Dg 1/10");
    expect(near(d.D_linear, stated_strength(diagonal, 2, 2, 2), 1e-15) && near(d.D_linear, 1.0, 1e-15),
           "S = (1,1,1): linear strength 1 = 2 dF + 2 dG + 2 dH");
    expect(near(d.D_circular, stated_strength(diagonal, -2, 3, 3), 1e-15) && near(d.D_circular, 0.0, 1e-15),
           "S = (1,1,1): circular strength 0 = -2 dF + 3 dG + 3 dH");
    // A symmetric S, as the code builds: the two conventions agree.
    const dq::Tensor2 symmetric{{{-0.0041, 0.0002, -2.93}, {0.0002, 0.0008, 1e-6}, {-2.93, 1e-6, -0.0086}}};
    const auto s = dq::two_photon_strengths(symmetric);
    expect(near(s.D_linear, stated_strength(symmetric, 2, 2, 2), 1e-14) &&
               near(s.D_circular, stated_strength(symmetric, -2, 3, 3), 1e-14),
           "symmetric S: strengths equal the stated convention");

    // The cross section as its definitions state it, with CODATA constants:
    // sigma = N pi^2 alpha a0^5 omega^2 D g / c, N = 8, omega = E/2, g = 1/Gamma,
    // Gamma = 0.1 eV (full width at half maximum); a0 in pm and c in cm/s give
    // GM. It agrees with the code's DALTON constants to their last digits.
    const double E = 0.5;
    const double gamma = 0.1 / madness::constants::hartree_electron_volt_relationship;
    const double a0_pm = madness::constants::atomic_unit_of_length * 1e12;
    const double c_cms = madness::constants::speed_of_light_in_vacuum * 100.0;
    const double stated = 8.0 * M_PI * M_PI * madness::constants::fine_structure_constant * std::pow(a0_pm, 5.0) *
                          (E / 2) * (E / 2) * d.D_linear / c_cms / gamma;
    const double code = dq::two_photon_cross_section_gm(d.D_linear, E);
    std::cout << "      S = (1,1,1), E = 0.5: sigma " << code << " GM; from the stated formula " << stated << " GM\n";
    expect(near(code, stated, stated * 1e-6), "cross section equals its stated formula (N = 8, g = 1/Gamma, FWHM)");
}

void test_results_layer() {
    using madness::iteration::add_derived_quantities;
    auto results = json::parse(R"json({"calculations": [
      {"calculation_id": "ground_state", "properties": {
         "dipole": {"units": "atomic_units", "status": "complete", "value": [0.0, 0.0, 0.7]},
         "orbital_energies": {"units": "hartree", "status": "complete", "value": {"alpha": [-20.5, -1.3, -0.49]},
                              "occupied": {"alpha": 2}, "orbitals": "canonical"}}},
      {"calculation_id": "response", "properties": {
         "polarizability": {"units": "bohr^3", "values": [
           {"frequency": 0.0, "status": "complete", "components": {"zz": 6.4}},
           {"frequency": 0.1, "status": "complete",
            "components": {"xx": 1, "yy": 2, "zz": 3, "xy": 0.1, "yx": 0.3, "xz": 0, "zx": 0, "yz": 0}},
           {"frequency": 0.2, "status": "partial", "components": {"xx": 1, "yy": 2, "zz": 3}}]},
         "hyperpolarizability": {"units": "atomic_units", "values": [
           {"frequency_b": 0.0, "frequency_c": 0.0, "status": "complete", "components": {"zzz": 1.0}}]},
         "excited_states": {"units": {"excitation_energy": "hartree"},
           "definitions": {"transition_quadrupole": "nontraceless second moment"},
           "values": [
           {"root": "es_root_0000", "status": "complete", "excitation_energy": 0.5},
           {"root": "es_root_0001", "status": "unavailable"}]}}}]})json");
    // The full hyperpolarizability of beta_zzz = 1.
    auto& beta = results["calculations"][1]["properties"]["hyperpolarizability"]["values"][0]["components"];
    for (const char i : {'x', 'y', 'z'})
        for (const char j : {'x', 'y', 'z'})
            for (const char k : {'x', 'y', 'z'})
                if (!beta.contains(std::string{i, j, k})) beta[std::string{i, j, k}] = 0.0;
    add_derived_quantities(results);
    const auto& ground = results["calculations"][0]["properties"];
    const auto& response = results["calculations"][1]["properties"];
    expect(near(ground["dipole"]["derived"]["dipole_magnitude"].get<double>(), 0.7, 1e-15) &&
               near(ground["dipole"]["derived"]["dipole_magnitude_debye"].get<double>(), 0.7 * 2.541746, 1e-6),
           "dipole: magnitude in atomic units and Debye");
    expect(ground["orbital_energies"]["derived"]["homo_energy"] == -1.3, "orbital energies: homo of 2 occupied");
    const auto& alpha = response["polarizability"]["values"];
    expect(!alpha[0].contains("derived"), "polarizability zz only: no derived quantities");
    expect(near(alpha[1]["derived"]["isotropic"].get<double>(), 2.0, 1e-15) &&
               near(alpha[1]["derived"]["anisotropy"].get<double>(), std::sqrt(3.0 + 3.0 * 0.04), 1e-14),
           "polarizability full: isotropic, anisotropy with xy, yx as their mean");
    expect(!alpha[2].contains("derived"), "polarizability partial: no derived quantities");
    expect(response["polarizability"]["derived_units"]["isotropic"] == "bohr^3" &&
               response["polarizability"]["definitions"].contains("anisotropy"),
           "polarizability: definitions and units");
    const auto& b = response["hyperpolarizability"]["values"][0]["derived"];
    expect(b["vector"]["z"] == 1.0 && near(b["hrs_depolarization"].get<double>(), 5.0, 1e-12) &&
               near(b["parallel"].get<double>(), 0.6, 1e-15),
           "hyperpolarizability: vector, depolarization, parallel along the dipole");
    const auto& es = response["excited_states"];
    expect(near(es["values"][0]["derived"]["energy_ev"].get<double>(), 0.5 * 27.211386245981, 1e-12) &&
               near(es["values"][0]["derived"]["wavelength_nm"].get<double>(), 1e7 / (0.5 * 219474.63136314), 1e-9) &&
               !es["values"][1].contains("derived"),
           "excited states: eV and nm for the complete root only");
    expect(es["definitions"].size() == 3 && es["definitions"].contains("transition_quadrupole"),
           "excited states: definitions added beside the existing one");

    // Written again, nothing changes; a value that is no longer complete, or
    // a dipole that is gone, loses what was derived from it.
    const auto once = results;
    add_derived_quantities(results);
    expect(results == once, "derived quantities are the same when written again");
    results["calculations"][1]["properties"]["polarizability"]["values"][1]["status"] = "partial";
    results["calculations"][0]["properties"].erase("dipole");
    add_derived_quantities(results);
    const auto& again = results["calculations"][1]["properties"];
    expect(!again["polarizability"]["values"][1].contains("derived") && !again["polarizability"].contains("definitions"),
           "a value no longer complete loses its derived quantities and definitions");
    expect(!again["hyperpolarizability"]["values"][0]["derived"].contains("parallel"),
           "hyperpolarizability: no parallel without a dipole");

    // A noise-level tensor: no depolarization ratio of two noise numbers.
    auto noise = json::parse(R"json({"calculations": [{"calculation_id": "response", "properties": {
      "hyperpolarizability": {"units": "atomic_units", "values": [
        {"frequency_b": 0.0, "frequency_c": 0.0, "status": "complete", "components": {"zzz": 1e-14, "xyy": -3e-15}}]}}}]})json");
    auto& small = noise["calculations"][0]["properties"]["hyperpolarizability"]["values"][0]["components"];
    for (const char i : {'x', 'y', 'z'})
        for (const char j : {'x', 'y', 'z'})
            for (const char k : {'x', 'y', 'z'})
                if (!small.contains(std::string{i, j, k})) small[std::string{i, j, k}] = 0.0;
    add_derived_quantities(noise);
    const auto& quiet = noise["calculations"][0]["properties"]["hyperpolarizability"]["values"][0]["derived"];
    expect(quiet.contains("hrs") && !quiet.contains("hrs_depolarization"),
           "noise-level hyperpolarizability (hrs < 1e-8): no hrs_depolarization");

    // Localized orbitals: no HOMO unless a spin has one occupied orbital;
    // open shell: one per spin.
    auto localized = json::parse(R"json({"calculations": [{"calculation_id": "ground_state", "properties": {
      "orbital_energies": {"units": "hartree", "status": "complete", "value": {"alpha": [-1.0, -0.5]},
                           "occupied": {"alpha": 2}, "orbitals": "localized"}}}]})json");
    add_derived_quantities(localized);
    expect(!localized["calculations"][0]["properties"]["orbital_energies"].contains("derived"),
           "localized orbitals, 2 occupied: no homo_energy");
    auto open = json::parse(R"json({"calculations": [{"calculation_id": "ground_state", "properties": {
      "orbital_energies": {"units": "hartree", "status": "complete",
                           "value": {"alpha": [-1.0, -0.6, 0.1], "beta": [-0.9]},
                           "occupied": {"alpha": 2, "beta": 1}, "orbitals": "localized"}}}]})json");
    add_derived_quantities(open);
    expect(open["calculations"][0]["properties"]["orbital_energies"].contains("derived") == false,
           "localized open shell with 2 alpha occupied: no homo_energy");
    open["calculations"][0]["properties"]["orbital_energies"]["orbitals"] = "canonical";
    add_derived_quantities(open);
    expect(open["calculations"][0]["properties"]["orbital_energies"]["derived"]["homo_energy"] ==
               json({{"alpha", -0.6}, {"beta", -0.9}}),
           "canonical open shell: homo per spin");
}

} // namespace

int main() {
    test_polarizability();
    test_hyperpolarizability();
    test_two_photon();
    test_results_layer();
    std::cout << (failures == 0 ? "all passed" : std::to_string(failures) + " failed") << "\n";
    return failures == 0 ? 0 : 1;
}
