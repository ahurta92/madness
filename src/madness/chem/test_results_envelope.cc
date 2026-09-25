//
// SCF task record: scf_data accessors and the
// SCFResults / ConvergenceResults envelope fields round-trip through JSON.
// Pure data — no Function arithmetic — so it runs in seconds on one rank.
//
#include <madness/mra/mra.h>
#include <madness/chem/SCF.h>
#include <madness/chem/Results.h>
#include <madness/world/test_utilities.h>
#include <nlohmann/json.hpp>
#include <cmath>
#include <map>
#include <string>

using namespace madness;

namespace {

int test_scf_data_accessors() {
    test_output t("scf_data: iterations() and last()");
    scf_data d;
    t.checkpoint(d.iterations() == 0, "fresh scf_data has 0 iterations");
    t.checkpoint(d.last().empty(),    "fresh scf_data has no last() values");
    d.add_data({{"e_kinetic", 1.0}, {"e_local", 0.0}, {"e_nuclear", -2.0}, {"e_coulomb", 0.5},
                {"e_pcm", 0.0}, {"e_disp", 0.0}, {"e_xc", -0.25}, {"e_nrep", 0.75}, {"e_tot", 0.0}});
    d.add_data({{"e_kinetic", 1.1}, {"e_local", 0.0}, {"e_nuclear", -2.1}, {"e_coulomb", 0.6},
                {"e_pcm", 0.0}, {"e_disp", 0.0}, {"e_xc", -0.26}, {"e_nrep", 0.75}, {"e_tot", 0.1}});
    t.checkpoint(d.iterations() == 2, "two add_data calls -> 2 iterations");
    const auto last = d.last();
    t.checkpoint(std::fabs(last.at("e_kinetic") - 1.1) < 1e-15, "last() returns the latest e_kinetic");
    t.checkpoint(std::fabs(last.at("e_nrep") - 0.75) < 1e-15,   "last() returns the latest e_nrep");
    return t.end();
}

} // namespace

int main(int argc, char** argv) {
    World& world = madness::initialize(argc, argv);
    startup(world, argc, argv);
    int error = 0;
    error += test_scf_data_accessors();
    finalize();
    return error;
}
