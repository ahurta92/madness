#include "GroundState.hpp"

#include <madness/chem/SCFOperators.h>
#include <madness/chem/xcfunctional.h>
#include <madness/tensor/tensor_json.hpp>
#include <madness/world/MADworld.h>

#include <apps/molresponse_v2/broadcast_json.hpp>

#include <filesystem>

namespace molresponse_v3 {

using namespace madness;

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

GroundState::GroundState(World& world, std::shared_ptr<SCF> scf)
    : scf_(std::move(scf)),
      original_k_(FunctionDefaults<3>::get_k()) {
}

GroundState GroundState::from_archive(World& world,
                                       const std::string& archive_path,
                                       const Molecule& molecule) {
    // Step 1: Read archive header to get L, xc, etc.
    auto header = read_archive_header(world, archive_path);

    // Step 2: Build CalculationParameters with values from the archive.
    // SCF::load_mos checks L == param.L(), so we must set L correctly.
    // SCF constructor calls set_derived_values(molecule) which computes
    // nalpha/nbeta from nuclear charge + nopen, so we must set nopen
    // for open-shell systems.
    CalculationParameters params;
    params.set_user_defined_value("l", header.L);
    params.set_user_defined_value("xc", header.xc);

    // Derive prefix from archive path: strip ".restartdata" suffix
    std::string prefix = archive_path;
    auto pos = prefix.rfind(".restartdata");
    if (pos != std::string::npos) {
        prefix = prefix.substr(0, pos);
    }
    params.set_user_defined_value("prefix", prefix);

    // For open-shell: infer nopen from nmo_alpha and total electrons.
    // Archive stores nmo_alpha directly; nopen = 2*nmo_alpha - nelec.
    // set_derived_values uses nopen to compute nalpha/nbeta.
    if (!header.spin_restricted) {
        int nelec = static_cast<int>(molecule.total_nuclear_charge());
        int nopen = 2 * static_cast<int>(header.nmo_alpha) - nelec;
        if (nopen < 0) nopen = 0;
        params.set_user_defined_value("nopen", nopen);
    }

    // Step 3: Construct SCF
    auto scf = std::make_shared<SCF>(world, params, molecule);

    // Step 4: Load orbitals via SCF's own archive reader
    scf->load_mos(world);

    // Step 5: Build GroundState
    GroundState gs(world, scf);
    gs.original_k_ = header.k;
    return gs;
}

GroundState::ArchiveHeader
GroundState::read_archive_header(World& world,
                                  const std::string& archive_path) {
    ArchiveHeader h;
    archive::ParallelInputArchive<archive::BinaryFstreamInputArchive>
        ar(world, archive_path.c_str());

    ar & h.version;
    MADNESS_CHECK(h.version == 4);
    ar & h.energy;
    ar & h.spin_restricted;
    ar & h.L;
    ar & h.k;

    // Read molecule (needed to advance archive position, but we discard it
    // since the caller provides the canonical molecule)
    Molecule mol_discard;
    ar & mol_discard;
    ar & h.xc;
    ar & h.localize_method;
    ar & h.converged_for_thresh;

    // Read nmo_alpha (needed to infer nopen for open-shell)
    ar & h.nmo_alpha;

    return h;
}

// ---------------------------------------------------------------------------
// Response-specific preparation
// ---------------------------------------------------------------------------

void GroundState::prepare(World& world, double vtol,
                           const poperatorT& coulop,
                           const std::string& fock_json_file) {
    auto current_k = FunctionDefaults<3>::get_k();
    auto thresh = FunctionDefaults<3>::get_thresh();

    // Re-project orbitals if k or thresh changed
    if (original_k_ != current_k) {
        if (original_k_ < current_k) {
            MADNESS_EXCEPTION(
                "Cannot project orbitals to higher k than archive", original_k_);
        }
        // Reload from archive at original k, then project down
        scf_->load_mos(world);
        reconstruct(world, scf_->amo);
        for (auto& orbital : scf_->amo) {
            orbital = project(orbital, current_k, thresh, true);
        }
        truncate(world, scf_->amo, thresh);

        if (!scf_->is_spin_restricted()) {
            reconstruct(world, scf_->bmo);
            for (auto& orbital : scf_->bmo) {
                orbital = project(orbital, current_k, thresh, true);
            }
            truncate(world, scf_->bmo, thresh);
        }
    } else {
        truncate(world, scf_->amo, thresh);
        if (!scf_->is_spin_restricted()) {
            truncate(world, scf_->bmo, thresh);
        }
    }

    // Build QProjector from current orbitals
    q_projector_ = QProjector<double, 3>(scf_->get_amo());

    // Compute density
    auto density = scf_->make_density(world, scf_->aocc, scf_->amo);

    // Nuclear potential (uses SCF's potentialmanager)
    scf_->make_nuclear_potential(world);
    world.gop.fence();

    auto V_nuc = scf_->potentialmanager->vnuclear();
    V_nuc.truncate(vtol);

    // Coulomb potential
    auto V_coul = 2.0 * apply(*coulop, density, true);
    V_coul.truncate(vtol);

    // Assemble local potential
    v_local_ = V_nuc + V_coul;

    // XC potential (for DFT with non-pure HF exchange)
    if (scf_->xc.is_dft() && scf_->xc.hf_exchange_coefficient() != 1.0) {
        XCOperator<double, 3> xc_op(world, scf_->param.xc(),
                                     scf_->is_spin_restricted(),
                                     density, density);
        v_local_ += xc_op.make_xc_potential();
    }

    V_nuc.clear();
    V_coul.clear();
    density.clear();

    // Hamiltonian / Fock matrix
    bool loaded_from_file = false;
    if (!fock_json_file.empty() && std::filesystem::exists(fock_json_file)) {
        auto fock_json = broadcast_json_file(world, fock_json_file);
        world.gop.fence();

        auto protocol_key = std::string("thresh: ") + std::to_string(thresh)
                          + std::string(" k: ") + std::to_string(current_k);
        if (fock_json.contains(protocol_key)) {
            hamiltonian_ = tensor_from_json<double>(fock_json[protocol_key]["focka"]);
            loaded_from_file = true;
        }
    }

    if (!loaded_from_file) {
        const auto& phi = scf_->get_amo();
        auto phi_copy = copy(world, phi, true);
        reconstruct(world, phi_copy);

        // Kinetic energy
        real_derivative_3d Dx(world, 0), Dy(world, 1), Dz(world, 2);
        auto fx = apply(world, Dx, phi_copy);
        auto fy = apply(world, Dy, phi_copy);
        auto fz = apply(world, Dz, phi_copy);
        compress(world, fx, true);
        compress(world, fy, true);
        compress(world, fz, true);
        world.gop.fence();
        tensorT T = 0.5 * (matrix_inner(world, fx, fx)
                          + matrix_inner(world, fy, fy)
                          + matrix_inner(world, fz, fz));
        fx.clear(); fy.clear(); fz.clear();

        // Potential energy
        auto V_local_phi = mul_sparse(world, v_local_, phi, vtol);
        vecfuncT V_hf_phi = zero_functions<double, 3>(world, phi.size());
        if (scf_->xc.hf_exchange_coefficient() > 0.0) {
            const double lo = 1.e-10;
            Exchange<double, 3> K(world, lo);
            auto phi_k = copy(world, phi, true);
            K.set_algorithm(Exchange<double, 3>::ExchangeAlgorithm::multiworld_efficient_row);
            K.set_bra_and_ket(phi_k, phi_k);
            V_hf_phi = -scf_->xc.hf_exchange_coefficient() * K(phi_k);
        }

        auto V_phi = gaxpy_oop(1.0, V_local_phi, 1.0, V_hf_phi);
        truncate(world, V_phi);
        auto phi_V_phi = matrix_inner(world, phi, V_phi);

        hamiltonian_ = T + phi_V_phi;
    }

    // Build off-diagonal Hamiltonian
    hamiltonian_no_diag_ = copy(hamiltonian_);
    for (long i = 0; i < num_orbitals(); i++) {
        hamiltonian_no_diag_(i, i) = 0.0;
    }

    prepared_ = true;

    if (world.rank() == 0) {
        print_info();
    }
}

// ---------------------------------------------------------------------------
// Cached operator accessors
// ---------------------------------------------------------------------------

const real_function_3d& GroundState::V_local() const {
    MADNESS_CHECK(prepared_);
    return v_local_;
}

const tensorT& GroundState::hamiltonian() const {
    MADNESS_CHECK(prepared_);
    return hamiltonian_;
}

const tensorT& GroundState::hamiltonian_no_diag() const {
    MADNESS_CHECK(prepared_);
    return hamiltonian_no_diag_;
}

const QProjector<double, 3>& GroundState::Q() const {
    MADNESS_CHECK(prepared_);
    return q_projector_;
}

// ---------------------------------------------------------------------------
// Info
// ---------------------------------------------------------------------------

void GroundState::print_info() const {
    print("GROUND_STATE_INFO (v3)");
    print("  xc =", scf_->param.xc());
    print("  spin_restricted =", scf_->is_spin_restricted());
    print("  num_alpha =", num_alpha());
    if (!scf_->is_spin_restricted()) {
        print("  num_beta =", num_beta());
    }
    print("  box_l =", L());
    print("  k =", k());
    print("  orbital_energies_alpha =", scf_->aeps);
    if (!scf_->is_spin_restricted() && scf_->beps.size() > 0) {
        print("  orbital_energies_beta =", scf_->beps);
    }
    print("  prepared =", prepared_);
    print("------------------------------------------------------------");
}

} // namespace molresponse_v3
