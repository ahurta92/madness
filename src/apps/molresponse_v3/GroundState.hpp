#ifndef MOLRESPONSE_V3_GROUNDSTATE_HPP
#define MOLRESPONSE_V3_GROUNDSTATE_HPP

#include <madness/chem/SCF.h>
#include <madness/chem/projector.h>
#include <madness/world/world.h>
#include <memory>
#include <string>

namespace molresponse_v3 {

using namespace madness;

/// Thin wrapper around SCF for response property calculations.
///
/// Delegates all ground-state data access to a shared SCF object.
/// Extends with response-specific cached operators (V_local, Hamiltonian,
/// QProjector) that are rebuilt per protocol step via prepare().
///
/// Supports both closed-shell (restricted) and open-shell (unrestricted)
/// through SCF's native alpha/beta orbital storage.
class GroundState {
public:
    /// Construct from an existing SCF calculation (shared ownership).
    /// Used by the madqc workflow where SCF is already available.
    explicit GroundState(World& world, std::shared_ptr<SCF> scf);

    /// Factory: load from a moldft restart archive.
    /// Reads the archive header to extract L, xc, etc., constructs an SCF
    /// with matching parameters, then delegates to SCF::load_mos().
    static GroundState from_archive(World& world,
                                     const std::string& archive_path,
                                     const Molecule& molecule);

    // -----------------------------------------------------------------
    // Delegated accessors — no local copies, all point to scf_->
    // -----------------------------------------------------------------

    const Molecule& molecule() const { return scf_->molecule; }

    /// Alpha orbitals (convenience alias for restricted case)
    const vecfuncT& orbitals() const { return scf_->get_amo(); }
    const vecfuncT& orbitals_alpha() const { return scf_->get_amo(); }
    const vecfuncT& orbitals_beta() const { return scf_->get_bmo(); }

    /// Orbital energies
    const tensorT& energies() const { return scf_->aeps; }
    const tensorT& energies_alpha() const { return scf_->aeps; }
    const tensorT& energies_beta() const { return scf_->beps; }

    /// Occupation numbers
    const tensorT& occupations_alpha() const { return scf_->get_aocc(); }
    const tensorT& occupations_beta() const { return scf_->get_bocc(); }

    /// Orbital counts
    long num_orbitals() const { return static_cast<long>(scf_->get_amo().size()); }
    long num_alpha() const { return static_cast<long>(scf_->get_amo().size()); }
    long num_beta() const { return static_cast<long>(scf_->get_bmo().size()); }

    bool is_spin_restricted() const { return scf_->is_spin_restricted(); }
    double L() const { return scf_->param.L(); }
    int k() const { return FunctionDefaults<3>::get_k(); }
    const CalculationParameters& params() const { return scf_->param; }
    double hf_exchange_coefficient() const { return scf_->xc.hf_exchange_coefficient(); }

    // -----------------------------------------------------------------
    // Response-specific cached operators (built by prepare())
    // -----------------------------------------------------------------

    /// Prepare ground-state operators for the current protocol.
    /// Must be called after FunctionDefaults<3> are set for the current
    /// protocol step. Reprojects orbitals if k changed, then computes
    /// V_local, Hamiltonian, and QProjector.
    void prepare(World& world, double vtol,
                 const poperatorT& coulop,
                 const std::string& fock_json_file = "");

    /// Local potential: V_nuc + 2*V_coul + V_xc
    const real_function_3d& V_local() const;

    /// Full Fock/Hamiltonian matrix
    const tensorT& hamiltonian() const;

    /// Fock matrix with diagonal zeroed (off-diagonal coupling)
    const tensorT& hamiltonian_no_diag() const;

    /// Projector onto virtual space: Q = 1 - |phi><phi|
    const QProjector<double, 3>& Q() const;

    // -----------------------------------------------------------------
    // Direct SCF access for advanced use
    // -----------------------------------------------------------------

    SCF& scf() { return *scf_; }
    const SCF& scf() const { return *scf_; }

    void print_info() const;

private:
    std::shared_ptr<SCF> scf_;
    int original_k_;

    // Cached response-specific data (rebuilt per protocol step)
    real_function_3d v_local_;
    tensorT hamiltonian_;
    tensorT hamiltonian_no_diag_;
    QProjector<double, 3> q_projector_;
    bool prepared_ = false;

    /// Read just the archive header to extract L, k, xc, etc.
    /// Used by from_archive to set up CalculationParameters before load_mos.
    struct ArchiveHeader {
        unsigned int version = 0;
        double energy = 0.0;
        bool spin_restricted = true;
        double L = 0.0;
        int k = 0;
        std::string xc;
        std::string localize_method;
        double converged_for_thresh = 0.0;
        unsigned int nmo_alpha = 0;
    };

    static ArchiveHeader read_archive_header(World& world,
                                              const std::string& archive_path);
};

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_GROUNDSTATE_HPP
