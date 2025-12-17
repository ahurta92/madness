#ifndef RESPONSE_STATE_HPP
#define RESPONSE_STATE_HPP
#include <SCF.h>
#include <madness/chem/projector.h>
#include <madness/mra/mra.h>

#include "ResponseVector.hpp"
#include <madness/external/nlohmann_json/json.hpp>
#include <sstream>
#include <string>
#include <variant>
#include <vector>

#include "GroundStateData.hpp"
#include "Perturbation.hpp"
#include "molecular_functors.h"
#include "vmra.h"

using json = nlohmann::json;
namespace fs = std::filesystem;

struct AbstractResponseDescriptor {
  [[nodiscard]] virtual bool is_spin_restricted() const = 0;
  [[nodiscard]] virtual std::string
  response_filename(const size_t &thresh_index,
                    const size_t &freq_index) const = 0;

  virtual ~AbstractResponseDescriptor() = default;
};

struct LinearResponseDescriptor : public AbstractResponseDescriptor {
  Perturbation perturbation;

  bool spin_restricted = false; // Is the system open shell?

  std::vector<double> frequencies;
  std::map<double, size_t> frequency_map; // Frequency to index map
  std::vector<double> thresholds;         // Accuracy levels to loop over

  LinearResponseDescriptor() = default;

  LinearResponseDescriptor(Perturbation pert, const std::vector<double> &freq,
                           const std::vector<double> &thresh,
                           bool spin_restricted);

  [[nodiscard]] size_t num_thresholds() const;
  [[nodiscard]] size_t num_frequencies() const;

  [[nodiscard]] double threshold(size_t ti) const;
  [[nodiscard]] double frequency(size_t fi) const;

  [[nodiscard]] bool is_static(size_t fi) const;

  [[nodiscard]] bool is_spin_restricted() const override;

  // helper that builds the core "<perturbation>_p<thresh>_f<freq>"
  [[nodiscard]] std::string make_key(double thresh, double freq) const;

  [[nodiscard]] std::string make_key(size_t ti, size_t fi) const;

  [[nodiscard]] std::string
  response_filename(const size_t &thresh_index,
                    const size_t &freq_index) const override;

  [[nodiscard]] std::string perturbationDescription() const;
};

struct LinearResponsePoint {
  const LinearResponseDescriptor &desc;
  size_t thresh_index;
  size_t freq_index;

  [[nodiscard]] double threshold() const;
  [[nodiscard]] double frequency() const;

  [[nodiscard]] bool is_static() const;
  [[nodiscard]] bool is_spin_restricted() const;

  [[nodiscard]] std::string response_filename() const;

  [[nodiscard]] std::string perturbationDescription() const;
};

struct SecondOrderResponseDescriptor : public AbstractResponseDescriptor {
  std::pair<Perturbation, Perturbation> perturbations_;
  std::pair<double, double> frequencies_;

  double thresh;
  bool spin_restricted_ = false; // Is the system open shell?

  SecondOrderResponseDescriptor(Perturbation p1, Perturbation p2, double f1,
                                double f2, double thresh, bool spin_restricted);

  [[nodiscard]] LinearResponseDescriptor B_state() const;

  [[nodiscard]] LinearResponseDescriptor C_state() const;

  [[nodiscard]] std::pair<LinearResponseDescriptor, LinearResponseDescriptor>
  get_states() const;

  [[nodiscard]] double current_threshold() const;

  [[nodiscard]] double current_frequency() const;

  [[nodiscard]] virtual const char *prefix() const = 0;

  [[nodiscard]] std::string perturbationDescription() const;

  [[nodiscard]] bool is_spin_restricted() const override;

  [[nodiscard]] bool is_static(size_t freq_index) const;

  // Build the core "<prefix><pert1>_<pert2>_f<f1>_<f2>_p<thresh>"
  [[nodiscard]] std::string make_key(double f1, double f2,
                                     double thresh) const;

  [[nodiscard]] std::string response_filename() const;

  [[nodiscard]] std::string
  response_filename(const size_t &thresh_index,
                    const size_t &freq_index) const override;

  [[nodiscard]] ResponseVector make_vector(int num_orbitals, size_t fi) const;
};

//-----------------------------------------------------------------------------
// Now two trivial subclasses for VBC vs. XBC
//-----------------------------------------------------------------------------

struct VBCResponseState : public SecondOrderResponseDescriptor {
  using SecondOrderResponseDescriptor::SecondOrderResponseDescriptor;
  [[nodiscard]] const char *prefix() const override { return "VBC"; }
};

struct XBCResponseState : public SecondOrderResponseDescriptor {
  using SecondOrderResponseDescriptor::SecondOrderResponseDescriptor;
  [[nodiscard]] const char *prefix() const override { return "XBC"; }
};

int dir_index(char c);

real_function_3d make_perturbation_operator(World &world,
                                            const GroundStateData &g_s,
                                            const DipolePerturbation &d);

real_function_3d
make_perturbation_operator(World &world, const GroundStateData &gs,
                           const NuclearDisplacementPerturbation &n);

real_function_3d make_perturbation_operator(World &world,
                                            const GroundStateData &gs,
                                            const MagneticPerturbation &m);

real_function_3d raw_perturbation_operator(World &world,
                                           const GroundStateData &gs,
                                           const Perturbation &p);

// -----------------------------------------------------------------------------
// 1) Raw operator in real space, before applying to orbitals:
//
//    e.g. for a dipole:  V(r) = x, y or z moment
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// 2) Apply it to the ground‐state orbitals to get your Vp basis functions.
//    (You already have this in ResponseState::perturbation_vector.)
// -----------------------------------------------------------------------------
vector_real_function_3d
project_perturbation_onto_orbitals(World &world, const GroundStateData &gs,
                                   const real_function_3d &raw_op);

// Linear (first-order) response
// For a specific (threshold, frequency) point, decide static/dynamic based
// on that frequency.
madness::vector_real_function_3d
perturbation_vector(madness::World &world, GroundStateData const &gs,
                    const LinearResponsePoint &pt);

// Descriptor-only overload: used in property code where we explicitly decide
// static/dynamic outside (e.g. duplicate for ω ≠ 0 in PropertyManager).
madness::vector_real_function_3d
perturbation_vector(madness::World &world, GroundStateData const &gs,
                    LinearResponseDescriptor const &state);

// Second-order (VBC) response
madness::vector_real_function_3d
perturbation_vector(madness::World &world, GroundStateData const &gs,
                    XBCResponseState const &sos);

#endif // RESPONSE_STATE_HPP
