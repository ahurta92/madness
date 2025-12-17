#include "ResponseState.hpp"

// ─────────────────────────────────────────────────────────────────────────────
// LinearResponseDescriptor
// ─────────────────────────────────────────────────────────────────────────────

LinearResponseDescriptor::LinearResponseDescriptor(
    Perturbation pert, const std::vector<double> &freq,
    const std::vector<double> &thresh, bool spin_restricted_in)
    : perturbation(pert), spin_restricted(spin_restricted_in),
      frequencies(freq), thresholds(thresh) {
  for (size_t i = 0; i < frequencies.size(); ++i) {
    frequency_map[frequencies[i]] = i;
  }
}

size_t LinearResponseDescriptor::num_thresholds() const {
  return thresholds.size();
}

size_t LinearResponseDescriptor::num_frequencies() const {
  return frequencies.size();
}

double LinearResponseDescriptor::threshold(size_t ti) const {
  return thresholds.at(ti);
}

double LinearResponseDescriptor::frequency(size_t fi) const {
  return frequencies.at(fi);
}

bool LinearResponseDescriptor::is_static(size_t fi) const {
  return std::abs(frequencies[fi]) < 1e-8;
}

bool LinearResponseDescriptor::is_spin_restricted() const {
  return spin_restricted;
}

std::string LinearResponseDescriptor::make_key(double thresh,
                                               double freq) const {
  std::ostringstream oss;
  double f = std::clamp(freq, 0.0, 100.0);

  oss << describe_perturbation(perturbation)
      << "_f" << std::fixed << std::setprecision(3) << f
      << "_p" << std::scientific << std::setprecision(2) << thresh;

  return oss.str();
}

std::string LinearResponseDescriptor::make_key(size_t ti, size_t fi) const {
  return make_key(thresholds[ti], frequencies[fi]);
}

std::string LinearResponseDescriptor::response_filename(
    const size_t &thresh_index, const size_t &freq_index) const {
  return make_key(thresh_index, freq_index);
}

std::string LinearResponseDescriptor::perturbationDescription() const {
  return describe_perturbation(perturbation);
}

// ─────────────────────────────────────────────────────────────────────────────
// LinearResponsePoint
// ─────────────────────────────────────────────────────────────────────────────

double LinearResponsePoint::threshold() const {
  return desc.threshold(thresh_index);
}

double LinearResponsePoint::frequency() const {
  return desc.frequency(freq_index);
}

bool LinearResponsePoint::is_static() const {
  return desc.is_static(freq_index);
}

bool LinearResponsePoint::is_spin_restricted() const {
  return desc.is_spin_restricted();
}

std::string LinearResponsePoint::response_filename() const {
  return desc.make_key(threshold(), frequency());
}

std::string LinearResponsePoint::perturbationDescription() const {
  return desc.perturbationDescription();
}

// ─────────────────────────────────────────────────────────────────────────────
// SecondOrderResponseDescriptor
// ─────────────────────────────────────────────────────────────────────────────

SecondOrderResponseDescriptor::SecondOrderResponseDescriptor(
    Perturbation p1, Perturbation p2, double f1, double f2, double thr,
    bool spin_restricted_in)
    : perturbations_(p1, p2), frequencies_(f1, f2), thresh(thr),
      spin_restricted_(spin_restricted_in) {}

LinearResponseDescriptor SecondOrderResponseDescriptor::B_state() const {
  return LinearResponseDescriptor(perturbations_.first, {frequencies_.first},
                                  {thresh}, spin_restricted_);
}

LinearResponseDescriptor SecondOrderResponseDescriptor::C_state() const {
  return LinearResponseDescriptor(perturbations_.second, {frequencies_.second},
                                  {thresh}, spin_restricted_);
}

std::pair<LinearResponseDescriptor, LinearResponseDescriptor>
SecondOrderResponseDescriptor::get_states() const {
  return {B_state(), C_state()};
}

double SecondOrderResponseDescriptor::current_threshold() const {
  return thresh;
}

double SecondOrderResponseDescriptor::current_frequency() const {
  return frequencies_.first + frequencies_.second;
}

std::string SecondOrderResponseDescriptor::perturbationDescription() const {
  return describe_perturbation(perturbations_.first) + "_" +
         describe_perturbation(perturbations_.second);
}

bool SecondOrderResponseDescriptor::is_spin_restricted() const {
  return spin_restricted_;
}

bool SecondOrderResponseDescriptor::is_static(size_t) const { return false; }

std::string SecondOrderResponseDescriptor::make_key(double f1, double f2,
                                                    double thr) const {
  std::ostringstream oss;

  f1 = std::clamp(f1, 0.0, 100.0);
  f2 = std::clamp(f2, 0.0, 100.0);

  oss << prefix() << describe_perturbation(perturbations_.first) << "_"
      << describe_perturbation(perturbations_.second) << "_f" << std::fixed
      << std::setprecision(3) << f1 << "_" << std::fixed
      << std::setprecision(3) << f2 << "_p" << std::scientific
      << std::setprecision(1) << thr;

  return oss.str();
}

std::string SecondOrderResponseDescriptor::response_filename() const {
  return make_key(frequencies_.first, frequencies_.second, thresh) +
         ".response";
}

std::string SecondOrderResponseDescriptor::response_filename(
    const size_t &, const size_t &) const {
  return response_filename();
}

ResponseVector
SecondOrderResponseDescriptor::make_vector(int num_orbitals, size_t) const {
  bool urstr = !spin_restricted_;
  return make_response_vector(num_orbitals, false, urstr);
}

// ─────────────────────────────────────────────────────────────────────────────
// Perturbation operators and helpers
// ─────────────────────────────────────────────────────────────────────────────

int dir_index(char c) {
  switch (std::tolower(static_cast<unsigned char>(c))) {
  case 'x':
    return 0;
  case 'y':
    return 1;
  case 'z':
    return 2;
  }
  throw std::invalid_argument("Direction must be x/y/z");
}

real_function_3d
make_perturbation_operator(World &world, const GroundStateData &g_s,
                           const DipolePerturbation &d) {
  std::map<char, int> dipole_map = {{'x', 0}, {'y', 1}, {'z', 2}};
  std::vector<int> dir(3, 0);
  dir[dipole_map.at(d.direction)] = 1;
  real_function_3d f =
      real_factory_3d(world).functor(real_functor_3d{new MomentFunctor(dir)});
  f.truncate(FunctionDefaults<3>::get_thresh());
  return f;
}

real_function_3d
make_perturbation_operator(World &world, const GroundStateData &gs,
                           const NuclearDisplacementPerturbation &n) {
  std::map<char, int> dipole_map = {{'x', 0}, {'y', 1}, {'z', 2}};

  madchem::MolecularDerivativeFunctor mdfunctor(gs.molecule, n.atom_index,
                                                dipole_map.at(n.direction));

  real_function_3d dvdx = real_factory_3d(world)
                              .functor(mdfunctor)
                              .truncate_on_project()
                              .truncate_mode(0);

  dvdx.get_impl()->set_truncate_mode(1);

  return dvdx;
}

real_function_3d
make_perturbation_operator(World &world, const GroundStateData &,
                           const MagneticPerturbation &) {
  throw std::runtime_error("Magnetic perturbation not implemented yet");
}

real_function_3d raw_perturbation_operator(World &world,
                                           const GroundStateData &gs,
                                           const Perturbation &p) {
  return std::visit(
      [&](const auto &pp) -> real_function_3d {
        return make_perturbation_operator(world, gs, pp);
      },
      p);
}

vector_real_function_3d
project_perturbation_onto_orbitals(World &world, const GroundStateData &gs,
                                   const real_function_3d &raw_op) {
  auto vp = mul(world, raw_op, gs.orbitals, /*fence=*/true);
  vp = gs.Qhat(vp);
  truncate(world, vp, FunctionDefaults<3>::get_thresh(), /*fence=*/true);
  return vp;
}

// Linear (first-order) response
vector_real_function_3d perturbation_vector(World &world,
                                            GroundStateData const &gs,
                                            const LinearResponsePoint &pt) {
  auto raw_op = raw_perturbation_operator(world, gs, pt.desc.perturbation);
  auto Vp = project_perturbation_onto_orbitals(world, gs, raw_op);
  if (!pt.is_static()) {
    Vp.insert(Vp.end(), Vp.begin(), Vp.end());
  }
  return Vp;
}

// Descriptor-only overload: used in property code where we explicitly decide
// static/dynamic outside (e.g. duplicate for ω ≠ 0 in PropertyManager).
vector_real_function_3d perturbation_vector(World &world,
                                            GroundStateData const &gs,
                                            LinearResponseDescriptor const &s) {
  auto raw_op = raw_perturbation_operator(world, gs, s.perturbation);
  return project_perturbation_onto_orbitals(world, gs, raw_op);
}

// Second-order (VBC) response (not yet implemented)
vector_real_function_3d perturbation_vector(World &, GroundStateData const &,
                                            XBCResponseState const &) {
  return {};
}

