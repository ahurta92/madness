#ifndef MOLRESPONSE_V3_SOLVERS_FREQUENCY_KEY_HPP
#define MOLRESPONSE_V3_SOLVERS_FREQUENCY_KEY_HPP

#include <cstdio>
#include <string>

namespace molresponse_v3 {

/// Canonical frequency key of a response state: f%.5f. It names driven-state
/// nodes (dipole_z-f0.05000), scheduler node ids and the snapshot's
/// fd_states entries, and is enough to distinguish any physically relevant
/// frequency we solve at; the exact frequency is in the node's settings.
inline std::string freq_key(double freq) {
  char buf[32];
  std::snprintf(buf, sizeof buf, "f%.5f", freq);
  return {buf};
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_FREQUENCY_KEY_HPP
