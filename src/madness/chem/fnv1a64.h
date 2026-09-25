#ifndef MADNESS_CHEM_FNV1A64_H
#define MADNESS_CHEM_FNV1A64_H

// FNV-1a, 64-bit: a small, stable, non-cryptographic hash. Used where a value
// must be reproducible across builds and hosts, because it is written to disk
// (the molresponse ground-state archive fingerprint, the DALTON geometry hash).
// Do not change the constants: stored hashes would stop matching.

#include <cstddef>
#include <cstdint>

namespace madness {

inline std::uint64_t fnv1a64_update(std::uint64_t h, const char *p,
                                    std::size_t n) {
  for (std::size_t i = 0; i < n; ++i) {
    h ^= static_cast<unsigned char>(p[i]);
    h *= 1099511628211ULL;
  }
  return h;
}

inline constexpr std::uint64_t kFnv1a64Basis = 14695981039346656037ULL;

} // namespace madness

#endif // MADNESS_CHEM_FNV1A64_H
