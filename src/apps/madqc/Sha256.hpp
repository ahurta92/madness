#pragma once

// SHA-256 (FIPS 180-4) of a file, for the Dalton seed's member kind (S5,
// skeleton spec §3.5): the seed block records each copied file by size and
// SHA-256, its own form (the storage core's crc64 is for native archives).
// Self-contained: MADNESS links no hash library.

#include <array>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace named_stages {

namespace sha256_detail {
inline constexpr std::array<std::uint32_t, 64> k
    = {0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1,
       0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
       0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786,
       0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
       0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147,
       0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
       0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
       0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
       0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a,
       0x5b9cca4f, 0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
       0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};

inline std::uint32_t rotr(std::uint32_t x, int n) {
  return (x >> n) | (x << (32 - n));
}

class Hash {
  std::array<std::uint32_t, 8> h_
      = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
         0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
  std::array<unsigned char, 64> block_{};
  std::size_t used_ = 0;
  std::uint64_t bits_ = 0;

  void compress() {
    std::array<std::uint32_t, 64> w{};
    for (int i = 0; i < 16; ++i)
      w[i] = (std::uint32_t(block_[4 * i]) << 24)
             | (std::uint32_t(block_[4 * i + 1]) << 16)
             | (std::uint32_t(block_[4 * i + 2]) << 8)
             | std::uint32_t(block_[4 * i + 3]);
    for (int i = 16; i < 64; ++i) {
      const auto s0
          = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
      const auto s1
          = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
      w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }
    auto [a, b, c, d, e, f, g, h] = h_;
    for (int i = 0; i < 64; ++i) {
      const auto t1 = h + (rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25))
                      + ((e & f) ^ (~e & g)) + k[i] + w[i];
      const auto t2 = (rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22))
                      + ((a & b) ^ (a & c) ^ (b & c));
      h = g;
      g = f;
      f = e;
      e = d + t1;
      d = c;
      c = b;
      b = a;
      a = t1 + t2;
    }
    const std::array<std::uint32_t, 8> add = {a, b, c, d, e, f, g, h};
    for (int i = 0; i < 8; ++i)
      h_[i] += add[i];
  }

public:
  void update(const unsigned char* data, std::size_t n) {
    for (std::size_t i = 0; i < n; ++i) {
      block_[used_++] = data[i];
      if (used_ == 64) {
        compress();
        used_ = 0;
      }
    }
    bits_ += 8 * static_cast<std::uint64_t>(n);
  }

  std::string hex() {
    const auto length = bits_;
    const unsigned char one = 0x80, zero = 0;
    update(&one, 1);
    while (used_ != 56)
      update(&zero, 1);
    for (int i = 7; i >= 0; --i) {
      const auto byte = static_cast<unsigned char>(length >> (8 * i));
      update(&byte, 1);
    }
    std::ostringstream out;
    for (const auto word : h_)
      out << std::hex << std::setw(8) << std::setfill('0') << word;
    return out.str();
  }
};
} // namespace sha256_detail

// The SHA-256 of a byte string, as 64 lowercase hex digits.
[[nodiscard]] inline std::string sha256_hex(const std::string& bytes) {
  sha256_detail::Hash hash;
  hash.update(reinterpret_cast<const unsigned char*>(bytes.data()),
              bytes.size());
  return hash.hex();
}

// The SHA-256 of a file's contents; throws when the file cannot be read.
[[nodiscard]] inline std::string sha256_file(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in)
    throw std::runtime_error("cannot read " + path);
  sha256_detail::Hash hash;
  std::vector<char> buffer(1 << 16);
  while (in) {
    in.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    hash.update(reinterpret_cast<const unsigned char*>(buffer.data()),
                static_cast<std::size_t>(in.gcount()));
  }
  if (in.bad())
    throw std::runtime_error("read failed: " + path);
  return hash.hex();
}

} // namespace named_stages
