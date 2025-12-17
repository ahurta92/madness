#pragma once
#include <filesystem>
#include <iostream>
#include <madness/external/nlohmann_json/json.hpp>

#include "ResponseState.hpp"
#include "ResponseVector.hpp"

using json = nlohmann::json;
namespace fs = std::filesystem;

template <typename Desc>
void save_response_vector(World &world, const Desc &desc,
                          const ResponseVector &response) {
  auto filename = desc.response_filename();
  if (world.rank() == 0) {
    std::cout << "✅ Saving response vector to: " << filename << std::endl;
  }
  auto ar = archive::ParallelOutputArchive(world, filename.c_str());
  auto current_k = FunctionDefaults<3>::get_k();
  ar & current_k;

  std::visit(
      [&](const auto &vec) {
        int i = 0;

        for (const auto &f : vec.flat) {
          if (world.rank() == 0) {
            print("Saving response vector ", i++);
          }
          ar & f;
        }
      },
      response);
}

template <typename Desc>
inline bool load_response_vector(World &world, const int &num_orbitals,
                                 const Desc &desc, const size_t &thresh_index,
                                 const size_t &freq_index,
                                 ResponseVector &load_response) {
  auto filename = desc.response_filename(thresh_index, freq_index);
  if (!fs::exists(filename + ".00000")) {
    if (world.rank() == 0) {
      std::cerr << "❌ Response file not found: " << filename << std::endl;
    }
    return false;
  }
  bool debug = false;
  if (world.rank() == 0 and debug) {
    std::cout << "Loading response vector from: " << filename << std::endl;
  }

  load_response = make_response_vector(num_orbitals, desc.is_static(freq_index),
                                       !desc.is_spin_restricted());

  archive::ParallelInputArchive ar(world, filename.c_str());
  auto current_k = FunctionDefaults<3>::get_k();
  int loaded_k;

  ar & loaded_k;
  FunctionDefaults<3>::set_k(loaded_k);

  std::visit(
      [&](auto &v) {
        int i = 0;
        for (auto &f : v.flat) {
          // if (world.rank() == 0) {
          //   print("Loading response vector ", i++);
          // }
          ar & f;
        }
        v.sync();
      },
      load_response);
  FunctionDefaults<3>::set_k(current_k);
  std::visit(
      [&](auto &v) {
        auto &flat_vec = v.flat;

        if (current_k != loaded_k) {
          if (world.rank() == 0) {
            print("Reconstructing response vector with k=", current_k);
          }
          reconstruct(world, flat_vec);
          for (auto &v : flat_vec) {
            v = project(v, current_k, FunctionDefaults<3>::get_thresh(), true);
          }

          truncate(world, flat_vec, FunctionDefaults<3>::get_thresh());
        } else {
          truncate(world, flat_vec, FunctionDefaults<3>::get_thresh());
        }
        v.sync();
      },
      load_response);

  return true;
}

// Convenience overloads for a single LinearResponsePoint
inline bool load_response_vector(World &world, const int &num_orbitals,
                                 const LinearResponsePoint &pt,
                                 ResponseVector &load_response) {
  struct PointAdapter {
    const LinearResponsePoint &pt;
    std::string response_filename(size_t, size_t) const {
      return pt.response_filename();
    }
    bool is_static(size_t) const { return pt.is_static(); }
    bool is_spin_restricted() const { return pt.is_spin_restricted(); }
  } adaptor{pt};

  return load_response_vector(world, num_orbitals, adaptor,
                              /*thresh_index=*/0, /*freq_index=*/0,
                              load_response);
}

inline void save_response_vector(World &world, const LinearResponsePoint &pt,
                                 const ResponseVector &response) {
  struct PointAdapter {
    const LinearResponsePoint &pt;
    std::string response_filename() const {
      return pt.response_filename();
    }
  } adaptor{pt};

  save_response_vector(world, adaptor, response);
}
