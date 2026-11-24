#pragma once
#include "response_metadata.hpp"
#include <madness/chem/RootCoordination.h>

namespace molresponse_v3 {
// Alpha's legacy analysis record remains available during shared adoption.
inline nlohmann::json read_response_property_metadata(madness::World &world,
                                                      const std::string &path) {
  const auto record = madness::checkpoint::detail::on_root(world, [&] {
    return ResponseMetadata::load_or_create(path).json().dump();
  });
  return nlohmann::json::parse(record);
}
inline void save_alpha_property(madness::World &world, const std::string &path,
                                const std::string &protocol,
                                const nlohmann::json &row) {
  madness::checkpoint::detail::on_root(world, [&] {
    auto metadata = ResponseMetadata::load_or_create(path);
    metadata.add_property("alpha", protocol, row);
    metadata.save();
    return std::string{};
  });
}
inline void replace_alpha_properties(madness::World& world,const std::string& path,
                                     const std::string& protocol,
                                     const nlohmann::json& rows) {
  madness::checkpoint::detail::on_root(world,[&] {
    auto metadata=ResponseMetadata::load_or_create(path);
    metadata.replace_property("alpha",protocol,rows);
    metadata.save();
    return std::string{};
  });
}
} // namespace molresponse_v3
