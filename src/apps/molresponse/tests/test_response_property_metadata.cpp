#include "../solvers/response_property_metadata.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <madness/world/MADworld.h>

int main(int argc, char **argv) {
  auto &world = madness::initialize(argc, argv);
  const bool write = argc > 1 && std::string(argv[1]) == "write";
  std::string directory;
  if (world.rank() == 0) {
    char pattern[] = "/tmp/response-property-XXXXXX";
    const auto *created = mkdtemp(pattern);
    MADNESS_CHECK_THROW(created, "temporary directory creation failed");
    directory = created;
  }
  world.gop.broadcast_serializable(directory, 0);
  const auto path = directory + "/response_metadata.json";
  nlohmann::json initial;
  if (world.rank() == 0) {
    initial = molresponse_v3::ResponseMetadata::load_or_create(path).json();
    if (write) {
      molresponse_v3::ResponseMetadata::load_or_create(path).save();
      std::filesystem::create_directory(path + ".tmp");
    } else {
      std::ofstream file(path);
      file << "{";
    }
  }
  world.gop.fence();
  std::string error;
  try {
    if (write)
      molresponse_v3::save_alpha_property(world, path, "1e-04_k6",
                                          {{"omega", 0.0}, {"alpha", {{1.0}}}});
    else
      (void)molresponse_v3::read_response_property_metadata(world, path);
  } catch (const std::exception &e) {
    error = e.what();
  }
  std::string root_error = error;
  world.gop.broadcast_serializable(root_error, 0);
  int failures = error.empty() || error != root_error;
  if (world.rank() == 0 && write)
    failures += molresponse_v3::ResponseMetadata::load_or_create(path).json() !=
                initial;
  world.gop.sum(failures);
  world.gop.fence();
  if (world.rank() == 0)
    std::filesystem::remove_all(directory);
  madness::finalize();
  return failures ? 1 : 0;
}
