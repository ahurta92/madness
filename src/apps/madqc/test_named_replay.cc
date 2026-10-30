// The from-store tier, level 1 (D-078, D-079): each stage's verdict
// recomputed from recorded files, no engine. A case is a stage of a run the
// madqc tests or the response qctest twins made: its deck, its final
// manifest (reduced to what the verdict reads) and the verdict the run
// recorded. The skeleton's shared_verdict, with the AppSpec rows' fields,
// must give the recorded verdict for every case.
//
// Usage: test_named_replay <fixture.json>
//   test_named_replay_fixture.json beside this file, written by
//   replay_verdict.py (pack); it changes only with a decided record change.
#include "MoldftCheckpoint.hpp"
#include "NamedStages.hpp"
#include "ResponseCheckpoint.hpp"
#include "RunNamedStages.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

using json = nlohmann::json;

int main(int argc, char** argv) {
  auto& world = madness::initialize(argc, argv);
  int failures = 0;
  std::size_t cases = 0;
  try {
    if (argc < 2)
      throw std::runtime_error("usage: test_named_replay <fixture.json>");
    std::ifstream in(argv[1]);
    if (!in)
      throw std::runtime_error(std::string("cannot read ") + argv[1]);
    const auto fixture = json::parse(in);
    const auto& decks = fixture.at("decks");
    for (const auto& c : fixture.at("cases")) {
      ++cases;
      const auto id = c.at("id").get<std::string>();
      try {
        const auto stages = named_stages::parse(
            world,
            decks.at(c.at("deck").get<std::string>()).get<std::string>());
        const auto stage
            = std::find_if(stages.begin(), stages.end(), [&](const auto& s) {
                return s.name == c.at("stage").get<std::string>();
              });
        if (stage == stages.end())
          throw std::runtime_error("no stage " + c.at("stage").dump());
        // The rows' fields alone, as test_named_stages tables them.
        const auto spec = stage->kind == "moldft"
                              ? named_stages::moldft_fields()
                              : named_stages::response_fields();
        const int natom = c.at("natom").get<int>();
        const auto manifest
            = madness::checkpoint::from_json_record(c.at("manifest"));
        const bool replayed
            = named_stages::shared_verdict(
                  spec, *stage, natom, spec.settings(*stage, natom), manifest)
                  .at("converged")
                  .get<bool>();
        if (replayed != c.at("recorded").get<bool>()) {
          ++failures;
          if (world.rank() == 0)
            std::cerr << id << ": recorded " << c.at("recorded")
                      << ", replayed " << replayed << '\n';
        }
      } catch (const std::exception& e) {
        ++failures;
        if (world.rank() == 0)
          std::cerr << id << ": " << e.what() << '\n';
      }
    }
  } catch (const std::exception& e) {
    ++failures;
    if (world.rank() == 0)
      std::cerr << e.what() << '\n';
  }
  if (cases == 0)
    ++failures; // an empty fixture proves nothing
  if (world.rank() == 0)
    std::cout << "replayed verdicts: " << cases << " cases, " << failures
              << " failures\n";
  madness::finalize();
  return failures ? 1 : 0;
}
