#pragma once

#include <algorithm>
#include <cctype>
#include <fstream>
#include <madness/chem/CheckpointStore.h>
#include <madness/mra/commandlineparser.h>
#include <set>
#include <sstream>
#include <vector>

namespace madqc {
// Validate only the outer group structure here. Existing typed readers own
// parameter syntax, values, precedence and engine-specific validation.
inline void validate_workflow_input(madness::World &world,
                                    const madness::commandlineparser &parser) {
  madness::checkpoint::detail::on_root(world, [&] {
    for (const auto *key : {"stage", "from-store", "from_store"})
      if (parser.key_exists(key))
        throw std::runtime_error(std::string("unknown input option: --") + key);
    std::ifstream file(parser.value_raw("input"));
    if (!file)
      return std::string{}; // CLI-only input is supported.
    std::ostringstream contents;
    contents << file.rdbuf();
    const auto text = contents.str();
    const std::set<std::string> groups{
        "dft",  "response", "nemo",     "optimization", "oep", "cc2",
        "tdhf", "pcm",      "molecule", "geometry",     "io"};
    if (nlohmann::json::accept(text)) {
      const auto value = nlohmann::json::parse(text);
      if (!value.is_object())
        throw std::runtime_error("unknown input: expected parameter groups");
      for (const auto &[name, unused] : value.items())
        if (!groups.contains(name))
          throw std::runtime_error("unknown input block: " + name);
      return std::string{};
    }
    std::istringstream input(text);
    std::string line;
    std::vector<std::string> open_groups;
    while (std::getline(input, line)) {
      line = line.substr(0, line.find('#'));
      std::istringstream words(line);
      std::string first, extra;
      if (!(words >> first))
        continue;
      auto canonical = first;
      std::transform(canonical.begin(), canonical.end(), canonical.begin(),
                     [](unsigned char c) { return std::tolower(c); });
      if (open_groups.empty()) {
        // Typed group readers locate canonical, case-sensitive headers.
        // Reject a spelling they would otherwise silently skip.
        if (!groups.contains(first) || words >> extra)
          throw std::runtime_error("unknown input block: " + line);
        open_groups.push_back(first);
      } else if (canonical == "end") {
        if (words >> extra)
          throw std::runtime_error("unknown input: " + line);
        open_groups.pop_back();
      } else if (open_groups.back() == "geometry" && first == "molecule") {
        // Preserve the existing nested geometry/molecule deck form.
        if (words >> extra)
          throw std::runtime_error("unknown input block: " + line);
        open_groups.push_back(first);
      } else if (groups.contains(canonical)) {
        throw std::runtime_error("unknown input: expected end before " + line);
      }
    }
    if (!open_groups.empty())
      throw std::runtime_error("unknown input: unterminated " +
                               open_groups.back() + " block");
    return std::string{};
  });
}
} // namespace madqc
