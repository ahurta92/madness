// Which saved ground-state Fock matrix a prepare may read (solvers/saved_fock.hpp).
// A node's file is read only at the level its description records, compared as
// numbers: 1e-7 and 1e-8 share the six-decimal entry name, and a node saved at
// 1e-7 must not be read at 1e-8.
#include "../solvers/saved_fock.hpp"

#include <cstdio>
#include <string>

using molresponse_v3::choose_saved_fock;
using molresponse_v3::fock_level_name;
using madness::iteration::Protocol;
using json = nlohmann::json;

namespace {
int failures = 0;
void check(bool ok, const std::string& what) {
  std::printf("%s  %s\n", ok ? "PASS" : "FAIL", what.c_str());
  if (!ok) ++failures;
}
json node_file(double thresh, int k) {
  return {{fock_level_name(thresh, k), {{"focka", json::object()}}}};
}
}  // namespace

int main() {
  const int k = 10;
  check(fock_level_name(1e-7, k) == fock_level_name(1e-8, k), "1e-7 and 1e-8 share an entry name");

  const auto saved_1e7 = node_file(1e-7, k);
  const auto at_1e8 = choose_saved_fock(saved_1e7, Protocol{1e-7, k}, 1e-8, k);
  check(at_1e8.entry.empty() && !at_1e8.reason.empty(), "a node saved at 1e-7 is not read at 1e-8: " + at_1e8.reason);
  const auto at_1e7 = choose_saved_fock(saved_1e7, Protocol{1e-7, k}, 1e-7, k);
  check(at_1e7.entry == fock_level_name(1e-7, k), "a node saved at 1e-7 is read at 1e-7");

  const auto other_k = choose_saved_fock(node_file(1e-4, 6), Protocol{1e-4, 6}, 1e-4, 8);
  check(other_k.entry.empty(), "a node saved at k 6 is not read at k 8");

  json two = node_file(1e-4, 6);
  two.update(node_file(1e-6, 8));
  check(choose_saved_fock(two, Protocol{1e-4, 6}, 1e-4, 6).entry.empty(),
        "a node file holding more than one level is not read");

  // A standalone file (no node level) is still matched by its entry names.
  check(choose_saved_fock(two, std::nullopt, 1e-6, 8).entry == fock_level_name(1e-6, 8),
        "a standalone file is read at a level it names");
  check(choose_saved_fock(two, std::nullopt, 1e-5, 7).entry.empty(),
        "a standalone file is not read at a level it does not name");

  std::printf("%s\n", failures == 0 ? "test_saved_fock passed" : "test_saved_fock FAILED");
  return failures == 0 ? 0 : 1;
}
