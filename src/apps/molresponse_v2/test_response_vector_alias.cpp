#include "ResponseState.hpp"
#include "ResponseVector.hpp"

#include <cassert>
#include <iostream>
#include <utility>

namespace {

template <typename ResponseType>
void check_x_only_alias_storage(const char *label, size_t initial_size) {
  ResponseType response(initial_size);
  assert(&response.x_alpha == &response.flat);
  assert(&response_x(response) == &response_all(response));
  assert(response.num_orbitals() == initial_size);

  response_x(response).resize(initial_size + 2);
  assert(response_all(response).size() == initial_size + 2);
  response.sync();
  response.flatten();
  assert(response_x(response).size() == initial_size + 2);
  assert(response.num_orbitals() == initial_size + 2);

  auto copied = response;
  assert(&copied.x_alpha == &copied.flat);
  assert(&response_x(copied) == &response_all(copied));
  assert(copied.num_orbitals() == response.num_orbitals());

  response_all(copied).resize(initial_size + 4);
  assert(response_x(copied).size() == initial_size + 4);
  assert(response.num_orbitals() == initial_size + 2);

  auto moved = std::move(copied);
  assert(&moved.x_alpha == &moved.flat);
  assert(&response_x(moved) == &response_all(moved));
  assert(moved.num_orbitals() == initial_size + 4);

  std::cout << "alias-storage-ok " << label << " size="
            << moved.num_orbitals() << "\n";
}

} // namespace

int main() {
  check_x_only_alias_storage<StaticRestrictedResponse>(
      "StaticRestrictedResponse", 3);
  check_x_only_alias_storage<TDARestrictedResponse>(
      "TDARestrictedResponse", 4);
  return 0;
}
