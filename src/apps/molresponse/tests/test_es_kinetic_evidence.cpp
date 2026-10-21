#include "../solvers/es_kinetic_evidence.hpp"
#include <iostream>

int main() {
  using nlohmann::json;
  const auto check = [](bool ok) { if (!ok) throw std::runtime_error("ES kinetic evidence test failed"); };
  const json blocks = {{"x_alpha", json::array({{{"function_index", 0}, {"value", 3.0}, {"availability", "measured"}}})},
                       {"y_alpha", json::array({{{"function_index", 0}, {"value", 1.0}, {"availability", "measured"}}})}};
  auto report = molresponse_v3::es_kinetic_evidence(blocks, 4, 1, 2, 7, "full", "closed_shell", 1e-6, 8, true);
  check(report.at("kinetic_energy").at("trace").at("value") == 4.0); // positive X+Y, never signed RPA X-Y
  check(report.at("delta_kinetic_energy").at("value").is_null());
  check(report.at("delta_kinetic_energy").at("reason") == "state_correspondence_unvalidated");
  check(report.at("identity").at("bookkeeping_index") == 7);
  check(report.at("kinetic_energy").at("sampled_this_iteration") == true);
  check(report.at("root_locked") == true);
  check(molresponse_v3::valid_es_kinetic_evidence(report));
  auto bad = report;
  bad["identity"]["slot"] = 2;
  check(!molresponse_v3::valid_es_kinetic_evidence(bad));
  bad = report;
  bad["identity"]["model"] = "tda"; // Y blocks cannot masquerade as TDA
  check(!molresponse_v3::valid_es_kinetic_evidence(bad));
  bad = report;
  bad["identity"]["shell"] = "open_shell"; // missing beta blocks
  check(!molresponse_v3::valid_es_kinetic_evidence(bad));
  report["kinetic_energy"]["trace"]["value"] = 2.0;
  check(!molresponse_v3::valid_es_kinetic_evidence(report));
  std::cout << "ES kinetic evidence tests passed\n";
}
