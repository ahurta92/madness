#include <madness/chem/ParameterManager.hpp>
#include <cmath>
#include <iostream>

int main() {
    int failures = 0;
    auto check = [&](bool ok, const char* why) {
        if (!ok) { ++failures; std::cerr << why << '\n'; }
    };
    OptimizationParameters p;
    for (double threshold : {1e-4, 1e-6, 1e-8}) {
        p.set_protocol_defaults(threshold);
        check(p.get_gtol() == 1e-4 && p.get_gradient_precision() == 1e-5,
              "optimizer gradient defaults are independent of electronic tolerance and protocol");
        check(p.get_etol() == std::max(1e-6, 2.0 * threshold) && p.get_xtol() == 1e-3 &&
              p.get_value_precision() == std::max(1e-6, threshold),
              "existing energy, geometry and value-precision defaults remain protocol based");
    }
    p.set_user_defined_value("gtol", 3e-6);
    p.set_user_defined_value("gradient_precision", 7e-7);
    p.set_user_defined_value("etol", 1e-9);
    p.set_user_defined_value("xtol", 2e-4);
    p.set_user_defined_value("value_precision", 2e-8);
    p.set_protocol_defaults(1e-3);
    check(p.get_gtol() == 3e-6 && p.get_gradient_precision() == 7e-7 && p.get_etol() == 1e-9 &&
          p.get_xtol() == 2e-4 && p.get_value_precision() == 2e-8,
          "explicit optimizer overrides retain precedence");
    return failures ? 1 : 0;
}
