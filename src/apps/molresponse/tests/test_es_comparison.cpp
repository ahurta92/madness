#include "../solvers/es_comparison.hpp"
#include <cmath>
#include <stdexcept>
using namespace molresponse_v3;
void require(bool ok) {if(!ok) throw std::runtime_error("excited-state comparison correspondence failed");}
int main() {
    const std::vector<double> norms{1,1};
    auto matches=match_prepared_roots({{0,-1},{1,0}},norms,norms);
    require(matches[0].previous==1 && matches[0].phase==-1);
    require(matches[1].previous==0 && matches[1].phase==1);
    const double c=std::sqrt(0.5);
    matches=match_prepared_roots({{c,c},{-c,c}},norms,norms);
    require(!matches[0].previous && !matches[1].previous);
    matches=match_prepared_roots({{0.98,0.02},{0.02,0.98}},norms,norms);
    require(matches[0].previous==0 && matches[1].previous==1);
    matches=match_prepared_roots({{0.7,0},{0,0.7}},norms,norms);
    require(!matches[0].previous && !matches[1].previous);
    matches=match_prepared_roots({{0,0},{0,1}},{0,1},norms);
    require(!matches[0].previous && matches[1].previous==1);
    // Competing roots cannot both claim the same previous root.
    matches=match_prepared_roots({{0.99,0.01},{0.98,0.02}},norms,norms);
    require(!matches[0].previous && !matches[1].previous);
    // Different raw amplitudes use normalized positive overlap, with no
    // RPA signed normalization or scalar-history rotation in this helper.
    matches=match_prepared_roots({{6,0},{0,-20}},{2,4},{3,5});
    require(matches[0].previous==0 && matches[1].previous==1 && matches[1].phase==-1);
}
