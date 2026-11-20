#include "../solvers/shared_response_checkpoint.hpp"
#include <madness/chem/LadderPlanner.h>
#include "../solvers/shared_response_input.hpp"
#include <stdexcept>
#include <source_location>

namespace it=madness::iteration;
void require(bool value,std::source_location loc=std::source_location::current()) {if(!value) throw std::runtime_error("FD protocol plan failed at line "+std::to_string(loc.line()));}
int main() {
  ResponseParameters parameters;
  const auto defaults=molresponse_v3::shared_response_request(parameters);
  require(defaults.has_value());
  const auto default_controls=it::resolve(*defaults,molresponse_v3::response_shared::model(),{1e-4,6},true);
  require(default_controls.criteria.size()==3);
  require(default_controls.criteria[0].absolute==5e-4);
  require(default_controls.criteria[1].absolute==5e-4);
  require(ResponseParameters::retired_input("dconv"));
  require(ResponseParameters::retired_input("excited.guess_max_iter"));
  it::Request request;
  request.criteria={it::Metric::bsh_residual,it::Metric::density_change,it::Metric::kinetic_change};
  for(auto metric:request.criteria) request.tolerances[metric]={1e-8,0.0};
  const std::vector<it::Protocol> ladder{{1e-4,6},{1e-6,8}};
  const auto model=molresponse_v3::response_shared::model(ladder.front());
  auto controls=it::resolve(request,model,ladder.front(),false);
  require(controls.criteria[0].absolute==5e-4 && controls.criteria[1].absolute==5e-4);
  require(controls.criteria[2].absolute==1e-8);
  controls=it::resolve(request,model,ladder.back(),true);
  require(controls.criteria[0].absolute==1e-8);
  it::SavedTask saved;
  saved.metrics.info={"fd","previous",ladder.front(),1,"same_ground"};
  it::StateMetrics metrics;metrics.id="response";
  for(auto metric:request.criteria) metrics.metrics[metric]=it::MetricValue::measured(0.0,1.0);
  saved.metrics.states={metrics};saved.conversion_supported=true;
  const it::ExpectedTask expected{"fd","same_ground",{"response"}};
  auto plan=it::plan_restart(ladder,request,model,expected,saved);
  require(plan.start_protocol==0 && plan.action==it::RestartAction::resume);
  saved.protocol_complete=true;
  plan=it::plan_restart(ladder,request,model,expected,saved);
  require(plan.start_protocol==1 && plan.action==it::RestartAction::resume);
  saved.metrics.info.protocol=ladder.back();
  require(it::plan_restart(ladder,request,model,expected,saved).action==it::RestartAction::resume);
  saved.metrics.states[0].comparison=it::Comparison{"fd","previous",ladder.back(),0,"response"};
  plan=it::plan_restart(ladder,request,model,expected,saved);
  require(plan.start_protocol==1 && plan.action==it::RestartAction::reuse);
  saved.metrics.states[0].metrics[it::Metric::bsh_residual]=it::MetricValue::measured(1.0,1.0);
  plan=it::plan_restart(ladder,request,model,expected,saved);
  require(plan.start_protocol==1 && plan.action==it::RestartAction::resume);
}
