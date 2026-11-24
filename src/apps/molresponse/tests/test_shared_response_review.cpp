#include "../orchestrator/response_workflow.hpp"
#include <iostream>
using namespace madness;
using namespace molresponse_v3;
namespace cp=madness::checkpoint;
namespace shared_review=molresponse_v3::response_shared;
using json=nlohmann::json;
void require(bool value) {if(!value) throw std::runtime_error("response review regression failed");}
struct NoSolve : ICalcExecutor {
  NodeResult run_protocol(const WorkItem&) override {throw std::runtime_error("unexpected numerical solve");}
};
void frequency_collision(World& world,const std::string& directory,bool explicit_leg,bool identical) {
  const auto key=protocol_key_at(1e-4);
  ResponsePlan plan;
  plan.es.push_back({true,explicit_leg ? 1 : 2,{1e-4}});
  plan.derived_fd.push_back({Perturbation::dipole(2),"*",{1e-4},0.5});
  if(explicit_leg) plan.fd.push_back({Perturbation::dipole(2),0.1,{1e-4}});
  const double second=identical ? 0.2 : 0.20000002;
  cp::detail::on_root(world,[&] {
    json roots=json::array({{{"root_id","es_root_0000"},{"omega",explicit_leg ? second : 0.2}}});
    if(!explicit_leg) roots.push_back({{"root_id","es_root_0001"},{"omega",second}});
    json metadata{{"schema_version",1},{"protocols",{{key,{{"thresh",1e-4},{"k",6},{"index",0}}}}},
      {"excited_states",{{key,{{"type","tda"},{"n_roots",explicit_leg ? 1 : 2},{"converged",true},{"roots",roots}}}}}};
    std::filesystem::create_directories(directory);
    std::ofstream(directory+"/response_metadata.json") << metadata.dump(); return std::string{};
  });
  CalcManager::Policy policy;
  policy.dependencies_ready=[](const CalcNode&){return false;};
  // The session plans every state; scheduling refuses a state without a plan.
  policy.state_plans[es_node_id(true,explicit_leg ? 1 : 2)]={StatePlan::Kind::load,1e-4,key,true};
  if(explicit_leg) policy.state_plans[fd_node_id(Perturbation::dipole(2),0.1)]={StatePlan::Kind::fresh,1e-4,""};
  policy.plan_derived=[](const CalcNode&){return StatePlan{StatePlan::Kind::fresh,1e-4,""};};
  CalcManager manager(plan,directory,policy);manager.build(0);NoSolve executor;
  bool collision=false;
  try {manager.run(world,executor);}
  catch(const std::exception& error) {
    collision=std::string(error.what()).find("distinct requested frequencies")!=std::string::npos;
    if(!collision) throw;
  }
  require(collision==!identical);
  if(identical) {
    int concrete=0;
    for(const auto& node:manager.dag()) if(node.kind==CalcKind::FD) ++concrete;
    require(concrete==1);
  }
}
void compatibility_output() {
  const std::string protocol="1e-04_k6";
  ResponseWorkflowOutput output;
  output.metadata={{"properties",json::object()},{"vibrational_analysis",{{"old",true}}},{"raman_spectra",{{"old",true}}}};
  json current=json::object();
  for(const auto& [shared,legacy]:std::vector<std::pair<std::string,std::string>>{
      {"polarizability","alpha"},{"hyperpolarizability","beta"},
      {"polarizability_gradient","raman"},{"two_photon_absorption","tpa"}}) {
    output.metadata["properties"][legacy]={{protocol,json::array({{{"current",true}}})},{"1e-03_k4",json::array({{{"old",true}}})}};
    current[shared]={{"values",json::array({{{"components",{{"zz",2.0}}}}})}};
  }
  detail_workflow::export_shared_properties(output,current,protocol);
  require(output.properties.size()==4);
  for(const auto& property:output.properties) require(property.size()==1 && property.contains(protocol));
  require(output.vibrational_analysis.is_null() && output.raman_spectra.is_null());
  // A later alpha-only request must not inherit older Raman or 2PA results.
  current.erase("hyperpolarizability");current.erase("polarizability_gradient");current.erase("two_photon_absorption");
  detail_workflow::export_shared_properties(output,current,protocol);
  require(output.properties.size()==1 && output.properties.contains("alpha"));
  current["polarizability"]["values"][0]["components"]=json::object();
  detail_workflow::export_shared_properties(output,current,protocol);
  require(output.properties.empty());
}
int main(int argc,char** argv) {
  initialize(argc,argv);int result=0;
  {
    World world(SafeMPI::COMM_WORLD);
    try {
      auto path=cp::detail::on_root(world,[] {
        char scratch[]="/tmp/shared-response-review-XXXXXX";
        auto* made=mkdtemp(scratch);if(!made) throw std::runtime_error("scratch creation failed");
        return std::string(made);
      });
      compatibility_output();
      for(bool explicit_leg:{false,true}) for(bool identical:{false,true})
        frequency_collision(world,path,explicit_leg,identical);
      cp::detail::on_root(world,[&]{std::filesystem::remove_all(path);return std::string{};});
      if(world.rank()==0) std::cout << "response review regressions passed\n";
    } catch(const std::exception& error) {if(world.rank()==0) std::cerr << error.what() << '\n';result=1;}
  }
  finalize();return result;
}
