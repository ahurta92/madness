#include "../solvers/shared_es_checkpoint.hpp"
#include <madness/chem/MetricHelpers.h>
#include <iostream>

using namespace madness;
using namespace molresponse_v3;
namespace it=madness::iteration;
double gaussian(const coord_3d& r) {return std::exp(-inner(r,r));}
double first(const coord_3d& r) {return r[0]*gaussian(r);}
double second(const coord_3d& r) {return r[2]*std::exp(-0.7*inner(r,r));}

template<class Type>
void verify_snapshot(World& world,const ResponseGroundState& gs,
                     const typename ESSolver<Type,ClosedShell>::State& state,
                     const it::IterationMetrics& metrics) {
  for(std::size_t i=0;i<state.roots.size();++i) {
    const auto& row=metrics.states[i];
    const auto& residual=row.metrics.at(it::Metric::bsh_residual);
    double maximum=0.,scale=0.;
    for(const auto& component:residual.evidence.components) {
      maximum=std::max(maximum,component.residual_norm);
      scale=std::max(scale,*component.function_norm);
    }
    MADNESS_CHECK(*residual.value==maximum && *residual.scale.value==scale);
    const auto& functions=state.roots[i];
    auto density=dot(world,gs.amo,functions.x_alpha);
    if constexpr(std::is_same_v<Type,Full>) density+=dot(world,gs.amo,functions.y_alpha);
    density.scale(2.);density.truncate();
    MADNESS_CHECK(std::abs(*row.metrics.at(it::Metric::density_change).evidence.current_norm-density.norm2())<1e-10);
    double kinetic=0.;
    for(const auto& [block,values]:measure_function_kinetics(world,functions))
      for(double value:values) kinetic+=value;
    MADNESS_CHECK(std::abs(*row.metrics.at(it::Metric::kinetic_change).evidence.current-kinetic)<1e-10);
    MADNESS_CHECK(*row.metrics.at(it::Metric::energy_change).evidence.current==state.omega(static_cast<long>(i)));
  }
}

template<class Type>
void exercise(World& world,bool stream,bool batch) {
  using Solver=ESSolver<Type,ClosedShell>;
  ResponseGroundState gs;
  real_function_3d phi=real_factory_3d(world).f(gaussian);
  phi.scale(1./phi.norm2()); gs.amo={phi};
  gs.aeps=Tensor<double>(1);gs.aeps(0L)=-1.;
  gs.focka=Tensor<double>(1,1);gs.focka(0,0)=-1.;
  gs.focka_no_diag=Tensor<double>(1,1);
  gs.Qa=QProjector<double,3>(gs.amo);
  gs.V_local_alpha=real_factory_3d(world);
  gs.c_xc=0.;gs.coulop=poperatorT(CoulombOperatorPtr(world,1e-6,1e-4));
  ConvergencePolicy policy; policy.stream_theta=stream;policy.kain=false;
  Solver solver(world,{gs,2},policy,PrintLevel::Silent);
  solver.set_batched(batch);
  typename Solver::State state;
  for (auto function:{first,second}) {
    typename Solver::Storage root;
    root.x_alpha={real_factory_3d(world).f(function)};
    if constexpr (std::is_same_v<Type,Full>) root.y_alpha={real_factory_3d(world)};
    state.roots.push_back(std::move(root));
  }
  {
    auto trials=state;
    trials.roots.clear();
    for(const auto& root:state.roots) trials.roots.push_back(root.copy(world));
    const auto prepared=solver.prepare_initial_guess(std::move(trials));
    MADNESS_CHECK(prepared.iter==0 && prepared.omega.size()==2);
    MADNESS_CHECK(prepared.omega(0L)<=prepared.omega(1L));
    bool rejected=false;
    try {solver.advance(state);} catch(...) {rejected=true;}
    MADNESS_CHECK(rejected); // preparation never creates an update proposal
  }
  it::IterationInfo info{"excited",{1e-4,6},0,"model"};
  it::Request request;request.max_updates=0;
  request.targets[it::Metric::bsh_residual]=1e6;
  it::ModelMetrics model{{it::Metric::bsh_residual,{it::Support::supported,false,{},it::LimitRule::plain,{}}}};
  namespace cp=madness::checkpoint;
  namespace shared=molresponse_v3::response_shared;
  const auto directory=cp::detail::on_root(world,[] {
    char path[]="/tmp/shared-es-native-XXXXXX";
    const auto made=mkdtemp(path);
    if(!made) throw std::runtime_error("cannot create ES native scratch");
    return std::string(made);
  });
  int observations=0;
  auto observer=[&](const auto& saved,const auto& metrics,const auto& decision,const auto& progress) {
    ++observations;
    MADNESS_CHECK(saved.iter==0 && metrics.states.size()==2 && decision.action==it::Action::complete);
    // The roots as ExcitedSet::publish describes them in a node.
    shared::save_excited_functions(world,saved,directory);
    world.gop.fence();
    cp::json roots=cp::json::array();
    for(std::size_t slot=0;slot<saved.roots.size();++slot)
      roots.push_back({{"slot",slot},{"root_id",molresponse_v3::make_root_id(saved.stable_index[slot])},
          {"stable_index",saved.stable_index[slot]},{"omega",saved.omega(static_cast<long>(slot))}});
    it::NodeDescription description;
    description.kind="excited_set";description.protocol=info.protocol;
    description.details={{"n_roots",roots.size()},{"roots",roots},{"writers",world.size()}};
    const auto loaded=shared::load_excited_functions<Type,ClosedShell>(world,directory,description,info.protocol);
    MADNESS_CHECK(loaded.stable_index==saved.stable_index && loaded.iter==0);
    for(std::size_t i=0;i<loaded.roots.size();++i) {
      MADNESS_CHECK(loaded.omega(static_cast<long>(i))==saved.omega(static_cast<long>(i)));
      MADNESS_CHECK(norm2(world,sub(world,loaded.roots[i].flatten(),saved.roots[i].flatten()))<1e-12);
    }
    verify_snapshot<Type>(world,gs,loaded,metrics);
    for (const auto& root:metrics.states) {
      MADNESS_CHECK(root.metrics.at(it::Metric::density_change).status==it::Availability::unavailable);
      MADNESS_CHECK(root.metrics.at(it::Metric::energy_change).status==it::Availability::unavailable);
      MADNESS_CHECK(root.metrics.at(it::Metric::kinetic_change).status==it::Availability::unavailable);
    }
  };
  solvers::iterate_shared(solver,state,it::resolve(request,model,info.protocol, true),true,
      {info.task_id,info.context_id,{"es_root_0000","es_root_0001"}},info,0,observer);
  MADNESS_CHECK(observations==1 && state.iter==0 && state.omega.size()==2);
  // The evaluated proposal is consumed once, with no second equation solve.
  solver.advance(state);
  MADNESS_CHECK(state.iter==1);
  bool rejected=false;
  try {solver.advance(state);} catch (...) {rejected=true;}
  MADNESS_CHECK(rejected);
  info.iteration=1;
  const auto measured=solver.evaluate(state,info);
  MADNESS_CHECK(measured.states.size()==2 && state.iter==1);
  verify_snapshot<Type>(world,gs,state,measured);
  for(const auto& root:measured.states) {
    for(auto id:{it::Metric::energy_change,it::Metric::kinetic_change}) {
      const auto& metric=root.metrics.at(id);
      if(root.comparison) {
        MADNESS_CHECK(metric.status==it::Availability::measured);
        MADNESS_CHECK(*metric.value==std::abs(*metric.evidence.current-*metric.evidence.previous));
        MADNESS_CHECK(*metric.scale.value==std::max(std::abs(*metric.evidence.current),std::abs(*metric.evidence.previous)));
      } else MADNESS_CHECK(metric.status==it::Availability::unavailable);
    }
  }
  solver.reset_shared();
  rejected=false;
  try {solver.advance(state);} catch (...) {rejected=true;}
  MADNESS_CHECK(rejected);
  cp::detail::on_root(world,[&] {std::filesystem::remove_all(directory);return std::string{};});
}
int main(int argc,char** argv) {
  initialize(argc,argv);int status=0;
  { World world(SafeMPI::COMM_WORLD);
  try {
    startup(world,argc,argv,true);
    FunctionDefaults<3>::set_cubic_cell(-8.,8.);
    FunctionDefaults<3>::set_k(6);FunctionDefaults<3>::set_thresh(1e-4);
    exercise<TDA>(world,false,false);exercise<TDA>(world,false,true);
    exercise<TDA>(world,true,false);exercise<Full>(world,false,false);
    exercise<Full>(world,true,false);
    if(world.rank()==0)std::cout<<"shared ES preparation/update interface passed\n";
  } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';status=1;}
  world.gop.max(status);world.gop.fence();
  }
  finalize();return status;
}
