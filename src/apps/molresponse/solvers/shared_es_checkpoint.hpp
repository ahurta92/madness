#pragma once
#include "shared_response_checkpoint.hpp"
#include "es_solver.hpp"

namespace molresponse_v3::response_shared {
inline std::map<it::Metric,std::string> excited_metric_set() {
  return {{it::Metric::bsh_residual,"es-bsh-v1"},{it::Metric::density_change,"es-density-v1"},
          {it::Metric::energy_change,"es-energy-v1"},{it::Metric::kinetic_change,"es-kinetic-v1"}};
}
inline std::map<it::Metric,std::string> excited_metric_labels() {
  auto labels=metric_labels();
  labels[it::Metric::energy_change]="Excitation-energy change (Ha)";
  return labels;
}
inline TaskDefinition excited_task(std::size_t roots) {
  auto ids=excited_metric_set();
  json definitions_es=json::object();
  const auto fd=definitions();
  for(const auto& [metric,id]:ids) {
    auto definition=fd.at(metric_set().at(metric));
    if(metric==it::Metric::energy_change) {
      definition["meaning"]="absolute change of the prepared root excitation energy";
      definition["units"]="hartree";definition["fields"]="prepared root Ritz energy";
    }
    if(metric!=it::Metric::bsh_residual)
      definition["comparison"]="preceding prepared root with unambiguous positive-overlap match; density phase aligned";
    definitions_es[id]=std::move(definition);
  }
  json subjects=json::array();
  for(std::size_t root=0;root<roots;++root) subjects.push_back(make_root_id(static_cast<int>(root)));
  return {{{"subjects",subjects},{"es_checkpoint_version",1}},ids,definitions_es,
      [](it::Protocol protocol) {
        // Excited states also judge the change in excitation energy, and keep
        // the plain rule for every metric.
        auto metrics=model(protocol);
        for(auto& [metric,definition]:metrics) definition.rule=it::LimitRule::plain;
        metrics[it::Metric::energy_change]={it::Support::supported,true,5e-4,it::LimitRule::plain,{}};
        return metrics;
      },
      [](const json& previous,const json& current) {
        cp::records::require(previous.at("layout")==current.at("layout"),
            "unsupported excited checkpoint block layout");
      },true,"excited states",excited_metric_labels()};
}

// Use the existing per-root native serializer. The shared manifest already
// carries root identity, omega and measured metrics, so no second sidecar can
// disagree with the archive used for convergence/restart decisions.
template<class State>
cp::Payload save_excited_functions(madness::World& world,const State& state,
                                  const std::filesystem::path& directory) {
  cp::Payload payload;
  cp::records::require(state.roots.size()==state.stable_index.size(),"excited bundle identity count differs");
  for(std::size_t slot=0;slot<state.roots.size();++slot) {
    const auto basename="root_"+std::to_string(slot);
    state.roots[slot].save(world,madness::native_archive_name((directory/basename).string()),IoBackend::Native);
  }
  // Finish every collective save before rank-zero inspection can throw.
  if(world.rank()==0) for(std::size_t slot=0;slot<state.roots.size();++slot) {
    const auto basename="root_"+std::to_string(slot);
    payload.archives.push_back(cp::inspect_archive(directory,"root",basename,{state.stable_index[slot]},1));
  }
  return payload;
}

template<class Type,class Shell>
typename ESSolver<Type,Shell>::State load_excited_functions(
    madness::World& world,const std::filesystem::path& directory,
    const cp::Manifest& manifest,it::Protocol target) {
  using Solver=ESSolver<Type,Shell>;
  const auto frame=cp::iteration_metrics_from_manifest(manifest);
  const auto stored=frame.info.protocol;
  cp::records::require(!frame.states.empty() && manifest.archives.size()==frame.states.size(),
      "incomplete excited native bundle");
  cp::records::require(it::Protocol{madness::FunctionDefaults<3>::get_thresh(),madness::FunctionDefaults<3>::get_k()}==target,
      "prepare excited representation before loading");
  typename Solver::State result;
  result.omega=madness::Tensor<double>(static_cast<long>(frame.states.size()));
  std::set<int> identities;
  for(std::size_t slot=0;slot<manifest.archives.size();++slot) {
    const auto& archive=manifest.archives[slot];
    cp::records::require(archive.role=="root" && archive.basename=="root_"+std::to_string(slot) &&
        archive.indices.size()==1 && archive.indices.front()>=0 &&
        archive.writers<=world.size() && identities.insert(archive.indices.front()).second,
        "unsupported excited native payload");
    const auto id=make_root_id(archive.indices.front());
    const auto metrics=std::find_if(frame.states.begin(),frame.states.end(),[&](const auto& root){return root.id==id;});
    cp::records::require(metrics!=frame.states.end(),"excited archive has no matching metrics");
    const auto omega=metrics->metrics.at(it::Metric::energy_change).evidence.current;
    cp::records::require(omega && std::isfinite(*omega),"excited checkpoint has no finite root energy");
    typename Solver::Storage functions;
    {
      struct Restore {
        it::Protocol target;
        ~Restore(){madness::FunctionDefaults<3>::set_k(target.k);madness::FunctionDefaults<3>::set_thresh(target.thresh);}
      } restore{target};
      madness::FunctionDefaults<3>::set_k(stored.k);madness::FunctionDefaults<3>::set_thresh(stored.thresh);
      functions=Solver::Storage::load(world,madness::native_archive_name((directory/archive.basename).string()),IoBackend::Native);
    }
    for(auto* block:functions.blocks()) {
      for(auto& function:*block) {
        cp::records::require(function.k()==stored.k,"excited archive polynomial order differs");
        if(stored.k!=target.k) function=madness::project(function,target.k,target.thresh,false);
      }
      madness::set_thresh(world,*block,target.thresh);
      if(stored!=target) madness::truncate(world,*block,target.thresh);
    }
    result.roots.push_back(std::move(functions));result.stable_index.push_back(archive.indices.front());
    result.omega(static_cast<long>(slot))=*omega;
  }
  world.gop.fence();
  return result;
}
} // namespace molresponse_v3::response_shared
