#pragma once
// A seed is a starting guess for one response node that came from somewhere
// else (a Dalton calculation, a converter tool), never a state this program
// solved. It lives beside the nodes, in the node layout without a history:
//   <run>/seeds/response/<name>/description.json   kind "seed"
//   <run>/seeds/response/<name>/archive/           the functions
// <name> is the name of the node the seed is for, and the archive has the form
// of that node's archive (archive/response for a driven state,
// archive/root_<slot> for an excited-state set), so the node loaders read both.
// A seed has no metrics and is never judged against convergence targets; a
// solve never writes it. It belongs to a calculation when its settings, the
// ground state's geometry and model plus the state's own physical settings,
// equal the calculation's. The ground state's save time is not compared: a
// seed is only a guess and stays usable when the ground state is refined.
#include "shared_es_checkpoint.hpp"
#include "es_root_identity.hpp"
#include <madness/chem/SCFScientificSettings.h>

namespace molresponse_v3::response_seed {
namespace fs = std::filesystem;
namespace it = madness::iteration;
namespace cp = madness::checkpoint;
using json = nlohmann::json;

// The seed folder for node `name`, relative to the run folder.
inline fs::path seed_folder(const std::string& name) {
  return fs::path("seeds")/response_shared::response_node_folder(name,"seed").lexically_relative("states");
}

// What a seed for a node with these scientific settings must hold as its
// settings: the ground state's model (geometry included) and the state's
// physical settings, without the ground state's save time or numerical options.
inline json seed_settings(const json& node_settings) {
  json state=node_settings;
  for(const auto* key:{"ground_settings","ground_saved","exchange_tensor","n_roots"}) state.erase(key);
  return {{"ground_model",madness::scf_model_settings(node_settings.at("ground_settings"))},{"state",state}};
}

// A seed found on disk: its folder relative to the run (for messages and the
// lineage), its absolute folder, and its description.
struct Seed {
  std::string name;
  fs::path folder;
  it::NodeDescription description;
  [[nodiscard]] std::string shown() const {return seed_folder(name).generic_string();}
  // The lineage entry of a state that starts from this seed.
  [[nodiscard]] json guess() const {
    json g=description.lineage;
    g["seed"]=shown();
    return g;
  }
};

// The seed description in `folder`, read without changing anything there.
// Nothing when it is missing, unreadable, not a seed, or has no archive.
inline std::optional<it::NodeDescription> inspect_seed(const fs::path& folder) {
  try {
    if(!fs::is_directory(it::Node(folder).archive())) return std::nullopt;
    std::ifstream file(folder/"description.json");
    auto d=it::read_description(json::parse(file));
    if(d.kind!="seed") return std::nullopt;
    return d;
  } catch(const std::exception&) {return std::nullopt;}
}

enum class SeedMatch {missing,matches,differs};
inline SeedMatch match_seed(const std::optional<it::NodeDescription>& d,const json& settings) {
  if(!d) return SeedMatch::missing;
  return d->settings==settings ? SeedMatch::matches : SeedMatch::differs;
}

// Whether node `name` has a seed under `root`, and whether its settings match.
// Collective: rank 0 reads, every process gets the same answer.
inline SeedMatch seed_status(madness::World& world,const fs::path& root,const std::string& name,const json& settings) {
  const auto text=cp::detail::on_root(world,[&] {
    switch(match_seed(inspect_seed(fs::absolute(root)/seed_folder(name)),settings)) {
      case SeedMatch::matches: return std::string("matches");
      case SeedMatch::differs: return std::string("differs");
      default: return std::string("missing");
    }
  });
  return text=="matches" ? SeedMatch::matches : text=="differs" ? SeedMatch::differs : SeedMatch::missing;
}

// The seed of node `name` whatever its settings. Collective.
inline std::optional<Seed> read_seed(madness::World& world,const fs::path& root,const std::string& name) {
  const auto folder=fs::absolute(root)/seed_folder(name);
  const auto text=cp::detail::on_root(world,[&] {
    const auto d=inspect_seed(folder);
    return d ? it::description_record(*d).dump() : std::string{};
  });
  if(text.empty()) return std::nullopt;
  return Seed{name,folder,it::read_description(json::parse(text))};
}

// The seed for node `name` when its settings match `settings`. A seed made for
// another geometry or model is ignored, with one line that names its folder.
// Collective.
inline std::optional<Seed> find_seed(madness::World& world,const fs::path& root,const std::string& name,
                                     const json& settings,it::Verbosity verbosity) {
  auto seed=read_seed(world,root,name);
  if(!seed) return std::nullopt;
  if(seed->description.settings==settings) return seed;
  if(world.rank()==0 && it::prints(verbosity,it::needed::seed_refused)) it::write_seed_refused(std::cout,seed->shown());
  return std::nullopt;
}

// Write the seed of node `name`: writer(folder) saves the archives collectively
// into an empty staging folder, then rank 0's one Node::save replaces the seed.
// `first_shard` names the archive file whose header holds the writer count.
template<class Writer>
void save_seed(madness::World& world,const fs::path& root,const std::string& name,it::Protocol protocol,
               json settings,json lineage,json details,const std::string& first_shard,const Writer& writer) {
  const auto folder=fs::absolute(root)/seed_folder(name);
  const auto staging=folder/"archive.staging";
  cp::detail::on_root(world,[&] {
    fs::remove_all(staging);
    fs::create_directories(staging);
    return std::string{};
  });
  world.gop.fence();
  writer(staging);
  world.gop.fence();
  cp::detail::on_root(world,[&] {
    madness::archive::BinaryFstreamInputArchive header((staging/first_shard).string());
    int writers=0;
    header & writers;
    it::NodeDescription d;
    d.kind="seed";
    d.protocol=protocol;
    d.settings=std::move(settings);
    d.metrics=json::array();
    d.lineage=std::move(lineage);
    d.code_version=madness::info::git_commit();
    d.level_finished=true;
    details["writers"]=writers;
    d.details=std::move(details);
    it::Node(folder).save(d,[&](const fs::path& archive) {
      for(const auto& file:fs::directory_iterator(staging)) fs::rename(file.path(),archive/file.path().filename());
    });
    fs::remove_all(staging);
    return std::string{};
  });
}

// A driven-state seed: one native archive "response", as a driven node holds.
template<class Storage>
void save_driven_seed(madness::World& world,const fs::path& root,const std::string& name,it::Protocol protocol,
                      const json& settings,const json& lineage,json details,const Storage& functions) {
  details["seed_for"]="driven";
  details["driven_node_version"]=1;
  save_seed(world,root,name,protocol,settings,lineage,std::move(details),"response.00000",[&](const fs::path& folder) {
    functions.save(world,madness::native_archive_name((folder/"response").string()),IoBackend::Native);
  });
}

// An excited-state seed: archive/root_<slot> per root and the roots' energies
// in its details, as an excited-state set node holds them. Roots keep their
// slot as identity.
template<class State>
void save_excited_seed(madness::World& world,const fs::path& root,const std::string& name,it::Protocol protocol,
                       const json& settings,const json& lineage,json details,const State& state) {
  json roots=json::array();
  for(std::size_t slot=0;slot<state.roots.size();++slot)
    roots.push_back({{"slot",slot},{"root_id",make_root_id(static_cast<int>(slot))},
                     {"stable_index",static_cast<int>(slot)},{"omega",state.omega(static_cast<long>(slot))}});
  details["seed_for"]="excited_set";
  details["excited_node_version"]=1;
  details["n_roots"]=roots.size();
  details["roots"]=roots;
  save_seed(world,root,name,protocol,settings,lineage,std::move(details),"root_0.00000",[&](const fs::path& folder) {
    response_shared::save_excited_functions(world,state,folder);
  });
}

// The functions of a driven-state seed at the active level: read at the seed's
// level, then projected to the active (k, thresh), as an imported guess always
// was. Collective.
template<class Storage>
Storage load_driven_seed(madness::World& world,const Seed& seed) {
  const auto& d=seed.description;
  cp::records::require(d.details.value("seed_for",std::string())=="driven" &&
      d.details.at("writers").get<int>()<=world.size(),"unsupported driven-state seed: "+seed.folder.string());
  const it::Protocol target{madness::FunctionDefaults<3>::get_thresh(),madness::FunctionDefaults<3>::get_k()};
  Storage functions;
  {
    struct Restore {
      it::Protocol target;
      ~Restore() {madness::FunctionDefaults<3>::set_k(target.k);madness::FunctionDefaults<3>::set_thresh(target.thresh);}
    } restore{target};
    madness::FunctionDefaults<3>::set_k(d.protocol.k);
    madness::FunctionDefaults<3>::set_thresh(d.protocol.thresh);
    functions=Storage::load(world,madness::native_archive_name((it::Node(seed.folder).archive()/"response").string()),
                            IoBackend::Native);
  }
  for(auto* block:functions.blocks())
    for(auto& function:*block) function=madness::project(function,target.k,target.thresh);
  world.gop.fence();
  return functions;
}

// An excited-state seed as a solver state at the active level, through the
// excited-set node loader. The residuals are unknown; the first update
// measures them. Collective.
template<class Type,class Shell>
typename ESSolver<Type,Shell>::State load_excited_seed(madness::World& world,const Seed& seed) {
  const auto& d=seed.description;
  cp::records::require(d.details.value("seed_for",std::string())=="excited_set",
      "unsupported excited-state seed: "+seed.folder.string());
  auto state=response_shared::load_excited_functions<Type,Shell>(world,it::Node(seed.folder).archive(),d,
      {madness::FunctionDefaults<3>::get_thresh(),madness::FunctionDefaults<3>::get_k()});
  state.last_bsh_residual.assign(state.roots.size(),1.0);
  state.last_density_residual.assign(state.roots.size(),1.0);
  state.iter=0;
  return state;
}
} // namespace molresponse_v3::response_seed
