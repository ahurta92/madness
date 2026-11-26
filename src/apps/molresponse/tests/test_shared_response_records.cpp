#include "../solvers/shared_response_session.hpp"
#include "../solvers/raman_results.hpp"
#include <madness/world/MADworld.h>
#include <iostream>
namespace rs = molresponse_v3::response_shared;
namespace cp = madness::checkpoint;
namespace it = madness::iteration;
using namespace madness;
void require(bool ok) { if (!ok) throw std::runtime_error("shared response record contract failed"); }
template<class F> void rejects(F&& f) {
    bool threw = false;
    try { f(); } catch (const std::exception&) { threw = true; }
    require(threw);
}
// Collective: one native archive "response" holding 123.
void payload(World& world, const std::filesystem::path& directory) {
    archive::ParallelOutputArchive out(world, (directory/"response").c_str(), 1);
    int value = 123;
    if (out.is_io_node()) out.local_archive() & value;
    cp::detail::collective_error(world, cp::close_archive(out));
    world.gop.fence();
}
cp::json read_root_json(World& world, const std::filesystem::path& path) {
    return cp::json::parse(cp::detail::on_root(world, [&] {return cp::detail::read_json(path).dump();}));
}
int main(int argc, char** argv) {
    initialize(argc, argv);
    int failed = 0;
    {
        World world(SafeMPI::COMM_WORLD);
        try {
            const auto root = cp::detail::on_root(world, [] {
                char path[] = "/tmp/shared-response-records-XXXXXX";
                const auto result = mkdtemp(path);
                if (!result) throw std::runtime_error("cannot make scratch directory");
                return std::string(result);
            });
            // The nearest saved frequency of the same ground state and model.
            {
                const cp::json g1{{"saved_at","T1"}},g2{{"saved_at","T2"}},m{{"perturbation","dipole_z"}};
                const auto c=[&](std::string name,double frequency,cp::json ground=cp::json()) {
                    return rs::StartCandidate{std::move(name),frequency,ground.is_null()?g1:ground,m};
                };
                // One candidate.
                require(rs::nearest_start(0.05,g1,m,{c("a",0.0)})==std::optional<std::string>("a"));
                require(!rs::nearest_start(0.05,g1,m,{}));
                // Several: the nearest, above or below.
                require(rs::nearest_start(0.05,g1,m,{c("a",0.0),c("b",0.04),c("c",0.2)})==std::optional<std::string>("b"));
                require(rs::nearest_start(0.03,g1,m,{c("a",0.0),c("c",0.04)})==std::optional<std::string>("c"));
                // A tie goes to the lower frequency, whatever the order.
                require(rs::nearest_start(0.04,g1,m,{c("hi",0.06),c("lo",0.02)})==std::optional<std::string>("lo"));
                require(rs::nearest_start(0.04,g1,m,{c("lo",0.02),c("hi",0.06)})==std::optional<std::string>("lo"));
                // Another ground-state identity, model or an unsettled save cannot serve.
                require(rs::nearest_start(0.05,g1,m,{c("a",0.0),c("b",0.04,g2)})==std::optional<std::string>("a"));
                require(!rs::nearest_start(0.05,g1,m,{c("b",0.04,g2)}));
                auto other_model=c("x",0.05);other_model.model={{"perturbation","dipole_x"}};
                auto unsettled=c("u",0.05);unsettled.settled=false;
                require(rs::nearest_start(0.05,g1,m,{c("a",0.0),other_model,unsettled})==std::optional<std::string>("a"));
            }
            // The two-photon states to remove: only derived nodes of the set
            // that no current root needs.
            {
                const auto derived=[](const std::string& bundle,const std::string& root) {
                    return cp::json{{"perturbation","dipole_x"},{"derived_from",{{"bundle",bundle},{"root",root},{"factor",0.5}}}};
                };
                const std::vector<std::pair<std::string,cp::json>> nodes{
                    {"dipole_x-f0.16139",derived("es:full","es_root_0000")},   // needed
                    {"dipole_x-f0.16205",derived("es:full","es_root_0000")},   // superseded
                    {"dipole_x-f0.20000",derived("es:tda","es_root_0000")},    // another set
                    {"dipole_x-f0.05000",{{"perturbation","dipole_x"},{"frequency",0.05}}}, // ordinary
                    {"dipole_x-f0.16000",{{"perturbation","dipole_x"},{"frequency",0.16},{"derived_from",nullptr}}}};
                const auto removed=rs::two_photon_removals("es:full",{"dipole_x-f0.16139"},nodes);
                require(removed==std::vector<std::string>{"dipole_x-f0.16205"});
                // Nothing needed: still only the set's derived nodes, never an ordinary state.
                require(rs::two_photon_removals("es:full",{},nodes)==
                        std::vector<std::string>{"dipole_x-f0.16139","dipole_x-f0.16205"});
            }
            const cp::json model{{"version",1},{"dft",{{"xc","HF"}}},
                {"molecule",{{"parameters",cp::json::object()},{"core_potentials",cp::json::array()}}},
                {"xc",{{"specification","HF"}}}};
            const rs::GroundSource ground{std::filesystem::path(root)/"ground.restartdata",{1e-4,6},
                {{"settings",model},{"protocol",{{"threshold",1e-4},{"k",6}}},{"saved_at","2026-11-23T00:00:00Z"}}};
            auto task = rs::driven_task();
            task.record = {{"subjects",{"root0","root1"}},{"binding_version",1}};
            task.validate_settings = [](const cp::json& previous,const cp::json& current) {
                cp::records::require(previous == current,"test model changed");
            };
            // These records carry only the residual; judge only it.
            task.model = [](it::Protocol) {
                return it::ModelMetrics{{it::Metric::bsh_residual,
                    {it::Support::supported,true,0.1,it::LimitRule::plain,{}}}};
            };
            it::Request request;
            request.targets[it::Metric::bsh_residual] = 0.1;
            request.max_updates = 0;
            const std::vector<it::Protocol> ladder{{1e-4,6}};
            // A driven node: saved, reused from its stored metrics, loaded.
            rs::DrivenState bundle(world,root,ground,"roots","roots",cp::json::object(),request,ladder,task);
            require(!bundle.has_saved() && bundle.folder()=="states/response/roots");
            it::IterationMetrics metrics; metrics.info=bundle.info();
            for (const auto* id : {"root0","root1"}) {
                it::StateMetrics state; state.id=id;
                state.metrics[it::Metric::bsh_residual]=it::MetricValue::measured(0.0,1.0);
                metrics.states.push_back(state);
            }
            const it::Progress progress{0,0,true};
            const auto decision=it::decide(bundle.controls(),bundle.expected(),progress,metrics,bundle.final_level());
            require(decision.action==it::Action::complete);
            bundle.publish(metrics,decision,progress,[&](const auto& dir){payload(world,dir);});
            require(bundle.has_node() && bundle.description().kind=="driven" && bundle.saved_metrics().states.size()==2);
            rs::DrivenState reused(world,root,ground,"roots","roots",cp::json::object(),request,ladder,task);
            require(reused.reusable() && reused.description()==bundle.description());
            const auto read_payload=[&](const auto& directory,const auto&) {
                archive::ParallelInputArchive in(world,(directory/"response").c_str());
                int value=0;
                if (in.is_io_node()) in.local_archive() & value;
                in.close();
                world.gop.broadcast_serializable(value,0);
                return value;
            };
            require(reused.load(read_payload,true)==123);
            require(reused.measured_values("root1").at("values").at("bsh_residual")==0.0);
            rejects([&]{reused.measured_values();});
            // A second-order node records its inputs' identities and is
            // reused only while its settings and every input are unchanged.
            rejects([&]{rs::second_order_folder("vbc:B");});
            rejects([&]{rs::second_order_folder("..");});
            const cp::json vbc_settings{{"kind","VBC"},{"protocol",cp::records::protocol_record(ladder.back())}};
            const auto lineage=cp::json::array({rs::lineage_entry("B",reused)});
            require(lineage[0].at("name")=="roots" &&
                    lineage[0].at("identity")==cp::json{{"saved_at",reused.description().saved_at}});
            const auto save_payload=[&](const auto& dir){payload(world,dir);};
            rs::SecondOrderNode built(world,root,"vbc:B","vbc-B",vbc_settings,lineage);
            require(!built.reusable());
            built.save({"response"},save_payload);
            require(built.folder()=="states/response/vbc-B" && built.description().kind=="second_order");
            rs::SecondOrderNode cached(world,root,"vbc:B","vbc-B",vbc_settings,lineage);
            require(cached.reusable());
            require(cached.load(read_payload)==123);
            require(cached.description()==built.description());
            auto other_settings=vbc_settings;other_settings["source"]="pq";
            rs::SecondOrderNode other(world,root,"vbc:B","vbc-B",other_settings,lineage);
            require(!other.reusable());
            rejects([&]{other.load(read_payload);});
            auto tighter=request; tighter.targets[it::Metric::bsh_residual]=0.0;
            rs::DrivenState tighter_bundle(world,root,ground,"roots","roots",cp::json::object(),tighter,ladder,task);
            require(tighter_bundle.reusable());
            require(rs::lineage_entry("B",tighter_bundle)==rs::lineage_entry("B",reused));
            metrics.states.back().metrics[it::Metric::bsh_residual]=it::MetricValue::measured(1.0,1.0);
            const auto capped=it::decide(bundle.controls(),bundle.expected(),progress,metrics,bundle.final_level());
            bundle.publish(metrics,capped,progress,[&](const auto& dir){payload(world,dir);});
            rejects([&]{(void)rs::lineage_entry("C",bundle);});
            // A changed input identity rebuilds the node over the old one.
            auto changed_lineage=lineage;
            changed_lineage[0]["identity"]["saved_at"]="2000-01-01T00:00:00Z";
            rs::SecondOrderNode changed_b(world,root,"vbc:B","vbc-B",vbc_settings,changed_lineage);
            require(!changed_b.reusable());
            rejects([&]{changed_b.load(read_payload);});
            changed_b.save({"response"},save_payload);
            require(rs::SecondOrderNode(world,root,"vbc:B","vbc-B",vbc_settings,changed_lineage).reusable());
            require(!rs::SecondOrderNode(world,root,"vbc:B","vbc-B",vbc_settings,lineage).reusable());
            cp::detail::on_root(world,[&] {
                const auto folder=std::filesystem::path(root)/built.folder();
                std::size_t entries=0;
                for(const auto& entry:std::filesystem::directory_iterator(folder)) {(void)entry;++entries;}
                require(entries==2 && std::filesystem::is_directory(folder/"archive"));
                return std::string{};
            });
            // An excited-state set is one node: saved whole, reused from its
            // stored metrics, resumed with added roots, refused with fewer.
            auto es_task=task;
            const auto es_record=[](int n) {
                cp::json subjects=cp::json::array();
                for(int r=0;r<n;++r) subjects.push_back(molresponse_v3::make_root_id(r));
                return cp::json{{"subjects",subjects}};
            };
            es_task.record=es_record(2);
            es_task.validate_settings=rs::excited_task(2).validate_settings;
            const cp::json es_settings{{"excited_model","tda"},{"layout","closed_shell_tda"}};
            rs::ExcitedSet es(world,root,ground,"es:tda",es_settings,request,ladder,es_task);
            require(!es.has_saved() && es.folder()=="states/response/excited-tda");
            it::IterationMetrics es_metrics; es_metrics.info=es.info();
            for(int r=0;r<2;++r) {
                it::StateMetrics state; state.id=molresponse_v3::make_root_id(r);
                state.metrics[it::Metric::bsh_residual]=it::MetricValue::measured(0.0,1.0);
                auto energy=it::MetricValue::not_applicable("first update");
                energy.evidence.current=0.5+r;
                state.metrics[it::Metric::energy_change]=energy;
                es_metrics.states.push_back(state);
            }
            const auto es_decision=it::decide(es.controls(),es.expected(),progress,es_metrics,es.final_level());
            es.publish(es_metrics,es_decision,progress,{0,1},[&](const auto& dir) {
                for(int r=0;r<2;++r) {
                    archive::ParallelOutputArchive out(world,(dir/("root_"+std::to_string(r))).c_str(),1);
                    if (out.is_io_node()) out.local_archive() & r;
                    require(cp::close_archive(out).empty());
                }
            });
            rs::ExcitedSet es_again(world,root,ground,"es:tda",es_settings,request,ladder,es_task);
            require(es_again.reusable() && es_again.description()==es.description());
            require(es_again.description().details.at("roots")[1].at("omega")==1.5);
            require(es_again.load([](const auto& dir,const auto&){return std::filesystem::exists(dir/"root_1.00000");}));
            auto three=es_task;three.record=es_record(3);
            rs::ExcitedSet es_more(world,root,ground,"es:tda",es_settings,request,ladder,three);
            require(es_more.has_saved() && !es_more.reusable() && es_more.plan().action==it::RestartAction::resume);
            auto one=es_task;one.record=es_record(1);
            rejects([&]{rs::ExcitedSet fewer(world,root,ground,"es:tda",es_settings,request,ladder,one);});
            auto moved_ground=ground;moved_ground.identity["saved_at"]="2026-11-24T00:00:00Z";
            rs::ExcitedSet es_guess(world,root,moved_ground,"es:tda",es_settings,request,ladder,es_task);
            require(es_guess.has_saved() && !es_guess.reusable());
            cp::detail::on_root(world,[&] {
                const auto folder=std::filesystem::path(root)/es.folder();
                std::size_t launches=0;std::string line;
                std::ifstream history(folder/"history.jsonl");
                while(std::getline(history,line)) launches+=cp::json::parse(line).contains("launch");
                require(launches==5 && std::filesystem::is_directory(folder/"archive"));
                return std::string{};
            });
            const auto derived_task=rs::derived_task();
            const cp::json origin{{"bundle","es:tda"},{"root","es_root_0000"},{"factor",0.5}};
            const cp::json old_settings{{"frequency",0.1},{"perturbation","dipole_z"},
                {"layout","closed_shell_full"},{"derived_from",origin}};
            rs::DrivenState old_leg(world,root,ground,"old-leg","dipole_z-f0.10000",old_settings,request,ladder,derived_task);
            it::IterationMetrics old_metrics;old_metrics.info=old_leg.info();old_metrics.info.iteration=1;
            it::StateMetrics leg;leg.id="response";
            leg.metrics[it::Metric::bsh_residual]=it::MetricValue::measured(0.0,1.0);
            // Driven response judges the density change by default, against a
            // comparable earlier iteration of the same task.
            const auto density=[](const it::IterationInfo& info) {
              auto value=it::MetricValue::measured(0.0,1.0);
              value.comparison=it::Comparison{info.task_id,info.protocol,0,"response"};
              return value;
            };
            leg.metrics[it::Metric::density_change]=density(old_leg.info());
            old_metrics.states={leg};
            // One update at the level: a stored comparison is the update before.
            const it::Progress leg_progress{1,1,true};
            const auto old_decision=it::decide(old_leg.controls(),old_leg.expected(),leg_progress,old_metrics,old_leg.final_level());
            old_leg.publish(old_metrics,old_decision,leg_progress,[&](const auto& dir){payload(world,dir);});
            const auto old_description=old_leg.description();
            // The moved root's state starts from the nearby frequency's node,
            // only reading it, and names it in its lineage.
            auto refined_settings=old_settings;refined_settings["frequency"]=0.11;
            rs::DrivenState refined(world,root,ground,"new-leg","dipole_z-f0.11000",refined_settings,request,ladder,
                                    derived_task,std::string("dipole_z-f0.10000"));
            require(refined.has_saved() && !refined.has_node() && !refined.reusable() && refined.plan().start_protocol==0);
            require(refined.load(read_payload)==123);
            rejects([&]{(void)rs::lineage_entry("B",refined);});
            auto new_metrics=old_metrics;new_metrics.info=refined.info();new_metrics.info.iteration=1;
            new_metrics.states[0].metrics[it::Metric::density_change]=density(refined.info());
            const auto new_decision=it::decide(refined.controls(),refined.expected(),leg_progress,new_metrics,refined.final_level());
            refined.publish(new_metrics,new_decision,leg_progress,[&](const auto& dir){payload(world,dir);});
            require(refined.has_node() && refined.description().lineage.at("started_from")=="dipole_z-f0.10000");
            require(rs::read_other_node(world,std::filesystem::path(root)/old_leg.folder())==old_description);
            // A later launch binds the derived state by its settings: its own
            // node is reused, and a further move starts from the nearest node.
            cp::detail::on_root(world,[&] {std::filesystem::create_directory(root+"/recovery");return std::string{};});
            it::TimingCollector<> recovery_timing;
            recovery_timing.begin("run"); recovery_timing.begin("response"); recovery_timing.begin("setup");
            rs::Session recovered(world,root,ground,request,ladder,{},root+"/recovery",std::move(recovery_timing),it::Verbosity::detail);
            auto& recovered_leg=recovered.bind_derived("new-leg",refined_settings);
            require(recovered_leg.reusable() && recovered_leg.description()==refined.description());
            auto moved_settings=refined_settings;moved_settings["frequency"]=0.12;
            auto& moved=recovered.bind_derived("moved-leg",moved_settings);
            require(moved.has_saved() && !moved.has_node() && moved.saved_metrics().states.size()==1 &&
                    moved.started_from()==std::optional<std::string>("dipole_z-f0.11000") && moved.previous_updates()==0);
            // A version 3 results file is replaced, not merged.
            const auto results_path=std::filesystem::path(root)/"results.json";
            cp::detail::on_root(world,[&] {
                cp::detail::replace_json(results_path,{{"version",3},{"calculations",cp::json::array({{{"calculation_id","old"}}})}});
                return std::string{};
            });
            // A property names the states its value used; the file carries
            // their worst measured values as its convergence instead.
            {
                cp::json single{{"units","bohr^3"},{"value_type","real"},{"value",1.0},{"status","complete"}};
                recovered.name_state(single,"new-leg");
                recovered.name_state(single,"not-bound"); // a task this request does not bind
                cp::json item{{"frequency",0.11},{"components",{{"zz",1.0}}},{"status","complete"},
                              {"availability",cp::json::object()}};
                recovered.name_state(item,"new-leg");
                recovered.set_properties({{"single",single},{"listed",{{"units","bohr^3"},{"values",cp::json::array({item})}}}});
            }
            recovered.report();
            // results.json is written directly: response's entry keys its
            // measured values by node folder name and names no launch; another
            // application's entry is left alone.
            {
                auto results=read_root_json(world,results_path);
                require(results.at("version")==4 && results.at("calculations").size()==1 && !results.contains("attempt_id"));
                const auto& provenance=results.at("provenance");
                require(provenance.at("code_version").is_string() && provenance.at("launched_at").is_string() &&
                        provenance.at("mpi_processes")==world.size() && provenance.at("archive_io")=="native");
                const auto& entry=results.at("calculations").at(0);
                require(entry.at("calculation_id")=="response" && !entry.contains("attempt_id") && entry.at("run").is_object());
                require(!entry.contains("driven") && !entry.contains("excited") && !entry.contains("second_order"));
                const auto& metrics=entry.at("state_metrics");
                require(metrics.contains(recovered_leg.name()) && !metrics.contains(moved.name()) &&
                        metrics.at(recovered_leg.name()).at("subjects").size()==1);
                require(recovered.node_name("new-leg")==recovered_leg.name());
                // Each value's convergence is its states' worst measured
                // values, to three significant figures; no states list and no
                // empty availability remain.
                const auto& measured=metrics.at(recovered_leg.name()).at("subjects").begin().value().at("values");
                const auto& properties=entry.at("properties");
                const auto& single=properties.at("single");
                const auto& listed=properties.at("listed").at("values").at(0);
                for(const auto* value:{&single,&listed}) {
                  require(!value->contains("sources") && !value->contains("states") && !value->contains("availability"));
                  require(value->at("convergence").size()==2);
                  for(const auto& metric:{"bsh_residual","density_change"})
                    require(value->at("convergence").at(metric)==it::round_significant(measured.at(metric).get<double>()));
                }
                require(single.at("value")==1.0 && listed.at("components").at("zz")==1.0);
                // Only driven states: one set of final-level targets.
                require(entry.at("targets").at("bsh_residual").is_number() &&
                        entry.at("targets").at("density_change").is_number() && !entry.at("targets").contains("energy_change"));
                cp::detail::on_root(world,[&] {
                    for(const auto& [name,value]:metrics.items())
                        require(std::filesystem::is_directory(std::filesystem::path(root)/"states"/"response"/name));
                    return std::string{};
                });
                // The ground entry takes its values from one state: one
                // convergence on the entry, none on its properties. It reports
                // every metric with a target (energy_change too, kinetic_change
                // not); targets are written to three significant figures.
                cp::json ground_entry{{"calculation_id","ground_state"},
                    {"properties",{{"total_energy",{{"value",-1.0},{"status","complete"}}}}},
                    {"targets",{{"bsh_residual",1.23456e-4},{"density_change",1e-4},{"energy_change",1e-5}}},
                    {"state_metrics",{{"ground_state",{{"subjects",{{"ground",{{"values",
                        {{"bsh_residual",2.04449e-5},{"density_change",3e-5},{"energy_change",2.73029e-8},{"kinetic_change",1.0}}}}}}}}}}}};
                it::add_results_source(ground_entry,"ground_state");
                it::write_results_entry(world,root,ground_entry);
                recovered.report();
                results=read_root_json(world,results_path);
                require(results.at("calculations").size()==2 &&
                        results.at("calculations").at(0).at("calculation_id")=="ground_state" &&
                        results.at("calculations").at(1).at("calculation_id")=="response");
                {
                    const auto& ground=results.at("calculations").at(0);
                    require(ground.at("convergence")==cp::json({{"bsh_residual",2.04e-5},{"density_change",3e-5},{"energy_change",2.73e-8}}) &&
                            !ground.contains("sources") && ground.at("targets").at("bsh_residual")==1.23e-4 &&
                            !ground.at("properties").at("total_energy").contains("convergence") &&
                            ground.at("properties").at("total_energy").at("value")==-1.0 &&
                            ground.at("state_metrics").at("ground_state").at("subjects").at("ground")
                                .at("values").at("bsh_residual")==2.04449e-5);
                }
                // Properties computed after an entry is written are added to
                // it; an entry with its own convergence gives them none.
                cp::json dipole{{"value",cp::json::array({0.0,0.0,1.0})},{"status","complete"}};
                it::add_results_source(dipole,"ground_state");
                it::add_results_properties(world,root,"ground_state",{{"dipole",dipole}});
                results=read_root_json(world,results_path);
                const auto& added=results.at("calculations").at(0).at("properties").at("dipole");
                require(added.at("value").size()==3 && !added.contains("sources") && !added.contains("convergence") &&
                        results.at("calculations").at(1).at("calculation_id")=="response");
                rejects([&]{it::add_results_properties(world,root,"missing",cp::json::object());});
            }
            // Convergence is the worst value of each metric over the states a
            // value used; an excited root counts through its own metrics, and
            // a node with no measured values adds nothing.
            {
                const cp::json table{
                    {"A",{{"subjects",{{"response",{{"values",{{"bsh_residual",1e-6},{"density_change",4e-5},{"kinetic_change",1.0}}}}}}}}},
                    {"set",{{"subjects",{{"root0",{{"values",{{"bsh_residual",5e-5},{"density_change",1e-6}}}}},
                                         {"root1",{{"values",{{"bsh_residual",9e-1},{"density_change",9e-1}}}}}}}}},
                    {"vbc",{{"subjects",cp::json::object()}}}};
                cp::json value=cp::json::object();
                it::add_results_source(value,"A");
                it::add_results_source(value,"set","root0");
                it::add_results_source(value,"vbc");
                it::add_results_source(value,"A"); // once
                require(value.at("sources").size()==3);
                require(it::worst_convergence(table,value.at("sources"))==cp::json({{"bsh_residual",5e-5},{"density_change",4e-5}}));
                require(it::worst_convergence(table,cp::json::array({{{"node","vbc"}},{{"node","absent"}}}))==cp::json::object());
                // A convergence reports every metric with a target, flat or
                // named by kind; with no targets, the residual and density.
                require(it::results_convergence_metrics({{"driven_states",{{"bsh_residual",1e-3},{"density_change",1e-3}}},
                            {"excited_states",{{"bsh_residual",1e-3},{"energy_change",1e-3}}}})==
                        std::vector<std::string>({"bsh_residual","density_change","energy_change"}));
                require(it::results_convergence_metrics(cp::json::object())==it::convergence_metrics());
                require(it::worst_convergence(table,value.at("sources"),{"bsh_residual","kinetic_change"})==
                        cp::json({{"bsh_residual",5e-5},{"kinetic_change",1.0}}));
                // Three significant figures, at any depth; integers, strings
                // and zero are left alone.
                require(it::round_significant(2.1578670402586745e-05)==2.16e-05 && it::round_significant(0.0005)==0.0005 &&
                        it::round_significant(0.13038009471965675)==0.13 && it::round_significant(-1.23456)==-1.23 &&
                        it::round_significant(0.0)==0.0);
                require(it::round_results_numbers({{"a",{{"b",1.23456e-4}}},{"n",3},{"s","x"}})==
                        cp::json({{"a",{{"b",1.23e-4}}},{"n",3},{"s","x"}}));
            }
            // Raman, which no check runs: one spectrum item per laser
            // frequency, a number, with the Hessian's states and that
            // frequency's polarizability-gradient states; the analysis is one
            // value.
            {
                const cp::json table{
                    {"nuc_0_x-f0.00000",{{"subjects",{{"response",{{"values",{{"bsh_residual",1e-6},{"density_change",2e-6}}}}}}}}},
                    {"dipole_x-f0.05000",{{"subjects",{{"response",{{"values",{{"bsh_residual",3e-6},{"density_change",1e-7}}}}}}}}},
                    {"dipole_x-f0.10000",{{"subjects",{{"response",{{"values",{{"bsh_residual",9e-1},{"density_change",9e-1}}}}}}}}}};
                cp::json hessian=cp::json::object();
                it::add_results_source(hessian,"nuc_0_x-f0.00000");
                cp::json at05{{"frequency",0.05},{"components",{{"nuc_0_x",{{"xx",1.0}}}}}},at10{{"frequency",0.10}};
                it::add_results_source(at05,"dipole_x-f0.05000");
                it::add_results_source(at10,"dipole_x-f0.10000");
                const cp::json gradient{{"values",cp::json::array({at05,at10})}};
                const cp::json whole{{"rows_found",9},{"rows_expected",9},{"all_rows_converged",true}};
                const cp::json short_rows{{"rows_found",3},{"rows_expected",9},{"all_rows_converged",true}};
                const cp::json rows=cp::json::array({{{"mode",1},{"freq_cm1",1000.0}}});
                const cp::json units{{"hessian","Eh/bohr^2"}};
                const cp::json analysis{{"hessian",cp::json::array({1.0})},{"frequencies",cp::json::array({0.01})}};
                const cp::json properties{
                    {"vibrational_analysis",molresponse_v3::vibrational_analysis_property(hessian,units,analysis)},
                    {"raman_spectrum",{{"values",cp::json::array({
                        molresponse_v3::raman_spectrum_item(0.05,hessian,gradient,whole,rows),
                        molresponse_v3::raman_spectrum_item(0.10,hessian,gradient,short_rows,cp::json()),
                        molresponse_v3::raman_spectrum_item(0.20,hessian,cp::json(),whole,cp::json())})}}}};
                const auto written=it::finish_results_properties(properties,table);
                const auto& va=written.at("vibrational_analysis");
                require(va.at("value")==analysis && va.at("status")=="complete" && va.at("units")==units &&
                        va.at("convergence")==cp::json({{"bsh_residual",1e-6},{"density_change",2e-6}}) &&
                        !va.contains("sources") && !va.contains("availability"));
                const auto& spectrum=written.at("raman_spectrum").at("values");
                require(spectrum.size()==3 && spectrum.at(0).at("frequency")==0.05 && spectrum.at(0).at("value")==rows &&
                        spectrum.at(0).at("status")=="complete" && !spectrum.at(0).contains("availability") &&
                        spectrum.at(0).at("convergence")==cp::json({{"bsh_residual",3e-6},{"density_change",2e-6}}));
                require(spectrum.at(1).at("status")=="unavailable" && !spectrum.at(1).contains("value") &&
                        spectrum.at(1).at("availability").contains("spectrum") &&
                        spectrum.at(1).at("convergence").at("bsh_residual")==9e-1);
                // No vibrational mode: a complete, empty spectrum.
                require(spectrum.at(2).at("status")=="complete" && spectrum.at(2).at("value")==cp::json::array() &&
                        spectrum.at(2).at("convergence").at("bsh_residual")==1e-6);
                const auto missing=molresponse_v3::vibrational_analysis_property(hessian,units,cp::json());
                require(missing.at("status")=="unavailable" && !missing.contains("value") &&
                        molresponse_v3::raman_spectrum_item(0.05,hessian,gradient,cp::json(),cp::json()).at("status")=="unavailable");
            }
            // The same collective metadata writer used by Raman must preserve
            // the old file and publish failure timing on every rank.
            const auto metadata_path=root+"/raman_metadata.json";
            rs::update_response_metadata(world,metadata_path,[](auto& metadata) {
                metadata.set_section("raman_spectra",{{"original",true}});
            });
            cp::detail::on_root(world,[&] {std::filesystem::create_directory(metadata_path+".tmp");return std::string{};});
            rejects([&]{recovered.with_failure_timing([&] {
                rs::update_response_metadata(world,metadata_path,[](auto& metadata) {
                    metadata.set_section("raman_spectra",{{"original",false}});
                });
            });});
            cp::detail::on_root(world,[&] {
                require(cp::detail::read_json(metadata_path).at("raman_spectra").at("original")==true);
                require(cp::detail::read_json(std::filesystem::path(root)/"timings.json").dump().find("cannot open for write")!=std::string::npos);
                return std::string{};
            });
            // A changed frequency inside the same rounded name also invalidates evidence.
            refined_settings["frequency"]=0.11000001;
            rs::DrivenState same_id(world,root,ground,"new-leg","dipole_z-f0.11000",refined_settings,request,ladder,derived_task);
            require(same_id.has_saved() && !same_id.reusable());
            require(same_id.load(read_payload)==123);
            auto wrong_origin=refined_settings;wrong_origin["derived_from"]["root"]="es_root_0001";
            rejects([&]{rs::DrivenState wrong(world,root,ground,"new-leg","dipole_z-f0.11000",wrong_origin,request,ladder,derived_task);});
            // timings.json: one entry per launch. A run's snapshot replaces its own
            // part; another run of the same launch adds a part; a new launch appends.
            {
                const auto entry=[&](const std::string& started,double seconds,const std::string& step) {
                    return it::timing_entry(started,{static_cast<std::size_t>(world.size()),2},
                        {{"name","run"},{"seconds",seconds},{"children",{{{"name",step},{"seconds",seconds}}}}});
                };
                const auto timings_path=std::filesystem::path(root)/"timed"/"timings.json";
                const auto timed=timings_path.parent_path();
                cp::detail::on_root(world,[&]{std::filesystem::create_directory(timed);return std::string{};});
                const auto first=it::write_timing_entry(world,timed,entry("T1",1.0,"scf"));
                const auto second=it::write_timing_entry(world,timed,entry("T1",2.0,"response"));
                require(first==0 && second==1);
                require(it::write_timing_entry(world,timed,entry("T1",3.0,"response"),second)==1);
                auto timings=read_root_json(world,timings_path);
                const auto& launch=timings.at("launches").at(0);
                require(timings.at("version")==2 && timings.at("launches").size()==1 && !launch.contains("attempt_id") &&
                        launch.at("started_at")=="T1" && launch.at("execution").at("threads_per_rank")==2);
                const auto& steps=launch.at("timing").at("children");
                require(steps.size()==2 && steps.at(0).at("name")=="scf" && steps.at(0).at("part")==0 &&
                        steps.at(1).at("seconds")==3.0 && steps.at(1).at("part")==1);
                (void)it::write_timing_entry(world,timed,entry("T2",1.0,"scf"));
                rejects([&]{(void)it::write_timing_entry(world,timed,entry("T3",1.0,"scf"),0);});
                require(read_root_json(world,timings_path).at("launches").size()==2);
                cp::detail::on_root(world,[&]{std::filesystem::remove(timings_path);
                    std::filesystem::create_directory(timings_path); return std::string{};});
                rejects([&]{(void)it::write_timing_entry(world,timed,entry("T2",1.0,"scf"));});
                // A run folder from an older version is refused, not read.
                cp::detail::on_root(world,[&]{std::filesystem::create_directory(timed/"checkpoints");return std::string{};});
                rejects([&]{it::require_current_run_folder(world,timed);});
            }
            // Each process group writes its own timing file; the main session
            // merges every one into its record and removes it.
            if(world.size()>1) {
              const auto comm=world.mpi.comm().Split(world.rank(),world.rank());
              {
                World sub(comm);
                const auto group=recovered.for_world(sub,{},world.rank());
                group->timing().end();group->timing().end();
                group->publish_timings();
                rejects([&]{group->report();});
              }
              world.gop.fence();
              cp::detail::on_root(world,[&] {
                for(int rank=0;rank<world.size();++rank)
                  require(std::filesystem::is_regular_file(rs::group_timing_file(root,rank)));
                return std::string{};
              });
              recovered.collect_subworlds();
              std::size_t merged=0;
              // NOLINTNEXTLINE(misc-no-recursion) -- a finite timing tree
              const auto count=[&](const auto& self,const cp::json& node)->void {
                merged+=node.at("name")=="subworld_work";
                for(const auto& child:node.value("children",cp::json::array())) self(self,child);
              };
              count(count,recovered.timing().record());
              require(merged==static_cast<std::size_t>(world.size()));
              cp::detail::on_root(world,[&] {
                for(int rank=0;rank<world.size();++rank)
                  require(!std::filesystem::exists(rs::group_timing_file(root,rank)));
                return std::string{};
              });
            }
            cp::detail::on_root(world,[&]{std::filesystem::remove_all(root);return std::string{};});
        } catch (const std::exception& e) { std::cerr<<e.what()<<'\n';failed=1; }
        world.gop.max(failed);
    }
    finalize();return failed;
}
