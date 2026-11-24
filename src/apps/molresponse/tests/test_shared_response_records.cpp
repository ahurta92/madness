#include "../solvers/shared_response_session.hpp"
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
            // node is reused, and a further move starts from the newest node.
            cp::detail::on_root(world,[&] {std::filesystem::create_directory(root+"/recovery");return std::string{};});
            it::TimingCollector<> recovery_timing;
            recovery_timing.begin("run"); recovery_timing.begin("response"); recovery_timing.begin("setup");
            rs::Session recovered(world,root,ground,request,ladder,{},root+"/recovery",std::move(recovery_timing),it::Verbosity::detail);
            auto& recovered_leg=recovered.bind_derived("new-leg",refined_settings);
            require(recovered_leg.reusable() && recovered_leg.description()==refined.description());
            auto moved_settings=refined_settings;moved_settings["frequency"]=0.12;
            auto& moved=recovered.bind_derived("moved-leg",moved_settings);
            require(moved.has_saved() && !moved.has_node() && moved.saved_metrics().states.size()==1);
            recovered.report();
            // results.json is written directly: response's entry names its nodes
            // and no launch; another application's entry is left alone.
            const auto results_path=std::filesystem::path(root)/"results.json";
            {
                auto results=read_root_json(world,results_path);
                require(results.at("version")==2 && results.at("calculations").size()==1 && !results.contains("attempt_id"));
                const auto& entry=results.at("calculations").at(0);
                require(entry.at("calculation_id")=="response" && !entry.contains("attempt_id") && entry.at("driven").is_array());
                it::write_results_entry(world,root,{{"calculation_id","ground_state"},{"node","states/ground_state"}});
                recovered.report();
                results=read_root_json(world,results_path);
                require(results.at("calculations").size()==2 &&
                        results.at("calculations").at(0).at("calculation_id")=="ground_state" &&
                        results.at("calculations").at(1).at("calculation_id")=="response");
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
