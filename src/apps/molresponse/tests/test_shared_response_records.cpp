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
cp::Payload payload(World& world, const std::filesystem::path& directory) {
    archive::ParallelOutputArchive out(world, (directory/"response").c_str(), 1);
    int value = 123;
    if (out.is_io_node()) out.local_archive() & value;
    cp::Payload result;
    result.local_error = cp::close_archive(out);
    world.gop.fence();
    if (world.rank() == 0 && result.local_error.empty())
        result.archives.push_back(cp::inspect_archive(directory, "functions", "response", {}, 1));
    return result;
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
            cp::CheckpointStore store(world, root);
            store.initialize("initial");
            store.register_task("ground_state", "scf", "scf", "direct", true);
            const cp::json model{{"version",1},{"dft",{{"xc","HF"}}},
                {"molecule",{{"parameters",cp::json::object()},{"core_potentials",cp::json::array()}}},
                {"xc",{{"specification","HF"}}}};
            const cp::json stage{{"application","scf"},{"settings",cp::json::object()},
                {"requested_properties",cp::json::array()},{"inputs",cp::json::array()},
                {"outcome",{{"complete",true},{"reason","computed"}}}};
            store.publish(cp::direct_manifest("ground_state","initial","ground",stage,"scf",
                {{"scientific_settings",model}},cp::json::array()),
                [&](const auto& dir){return payload(world,dir);});
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
            rs::DrivenState bundle(world,root,ground,"roots","roots",cp::json::object(),request,ladder,"solve",task);
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
            bundle.publish(metrics,decision,progress,[&](const auto& dir){return payload(world,dir);});
            require(bundle.has_node() && bundle.description().kind=="driven" && bundle.saved_metrics().states.size()==2);
            rs::DrivenState reused(world,root,ground,"roots","roots",cp::json::object(),request,ladder,"restart",task);
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
            const auto save_payload=[&](const auto& dir){(void)payload(world,dir);};
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
            rs::DrivenState tighter_bundle(world,root,ground,"roots","roots",cp::json::object(),tighter,ladder,"tighter",task);
            require(tighter_bundle.reusable());
            require(rs::lineage_entry("B",tighter_bundle)==rs::lineage_entry("B",reused));
            metrics.states.back().metrics[it::Metric::bsh_residual]=it::MetricValue::measured(1.0,1.0);
            const auto capped=it::decide(bundle.controls(),bundle.expected(),progress,metrics,bundle.final_level());
            bundle.publish(metrics,capped,progress,[&](const auto& dir){return payload(world,dir);});
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
            rs::ExcitedSet es(world,root,ground,"es:tda",es_settings,request,ladder,"es-solve",es_task);
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
            rs::ExcitedSet es_again(world,root,ground,"es:tda",es_settings,request,ladder,"es-again",es_task);
            require(es_again.reusable() && es_again.description()==es.description());
            require(es_again.description().details.at("roots")[1].at("omega")==1.5);
            require(es_again.load([](const auto& dir,const auto&){return std::filesystem::exists(dir/"root_1.00000");}));
            auto three=es_task;three.record=es_record(3);
            rs::ExcitedSet es_more(world,root,ground,"es:tda",es_settings,request,ladder,"es-more",three);
            require(es_more.has_saved() && !es_more.reusable() && es_more.plan().action==it::RestartAction::resume);
            auto one=es_task;one.record=es_record(1);
            rejects([&]{rs::ExcitedSet fewer(world,root,ground,"es:tda",es_settings,request,ladder,"es-fewer",one);});
            auto moved_ground=ground;moved_ground.identity["saved_at"]="2026-11-24T00:00:00Z";
            rs::ExcitedSet es_guess(world,root,moved_ground,"es:tda",es_settings,request,ladder,"es-guess",es_task);
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
            rs::DrivenState old_leg(world,root,ground,"old-leg","dipole_z-f0.10000",old_settings,request,ladder,"derived-old",derived_task);
            it::IterationMetrics old_metrics;old_metrics.info=old_leg.info();
            it::StateMetrics leg;leg.id="response";
            leg.metrics[it::Metric::bsh_residual]=it::MetricValue::measured(0.0,1.0);
            // Driven response judges the density change by default, against a
            // comparable earlier iteration of the same task.
            const auto density=[](const it::IterationInfo& info) {
              auto value=it::MetricValue::measured(0.0,1.0);
              value.comparison=it::Comparison{info.task_id,"earlier",info.protocol,0,"response"};
              return value;
            };
            leg.metrics[it::Metric::density_change]=density(old_leg.info());
            old_metrics.states={leg};
            // One update at the level: a stored comparison is the update before.
            const it::Progress leg_progress{0,1,true};
            const auto old_decision=it::decide(old_leg.controls(),old_leg.expected(),leg_progress,old_metrics,old_leg.final_level());
            old_leg.publish(old_metrics,old_decision,leg_progress,[&](const auto& dir){(void)payload(world,dir);});
            const auto old_description=old_leg.description();
            // The moved root's state starts from the nearby frequency's node,
            // only reading it, and names it in its lineage.
            auto refined_settings=old_settings;refined_settings["frequency"]=0.11;
            rs::DrivenState refined(world,root,ground,"new-leg","dipole_z-f0.11000",refined_settings,request,ladder,
                                    "derived-new",derived_task,std::string("dipole_z-f0.10000"));
            require(refined.has_saved() && !refined.has_node() && !refined.reusable() && refined.plan().start_protocol==0);
            require(refined.load(read_payload)==123);
            rejects([&]{(void)rs::lineage_entry("B",refined);});
            auto new_metrics=old_metrics;new_metrics.info=refined.info();
            new_metrics.states[0].metrics[it::Metric::density_change]=density(refined.info());
            const auto new_decision=it::decide(refined.controls(),refined.expected(),leg_progress,new_metrics,refined.final_level());
            refined.publish(new_metrics,new_decision,leg_progress,[&](const auto& dir){(void)payload(world,dir);});
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
            rs::DrivenState same_id(world,root,ground,"new-leg","dipole_z-f0.11000",refined_settings,request,ladder,"derived-bin",derived_task);
            require(same_id.has_saved() && !same_id.reusable());
            require(same_id.load(read_payload)==123);
            auto wrong_origin=refined_settings;wrong_origin["derived_from"]["root"]="es_root_0001";
            rejects([&]{rs::DrivenState wrong(world,root,ground,"new-leg","dipole_z-f0.11000",wrong_origin,request,ladder,"wrong-root",derived_task);});
            require(store.task_inventory("response").is_array());
            const auto timing=cp::json{{"attempt_id","timed"},{"execution",{{"mpi_ranks",world.size()},{"threads_per_rank",2}}},
                {"timing",{{"name","run"},{"seconds",0.0},{"clock_scope","driver_rank"}}}};
            store.publish_timings(timing);
            cp::detail::on_root(world,[&]{std::filesystem::remove(std::filesystem::path(root)/"timings.json");
                std::filesystem::create_directory(std::filesystem::path(root)/"timings.json"); return std::string{};});
            rejects([&]{store.publish_timings(timing);});
            if(world.size()>1) {
              // Each rank is an independent task World. The shared coordinator
              // serializes publication; no universe collective enters a writer.
              const auto comm=world.mpi.comm().Split(world.rank(),world.rank());
              {
                World sub(comm);
                const auto id="parallel-"+std::to_string(world.rank());
                cp::CheckpointStore concurrent(sub,root,{},"same-invocation");
                concurrent.initialize(id);
                concurrent.register_task("response","molresponse",id,"direct",true);
                auto direct_stage=stage;direct_stage["application"]="molresponse";
                for(int iteration=0;iteration<3;++iteration) {
                  const auto manifest=cp::direct_manifest("response",id,id+"-"+std::to_string(iteration),
                      direct_stage,id,cp::json::object(),cp::json::array());
                  concurrent.publish(manifest,[&](const auto& dir){return payload(sub,dir);});
                }
                if(world.rank()==0) {
                  const auto broken=cp::direct_manifest("response",id,id+"-broken",direct_stage,id,cp::json::object(),cp::json::array());
                  rejects([&]{concurrent.publish(broken,[&](const auto&){cp::Payload bad;bad.local_error="injected subgroup write failure";return bad;});});
                }
              }
              world.gop.fence();
              for(int rank=0;rank<world.size();++rank)
                require(store.read_latest("response","parallel-"+std::to_string(rank))->generation=="parallel-"+std::to_string(rank)+"-2");
            }
            cp::detail::on_root(world,[&]{std::filesystem::remove_all(root);return std::string{};});
        } catch (const std::exception& e) { std::cerr<<e.what()<<'\n';failed=1; }
        world.gop.max(failed);
    }
    finalize();return failed;
}
