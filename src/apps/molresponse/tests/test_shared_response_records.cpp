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
            const rs::GroundSource ground{root,{"ground_state","scf","ground"}};
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
            rs::Checkpoint bundle(world,ground,"roots",cp::json::object(),request,ladder,"solve",task);
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
            require(cp::iteration_metrics_from_manifest(bundle.manifest()).states.size()==2);
            rs::Checkpoint reused(world,ground,"roots",cp::json::object(),request,ladder,"restart",task);
            require(reused.reusable());
            const auto read_payload=[&](const auto& directory,const cp::Manifest&) {
                archive::ParallelInputArchive in(world,(directory/"response").c_str());
                int value=0;
                if (in.is_io_node()) in.local_archive() & value;
                in.close();
                world.gop.broadcast_serializable(value,0);
                return value;
            };
            require(reused.load(read_payload,true)==123);
            require(reused.source_record("root1").at("subject")=="root1");
            require(reused.measured_values("root1").at("values").at("bsh_residual")==0.0);
            rejects([&]{reused.measured_values();});
            const auto inputs=cp::json::array({reused.qualified_input("B")});
            rs::DirectCheckpoint direct(world,root,"VBC",{{"kind","VBC"}},inputs,"build");
            require(!direct.reusable());
            direct.publish([&](const auto& dir){return payload(world,dir);});
            require(!direct.manifest().record.at("task").contains("iteration"));
            rs::DirectCheckpoint cached(world,root,"VBC",{{"kind","VBC"}},inputs,"again");
            require(cached.reusable());
            require(cached.load(read_payload)==123);
            require(cached.qualified_input("VBC").at("generation")==direct.manifest().generation);
            require(cached.source_record().at("attempt_id")=="build");
            // A is a property input, not part of the direct B/C definition.
            rs::DirectCheckpoint new_a(world,root,"VBC",{{"kind","VBC"}},inputs,"new-A");
            require(new_a.reusable());
            auto tighter=request; tighter.targets[it::Metric::bsh_residual]=0.0;
            rs::Checkpoint tighter_bundle(world,ground,"roots",cp::json::object(),tighter,ladder,"tighter",task);
            require(tighter_bundle.reusable());
            require(tighter_bundle.qualified_input("B")==reused.qualified_input("B"));
            auto bad=inputs; bad[0]["generation"]="missing";
            rejects([&]{rs::DirectCheckpoint invalid(world,root,"VBC",{{"kind","VBC"}},bad,"bad");});
            metrics.states.back().metrics[it::Metric::bsh_residual]=it::MetricValue::measured(1.0,1.0);
            const auto capped=it::decide(bundle.controls(),bundle.expected(),progress,metrics,bundle.final_level());
            bundle.publish(metrics,capped,progress,[&](const auto& dir){return payload(world,dir);});
            rejects([&]{bundle.qualified_input("C");});
            auto changed_inputs=inputs;
            changed_inputs[0]["generation"]=bundle.manifest().generation;
            rs::DirectCheckpoint changed_b(world,root,"VBC",{{"kind","VBC"}},changed_inputs,"new-B");
            require(!changed_b.reusable());
            rejects([&]{changed_b.qualified_input("VBC");});
            rejects([&]{changed_b.load(read_payload);});
            auto expanded_task=task;
            expanded_task.allow_subject_expansion=true;
            expanded_task.record["subjects"]={"root0","root1","root2"};
            rs::Checkpoint expanded(world,ground,"roots",cp::json::object(),request,ladder,"expanded",expanded_task);
            require(expanded.has_saved() && !expanded.reusable());
            require(expanded.plan().action==it::RestartAction::resume && expanded.plan().start_protocol==0);
            require(expanded.load(read_payload)==123);
            rejects([&]{expanded.qualified_input("root");});
            auto replacement=expanded_task;replacement.record["subjects"]={"root0","different","root2"};
            rejects([&]{rs::Checkpoint invalid(world,ground,"roots",cp::json::object(),request,ladder,"replaced",replacement);});
            const auto derived_task=rs::derived_task();
            const cp::json origin{{"bundle","es:tda"},{"root","es_root_0000"},{"factor",0.5}};
            const cp::json old_settings{{"frequency",0.1},{"perturbation","dipole_z"},
                {"layout","closed_shell_full"},{"derived_from",origin}};
            rs::Checkpoint old_leg(world,ground,"old-leg",old_settings,request,ladder,"derived-old",derived_task);
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
            const auto old_decision=it::decide(old_leg.controls(),old_leg.expected(),progress,old_metrics,old_leg.final_level());
            old_leg.publish(old_metrics,old_decision,progress,[&](const auto& dir){return payload(world,dir);});
            const cp::GenerationRef previous{"response","old-leg",old_leg.manifest().generation};
            auto refined_settings=old_settings;refined_settings["frequency"]=0.11;
            rs::Checkpoint refined(world,ground,"new-leg",refined_settings,request,ladder,"derived-new",derived_task,previous);
            require(refined.has_saved() && !refined.reusable() && refined.plan().start_protocol==0);
            require(refined.load(read_payload)==123);
            rejects([&]{refined.qualified_input("B");});
            auto new_metrics=old_metrics;new_metrics.info=refined.info();
            new_metrics.states[0].metrics[it::Metric::density_change]=density(refined.info());
            const auto new_decision=it::decide(refined.controls(),refined.expected(),progress,new_metrics,refined.final_level());
            refined.publish(new_metrics,new_decision,progress,[&](const auto& dir){return payload(world,dir);});
            require(refined.manifest().record.at("task").at("definition").at("replaces").at("task_id")=="old-leg");
            // Simulate a crash after replacement publication and before report/retirement.
            // Both origins remain active; the exact current task must win.
            cp::detail::on_root(world,[&] {std::filesystem::create_directory(root+"/recovery");return std::string{};});
            it::TimingCollector<> recovery_timing;
            recovery_timing.begin("run"); recovery_timing.begin("response"); recovery_timing.begin("setup");
            rs::Session recovered(world,ground,request,ladder,{},root+"/recovery",std::move(recovery_timing));
            auto& recovered_leg=recovered.bind_derived("new-leg",refined_settings);
            require(recovered_leg.reusable());
            require(recovered_leg.manifest().generation==refined.manifest().generation);
            // Retiring the predecessor must preserve any exact live dependency.
            rs::DirectCheckpoint pinned(world,root,"old-dependent",{{"kind","fixture"}},
                cp::json::array({old_leg.qualified_input("old")}),"pin-old");
            pinned.publish([&](const auto& dir){return payload(world,dir);});
            recovered.report();
            require(store.read(previous).generation==previous.generation);
            rejects([&]{store.read_latest("response","old-leg");});
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
            require(store.read_latest("response","new-leg")->generation==refined.manifest().generation);
            // A changed frequency inside the same rounded ID also invalidates evidence.
            refined_settings["frequency"]=0.11000001;
            rs::Checkpoint same_id(world,ground,"new-leg",refined_settings,request,ladder,"derived-bin",derived_task);
            require(same_id.has_saved() && !same_id.reusable());
            require(same_id.load(read_payload)==123);
            auto wrong_origin=refined_settings;wrong_origin["derived_from"]["root"]="es_root_0001";
            rejects([&]{rs::Checkpoint wrong(world,ground,"new-leg",wrong_origin,request,ladder,"wrong-root",derived_task);});
            rejects([&]{store.read_latest("response","old-leg");});
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
