// Native archives, collective failures and live-reference retention; no
// chemistry solve.
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <madness/chem/CheckpointStore.h>
#include <madness/world/MADworld.h>

namespace cp = madness::checkpoint;
namespace it = madness::iteration;
namespace fs = std::filesystem;
using madness::World;
namespace {
int failures = 0;
void check(bool ok, const char *why) {
  if (!ok) {
    ++failures;
    std::cerr << "FAILED: " << why << '\n';
  }
}
template <class F> void rejects(F &&f, const char *why) {
  bool threw = false;
  try {
    f();
  } catch (const std::exception &) {
    threw = true;
  }
  check(threw, why);
}
void text_file(const fs::path &path, const std::string &text) {
  std::ofstream out(path);
  out.exceptions(std::ios::failbit | std::ios::badbit);
  out << text;
  out.close();
}
cp::json stage() {
  return {{"application", "scf"},
          {"settings", {{"model", "HF"}}},
          {"requested_properties", {"energy"}},
          {"inputs", cp::json::array()},
          {"outcome", {{"complete", false}, {"reason", "in_progress"}}}};
}
cp::json definitions() {
  return {{"bsh.v1",
           {{"metric", "bsh_residual"},
            {"meaning", "raw fixed-point difference"},
            {"units", "dimensionless"},
            {"fields", {"orbitals"}},
            {"weights", "none"},
            {"normalization", "none"},
            {"reduction", "maximum"},
            {"comparison", "none"},
            {"scale", "maximum current function norm"}}}};
}
cp::Manifest evaluated(const std::string &calculation,
                       const std::string &generation,
                       const std::string &attempt = "a1",
                       const cp::json &inputs = cp::json::array()) {
  it::IterationMetrics e;
  e.info = {"state", attempt, {1e-6, 8}, 2, "model-context"};
  it::StateMetrics subject;
  subject.id = "ground";
  subject.metrics[it::Metric::bsh_residual] =
      it::MetricValue::measured(1e-7, 1.0);
  e.states = {subject};
  it::Decision decision{
      e.info, it::Action::complete, "final_criteria_met", true, {}};
  return cp::iterative_manifest(calculation, generation, stage(),
                                {{"subject", "ground"}}, inputs, definitions(),
                                {{it::Metric::bsh_residual, "bsh.v1"}}, e,
                                decision, {2, 2, true});
}
cp::GenerationRef ref(const std::string &calc, const std::string &generation) {
  return {calc, "state", generation};
}
cp::json inputs(const std::string &calc, const std::string &generation) {
  return {{{"role", "reference"},
           {"calculation_id", calc},
           {"task_id", "state"},
           {"generation", generation}}};
}
cp::json source_map(const std::string &calc, const std::string &owner,
                    const std::string &generation,
                    const std::string &producer = "a1") {
  return {{"version", 2},
          {"calculation_id", calc},
          {"attempt_id", owner},
          {"states",
           {{"ground",
             {{"task_id", "state"},
              {"subject", "ground"},
              {"attempt_id", producer},
              {"protocol", {{"threshold", 1e-6}, {"k", 8}}},
              {"iteration", 2},
              {"generation", generation}}}}},
          {"tasks", cp::json::array()}};
}
cp::json result(const std::string &calc, const std::string &owner) {
  return {{"version", 2},
          {"attempt_id", owner},
          {"calculations", {{{"calculation_id", calc}}}}};
}
cp::Payload write_payload(World &world, const fs::path &dir, int value) {
  madness::archive::ParallelOutputArchive out(world, (dir / "state").c_str(),
                                              1);
  if (out.is_io_node())
    out.local_archive() & value;
  cp::Payload payload;
  payload.local_error = cp::close_archive(out);
  world.gop.fence();
  if (world.rank() == 0 && payload.local_error.empty()) {
    try {
      payload.archives.push_back(
          cp::inspect_archive(dir, "state", "state", {}, 1));
    } catch (const std::exception &e) {
      payload.local_error = e.what();
    }
  }
  return payload;
}
void load(World &world, cp::CheckpointStore &store,
          const cp::GenerationRef &source, int expected) {
  auto m = store.read(source);
  madness::archive::ParallelInputArchive in(
      world, (store.task_directory(source.calculation_id, source.task_id) /
              source.generation / "state")
                 .c_str());
  if (in.is_io_node()) {
    int value = 0;
    in.local_archive() & value;
    check(value == expected, "exact original native payload loads");
  }
  in.close();
  check(m.version == 2, "shared version-2 record loaded");
}
void iteration_record_cases() {
  it::IterationMetrics metrics;
  metrics.info = {"scf", "a2", {1e-6, 8}, 0, "same-model"};
  it::StateMetrics ground;
  ground.id = "ground";
  auto value = it::MetricValue::measured(1e-8, 1.0);
  value.comparison = it::Comparison{"scf", "a1", {1e-6, 8}, 7, "ground"};
  ground.metrics[it::Metric::energy_change] = value;
  metrics.states = {ground};
  auto meaning = definitions().at("bsh.v1");
  meaning["metric"] = "energy_change";
  meaning["units"] = "hartree";
  auto decision = it::Decision{
      metrics.info, it::Action::complete, "final_criteria_met", true, {}};
  const auto make = [&](const it::Progress &progress) {
    return cp::iterative_manifest("ground_state", "g1", stage(),
                                  {{"subject", "ground"}}, cp::json::array(),
                                  {{"energy.v1", meaning}},
                                  {{it::Metric::energy_change, "energy.v1"}},
                                  metrics, decision, progress);
  };
  const auto manifest = make({0, 7, true});
  const auto restored = cp::iteration_metrics_from_manifest(manifest);
  const auto comparison =
      restored.states[0].metrics.at(it::Metric::energy_change).comparison;
  check(restored.info == metrics.info && comparison &&
            comparison->task_id == "scf" && comparison->attempt_id == "a1" &&
            comparison->iteration == 7 &&
            comparison->protocol == it::Protocol{1e-6, 8},
        "zero-update restart record and prior iteration identity round-trip "
        "exactly");
  rejects([&] { make({1, 8, true}); },
          "stale iteration record cannot publish at a later update count");
  decision.info.iteration = 1;
  rejects([&] { make({0, 7, true}); },
          "decision must bind the same full iteration identity");
  auto invalid = manifest;
  auto &prior = invalid.record["task"]["iteration"]["metrics"][0]["comparison"];
  prior["attempt_id"] = "a2";
  prior["iteration"] = 0;
  rejects([&] { cp::records::validate(invalid.record); },
          "saved change cannot compare an iteration with itself");
  prior["iteration"] = 1;
  rejects(
      [&] { cp::records::validate(invalid.record); },
      "saved comparison cannot refer to a future iteration in the same solve");
  prior["iteration"] = 0;
  prior["protocol"]["threshold"] = 1e-4;
  cp::records::validate(invalid.record);
  prior["task_id"] = "other_task";
  rejects([&] { cp::records::validate(invalid.record); },
          "saved comparison must belong to its task");
  invalid = manifest;
  auto &unavailable = invalid.record["task"]["iteration"]["metrics"][0];
  unavailable["status"] = "unavailable";
  unavailable["reason"] = "not_computed";
  unavailable.erase("value");
  cp::records::validate(invalid.record);
  unavailable["comparison"]["task_id"] = "other_task";
  rejects([&] { cp::records::validate(invalid.record); },
          "unavailable metric cannot hide a cross-task comparison");
  unavailable["comparison"]["task_id"] = "scf";
  unavailable["comparison"]["attempt_id"] = "a2";
  unavailable["comparison"]["iteration"] = 0;
  rejects([&] { cp::records::validate(invalid.record); },
          "unavailable metric cannot hide a self comparison");
  unavailable["comparison"] = "malformed";
  rejects([&] { cp::records::validate(invalid.record); },
          "unavailable malformed comparison fails manifest validation");
  invalid = manifest;
  invalid.record["task"]["iteration"]["updates_this_attempt"] = -1;
  rejects([&] { cp::records::validate(invalid.record); },
          "negative iteration is not a valid identity");
}
void exercise(World &world, const fs::path &run) {
  cp::CheckpointStore store(world, run);
  store.initialize("a1");
  store.register_task("ground_state", "scf", "state", "iterative", true);
  check(!store.read_latest("ground_state", "state"),
        "registered null pointer means no numerical publication");
  int writes = 0;
  auto publish = [&](const cp::Manifest &manifest, int value) {
    store.publish(manifest, [&](const fs::path &dir) {
      ++writes;
      return write_payload(world, dir, value);
    });
  };
  auto first = evaluated("ground_state", "g1");
  publish(first, 11);
  auto loaded = store.read_latest("ground_state", "state");
  check(loaded && loaded->generation == "g1",
        "publication discovers exact generation without inventory update");
  const auto evidence = cp::iteration_metrics_from_manifest(*loaded);
  check(evidence.info.iteration == 2 && evidence.info.attempt_id == "a1" &&
            evidence.info.protocol == it::Protocol{1e-6, 8} &&
            evidence.states[0].metrics.at(it::Metric::bsh_residual).value ==
                1e-7,
        "metrics retain the evaluated sample and measured values");
  check(!cp::to_json_record(*loaded)
             .at("task")
             .at("iteration")
             .at("metrics")[0]
             .contains("absolute"),
        "checkpoint measurements do not contain stopping tolerances");
  check(!loaded->record.at("task").at("iteration").contains("sample_id"),
        "checkpoint identity uses the update counter without a separate sample "
        "label");
  store.publish_source_map(source_map("ground_state", "a1", "g1"));
  store.publish_results(result("ground_state", "a1"));
  publish(evaluated("ground_state", "g2", "a2"), 22);
  store.cleanup();
  load(world, store, ref("ground_state", "g1"), 11);
  load(world, store, ref("ground_state", "g2"), 22);
  check(writes == 2, "cleanup and reads do not publish numerical archives");

  // An active pending task protects A transitively, even after result B
  // replaces A.
  store.register_task("response", "molresponse", "state", "direct", true,
                      inputs("ground_state", "g1"));
  auto direct_stage = stage();
  direct_stage["application"] = "molresponse";
  auto direct =
      cp::direct_manifest("response", "a1", "vbc1", direct_stage, "state",
                          {{"meaning", "VBC"}}, inputs("ground_state", "g1"));
  publish(direct, 33);
  check(!cp::to_json_record(store.read(ref("response", "vbc1")))
             .at("task")
             .contains("iteration"),
        "direct checkpoint has completion and sources without fabricated "
        "iteration metrics");
  store.publish_source_map(source_map("ground_state", "a2", "g2", "a2"));
  store.publish_results(result("ground_state", "a2"));
  store.cleanup();
  load(world, store, ref("ground_state", "g1"), 11);
  rejects([&] { store.retire_calculation("ground_state"); },
          "current results and live dependencies prohibit retirement");
  store.retire_calculation("response");
  if (world.rank() == 0)
    fs::rename(store.task_directory("response", "state") / "vbc1",
               store.task_directory("response", "state") / ".deleting-vbc1");
  world.gop.fence();
  store.cleanup();
  if (world.rank() == 0) {
    check(!fs::exists(store.task_directory("ground_state", "state") /
                      "g1/state.00000"),
          "last live dependency ending permits payload reclamation");
    check(fs::exists(store.task_directory("ground_state", "state") /
                     "g1/manifest.json"),
          "historical metadata retains original identities and reserves its "
          "generation");
    check(!fs::exists(store.task_directory("response", "state") /
                      "vbc1/state.00000"),
          "retired latest payload can be reclaimed");
    check(!fs::exists(store.task_directory("response", "state") /
                      ".deleting-vbc1"),
          "interrupted deletion finishes on repeat cleanup");
    check(fs::exists(
              run / "checkpoints/ground_state/attempts/a1/result_sources.json"),
          "historical source map survives payload reclamation");
  }
  // Interrupted remove_all after its manifest was removed: the retained
  // metadata must still reserve the original generation identity.
  if (world.rank() == 0) {
    const auto doomed =
        store.task_directory("response", "state") / ".deleting-vbc1";
    fs::create_directory(doomed);
    text_file(doomed / "remaining-shard", "partial deletion");
  }
  world.gop.fence();
  check(store.cleanup().empty(), "partial payload deletion resumes");
  if (world.rank() == 0) {
    check(fs::exists(store.task_directory("response", "state") /
                     "vbc1/manifest.json"),
          "partial deletion preserves original metadata");
    check(!fs::exists(store.task_directory("response", "state") /
                      ".deleting-vbc1"),
          "partial deletion removes the leftover shard");
  }
  rejects([&] { store.read(ref("ground_state", "g1")); },
          "missing historical payload never substitutes newer state");
  rejects([&] { publish(evaluated("ground_state", "g1"), 99); },
          "reclaimed generation identity cannot be reused");
  store.cleanup(); // repeatable after retirement/deletion
  store.initialize("a3");
  check(store.read_latest("ground_state", "state")->generation == "g2",
        "new attempt does not discard omitted saved tasks");

  // Failed first publication leaves a null pointer, collectively.
  store.register_task("new_state", "scf", "state", "iterative", true);
  auto bad = evaluated("new_state", "bad");
  rejects(
      [&] {
        store.publish(bad, [&](const fs::path &dir) {
          auto p = write_payload(world, dir, 44);
          if (world.rank() == world.size() - 1)
            p.local_error = "injected close error";
          return p;
        });
      },
      "error on any writer prevents publication on every rank");
  check(!store.read_latest("new_state", "state"),
        "failed first write leaves null");
  auto unregistered = evaluated("absent", "g1");
  const auto before = writes;
  rejects([&] { publish(unregistered, 55); },
          "unregistered task rejected before writer");
  check(writes == before,
        "registration validation precedes numerical callback");

  // Simulate interruption after initializing a null pointer but before
  // registration.
  if (world.rank() == 0) {
    auto path = store.task_directory("interrupted", "state");
    fs::create_directories(path);
    text_file(path / "latest", "{\"generation\":null}");
  }
  world.gop.fence();
  store.register_task("interrupted", "scf", "state", "iterative", true);
  check(!store.read_latest("interrupted", "state"),
        "valid orphan null pointer can finish registration");
  if (world.rank() == 0)
    fs::remove(store.task_directory("interrupted", "state") / "latest");
  world.gop.fence();
  rejects(
      [&] {
        store.register_task("interrupted", "scf", "state", "iterative", true);
      },
      "registered missing pointer is not reinitialized");
  rejects([&] { store.cleanup(); },
          "unreadable live record prevents every deletion");
  load(world, store, ref("ground_state", "g2"), 22);
  if (world.rank() == 0)
    text_file(store.task_directory("interrupted", "state") / "latest",
              "{\"generation\":null}");
  world.gop.fence();

  // Publication failure after B's manifest but before pointer replacement keeps
  // A.
  if (world.rank() == 0)
    fs::create_directory(store.task_directory("ground_state", "state") /
                         ".latest-g3.tmp");
  world.gop.fence();
  rejects([&] { publish(evaluated("ground_state", "g3", "a3"), 66); },
          "pointer failure leaves previous checkpoint selected");
  check(store.read_latest("ground_state", "state")->generation == "g2",
        "replacement failure preserves current numerical progress");
  load(world, store, ref("ground_state", "g2"), 22);
  rejects([&] { store.read(ref("ground_state", "g3")); },
          "unpublished manifest cannot be read as a saved source");
  rejects(
      [&] {
        store.publish_source_map(source_map("ground_state", "a3", "g3", "a3"));
      },
      "unpublished manifest cannot become a result source");
  rejects(
      [&] {
        store.register_task("pending", "molresponse", "state", "direct", true,
                            inputs("ground_state", "g3"));
      },
      "unpublished manifest cannot become a dependency");
  auto wrong = evaluated("ground_state", "invalid");
  wrong.record["task"]["iteration"]["metrics"].push_back(
      wrong.record["task"]["iteration"]["metrics"][0]);
  rejects([&] { publish(wrong, 1); },
          "duplicate metric/subject pair rejected before archive writer");
  wrong = evaluated("ground_state", "invalid");
  wrong.record["task"]["iteration"]["metrics"][0]["absolute"] = 1e-5;
  rejects([&] { publish(wrong, 1); },
          "checkpoint metric cannot carry request-specific checks");
  auto mixed = result("ground_state", "a3");
  mixed["calculations"][0]["attempt_id"] = "a2";
  store.publish_results(mixed);
  mixed["calculations"][0]["attempt_id"] = "missing";
  rejects([&] { store.publish_results(mixed); },
          "invalid explicit source-map binding cannot inherit root or latest");
  auto rebind = source_map("ground_state", "iteration_mismatch", "g2", "a2");
  rebind["states"]["ground"]["iteration"] = 3;
  rejects([&] { store.publish_source_map(rebind); },
          "source map must bind the correct iteration");
  auto wrong_protocol =
      source_map("ground_state", "protocol_mismatch", "g2", "a2");
  wrong_protocol["states"]["ground"]["protocol"]["k"] = 6;
  rejects([&] { store.publish_source_map(wrong_protocol); },
          "result source must bind the exact protocol as well as iteration");
  auto changed_map = source_map("ground_state", "a2", "g2", "a2");
  changed_map["states"] = cp::json::object();
  rejects([&] { store.publish_source_map(changed_map); },
          "valid new mapping cannot rebind an immutable map");
  store.register_task("property", "response", "report", "direct", false);
  if (world.rank() == 0)
    check(!fs::exists(store.task_directory("property", "report") / "latest"),
          "payload-free task gets no fake checkpoint pointer");
  // Rank-zero publication metadata selects every writer's destination.
  auto fourth = evaluated("ground_state", "g4", "a4");
  if (world.rank() != 0)
    fourth.generation = "rank-local-ignored";
  store.publish(fourth, [&](const fs::path &dir) {
    check(dir.filename() == "g4", "all writers use canonical sample directory");
    return write_payload(world, dir, 77);
  });
  if (world.rank() == 0)
    check(!fs::exists(store.task_directory("new_state", "state") / "bad"),
          "successful cleanup removes abandoned unpublished payloads");
  if (world.rank() == 0)
    fs::rename(run / "results.json", run / "saved-results.json");
  world.gop.fence();
  rejects([&] { store.cleanup(); },
          "missing current results cannot be mistaken for no result roots");
  bool resumed = true;
  try {
    store.initialize("a5");
  } catch (const std::exception &) {
    resumed = false;
  }
  check(resumed,
        "missing scientific output does not control numerical restart");
  if (world.rank() == 0)
    check(!fs::exists(run / "results.json"),
          "restart does not recreate missing results");
  const auto deferred = store.publish(
      evaluated("ground_state", "g5", "a5"),
      [&](const fs::path &dir) { return write_payload(world, dir, 88); });
  check(!deferred.empty(),
        "missing results defers cleanup after successful numerical save");
  load(world, store, ref("ground_state", "g5"), 88);
  load(world, store, ref("ground_state", "g2"), 22);
  store.publish_source_map(source_map("ground_state", "a5", "g5", "a5"));
  check(store.publish_results(result("ground_state", "a5")).empty(),
        "validated new results replace missing output and resume cleanup");
  if (world.rank() == 0)
    check(!fs::exists(store.task_directory("ground_state", "state") /
                      "g2/state.00000"),
          "only successful result replacement releases the old source");
  // Same-size corruption is caught before constructing the native input
  // archive.
  std::string original;
  const auto shard =
      store.task_directory("ground_state", "state") / "g5/state.00000";
  if (world.rank() == 0) {
    std::ifstream in(shard, std::ios::binary);
    original.assign(std::istreambuf_iterator<char>(in),
                    std::istreambuf_iterator<char>());
    auto corrupt = original;
    corrupt.back() ^= 1;
    text_file(shard, corrupt);
  }
  world.gop.fence();
  rejects([&] { store.read_latest("ground_state", "state"); },
          "corrupt selected payload fails collectively without fallback");
  if (world.rank() == 0) {
    const auto abandoned =
        store.task_directory("ground_state", "state") / "abandoned";
    fs::create_directory(abandoned);
    text_file(abandoned / "payload", "do not delete with corrupt live source");
  }
  world.gop.fence();
  rejects([&] { store.cleanup(); }, "corrupt live payload prevents cleanup");
  if (world.rank() == 0)
    check(fs::exists(store.task_directory("ground_state", "state") /
                     "abandoned/payload"),
          "failed live validation preserves all candidates");
  rejects([&] { store.initialize("a6"); },
          "corrupt selected payload fails restart preflight");
  rejects(
      [&] { publish(evaluated("ground_state", "g6", "a6"), 99); },
      "corrupt selected payload cannot be replaced without a valid restart");

  if (world.rank() == 0)
    text_file(shard, original);
  world.gop.fence();
}
void workflow_cases(World &world, const fs::path &run) {
  const auto provider = [](const fs::path &root, const cp::json &) {
    const auto plan =
        cp::detail::read_json(root / "finite_difference_plan.json");
    cp::records::require(plan.at("version") == 2,
                         "invalid workflow plan fixture");
    cp::WorkflowRoots roots;
    for (const auto &value : plan.at("sources"))
      roots.generations.insert(cp::records::reference(value));
    for (const auto &value : plan.at("calculations"))
      roots.calculations.insert(value.get<std::string>());
    return roots;
  };
  cp::CheckpointStore store(world, run, provider);
  // Workflow-owned fixture: the provider, not shared storage, interprets this
  // plan.
  if (world.rank() == 0) {
    fs::create_directories(run / "checkpoints");
    text_file(run / "checkpoints/finite_difference_plan.json",
              cp::json{{"version", 2},
                       {"sources", cp::json::array()},
                       {"calculations", cp::json::array()}}
                  .dump());
  }
  world.gop.fence();
  store.initialize("a1");
  store.register_task("ground_state", "scf", "state", "iterative", true);
  store.publish(evaluated("ground_state", "g1"), [&](const fs::path &dir) {
    return write_payload(world, dir, 1);
  });
  if (world.rank() == 0) {
    auto plan = cp::json{{"version", 2},
                         {"sources", inputs("ground_state", "g1")},
                         {"calculations", {"ground_state"}}};
    text_file(run / "checkpoints/finite_difference_plan.json", plan.dump());
    auto inventory = cp::detail::read_json(run / "checkpoints/run.json");
    inventory["finite_difference_plan"] = true;
    text_file(run / "checkpoints/run.json", inventory.dump());
  }
  world.gop.fence();
  store.publish(
      evaluated("ground_state", "g2", "a2"),
      [&](const fs::path &dir) { return write_payload(world, dir, 2); });
  load(world, store, ref("ground_state", "g1"), 1);
  cp::CheckpointStore no_provider(world, run);
  rejects([&] { no_provider.cleanup(); },
          "registered plan without validated provider prevents cleanup");
  rejects([&] { store.retire_calculation("ground_state"); },
          "active plan prevents retiring a required calculation");
  if (world.rank() == 0)
    text_file(run / "checkpoints/finite_difference_plan.json", "broken");
  world.gop.fence();
  rejects([&] { store.cleanup(); },
          "unreadable required workflow plan prevents cleanup");
  load(world, store, ref("ground_state", "g1"), 1);
  if (world.rank() == 0)
    text_file(run / "checkpoints/finite_difference_plan.json",
              cp::json{{"version", 2},
                       {"sources", cp::json::array()},
                       {"calculations", cp::json::array()}}
                  .dump());
  world.gop.fence();
  store.cleanup();
  if (world.rank() == 0)
    check(!fs::exists(store.task_directory("ground_state", "state") /
                      "g1/state.00000"),
          "published plan release permits old payload reclamation");
  store.publish_source_map(source_map("ground_state", "a2", "g2", "a2"));
  if (world.rank() == 0)
    fs::create_directory(run / ".results.json.tmp");
  world.gop.fence();
  rejects([&] { store.publish_results(result("ground_state", "a2")); },
          "failed result publication preserves existing snapshot");
  if (world.rank() == 0)
    check(
        cp::detail::read_json(run / "results.json").at("calculations").empty(),
        "old result is intact after replacement failure");
  if (world.rank() == 0)
    fs::remove(run / ".results.json.tmp");
  world.gop.fence();
  store.publish_results(result("ground_state", "a2"));
}

void corrupt_replacement_case(World &world, const fs::path &run) {
  cp::CheckpointStore store(world, run);
  store.initialize("a1");
  store.register_task("ground_state", "scf", "state", "iterative", true);
  store.publish(evaluated("ground_state", "g1"), [&](const fs::path &dir) {
    return write_payload(world, dir, 1);
  });
  store.publish_source_map(source_map("ground_state", "a1", "g1"));
  store.publish_results(result("ground_state", "a1"));
  store.publish(
      evaluated("ground_state", "g2", "a2"),
      [&](const fs::path &dir) { return write_payload(world, dir, 2); });
  store.publish_source_map(source_map("ground_state", "a2", "g2", "a2"));
  if (world.rank() == 0) {
    const auto shard =
        store.task_directory("ground_state", "state") / "g2/state.00000";
    std::ifstream in(shard, std::ios::binary);
    std::string contents{std::istreambuf_iterator<char>(in),
                         std::istreambuf_iterator<char>()};
    contents.back() ^= 1;
    text_file(shard, contents);
  }
  world.gop.fence();
  rejects(
      [&] {
        store.publish_source_map(source_map("ground_state", "a3", "g2", "a2"));
      },
      "corrupt payload cannot become new source-map evidence");
  rejects([&] { store.publish_results(result("ground_state", "a2")); },
          "corrupt replacement cannot replace current result");
  load(world, store, ref("ground_state", "g1"), 1);
  if (world.rank() == 0)
    check(cp::detail::read_json(run / "results.json").at("attempt_id") == "a1",
          "failed corrupt replacement preserves old results and valid source");
}

} // namespace
int main(int argc, char **argv) {
  World &world = madness::initialize(argc, argv);
  std::string root;
  if (world.rank() == 0) {
    char pattern[] = "/tmp/madness-shared-store-XXXXXX";
    const auto *p = mkdtemp(pattern);
    if (p)
      root = p;
  }
  world.gop.broadcast_serializable(root, 0);
  if (root.empty()) {
    madness::finalize();
    return 1;
  }
  try {
    iteration_record_cases();
    exercise(world, root);
    workflow_cases(world, fs::path(root) / "workflow");
    corrupt_replacement_case(world, fs::path(root) / "corruption");
  } catch (const std::exception &e) {
    check(false, e.what());
  }
  world.gop.sum(failures);
  if (world.rank() == 0) {
    std::cout << "shared checkpoint failures: " << failures << '\n';
    fs::remove_all(root);
  }
  madness::finalize();
  return failures ? 1 : 0;
}
