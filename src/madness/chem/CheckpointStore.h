#pragma once

#include <functional>
#include <cerrno>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <madness/chem/StageCheckpoint.h>
#include <memory>
#include <optional>

namespace madness::checkpoint {

inline Manifest iterative_manifest(
    const std::string &calculation_id, const std::string &generation,
    const json &stage, json definition, const json &inputs,
    const json &metric_definitions,
    const std::map<iteration::Metric, std::string> &metric_set,
    const iteration::IterationMetrics &metrics,
    const iteration::Decision &decision, const iteration::Progress &progress) {
  const auto &info = metrics.info;
  const auto &bound = decision.info;
  records::require(progress.state_valid && decision.checkpoint,
                   "sample is not eligible for publication");
  records::require(info == bound,
                   "decision belongs to another evaluated iteration");
  records::require(info.iteration == progress.main_updates,
                   "metrics do not describe the current update count");
  definition["protocol"] = {{"threshold", info.protocol.thresh},
                            {"k", info.protocol.k}};
  definition["context_id"] = info.context_id;
  Manifest manifest;
  manifest.version = 2;
  manifest.generation = generation;
  manifest.record = {
      {"version", 2},
      {"calculation_id", calculation_id},
      {"attempt_id", info.attempt_id},
      {"generation", generation},
      {"stage", stage},
      {"metric_definitions", metric_definitions},
      {"task",
       {{"task_id", info.task_id},
        {"kind", "iterative"},
        {"definition", definition},
        {"inputs", inputs},
        {"outcome",
         {{"complete", decision.action == iteration::Action::complete},
          {"reason", decision.reason}}},
        {"iteration",
         {{"updates_this_attempt", info.iteration},
          {"updates_at_protocol", progress.updates_at_protocol},
          {"protocol_complete",
           decision.action == iteration::Action::finish_protocol ||
               decision.action == iteration::Action::complete},
          {"metrics", records::measurements(metrics, metric_set)}}}}}};
  records::validate(manifest.record);
  return manifest;
}
inline Manifest direct_manifest(const std::string &calculation_id,
                                const std::string &attempt_id,
                                const std::string &generation,
                                const json &stage, const std::string &task_id,
                                const json &definition, const json &inputs) {
  Manifest manifest;
  manifest.version = 2;
  manifest.generation = generation;
  manifest.record = {
      {"version", 2},
      {"calculation_id", calculation_id},
      {"attempt_id", attempt_id},
      {"generation", generation},
      {"stage", stage},
      {"metric_definitions", json::object()},
      {"task",
       {{"task_id", task_id},
        {"kind", "direct"},
        {"definition", definition},
        {"inputs", inputs},
        {"outcome", {{"complete", true}, {"reason", "computed"}}}}}};
  records::validate(manifest.record);
  return manifest;
}
inline iteration::IterationMetrics
iteration_metrics_from_manifest(const Manifest &manifest) {
  records::require(manifest.version == 2, "not a shared checkpoint");
  records::validate(manifest.record);
  const auto &record = manifest.record;
  const auto &task = record.at("task");
  records::require(task.at("kind") == "iterative",
                   "direct tasks have no iteration measurements");
  const auto &definition = task.at("definition");
  const auto &values = task.at("iteration");
  iteration::IterationMetrics result;
  result.info = {task.at("task_id"), record.at("attempt_id"),
                 records::read_protocol(definition.at("protocol")),
                 records::count(values.at("updates_this_attempt")),
                 definition.at("context_id")};
  for (const auto &id : records::subjects(definition)) {
    iteration::StateMetrics subject;
    subject.id = id;
    for (const auto &metric : values.at("metrics")) {
      if (metric.value("subject", id) != id)
        continue;
      iteration::MetricValue measurement;
      const auto status = metric.at("status").get<std::string>();
      if (status == "measured")
        measurement =
            iteration::MetricValue::measured(metric.at("value").get<double>());
      else if (status == "invalid")
        measurement = iteration::MetricValue::invalid(metric.at("reason"));
      else if (status == "not_applicable")
        measurement =
            iteration::MetricValue::not_applicable(metric.at("reason"));
      else
        measurement = iteration::MetricValue::unavailable(metric.at("reason"));
      if (metric.contains("scale")) {
        const auto &scale = metric.at("scale");
        if (scale.at("status") == "measured")
          measurement.scale = iteration::Scale::measured(scale.at("value"));
        else
          measurement.scale.reason = scale.at("reason");
      }
      if (metric.contains("comparison")) {
        const auto &comparison = metric.at("comparison");
        measurement.comparison = records::read_comparison(comparison);
      }
      if (metric.contains("evidence")) {
        const auto &evidence = metric.at("evidence");
        auto &out = measurement.evidence;
        if (evidence.contains("current"))
          out.current = evidence.at("current").get<double>();
        if (evidence.contains("previous"))
          out.previous = evidence.at("previous").get<double>();
        if (evidence.contains("current_norm"))
          out.current_norm = evidence.at("current_norm").get<double>();
        if (evidence.contains("previous_norm"))
          out.previous_norm = evidence.at("previous_norm").get<double>();
        if (evidence.contains("limiting_id"))
          out.limiting_id = evidence.at("limiting_id").get<std::string>();
        if (evidence.contains("components"))
          for (const auto &component : evidence.at("components")) {
            iteration::ComponentEvidence entry{component.at("id"),
                                               component.at("residual_norm"),
                                               std::nullopt};
            if (component.contains("function_norm"))
              entry.function_norm = component.at("function_norm").get<double>();
            out.components.push_back(std::move(entry));
          }
      }
      subject.metrics.emplace(iteration::metric_from_name(metric.at("metric")),
                              std::move(measurement));
    }
    result.states.push_back(std::move(subject));
  }
  return result;
}

struct Payload {
  std::vector<Archive> archives;
  std::string local_error;
};
struct WorkflowRoots {
  std::set<GenerationRef> generations;
  std::set<std::string>
      calculations; // required calculations, including not-yet-published work
};

/// One run's storage coordinator. All calls are collective in the same order.
/// Rank zero serializes mutations and cleanup with one run lock. Concurrent
/// Independent invocations and process-loss recovery are unsupported.
/// Participating task Worlds use an invocation-scoped writer group, as in
/// StageCheckpoint. A stale lock requires manual recovery, never automatic
/// theft. Workflow callbacks run on rank zero under that lock, validate their
/// current plan and return its exact roots. They must not perform MPI or mutate
/// the store.
class CheckpointStore {
  World &world_;
  fs::path run_, root_;
  using RootProvider =
      std::function<WorkflowRoots(const fs::path &, const json &)>;
  RootProvider workflow_roots_;
  std::string writer_group_;
  // Valid only while this process holds the exclusive store lock. No stat-only
  // cache survives an operation, so same-size corruption is checked next time.
  mutable std::map<GenerationRef, Manifest> validated_;
  // Only task Worlds from one invocation share this token. Serialize their
  // ordinary store locks without changing the rejection of outside writers or
  // stale .store-lock directories. The kernel releases this coordination lock
  // if a process exits; the existing store recovery policy remains unchanged.
  struct GroupLock {
    int fd=-1;
    explicit GroupLock(const fs::path& path) {
      fd=::open(path.c_str(),O_CREAT|O_RDWR|O_CLOEXEC|O_NOFOLLOW,0600);
      detail::require(fd>=0,"cannot open task-world publication coordinator");
      struct stat info{};
      if(::fstat(fd,&info)!=0 || !S_ISREG(info.st_mode)) {
        ::close(fd);fd=-1;
        throw std::runtime_error("checkpoint: invalid publication coordinator");
      }
      struct flock request{};request.l_type=F_WRLCK;request.l_whence=SEEK_SET;
      int status;
      do {status=::fcntl(fd,F_SETLKW,&request);} while(status<0 && errno==EINTR);
      if(status<0) {
        ::close(fd);fd=-1;
        throw std::runtime_error("checkpoint: cannot coordinate task-world publication");
      }
    }
    ~GroupLock() {if(fd>=0) ::close(fd);}
  };
  struct Lock {
    fs::path path;
    std::unique_ptr<GroupLock> group;
    Lock(fs::path path,const std::string& writer_group):path(std::move(path)) {
      if(!writer_group.empty()) {
        records::identifier(writer_group);
        group=std::make_unique<GroupLock>(this->path.parent_path()/(".writers-"+writer_group));
      }
      detail::require(fs::create_directory(this->path),
                      "another store writer or stale .store-lock");
    }
    ~Lock() {std::error_code error;fs::remove(path,error);}
  };
  std::unique_ptr<Lock> acquire_lock() {
    auto lock=std::make_unique<Lock>(root_/".store-lock",writer_group_);
    validated_.clear();return lock;
  }
  static json receipt(const Manifest &manifest) {
    return {{"version", 2},
            {"calculation_id", manifest.record.at("calculation_id")},
            {"task_id", manifest.record.at("task").at("task_id")},
            {"generation", manifest.generation},
            {"attempt_id", manifest.record.at("attempt_id")}};
  }
  [[nodiscard]] fs::path receipt_path(const GenerationRef &ref) const {
    records::identifier(ref.generation);
    return detail::member_path(checked_task(ref.calculation_id, ref.task_id),
                               ".published-" + ref.generation + ".json");
  }
  void remember_selected(const GenerationRef &ref, const Manifest &manifest) {
    const auto path = receipt_path(ref);
    const auto expected = receipt(manifest);
    if (fs::exists(fs::symlink_status(path)))
      records::require(detail::read_json(path) == expected,
                       "publication receipt mismatch");
    else
      replace_json(path, expected);
  }
  static void directory(const fs::path &root, const std::string &relative) {
    detail::require(detail::safe_relative(relative), "unsafe directory");
    auto path = root;
    for (const auto &part : fs::path(relative)) {
      path /= part;
      if (!fs::exists(fs::symlink_status(path)))
        fs::create_directory(path);
      detail::require(fs::is_directory(fs::symlink_status(path)),
                      "not a real directory: " + path.string());
    }
  }
  static void replace_json(const fs::path &path, const json &record) {
    records::finite_numbers(record);
    const auto temp =
        path.parent_path() / ("." + path.filename().string() + ".tmp");
    if (fs::is_regular_file(fs::symlink_status(temp)))
      fs::remove(temp); // abandoned, never a recovery source
    detail::write_exclusive(temp, record.dump(2) + "\n");
    detail::require(!fs::is_symlink(fs::symlink_status(path)),
                    "symlink record destination");
    fs::rename(temp, path);
  }
  [[nodiscard]] json inventory() const {
    auto run = detail::read_json(detail::member_path(root_, "run.json"));
    records::require(run.at("version") == 2, "unsupported run inventory");
    records::id(run.at("attempt_id"));
    records::require(run.at("calculations").is_array(),
                     "invalid calculation inventory");
    std::set<std::string> calculations;
    for (const auto &calculation : run.at("calculations")) {
      records::require(
          calculations.insert(records::id(calculation.at("calculation_id")))
              .second,
          "duplicate calculation");
      records::require(
          !calculation.at("application").get<std::string>().empty(),
          "missing application");
      if (calculation.contains("retired"))
        records::require(calculation.at("retired").is_boolean(),
                         "invalid retirement");
      std::set<std::string> tasks;
      records::require(calculation.at("tasks").is_array(),
                       "invalid task inventory");
      for (const auto &task : calculation.at("tasks")) {
        records::require(tasks.insert(records::id(task.at("task_id"))).second,
                         "duplicate task");
        records::require(task.at("kind") == "iterative" ||
                             task.at("kind") == "direct",
                         "invalid registered kind");
        records::require(task.at("checkpoint").is_boolean(),
                         "invalid checkpoint capability");
        records::require(task.at("kind") == "direct" ||
                             task.at("checkpoint") == true,
                         "iterative task needs a checkpoint");
        if(task.contains("retired")) records::require(task.at("retired").is_boolean(),"invalid task retirement");
        records::inputs(task.at("inputs"));
      }
    }
    return run;
  }
  static json &calculation(json &run, const std::string &id) {
    for (auto &calculation : run.at("calculations"))
      if (calculation.at("calculation_id") == id)
        return calculation;
    throw std::runtime_error("checkpoint: calculation is not registered: " +
                             id);
  }
  static const json &calculation(const json &run, const std::string &id) {
    for (const auto &calculation : run.at("calculations"))
      if (calculation.at("calculation_id") == id)
        return calculation;
    throw std::runtime_error("checkpoint: calculation is not registered: " +
                             id);
  }
  static json &task(json &calculation,const std::string &id) {
    for(auto& entry:calculation.at("tasks"))
      if(entry.at("task_id")==id) return entry;
    throw std::runtime_error("checkpoint: task is not registered: "+id);
  }
  static const json &task(const json &calculation, const std::string &id) {
    for (const auto &task : calculation.at("tasks"))
      if (task.at("task_id") == id)
        return task;
    throw std::runtime_error("checkpoint: task is not registered: " + id);
  }
  [[nodiscard]] fs::path checked_task(const std::string &calculation,
                                      const std::string &task) const {
    records::identifier(calculation);
    records::identifier(task);
    return detail::member_path(root_, calculation + "/tasks/" + task);
  }
  [[nodiscard]] std::optional<std::string>
  selected(const std::string &calculation, const std::string &task) const {
    auto pointer = detail::read_json(
        detail::member_path(checked_task(calculation, task), "latest"));
    records::require(pointer.is_object() && pointer.size() == 1 &&
                         pointer.contains("generation"),
                     "invalid latest pointer");
    if (pointer.at("generation").is_null())
      return std::nullopt;
    return records::id(pointer.at("generation"));
  }
  [[nodiscard]] Manifest source(const json &run, const GenerationRef &ref,
                                bool for_read = false) const {
    const auto &calc = calculation(run, ref.calculation_id);
    const auto &registered = task(calc, ref.task_id);
    records::require(registered.at("checkpoint") == true,
                     "source has no numerical checkpoint");
    records::identifier(ref.generation);
    auto cached = validated_.find(ref);
    if (cached == validated_.end()) {
      const auto dir = detail::member_path(
          checked_task(ref.calculation_id, ref.task_id), ref.generation);
      auto m = read_manifest(dir);
      records::require(
          m.version == 2 &&
              m.record.at("calculation_id") == ref.calculation_id &&
              m.record.at("task").at("task_id") == ref.task_id &&
              m.record.at("stage").at("application") ==
                  calc.at("application") &&
              m.record.at("task").at("kind") == registered.at("kind"),
          "source identity mismatch");
      if (selected(ref.calculation_id, ref.task_id) != ref.generation) {
        const auto path = receipt_path(ref);
        records::require(fs::exists(fs::symlink_status(path)),
                         "source generation was never published");
        records::require(detail::read_json(path) == receipt(m),
                         "publication receipt mismatch");
      }
      records::require(
          !m.archives.empty(),
          "registered numerical checkpoint has no archive inventory");
      const auto valid =
          validate_files(dir, m); // storage validity, not reader rank count
      detail::require(valid.valid(), valid.message());
      cached = validated_.emplace(ref, std::move(m)).first;
    }
    if (for_read)
      for (const auto &archive : cached->second.archives)
        records::require(archive.writers <= world_.size(),
                         "reader count is smaller than archive writer count");
    return cached->second;
  }
  [[nodiscard]] Manifest dependency(const json &run,
                                    const GenerationRef &ref) const {
    records::require(
        !calculation(run, ref.calculation_id).value("retired", false),
        "dependency belongs to a retired calculation");
    return source(run, ref);
  }
  [[nodiscard]] fs::path map_path(const std::string &calculation,
                                  const std::string &attempt) const {
    records::identifier(calculation);
    records::identifier(attempt);
    return detail::member_path(root_, calculation + "/attempts/" + attempt +
                                          "/result_sources.json");
  }
  [[nodiscard]] std::set<GenerationRef>
  map_sources(const json &run, const json &map, const std::string &calc,
              const std::string &attempt) const {
    records::require(map.at("version") == 2 &&
                         map.at("calculation_id") == calc &&
                         map.at("attempt_id") == attempt,
                     "source-map identity mismatch");
    records::require(map.at("states").is_object() && map.at("tasks").is_array(),
                     "invalid source-map collections");
    std::set<GenerationRef> sources;
    for (const auto &state : map.at("states")) {
      GenerationRef ref{calc, records::id(state.at("task_id")),
                        records::id(state.at("generation"))};
      const auto m = dependency(run, ref);
      const auto &task = m.record.at("task");
      records::require(
          m.record.at("attempt_id") == state.at("attempt_id") &&
              task.at("kind") == "iterative" &&
              records::count(task.at("iteration").at("updates_this_attempt")) ==
                  records::count(state.at("iteration")) &&
              records::read_protocol(task.at("definition").at("protocol")) ==
                  records::read_protocol(state.at("protocol")) &&
              records::subjects(task.at("definition"))
                  .contains(records::id(state.at("subject"))),
          "result sample binding mismatch");
      sources.insert(ref);
    }
    for (const auto &entry : map.at("tasks")) {
      GenerationRef ref{entry.value("calculation_id", calc),
                        records::id(entry.at("task_id")),
                        records::id(entry.at("generation"))};
      const auto m = dependency(run, ref);
      records::require(m.record.at("task").at("kind") == "direct",
                       "direct result source is iterative");
      if (entry.contains("attempt_id"))
        records::require(m.record.at("attempt_id") == entry.at("attempt_id"),
                         "direct source attempt mismatch");
      sources.insert(ref);
    }
    return sources;
  }
  [[nodiscard]] std::set<GenerationRef>
  result_sources(const json &run, const json &result) const {
    records::require(result.at("version") == 2 &&
                         result.at("calculations").is_array(),
                     "invalid results envelope");
    const auto root_attempt = records::id(result.at("attempt_id"));
    std::set<GenerationRef> roots;
    std::set<std::string> seen;
    for (const auto &entry : result.at("calculations")) {
      const auto calc = records::id(entry.at("calculation_id"));
      records::require(seen.insert(calc).second,
                       "duplicate result calculation");
      const auto attempt = entry.contains("attempt_id")
                               ? records::id(entry.at("attempt_id"))
                               : root_attempt;
      const auto sources = map_sources(
          run, detail::read_json(map_path(calc, attempt)), calc, attempt);
      roots.insert(sources.begin(), sources.end());
    }
    return roots;
  }
  [[nodiscard]] std::set<GenerationRef> live(const json &run) const {
    std::set<GenerationRef> roots;
    for (const auto &calc : run.at("calculations")) {
      if (calc.value("retired", false))
        continue;
      for (const auto &task : calc.at("tasks")) {
        if(task.value("retired",false)) continue;
        if (task.at("checkpoint") == true) {
          const auto generation =
              selected(calc.at("calculation_id"), task.at("task_id"));
          if (generation)
            roots.insert(
                {calc.at("calculation_id"), task.at("task_id"), *generation});
        }
        const auto refs = records::inputs(task.at("inputs"));
        roots.insert(refs.begin(), refs.end());
      }
    }
    // initialize publishes an empty snapshot before run.json. Thereafter
    // absence is an error, never permission to forget current result roots.
    auto result_refs = result_sources(
        run, detail::read_json(detail::member_path(run_, "results.json")));
    roots.insert(result_refs.begin(), result_refs.end());
    if (run.contains("optimizer_electronic") &&
        !run.at("optimizer_electronic").is_null()) {
      const auto id =
          records::id(run.at("optimizer_electronic").at("calculation_id"));
      records::id(run.at("optimizer_electronic").at("attempt_id"));
      records::require(!calculation(run, id).value("retired", false),
                       "selected electronic calculation is retired");
    }
    const bool needs_workflow = run.contains("finite_difference_plan");
    if (needs_workflow) {
      records::require(run.at("finite_difference_plan") == true &&
                           static_cast<bool>(workflow_roots_),
                       "required workflow plan has no validated root provider");
    }
    if (workflow_roots_) {
      const auto extra = workflow_roots_(root_, run);
      roots.insert(extra.generations.begin(), extra.generations.end());
      for (const auto &id : extra.calculations)
        records::require(!calculation(run, id).value("retired", false),
                         "workflow requires a retired calculation");
    }
    DependencyGraph graph;
    // Finite generation graph, visited once per node.
    // NOLINTNEXTLINE(misc-no-recursion)
    const auto visit = [&](const auto &self, const GenerationRef &ref) -> void {
      if (graph.contains(ref))
        return;
      const auto &calc = calculation(run, ref.calculation_id);
      records::require(!calc.value("retired", false),
                       "live source belongs to retired calculation");
      const auto m = dependency(run, ref);
      auto deps = records::inputs(m.record.at("task").at("inputs"));
      auto stage_deps = records::inputs(m.record.at("stage").at("inputs"));
      deps.insert(deps.end(), stage_deps.begin(), stage_deps.end());
      graph.emplace(ref, deps);
      for (const auto &dependency : deps)
        self(self, dependency);
    };
    for (const auto &root : roots)
      visit(visit, root);
    return retained_generations(roots, graph);
  }
  std::string cleanup_root(const json &run) {
    const auto retained =
        live(run); // all reads/validation finish before any deletion
    struct Batch {
      fs::path dir;
      std::vector<std::string> names, abandoned;
      std::set<std::string> protected_names;
    };
    std::vector<Batch> candidates;
    for (const auto &calc : run.at("calculations"))
      for (const auto &task : calc.at("tasks")) {
        if (task.at("checkpoint") != true)
          continue;
        Batch batch;
        batch.dir = checked_task(calc.at("calculation_id"), task.at("task_id"));
        for (const auto &ref : retained)
          if (ref.calculation_id ==
                  calc.at("calculation_id").get<std::string>() &&
              ref.task_id == task.at("task_id").get<std::string>())
            batch.protected_names.insert(ref.generation);
        std::set<std::string> names;
        for (const auto &entry : fs::directory_iterator(batch.dir)) {
          auto name = entry.path().filename().string();
          if (name.starts_with(".deleting-"))
            name = name.substr(10);
          if (detail::generation_name(name) &&
              fs::is_directory(entry.symlink_status()))
            names.insert(name);
        }
        for (const auto &generation : names) {
          if (batch.protected_names.contains(generation))
            continue;
          bool metadata = false, valid = true;
          // Inspect both locations before scheduling deletion. A previous
          // remove_all may have deleted the doomed manifest first, while the
          // retained copy still reserves the published generation identity.
          for (const auto &name : {generation, ".deleting-" + generation}) {
            try {
              const auto path =
                  detail::member_path(batch.dir, name + "/manifest.json");
              if (!fs::exists(fs::symlink_status(path)))
                continue;
              const auto m = from_json_record(detail::read_json(path));
              if (m.version != 2 || m.generation != generation ||
                  m.record.at("calculation_id") != calc.at("calculation_id") ||
                  m.record.at("task").at("task_id") != task.at("task_id")) {
                valid = false;
                break;
              }
              metadata = true;
            } catch (const std::exception &) {
              valid = false;
              break;
            }
          }
          // Unknown/corrupt metadata is preserved for inspection. No manifest
          // in either location means an unpublished write, with no live refs.
          if (valid)
            (metadata ? batch.names : batch.abandoned).push_back(generation);
        }
        candidates.push_back(std::move(batch));
      }
    std::string warnings;
    for (const auto &batch : candidates) {
      const auto warning = detail::remove_generation_directories(
          batch.dir, batch.names, batch.protected_names, true);
      if (!warning.empty())
        warnings += warning + "\n";
      const auto abandoned = detail::remove_generation_directories(
          batch.dir, batch.abandoned, batch.protected_names);
      if (!abandoned.empty())
        warnings += abandoned + "\n";
    }
    return warnings;
  }

public:
  CheckpointStore(World &world, const fs::path &run,
                  RootProvider workflow_roots = {}, std::string writer_group = {})
      : world_(world), run_(fs::absolute(run)), root_(run_ / "checkpoints"),
        workflow_roots_(std::move(workflow_roots)),writer_group_(std::move(writer_group)) {}
  [[nodiscard]] fs::path task_directory(const std::string &calculation,
                                        const std::string &task) const {
    records::identifier(calculation);
    records::identifier(task);
    return root_ / calculation / "tasks" / task;
  }
  void initialize(const std::string &attempt) {
    detail::on_root(world_, [&] {
      records::identifier(attempt);
      fs::create_directories(run_);
      detail::require(fs::is_directory(fs::symlink_status(run_)),
                      "run is not a real directory");
      directory(run_, "checkpoints");
      const auto lock = acquire_lock();
      const bool registered =
          fs::exists(fs::symlink_status(root_ / "run.json"));
      auto run = registered ? inventory()
                            : json{{"version", 2},
                                   {"attempt_id", attempt},
                                   {"calculations", json::array()}};
      if (!registered) {
        const auto results = detail::member_path(run_, "results.json");
        if (fs::exists(fs::symlink_status(results))) {
          const auto initial = detail::read_json(results);
          records::require(initial.at("version") == 2 &&
                               initial.at("calculations").is_array() &&
                               initial.at("calculations").empty(),
                           "unexpected results without run inventory");
        } else
          replace_json(results, {{"version", 2},
                                 {"attempt_id", attempt},
                                 {"calculations", json::array()}});
      }
      // Numerical restart does not read scientific output. Missing output
      // disables cleanup until a valid new result snapshot is published.
      for (const auto &calc : run.at("calculations")) {
        if (calc.value("retired", false))
          continue;
        for (const auto &task : calc.at("tasks")) {
          if (task.at("checkpoint") != true)
            continue;
          const auto generation =
              selected(calc.at("calculation_id"), task.at("task_id"));
          if (generation)
            (void)source(run, {calc.at("calculation_id"), task.at("task_id"),
                               *generation});
        }
      }
      run["attempt_id"] = attempt;
      replace_json(root_ / "run.json", run);
      return std::string{};
    });
  }
  /// Atomic current-attempt timing updates, independent of numerical archives.
  void publish_timings(const json& attempt) {
    detail::on_root(world_, [&] {
      const auto lock = acquire_lock();
      records::id(attempt.at("attempt_id"));
      records::finite_numbers(attempt);
      records::require(attempt.at("execution").is_object() && attempt.at("timing").is_object(),
                       "invalid timing attempt");
      const auto path = detail::member_path(run_, "timings.json");
      auto timings = fs::exists(fs::symlink_status(path)) ? detail::read_json(path)
          : json{{"version",1},{"attempts",json::array()}};
      records::require(timings.at("version")==1 && timings.at("attempts").is_array(),
                       "unsupported timing record");
      auto& attempts = timings.at("attempts");
      auto found = std::find_if(attempts.begin(), attempts.end(), [&](const json& previous) {
        return previous.at("attempt_id")==attempt.at("attempt_id");
      });
      if (found==attempts.end()) attempts.push_back(attempt); else *found=attempt;
      replace_json(path,timings);
      return std::string{};
    });
  }
  void register_task(const std::string &calc, const std::string &application,
                     const std::string &task_id, const std::string &kind,
                     bool checkpoint, const json &inputs = json::array(),
                     bool require_same_inputs = false) {
    detail::on_root(world_, [&] {
      const auto lock = acquire_lock();
      auto run = inventory();
      records::identifier(calc);
      records::identifier(task_id);
      records::require(!application.empty() &&
                           (kind == "iterative" || kind == "direct") &&
                           (checkpoint || kind == "direct"),
                       "invalid task registration");
      for (const auto &ref : records::inputs(inputs))
        (void)dependency(run, ref);
      auto found =
          std::find_if(run["calculations"].begin(), run["calculations"].end(),
                       [&](const json &value) {
                         return value.at("calculation_id") == calc;
                       });
      if (found == run["calculations"].end()) {
        run["calculations"].push_back({{"calculation_id", calc},
                                       {"application", application},
                                       {"tasks", json::array()}});
        found = std::prev(run["calculations"].end());
      }
      auto &calculation = *found;
      records::require(!calculation.value("retired", false) &&
                           calculation.at("application") == application,
                       "incompatible or retired calculation");
      auto &tasks = calculation["tasks"];
      auto existing =
          std::find_if(tasks.begin(), tasks.end(), [&](const json &value) {
            return value.at("task_id") == task_id;
          });
      if (existing != tasks.end()) {
        records::require(existing->at("kind") == kind &&
                             existing->at("checkpoint") == checkpoint,
                         "task registration changed kind");
        // Adapters without a source-conversion path can refuse replacement
        // under the registration lock, before changing provenance or retention.
        records::require(!require_same_inputs || existing->at("inputs") == inputs,
                         "task registration changed inputs; source conversion is unsupported");
        if (checkpoint) {
          const auto generation = selected(calc, task_id); // never recreate
          if (generation)
            (void)source(run, {calc, task_id, *generation});
        }
        (*existing)["inputs"] = inputs;
        existing->erase("retired");
      } else {
        if (checkpoint) {
          directory(root_, calc + "/tasks/" + task_id);
          const auto pointer = checked_task(calc, task_id) / "latest";
          if (fs::exists(fs::symlink_status(pointer)))
            records::require(!selected(calc, task_id),
                             "unexpected published pointer for new task");
          else
            replace_json(pointer, {{"generation", nullptr}});
        }
        tasks.push_back({{"task_id", task_id},
                         {"kind", kind},
                         {"checkpoint", checkpoint},
                         {"inputs", inputs}});
      }
      replace_json(root_ / "run.json", run);
      return std::string{};
    });
  }
  std::optional<Manifest> read_latest(const std::string &calc,
                                      const std::string &task_id,
                                      bool allow_unregistered = false) {
    const auto record = detail::on_root(world_, [&] {
      const auto lock = acquire_lock();
      const auto run = inventory();
      if (allow_unregistered) {
        const auto& calculations = run.at("calculations");
        const auto found = std::find_if(calculations.begin(), calculations.end(), [&](const json& c) {
          return c.at("calculation_id") == calc;
        });
        if (found == calculations.end()) return std::string("null");
        const auto& tasks = found->at("tasks");
        if (std::none_of(tasks.begin(), tasks.end(), [&](const json& t) { return t.at("task_id") == task_id; }))
          return std::string("null");
      }
      const auto &c = calculation(run, calc);
      records::require(!c.value("retired", false) &&
                           task(c, task_id).at("checkpoint") == true && !task(c,task_id).value("retired",false),
                       "task is retired or has no checkpoint");
      const auto generation = selected(calc, task_id);
      return generation ? to_json_record(
                              source(run, {calc, task_id, *generation}, true))
                              .dump()
                        : std::string("null");
    });
    const auto parsed = json::parse(record);
    if (parsed.is_null())
      return std::nullopt;
    return from_json_record(parsed);
  }
  Manifest read(const GenerationRef &ref) {
    return from_json_record(json::parse(detail::on_root(world_, [&] {
      const auto lock = acquire_lock();
      return to_json_record(source(inventory(), ref, true)).dump();
    })));
  }
  template <class Writer>
  std::string publish(Manifest manifest, Writer &&writer) {
    std::unique_ptr<Lock> lock;
    const json publication = json::parse(detail::on_root(world_, [&] {
      lock = acquire_lock();
      const auto run = inventory();
      records::require(manifest.version == 2,
                       "shared store requires version-2 evidence");
      records::validate(manifest.record);
      const auto calc = records::id(manifest.record.at("calculation_id"));
      const auto task_id =
          records::id(manifest.record.at("task").at("task_id"));
      const auto &c = calculation(run, calc);
      const auto &t = task(c, task_id);
      records::require(
          !c.value("retired", false) && t.at("checkpoint") == true &&
              t.at("kind") == manifest.record.at("task").at("kind") &&
              c.at("application") ==
                  manifest.record.at("stage").at("application"),
          "checkpoint does not match registration");
      const auto generation = selected(calc, task_id);
      if (generation) {
        const GenerationRef previous{calc, task_id, *generation};
        // The current pointer proves publication. Preserve that proof before
        // advancing it; a manifest from a failed save proves nothing alone.
        remember_selected(previous, source(run, previous));
      }
      detail::structure(manifest);
      for (const auto *list : {&manifest.record.at("task").at("inputs"),
                               &manifest.record.at("stage").at("inputs")})
        for (const auto &ref : records::inputs(*list))
          (void)dependency(run, ref);
      const auto dir = checked_task(calc, task_id);
      detail::require(!fs::exists(fs::symlink_status(
                          dir / (".deleting-" + manifest.generation))),
                      "generation is being reclaimed");
      detail::require(fs::create_directory(dir / manifest.generation),
                      "generation already exists");
      return json{{"directory", dir.string()},
                  {"manifest", to_json_record(manifest)}}
          .dump();
    }));
    const auto location = publication.at("directory").get<std::string>();
    manifest = from_json_record(publication.at("manifest"));
    Payload payload;
    try {
      payload = writer(fs::path(location) / manifest.generation);
    } catch (const std::exception &error) {
      payload.local_error = error.what();
    } catch (...) {
      payload.local_error = "checkpoint: unknown payload writer error";
    }
    detail::collective_error(world_, payload.local_error);
    if (world_.rank() == 0)
      manifest.archives = std::move(payload.archives);
    detail::on_root(world_, [&] {
      detail::require(!manifest.archives.empty(),
                      "numerical checkpoint has no archives");
      return std::string{};
    });
    checkpoint::publish(world_, location, manifest);
    // A cleanup failure after publication must not turn a successful save
    // into a false failure; report it and retain extra archives instead.
    return detail::on_root(world_, [&] {
      try {
        return cleanup_root(inventory());
      } catch (const std::exception &error) {
        return std::string("cleanup deferred: ") + error.what();
      }
    });
  }
  void publish_source_map(const json &map) {
    detail::on_root(world_, [&] {
      const auto lock = acquire_lock();
      const auto run = inventory();
      const auto calc = records::id(map.at("calculation_id"));
      const auto attempt = records::id(map.at("attempt_id"));
      (void)map_sources(run, map, calc, attempt);
      directory(root_, calc + "/attempts/" + attempt);
      const auto path = map_path(calc, attempt);
      if (fs::exists(fs::symlink_status(path)))
        records::require(detail::read_json(path) == map,
                         "immutable source map cannot be rebound");
      else
        replace_json(path, map);
      return std::string{};
    });
  }
  std::string publish_results(const json &result) {
    return detail::on_root(world_, [&] {
      const auto lock = acquire_lock();
      const auto run = inventory();
      (void)result_sources(run, result);
      replace_json(run_ / "results.json", result);
      try {
        return cleanup_root(run);
      } catch (const std::exception &error) {
        return std::string("cleanup deferred: ") + error.what();
      }
    });
  }
  /// Small manifest inventory for application-owned explicit task lineage.
  json task_inventory(const std::string& calc) {
    return json::parse(detail::on_root(world_,[&] {
      const auto run=inventory();
      for(const auto& entry:run.at("calculations"))
        if(entry.at("calculation_id")==calc) return entry.at("tasks").dump();
      return std::string("[]");
    }));
  }
  /// Release an obsolete task's current selection. Exact result/dependency
  /// references still retain its payload through the normal retention graph.
  void retire_task(const std::string& calc,const std::string& id) {
    detail::on_root(world_,[&] {
      const auto lock=acquire_lock();auto run=inventory();(void)live(run);
      auto& entry=task(calculation(run,calc),id);
      if(entry.value("retired",false)) return std::string{};
      if(entry.at("checkpoint")==true) {
        if(const auto generation=selected(calc,id))
          remember_selected({calc,id,*generation},source(run,{calc,id,*generation}));
        replace_json(checked_task(calc,id)/"latest",{{"generation",nullptr}});
      }
      entry["retired"]=true;replace_json(root_/"run.json",run);
      return std::string{};
    });
  }
  void retire_calculation(const std::string &id) {
    detail::on_root(world_, [&] {
      const auto lock = acquire_lock();
      auto run = inventory();
      (void)live(run);
      auto &c = calculation(run, id);
      c["retired"] = true;
      (void)live(run);
      replace_json(root_ / "run.json", run);
      return std::string{};
    });
  }
  std::string cleanup() {
    return detail::on_root(world_, [&] {
      const auto lock = acquire_lock();
      return cleanup_root(inventory());
    });
  }
};
} // namespace madness::checkpoint
