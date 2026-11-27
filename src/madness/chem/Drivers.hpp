
#pragma once

// Workflow.hpp
// Defines the Driver interface, a generic SinglePointDriver, and the Workflow
// class

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
#include <utility>
#include <nlohmann/json.hpp>
#include <memory>
#include <stdexcept>
#include <system_error>
#include <vector>

#include <madness/chem/Applications.hpp>  // Interface for SCFApplication / ResponseApplication
#include <madness/chem/RunFiles.h>        // write_results_entry
#include <madness/chem/SCFTargetAdapter.hpp>  // SCFTarget

namespace qcapp {

/// A controlled incomplete execution, reported without losing its results.
struct WorkflowStopped : std::runtime_error {
  using std::runtime_error::runtime_error;
};

/**
 * @brief Abstract base class for all drivers that encapsulate one or more
 * Applications.
 */
class Driver {
 public:
  virtual ~Driver() = default;

  virtual void print_parameters(World& world) const =0;

  /**
   * @brief Execute the driver, writing outputs under the given directory.
   * @param workdir Base directory for this driver's outputs.
   * @param ctx     Shared StepContext threaded task-to-task; a driver reads
   *                artifacts published by upstream steps and publishes its own.
   */
  virtual void execute(const std::filesystem::path& workdir,
                       madness::StepContext& ctx) = 0;

  /**
   * @brief Return a JSON summary of results produced by this driver.
   */
  [[nodiscard]] virtual nlohmann::json summary() const = 0;

  /// The results.json entry this driver's calculation writes in a shared run
  /// ({calculation_id, application}), if any.
  [[nodiscard]] virtual std::optional<std::pair<std::string, std::string>> results_entry() const {
    return std::nullopt;
  }
};

/**
 * @brief Runs a single Application (e.g. SCF or Response) in its own
 * subdirectory.
 */
class SinglePointDriver : public Driver {
 public:
  explicit SinglePointDriver(std::shared_ptr<Application> app) : app_(app) {}

  void print_parameters(World& world) const override {
    app_->print_parameters(world);
  }

  void execute(const std::filesystem::path& workdir,
               madness::StepContext& ctx) override {
    // Create workdir for this application
    std::filesystem::create_directories(workdir);

    // Read upstream artifacts (no-op unless the app overrides), run, then
    // publish this step's artifacts for downstream steps.
    app_->set_shared_run_root(ctx.shared_run_root);
    app_->consume_context(ctx);
    app_->run(workdir);
    result_ = app_->results();
    if (const auto outcome = app_->outcome()) {
      result_["outcome"] = {{"complete", outcome->complete}, {"reason", outcome->reason}};
      if (!outcome->complete)
        throw WorkflowStopped("Calculation stopped incomplete: " + outcome->reason);
    }
    app_->publish_to_context(ctx);
  }

  nlohmann::json summary() const override {
    return result_;
  }

  [[nodiscard]] std::optional<std::pair<std::string, std::string>> results_entry() const override {
    return app_->results_entry();
  }

 private:
  std::shared_ptr<Application> app_;
  nlohmann::json result_;
};

/**
 * @brief Orchestrates multiple drivers in sequence and writes a global
 * output.json.
 */
class Workflow {
public:
  Workflow() = default;

  /**
   * @brief Add a driver to the workflow.
   * @param driver Unique pointer to a Driver instance.
   */
  void addDriver(std::unique_ptr<Driver> driver) { drivers_.push_back(std::move(driver)); }

  /// Run-level provenance stamped into EVERY calc_info write (e.g. the
  /// effective restart-I/O configuration — roadmap change 5: io provenance
  /// for all tasks, not just response). Recorded here and re-applied at the
  /// top of run() (which resets the aggregate), so it is present even in the
  /// partial calc_info a failed task leaves behind.
  void set_provenance(const std::string &key, nlohmann::json value) {
    provenance_[key] = std::move(value);
  }

  void print_parameters(World &world) const {
    for (const auto &d : drivers_)
      d->print_parameters(world);
  }

  /**
   * @brief Run all added drivers under the top-level directory, then emit
   * output.json.
   * @param topDir Root directory for the entire workflow.
   * @param outputfile Name of the output file to write the aggregated results.
   */
  void run(const std::string prefix, bool shared_output = false) {
    std::filesystem::path topDir = prefix;
    std::filesystem::create_directories(topDir);
    all_ = nlohmann::json::object();
    for (const auto &kv : provenance_.items()) all_[kv.key()] = kv.value();
    all_["tasks"] = nlohmann::json::array();

    // One StepContext threaded through the whole chain: each driver reads what
    // upstream steps published and publishes its own artifacts (roadmap 1).
    madness::StepContext ctx;
    if (shared_output) ctx.shared_run_root = std::filesystem::absolute(topDir);

    for (size_t i = 0; i < drivers_.size(); ++i) {
      auto taskDir = shared_output ? topDir / "analysis"
                                  : topDir / ("task_" + std::to_string(i));
      try {
        drivers_[i]->execute(taskDir, ctx);
      } catch (const WorkflowStopped &) {
        append_summary(drivers_[i]->summary());
        write_calc_info(prefix);
        throw;
      } catch (...) {
        // A failing task must still leave a complete, diagnosable
        // calc_info.json: record the failure as a task entry, persist, and
        // rethrow so the app-level handler reports the error (madqc.cpp).
        nlohmann::json failed;
        failed["type"] = "task_failed";
        failed["task_index"] = i;
        try {
          throw;
        } catch (const std::exception &e) {
          failed["error"] = e.what();
        } catch (...) {
          failed["error"] = "unknown exception";
        }
        all_["tasks"].push_back(failed);
        write_calc_info(prefix);
        record_failure(ctx, *drivers_[i], failed);
        throw;
      }
      append_summary(drivers_[i]->summary());
      write_calc_info(prefix);
    }
  }

  /// Aggregated results of all drivers (the calc_info JSON); valid after run().
  const nlohmann::json &results() const { return all_; }

private:
  /// A failed calculation's entry in the shared run's results.json: status
  /// "failed" and what calc_info.json records for it. It replaces any entry an
  /// earlier launch left under the same id. A failure to write it must not
  /// hide the error being reported, so it is only printed.
  static void record_failure(const madness::StepContext &ctx, const Driver &driver,
                             const nlohmann::json &failed) {
    const auto id = driver.results_entry();
    if (ctx.shared_run_root.empty() || !id) return;
    try {
      madness::iteration::write_results_entry(
          madness::World::get_default(), ctx.shared_run_root,
          {{"calculation_id", id->first}, {"application", id->second}, {"status", "failed"},
           {"error", {{"task_index", failed.at("task_index")}, {"message", failed.at("error")}}}});
    } catch (const std::exception &e) {
      if (madness::World::get_default().rank() == 0)
        std::cerr << "could not record the failure in results.json: " << e.what() << std::endl;
    }
  }

  void append_summary(const nlohmann::json &output) {
    if (output.is_array()) {
      for (const auto &item : output) all_["tasks"].push_back(item);
    } else {
      all_["tasks"].push_back(output);
    }
  }

  /// Rank-0-only, ATOMIC (tmp+rename) aggregate write. Every rank used to
  /// stream the file directly — N concurrent writers to one path on a shared
  /// FS, and a crash mid-write truncated it.
  ///
  /// Two additional hazards handled here (raman thread brief, defect 2):
  ///  (a) an unchecked ofstream means an ENOSPC-truncated tmp still SUCCEEDS
  ///      the rename and replaces a good calc_info.json — so we verify the
  ///      stream and the rename before committing, and drop the tmp on failure;
  ///  (b) a rank-0-only filesystem failure would unwind only rank 0 while the
  ///      other ranks march into the next driver's collectives → a hang that
  ///      burns the whole allocation. We broadcast an ok/fail flag so every
  ///      rank throws together (or none does).
  void write_calc_info(const std::string &prefix) const {
    auto &world = madness::World::get_default();
    std::string error;
    if (world.rank() == 0) {
      const std::string outputfile = prefix + ".calc_info.json";
      const std::string tmpfile = outputfile + ".tmp";
      try {
        std::ofstream ofs(tmpfile);
        ofs.exceptions(std::ios::failbit | std::ios::badbit);
        ofs << std::setw(4) << all_;
        ofs.close();
        std::filesystem::rename(tmpfile, outputfile);
      } catch (const std::exception &exception) {
        error = exception.what();
        std::error_code ec;
        std::filesystem::remove(tmpfile, ec);
      }
    }
    world.gop.broadcast_serializable(error, 0);
    if (!error.empty())
      throw std::runtime_error(
          "write_calc_info: failed to persist " + prefix +
          ".calc_info.json: " + error);
  }

  std::vector<std::unique_ptr<Driver>> drivers_;
  nlohmann::json all_;
  nlohmann::json provenance_ = nlohmann::json::object();
};

/**
 * @brief Geometry optimization as its own task (madqc ARCHITECTURE_ROADMAP change 2).
 *
 * Drives MolOpt over the reference engine the Library policy supplies, so the same
 * code optimizes on moldft (`Calc = SCF`) and nemo (`Calc = Nemo`). Selected by
 * `--optimize` together with `--wf=<scf|nemo>`, which names the reference method.
 *
 * A Driver rather than an Application, deliberately: numerical gradients are to be
 * computed from displaced sub-runs, each with its own directory and calc_info (see
 * GeometryTarget), and owning sub-runs is a Driver's job. The same seam serves
 * roadmap changes 3 and 4.
 *
 * Not restartable as a step: the optimizer keeps no checkpoint of its own (hessian
 * and step history), so an interrupted optimization restarts from the input
 * geometry. The SCF underneath still restarts per geometry as usual.
 */
template <typename Library> class OptimizeDriver : public Driver {
public:
  using Calc = typename Library::Calc;

  OptimizeDriver(World &world, const Params &params)
      : world_(world), params_(params) {}

  void print_parameters(World &world) const override {
    if (world.rank() != 0)
      return;
    params_.get<OptimizationParameters>().print(OptimizationParameters::tag,
                                                "end");
    if constexpr (std::is_same_v<Calc, madness::SCF>) {
      params_.get<CalculationParameters>().print(CalculationParameters::tag,
                                                 "end");
    } else {
      params_.get<CalculationParameters>().print(CalculationParameters::tag);
      params_.get<madness::Nemo::NemoCalculationParameters>().print();
      madness::print("end");
    }
  }

  void execute(const std::filesystem::path &workdir,
               madness::StepContext &ctx) override {
    // An upstream step may have moved the molecule; honour it before the engine is
    // built, since construction freezes the geometry.
    if (ctx.molecule && ctx.molecule->natom() > 0)
      params_.get<madness::Molecule>() = *ctx.molecule;

    madness::PathManager pm(workdir, Library::label());
    pm.create();
    world_.gop.fence();
    {
      madness::ScopedCWD scwd(pm.dir());
      if (world_.rank() == 0)
        madness::print("Running geometry optimization on", Library::label(), "in",
                       pm.dir().string());

      // Keep the existing protocol-derived energy and step defaults. Gradient
      // defaults are independent of SCF's density tolerance; user overrides
      // still win and MolOpt retains its precision clamp.
      auto &op = params_.get<OptimizationParameters>();
      const auto &cp = params_.get<CalculationParameters>();
      op.set_protocol_defaults(cp.protocol().back());

      MADNESS_CHECK_THROW(
          !op.get_initial_hessian(),
          "optimization: initial_hessian is not implemented for the optimizer -- "
          "MolOpt starts from a mass-weighted diagonal guess. Remove the key "
          "rather than have it silently ignored");

      auto engine = lib_.calc(world_, params_);
      engine->work_dir = pm.dir();

      if constexpr (std::is_same_v<Calc, madness::SCF>) {
        // The same preparation moldft_lib::run does before handing the engine to
        // MolOpt. It is not optional: MolecularEnergy::value only calls
        // set_protocol when FunctionDefaults' thresh differs from protocol[0], so
        // when they already agree nothing would build the nuclear potential / data
        // map and the first geometry segfaults in the initial guess. Nemo needs no
        // equivalent -- Nemo::value sets its own protocol on every call.
        // `template` disambiguator: Calc is a template parameter here.
        if (world_.size() > 1) {
          engine->template set_protocol<3>(world_, 1e-4);
          engine->make_nuclear_potential(world_);
          engine->initial_load_bal(world_);
        }
        engine->template set_protocol<3>(world_,
                                         engine->param.protocol()[0]);
      } else {
        // Nemo carries work_dir both on itself and on its inner SCF; downstream
        // consumers read one or the other, so set both.
        engine->get_calc()->work_dir = pm.dir();
      }

      // The engine as an OptimizationTargetInterface: MolecularEnergy wraps an SCF,
      // Nemo is one itself. `scf_target` must outlive the optimization -- it holds
      // an SCF&.
      std::unique_ptr<madness::MolecularEnergy> scf_target;
      madness::OptimizationTargetInterface *engine_target = nullptr;
      if constexpr (std::is_same_v<Calc, madness::SCF>) {
        scf_target = std::make_unique<madness::MolecularEnergy>(world_, *engine);
        engine_target = scf_target.get();
      } else {
        engine_target = engine.get();
      }
      MADNESS_CHECK_THROW(engine_target->provides_gradient(),
                          "the reference engine provides no analytic gradient; "
                          "numerical gradients from displaced sub-runs are not "
                          "implemented yet (see GeometryTarget)");

      // THE SEAM: swap in a displaced-sub-run target here to get numerical
      // gradients; everything below is unchanged by that choice.
      madness::AnalyticTarget target(*engine_target);

      // Keep the raw derivative table out of the log: MolOpt now prints the
      // projected gradient it actually uses, and the two side by side (differing
      // by ~100x) are what made the projection look like it was not happening.
      if constexpr (std::is_same_v<Calc, madness::SCF>)
        engine->suppress_raw_gradient_print = true;
      else
        engine->get_calc()->suppress_raw_gradient_print = true;

      // The deck's print_level drives the optimizer too, rather than a hardcoded
      // 1 that no input could raise. The mapping keeps the default (print_level
      // 3 -> 1) exactly as it was and still shows the iteration table at 2;
      // only 4 and above turn on the projector/hessian dumps.
      const int molopt_print =
          (cp.print_level() >= 2) ? std::max(1, cp.print_level() - 2) : 0;
      madness::MolOpt opt(op.get_maxiter(), op.get_maxstep(), op.get_etol(),
                          op.get_gtol(), op.get_xtol(), op.get_value_precision(),
                          op.get_gradient_precision(),
                          (world_.rank() == 0) ? molopt_print : 0,
                          op.get_algopt());

      madness::OptimizationResults opt_res;
      if constexpr (std::is_same_v<Calc, madness::SCF>) {
        opt_res = opt.optimize_app(engine->molecule, target);
      } else {
        opt_res = opt.optimize_app(engine->molecule(), target);
      }

      // Leave the engine AT the optimized geometry and pick up the final energy and
      // gradient there (mirrors the in-SCF gopt path).
      double energy = 0.0;
      madness::Tensor<double> gradient;
      target.energy_and_gradient(opt_res.final_geometry, energy, gradient);
      opt_res.final_energy = energy;

      madness::PropertyResults prop_res;
      prop_res.energy = energy;
      prop_res.gradient = gradient;

      summary_["model"] = "optimize";
      summary_["optimization_results"] = opt_res.to_json();
      summary_["molecule"] = opt_res.final_geometry.to_json();
      summary_["properties"] = prop_res.to_json();
      summary_["metadata"] = {{"mpi_size", world_.size()},
                              {"method", Library::label()}};

      if (world_.rank() == 0) {
        const std::string geomfile =
            params_.get<CalculationParameters>().prefix() + "_opt.xyz";
        std::ofstream ofs(geomfile);
        opt_res.final_geometry.print(ofs);
        ofs.close();
        madness::print("optimized geometry written to", geomfile);
        opt_res.final_geometry.print();
      }
      final_geometry_ = opt_res.final_geometry;
    }

    // Hand the optimized geometry to whatever runs next -- the point of making this
    // a first-class step.
    if (final_geometry_.natom() > 0)
      ctx.molecule = final_geometry_;
    if constexpr (std::is_same_v<Calc, madness::SCF>)
      ctx.reference = lib_.calc(world_, params_);
    try {
      const auto &cp = params_.get<CalculationParameters>();
      ctx.archives["restartdata"] = pm.dir() / (cp.prefix() + ".restartdata");
    } catch (...) {
      // best effort, as in SCFApplication
    }
  }

  nlohmann::json summary() const override { return summary_; }

private:
  World &world_;
  Params params_;
  Library lib_;
  nlohmann::json summary_;
  madness::Molecule final_geometry_;
};

} // namespace qcapp
