/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680


  $Id$
 */
#pragma once

#include "ResponseParameters.hpp"
#include "madness_exception.h"
#include <apps/molresponse_v2/DerivedStatePlanner.hpp>
#include <apps/molresponse_v2/FrequencyLoop.hpp>
#include <apps/molresponse_v2/GroundStateData.hpp>
#include <apps/molresponse_v2/PropertyManager.hpp>
#include <apps/molresponse_v2/ResponseDebugLogger.hpp>
#include <apps/molresponse_v2/ResponseManager.hpp>
#include <apps/molresponse_v2/ResponseRecord.hpp>
#include <apps/molresponse_v2/StateParallelPlanner.hpp>
#include <apps/molresponse_v2/StateGenerator.hpp>
#include <filesystem>
#include <madness/chem/InputWriter.hpp>
#include <madness/chem/ParameterManager.hpp>
#include <madness/chem/Results.h>

struct molresponse_lib {
  struct Results {
    nlohmann::json metadata;
    nlohmann::json properties;
    nlohmann::json vibrational_analysis;
    nlohmann::json raman_spectra;
    nlohmann::json debug_log;
  };

  static constexpr const char *label() { return "molresponse"; }

private:
  struct GroundContext {
    Molecule molecule;
    GroundStateData ground;
    ResponseManager response_manager;
    std::string fock_json_file;
  };

  struct PlannedStates {
    GeneratedStateData generated_states;
    DerivedStatePlan derived_state_plan;
    StateParallelPlan state_parallel_plan;
  };

  struct SolvedStates {
    PlannedStates planned_states;
    nlohmann::json metadata;
    nlohmann::json debug_log;
  };

  struct PropertyStageOutput {
    nlohmann::json properties;
    VibrationalResults vibrational_analysis;
    RamanResults raman_spectra;
  };

  class JsonStateSolvePersistence final : public StateSolvePersistence {
  public:
    JsonStateSolvePersistence(World &world, const std::string &meta_file,
                              const std::string &debug_file)
        : response_record_(world, meta_file), debug_logger_(debug_file) {}

    void initialize_states(const std::vector<LinearResponseDescriptor> &states) {
      response_record_.initialize_states(states);
    }

    void print_summary() const { response_record_.print_summary(); }

    [[nodiscard]] bool is_saved(const LinearResponsePoint &pt) const override {
      return response_record_.is_saved(pt.perturbationDescription(),
                                       pt.threshold(), pt.frequency());
    }

    [[nodiscard]] bool
    is_converged(const LinearResponsePoint &pt) const override {
      return response_record_.is_converged(pt.perturbationDescription(),
                                           pt.threshold(), pt.frequency());
    }

    void record_status(const LinearResponsePoint &pt, bool c) override {
      response_record_.record_status(pt, c);
    }

    ResponseDebugLogger &logger() override { return debug_logger_; }

    void flush_debug_log(World &world) override {
      if (debug_logger_.enabled() && world.rank() == 0) {
        debug_logger_.write_to_disk();
      }
    }

    [[nodiscard]] nlohmann::json metadata_json() const {
      return response_record_.to_json();
    }

    [[nodiscard]] nlohmann::json debug_log_json() const {
      return debug_logger_.to_json();
    }

  private:
    ResponseRecord2 response_record_;
    ResponseDebugLogger debug_logger_;
  };

  static GroundContext
  make_ground_context(World &world, const CalculationParameters &calc_params,
                      const std::shared_ptr<SCF> &scf_calc,
                      const std::filesystem::path &outdir) {
    auto indir = scf_calc->work_dir;
    auto rel = std::filesystem::relative(indir, outdir);
    auto prox = std::filesystem::proximate(indir, outdir);

    if (world.rank() == 0) {
      print("Running MolresponseLib::run_response() in directory: ", outdir);
      print("Ground state archive: ", indir);
      print("Relative path: ", rel);
      print("Proximate path: ", prox);
    }

    const auto &prefix = calc_params.prefix();
    std::string archive_name = prefix + ".restartdata";
    std::string fock_json_file = (prox / (prefix + ".fock.json")).string();
    std::string moldft_checkpt = prox / "moldft.calc_info.json";
    auto relative_archive = prox / archive_name;

    if (!std::filesystem::exists(moldft_checkpt)) {
      if (world.rank() == 0) {
        print("Error: Missing ground-state checkpoint file: ", moldft_checkpt);
      }
      throw std::runtime_error("Missing ground-state checkpoint file");
    }

    auto read_molecule = [](const std::string &scf_ckpt) -> Molecule {
      std::ifstream ifs(scf_ckpt);
      nlohmann::json j;
      ifs >> j;
      ifs.close();
      Molecule mol;
      if (j.contains("molecule")) {
        mol.from_json(j["molecule"]);
      } else {
        throw std::runtime_error(
            "Molecule information missing from checkpoint JSON.");
      }
      return mol;
    };

    Molecule molecule = read_molecule(moldft_checkpt);
    print("Read molecule with ", molecule.natom(), " atoms.");

    GroundStateData ground(world, relative_archive.string(), molecule);
    ResponseManager response_manager(world, calc_params);

    return GroundContext{std::move(molecule), std::move(ground),
                         std::move(response_manager),
                         std::move(fock_json_file)};
  }

  static PlannedStates
  plan_required_states(World &world, const CalculationParameters &calc_params,
                       const GroundContext &ctx,
                       const ResponseParameters &response_params) {
    StateGenerator state_generator(ctx.molecule, calc_params.protocol(),
                                   ctx.ground.isSpinRestricted(),
                                   response_params);
    auto generated_states = state_generator.generateStates();

    if (world.rank() == 0) {
      GeneratedStateData::print_generated_state_map(generated_states.state_map);
    }
    DerivedStatePlan derived_state_plan =
        DerivedStatePlanner::build_vbc_driven_quadratic_plan(
            response_params, ctx.molecule, ctx.ground.isSpinRestricted(),
            calc_params.protocol());
    StateParallelPlan state_parallel_plan = StateParallelPlanner::build(
        response_params, world.size(), generated_states.states);
    if (world.rank() == 0 && !derived_state_plan.requests.empty()) {
      print("🧩 Planned ", derived_state_plan.requests.size(),
            " VBC-derived quadratic-state requests (scaffolding only).");
    }
    if (world.rank() == 0 && response_params.state_parallel() != "off") {
      print("State-parallel plan: mode=", state_parallel_plan.effective_mode,
            " requested_groups=", state_parallel_plan.requested_groups,
            " mapping_groups=", state_parallel_plan.mapping_groups,
            " reason=", state_parallel_plan.reason);
    }
    world.gop.fence();
    return PlannedStates{std::move(generated_states),
                         std::move(derived_state_plan),
                         std::move(state_parallel_plan)};
  }

  static SolvedStates
  solve_all_states(World &world, const CalculationParameters &calc_params,
                   GroundContext &ctx,
                   const ResponseParameters &response_params,
                   PlannedStates planned_states) {
    (void)response_params;

    const auto &state_parallel_plan = planned_states.state_parallel_plan;
    const bool owner_group_schedule =
        state_parallel_plan.effective_mode == "owner_group_serial" &&
        state_parallel_plan.mapping_groups > 1;
    if (world.rank() == 0 && owner_group_schedule) {
      print("State ownership mapping is active across ",
            state_parallel_plan.mapping_groups,
            " groups; executing deterministic owner-group solve passes "
            "serially on the universe communicator.");
    } else if (world.rank() == 0 &&
               state_parallel_plan.mapping_groups > 1) {
      print("State ownership mapping is active across ",
            state_parallel_plan.mapping_groups,
            " groups; falling back to plain serial state loop.");
    }

    std::vector<size_t> owner_by_state_index(
        planned_states.generated_states.states.size(), 0);
    for (const auto &assignment : state_parallel_plan.assignments) {
      if (assignment.state_index < owner_by_state_index.size()) {
        owner_by_state_index[assignment.state_index] = assignment.owner_group;
      }
    }

    JsonStateSolvePersistence persistence(world, "response_metadata.json",
                                          "response_log.json");
    persistence.initialize_states(planned_states.generated_states.states);

    if (world.rank() == 0) {
      persistence.print_summary();
    }
    world.gop.fence();

    auto needs_solving_at_protocol = [&](double protocol_thresh,
                                         size_t thresh_index) {
      bool at_final_protocol =
          protocol_thresh == calc_params.protocol().back();

      for (const auto &state : planned_states.generated_states.states) {
        for (size_t freq_idx = 0; freq_idx < state.num_frequencies();
             ++freq_idx) {
          LinearResponsePoint pt{state, thresh_index, freq_idx};
          bool is_saved = persistence.is_saved(pt);
          bool should_solve =
              !is_saved || (at_final_protocol && !persistence.is_converged(pt));

          if (world.rank() == 0) {
            print("Checking state ", pt.perturbationDescription(),
                  " at thresh ", protocol_thresh, " freq ", pt.frequency(),
                  " is_saved=", is_saved,
                  " at_final_protocol=", at_final_protocol,
                  " should_solve=", should_solve);
          }

          if (should_solve) {
            return true;
          }
        }
      }

      return false;
    };

    const auto &protocol = calc_params.protocol();
    for (size_t ti = 0; ti < protocol.size(); ++ti) {
      double thresh = protocol[ti];

      if (!needs_solving_at_protocol(thresh, ti)) {
        if (world.rank() == 0) {
          madness::print("✓ All states converged at thresh", thresh,
                         "skipping to next protocol.");
        }
        continue;
      }

      ctx.response_manager.setProtocol(world, ctx.ground.getL(), thresh);
      ctx.ground.prepareOrbitals(world, FunctionDefaults<3>::get_k(), thresh);
      ctx.ground.computePreliminaries(
          world, *ctx.response_manager.getCoulombOp(),
          ctx.response_manager.getVtol(), ctx.fock_json_file);

      const bool at_final_protocol = (ti + 1 == protocol.size());
      auto solve_state_index = [&](size_t state_index) {
        auto &state = planned_states.generated_states.states[state_index];
        computeFrequencyLoop(world, ctx.response_manager, state, ti, ctx.ground,
                             persistence, at_final_protocol);
        persistence.flush_debug_log(world);
      };

      if (owner_group_schedule) {
        for (size_t gid = 0; gid < state_parallel_plan.mapping_groups; ++gid) {
          if (world.rank() == 0) {
            print("State solve lane ", gid, "/",
                  state_parallel_plan.mapping_groups - 1,
                  " at protocol thresh ", thresh);
          }
          for (size_t state_index = 0;
               state_index < planned_states.generated_states.states.size();
               ++state_index) {
            if (owner_by_state_index[state_index] == gid) {
              solve_state_index(state_index);
            }
          }
        }
      } else {
        for (size_t state_index = 0;
             state_index < planned_states.generated_states.states.size();
             ++state_index) {
          solve_state_index(state_index);
        }
      }
    }

    double final_thresh = calc_params.protocol().back();
    ctx.response_manager.setProtocol(world, ctx.ground.getL(), final_thresh);
    ctx.ground.prepareOrbitals(world, FunctionDefaults<3>::get_k(),
                               final_thresh);
    ctx.ground.computePreliminaries(world, *ctx.response_manager.getCoulombOp(),
                                    ctx.response_manager.getVtol(),
                                    ctx.fock_json_file);

    // Verify that all states are converged at the final protocol
    const size_t final_ti = calc_params.protocol().size() - 1;
    bool all_are_converged = true;
    for (const auto &state : planned_states.generated_states.states) {
      for (size_t fi = 0; fi < state.num_frequencies(); ++fi) {
        LinearResponsePoint pt{state, final_ti, fi};
        if (!persistence.is_converged(pt)) {
          all_are_converged = false;
          break;
        }
      }
      if (!all_are_converged)
        break;
    }

    MADNESS_ASSERT(all_are_converged);

    DerivedStateGateReport derived_gate =
        DerivedStatePlanner::evaluate_dependency_gate(
            planned_states.derived_state_plan, planned_states.generated_states,
            final_ti, [&](const LinearResponsePoint &pt) {
              return persistence.is_saved(pt) && persistence.is_converged(pt);
            });
    if (world.rank() == 0 && derived_gate.total_requests > 0) {
      print("Derived-state dependency gate: ready ",
            derived_gate.ready_requests, "/", derived_gate.total_requests,
            " blocked=", derived_gate.blocked_requests);
    }

    auto metadata = persistence.metadata_json();
    metadata["state_parallel_planner"] =
        planned_states.state_parallel_plan.to_json();
    metadata["derived_state_planner"] = {
        {"note",
         "Stage 2b scaffolding only: derived quadratic-state solves are not "
         "executed yet."},
        {"plan", planned_states.derived_state_plan.to_json()},
        {"dependency_gate", derived_gate.to_json()}};

    return SolvedStates{std::move(planned_states), std::move(metadata),
                        persistence.debug_log_json()};
  }

  struct PropertyContext {
    World &world;
    const ResponseParameters &response_params;
    const GroundContext &ground_ctx;
    const SolvedStates &solved_states;
    PropertyManager &properties;
    std::shared_ptr<SCF> scf_calc;
  };

  enum class PropertyType { Alpha, Beta, Raman };

  inline static PropertyType parse_property_name(const std::string &raw) {
    auto key = raw;
    if (key.size() >= 2 &&
        ((key.front() == '"' && key.back() == '"') ||
         (key.front() == '\'' && key.back() == '\''))) {
      key = key.substr(1, key.size() - 2);
    }
    if (key == "polarizability")
      return PropertyType::Alpha;
    if (key == "hyperpolarizability")
      return PropertyType::Beta;
    if (key == "raman")
      return PropertyType::Raman;
    MADNESS_EXCEPTION(std::string("Unknown property: " + key).c_str(), 0);
  }

  inline static void compute_polarizability(PropertyContext &ctx) {
    if (ctx.world.rank() == 0)
      madness::print("▶️ Computing polarizability α...");

    compute_alpha(ctx.world,
                  ctx.solved_states.planned_states.generated_states.state_map,
                  ctx.ground_ctx.ground,
                  ctx.response_params.dipole_frequencies(),
                  ctx.response_params.dipole_directions(), ctx.properties);

    ctx.properties.save();
  }

  inline static void compute_hyperpolarizability(PropertyContext &ctx) {
    if (ctx.world.rank() == 0)
      madness::print("▶️ Computing hyperpolarizability β...");

    auto dip_dirs = ctx.response_params.dipole_directions();
    ::compute_hyperpolarizability(ctx.world, ctx.ground_ctx.ground,
                                  ctx.response_params.dipole_frequencies(),
                                  dip_dirs, ctx.properties);

    ctx.properties.save();
  }

  inline static void compute_raman(PropertyContext &ctx,
                                   VibrationalResults &vib,
                                   RamanResults &raman) {
    if (ctx.world.rank() == 0) {
      madness::print("------------------------- Raman Computation "
                     "-------------------------");
    }

    vib = compute_hessian(
        ctx.world, ctx.solved_states.planned_states.generated_states.state_map,
        ctx.ground_ctx.ground, ctx.response_params.dipole_directions(),
        ctx.scf_calc);

    const double csg_factor = 142.9435756;
    Tensor<double> normal_modes = *vib.normalmodes_atomic;
    if (ctx.world.rank() == 0) {
      print("normal modes in atomic coordinates (au): \n", normal_modes);
    }
    auto vib_freq = *vib.frequencies * constants::au2invcm;

    std::vector<int> mode;
    for (int i = 0; i < vib_freq.size(); ++i) {

      if (abs(vib_freq(i)) > 1e-2) {
        mode.push_back(i);
        raman.vibrational_frequencies.push_back(vib_freq(i));
      }
    }
    auto nnmodes = normal_modes(_, Slice(mode[0], -1, 1));

    raman.normal_modes = nnmodes;

    if (ctx.world.rank() == 0) {
      print(mode);
      print(mode[0], mode[mode.size() - 1]);
    }
    auto alpha_derivatives =
        compute_Raman(ctx.world, ctx.ground_ctx.ground,
                      ctx.response_params.dipole_frequencies(),
                      ctx.response_params.dipole_directions(),
                      ctx.response_params.nuclear_directions(), ctx.properties);
    raman.polarization_frequencies = ctx.response_params.dipole_frequencies();
    ctx.properties.save();
    ctx.world.gop.fence();
    auto compute_alpha2 = [](const Tensor<double> &alpha) {
      auto alpha_mean = 0.0;
      for (int i = 0; i < 3; ++i) {
        alpha_mean += alpha(i, i);
      }
      alpha_mean *= (1.0 / 3.0);
      return alpha_mean * alpha_mean;
    };
    auto compute_beta2 = [](const Tensor<double> &alpha) {
      auto beta2 = 0.0;
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          beta2 +=
              0.5 * (3 * alpha(i, j) * alpha(i, j) - alpha(i, i) * alpha(j, j));
        }
      }
      return beta2;
    };

    // Raman Intensity linearly polarized light
    auto RamanIntensityL = [](const double alpha2, const double beta2) {
      return 45 * alpha2 + 4 * beta2; // in a.u.
    };
    auto DepolarizationRatio = [](const double alpha2, const double beta2) {
      return 3 * beta2 / (45 * alpha2 + 4 * beta2);
    };
    bool debug = ctx.world.rank() == 0 && ctx.response_params.print_level() > 1;
    if (debug) {
      for (size_t freq_idx = 0;
           freq_idx < raman.polarization_frequencies.size(); ++freq_idx) {
        print("Alpha derivatives for frequency ",
              raman.polarization_frequencies[freq_idx], " (a.u.): \n",
              alpha_derivatives[freq_idx]);
      };
    }
    // Compute Normal Mode
    for (size_t freq_idx = 0;
         freq_idx < raman.polarization_frequencies.size(); ++freq_idx) {

      const auto &alpha_dxyz = alpha_derivatives[freq_idx];
      double pol_freq = raman.polarization_frequencies[freq_idx];
      // Save polarizability derivatives to Results
      raman.polarizability_derivatives.push_back(alpha_dxyz);
      auto alpha_qi = inner(alpha_dxyz, nnmodes);
      raman.polarizability_derivatives_normal_modes.push_back(alpha_qi);
      if (debug) {
        print("Alpha derivative projected onto normal modes for frequency ",
              pol_freq, " (a.u.): \n", alpha_qi);
      }
    }
    // 1 a.u. of photon energy → wavelength (nm)
    auto wavelength_nm_from_au = [](double omega_au) {
      // lambda(nm) = 1239.841973 eV*nm / (27.211386 eV * omega_au)
      //            ≈ 45.56335 / omega_au
      const double AU_TO_NM_FACTOR = 45.5633525316;
      return (omega_au > 0.0) ? (AU_TO_NM_FACTOR / omega_au) : 0.0;
    };
    struct ColumnSpec {
      int width;
      int precision;
    };

    static const ColumnSpec COL_FREQ{10, 2};  // "Freq." like Dalton (xx.xx)
    static const ColumnSpec COL_FLOAT{11, 6}; // six decimals for the rest
                                              //
    auto print_header = [&](std::ostream &os, double omega_au) {
      const double lambda_nm = wavelength_nm_from_au(omega_au);
      os << "     Raman related properties for freq.  " << std::fixed
         << std::setprecision(6) << std::setw(9) << omega_au
         << " au  = " << std::setw(9) << std::setprecision(2) << lambda_nm
         << " nm\n";
      os << "     "
            "----------------------------------------------------------"
            "----"
            "-\n\n";
      os << " Mode    Freq.     Alpha**2   Beta(a)**2   Pol.Int.   "
            "Depol.Int.  Dep. Ratio \n\n";
    };

    auto print_row = [&](std::ostream &os,
                         const RamanResults::RamanModeRow &r) {
      // Compute missing fields if not provided
      double pol = r.pol_int ? *r.pol_int : (45.0 * r.alpha2 + 4.0 * r.beta2);
      double depi = r.depol_int ? *r.depol_int : (3.0 * r.beta2);
      double rho;
      if (r.dep_ratio) {
        rho = *r.dep_ratio;
      } else {
        rho = (pol != 0.0) ? (depi / pol) : 0.0; // guard div-by-zero
      }

      // Mode index
      os << std::setw(5) << r.mode
         << " "
         // Freq.
         << std::fixed << std::setw(COL_FREQ.width)
         << std::setprecision(COL_FREQ.precision) << r.freq_cm1
         << "  "
         // Alpha**2
         << std::setw(COL_FLOAT.width) << std::setprecision(COL_FLOAT.precision)
         << r.alpha2
         // Beta(a)**2
         << std::setw(COL_FLOAT.width) << std::setprecision(COL_FLOAT.precision)
         << r.beta2
         // Pol.Int.
         << std::setw(COL_FLOAT.width) << std::setprecision(COL_FLOAT.precision)
         << pol
         // Depol.Int.
         << std::setw(COL_FLOAT.width) << std::setprecision(COL_FLOAT.precision)
         << depi
         // Dep. Ratio
         << std::setw(COL_FLOAT.width) << std::setprecision(COL_FLOAT.precision)
         << rho << "\n";
    };

    auto print_block =
        [&](std::ostream &os, double omega_au,
            const std::vector<RamanResults::RamanModeRow> &rows) {
          print_header(os, omega_au);
          // print in reverse order to match Dalton output
          // for (auto const& x : range | std::views::reverse)
          for (auto r_it = rows.rbegin(); r_it != rows.rend(); ++r_it)
            print_row(os, *r_it);
          os << "\n";
        };
    using RamanModeRow = RamanResults::RamanModeRow;

    for (size_t freq_idx = 0;
         freq_idx < raman.polarization_frequencies.size(); ++freq_idx) {
      auto alpha_qi = raman.polarizability_derivatives_normal_modes[freq_idx];
      vector<RamanModeRow> raman_rows;
      for (size_t i = 0; i < mode.size(); ++i) {
        RamanModeRow row;
        double vib_freq_i = vib_freq(mode[i]);
        row.mode = static_cast<int>(i) + 1;
        row.freq_cm1 = vib_freq_i;
        auto alpha_i = copy(alpha_qi(_, i));
        auto alpha = alpha_i.reshape(3, 3);
        auto alpha2_au = compute_alpha2(alpha);
        row.alpha2 = alpha2_au * csg_factor;
        auto beta2_au = compute_beta2(alpha);
        row.beta2 = beta2_au * csg_factor;
        auto depol_int_au = RamanIntensityL(alpha2_au, beta2_au);
        row.pol_int = depol_int_au * csg_factor;
        row.dep_ratio = DepolarizationRatio(alpha2_au, beta2_au);
        row.depol_int = (*row.dep_ratio) * (*row.pol_int);
        ctx.world.gop.fence();
        raman_rows.push_back(row);
      }
      raman.raman_spectra[raman.polarization_frequencies[freq_idx]] =
          raman_rows;
      ctx.world.gop.fence();
      if (ctx.world.rank() == 0) {
        print_block(std::cout, raman.polarization_frequencies[freq_idx],
                    raman_rows);
      }

      ctx.world.gop.fence();
    }
  }

  static PropertyStageOutput
  compute_requested_properties(World &world,
                               const ResponseParameters &response_params,
                               const GroundContext &ground_ctx,
                               const SolvedStates &solved_states,
                               const std::shared_ptr<SCF> &scf_calc) {
    VibrationalResults vib;
    RamanResults raman;
    PropertyManager properties(world, "properties.json");

    PropertyContext prop_ctx{world,         response_params, ground_ctx,
                             solved_states, properties,      scf_calc};

    for (const std::string &prop : response_params.requested_properties()) {
      auto prop_type = parse_property_name(prop);
      switch (prop_type) {
      case PropertyType::Alpha:
        compute_polarizability(prop_ctx);
        break;
      case PropertyType::Beta:
        compute_hyperpolarizability(prop_ctx);
        break;
      case PropertyType::Raman:
        compute_raman(prop_ctx, vib, raman);
        break;
      }
    }

    return PropertyStageOutput{properties.to_json(), std::move(vib),
                               std::move(raman)};
  }

  public:

  /**
   * @brief Run the full molecular response & property workflow.
   *
   * @param world      The MADNESS world communicator
   * @param params     Unified parameters containing response and molecule info
   * @param outdir     Directory where all outputs will be written
   * @return Results   Structured JSON fragments: metadata + properties
   */
  inline static Results run_response(World &world, const Params &params,
                                     const std::shared_ptr<SCF> &scf_calc,
                                     const std::filesystem::path &outdir) {
    const auto &calc_params = params.get<CalculationParameters>();
    const auto &rp_copy = params.get<ResponseParameters>();

    if (world.rank() == 0) {
      json response_input_json = {};
      response_input_json["response"] =
          rp_copy.to_json_if_precedence("defined");
      print("response_input_json: ", response_input_json.dump(4));
      std::ofstream ofs("response.in");
      write_json_to_input_file(response_input_json, {"response"}, ofs);
      ofs.close();
    }
    world.gop.fence();
    commandlineparser parser;
    parser.set_keyval("input", "response.in");
    if (world.rank() == 0) {
      ::print("input filename: ", parser.value("input"));
    }

    auto response_params = ResponseParameters(world, parser);
    GroundContext ground_ctx =
        make_ground_context(world, calc_params, scf_calc, outdir);

    PlannedStates planned_states =
        plan_required_states(world, calc_params, ground_ctx, response_params);
    SolvedStates solved_states = solve_all_states(
        world, calc_params, ground_ctx, response_params, std::move(planned_states));
    PropertyStageOutput property_stage =
        compute_requested_properties(world, response_params, ground_ctx,
                                     solved_states, scf_calc);

    // finalize & stats
    if (world.rank() == 0) {
      madness::print(
          "\n✅ Molecular response & property calculation complete.");
    }
    world.gop.fence();
    world.gop.fence();
    print_stats(world);

    // aggregate JSON results
    Results results;
    results.metadata = std::move(solved_states.metadata);
    results.properties = std::move(property_stage.properties);
    results.debug_log = std::move(solved_states.debug_log);
    results.vibrational_analysis = property_stage.vibrational_analysis.to_json();
    results.raman_spectra = property_stage.raman_spectra.to_json();
    return results;
  }
}; // namespace molresponse_lib
