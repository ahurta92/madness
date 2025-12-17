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
#include <apps/molresponse_v2/GroundStateData.hpp>
#include <apps/molresponse_v2/PropertyManager.hpp>
#include <apps/molresponse_v2/ResponseDebugLogger.hpp>
#include <apps/molresponse_v2/ResponseManager.hpp>
#include <apps/molresponse_v2/ResponseRecord.hpp>
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

  struct StateSolution {
    GeneratedStateData generated_states;
    ResponseRecord2 response_record;
    ResponseDebugLogger debug_logger;
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

    std::string archive_name = "mad.restartdata";
    std::string fock_json_file = prox / "moldft.fock.json";
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

  static StateSolution
  solve_all_states(World &world, const CalculationParameters &calc_params,
                   GroundContext &ctx,
                   const ResponseParameters &response_params) {
    StateGenerator state_generator(ctx.molecule, calc_params.protocol(),
                                   ctx.ground.isSpinRestricted(),
                                   response_params);
    auto generated_states = state_generator.generateStates();

    if (world.rank() == 0) {
      GeneratedStateData::print_generated_state_map(generated_states.state_map);
    }
    world.gop.fence();

    std::string meta_file = "response_metadata.json";
    ResponseRecord2 response_record(world, meta_file);
    response_record.initialize_states(generated_states.states);

    if (world.rank() == 0) {
      response_record.print_summary();
    }
    world.gop.fence();

    ResponseDebugLogger debug_logger("response_log.json", true);

    auto needs_solving_at_protocol = [&](double protocol_thresh,
                                         size_t thresh_index) {
      bool at_final_protocol =
          protocol_thresh == calc_params.protocol().back();

      for (const auto &state : generated_states.states) {
        for (size_t freq_idx = 0; freq_idx < state.num_frequencies();
             ++freq_idx) {
          LinearResponsePoint pt{state, thresh_index, freq_idx};
          bool is_saved =
              response_record.is_saved(pt.perturbationDescription(),
                                       pt.threshold(), pt.frequency());
          bool should_solve =
              !is_saved ||
              (at_final_protocol &&
               !response_record.is_converged(pt.perturbationDescription(),
                                             pt.threshold(), pt.frequency()));

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

      for (auto &state : generated_states.states) {
        bool at_final_protocol = (ti + 1 == protocol.size());
        computeFrequencyLoop(world, ctx.response_manager, state, ti,
                             ctx.ground, response_record, debug_logger,
                             at_final_protocol);

        if (debug_logger.enabled() && world.rank() == 0) {
          debug_logger.write_to_disk();
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
    bool all_are_converged = true;
    for (const auto &state : generated_states.states) {
      for (double f : state.frequencies) {
        if (!response_record.is_converged(state.perturbationDescription(),
                                          final_thresh, f)) {
          all_are_converged = false;
          break;
        }
      }
      if (!all_are_converged)
        break;
    }

    MADNESS_ASSERT(all_are_converged);

    return StateSolution{std::move(generated_states),
                         std::move(response_record), std::move(debug_logger)};
  }

  struct PropertyContext {
    World &world;
    const ResponseParameters &response_params;
    const GroundContext &ground_ctx;
    const StateSolution &state_solution;
    PropertyManager &properties;
    std::shared_ptr<SCF> scf_calc;
  };

  enum class PropertyType { Alpha, Beta, Raman };

  inline static PropertyType parse_property_name(const std::string &raw) {
    // strip leading/trailing quotes as you already do
    auto key = raw.substr(1, raw.size() - 2);
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

    compute_alpha(ctx.world, ctx.state_solution.generated_states.state_map,
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
        ctx.world, ctx.state_solution.generated_states.state_map,
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
      for (int freq_idx = 0; freq_idx < raman.polarization_frequencies.size();
           ++freq_idx) {
        print("Alpha derivatives for frequency ",
              raman.polarization_frequencies[freq_idx], " (a.u.): \n",
              alpha_derivatives[freq_idx]);
      };
    }
    // Compute Normal Mode
    for (int freq_idx = 0; freq_idx < raman.polarization_frequencies.size();
         ++freq_idx) {

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

    for (int freq_idx = 0; freq_idx < raman.polarization_frequencies.size();
         ++freq_idx) {
      auto alpha_qi = raman.polarizability_derivatives_normal_modes[freq_idx];
      vector<RamanModeRow> raman_rows;
      for (int i = 0; i < mode.size(); ++i) {
        RamanModeRow row;
        double vib_freq_i = vib_freq(mode[i]);
        row.mode = i + 1;
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

    StateSolution state_solution =
        solve_all_states(world, calc_params, ground_ctx, response_params);

    VibrationalResults vib;
    RamanResults raman;

    // compute requested properties
    PropertyManager properties(world, "properties.json");

    PropertyContext prop_ctx{world,          response_params, ground_ctx,
                             state_solution, properties,      scf_calc};

    std::string dip_dirs = response_params.dipole_directions();

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
    results.metadata = state_solution.response_record.to_json();
    results.properties = properties.to_json();
    results.debug_log = state_solution.debug_logger.to_json();
    results.vibrational_analysis = vib.to_json();
    results.raman_spectra = raman.to_json();
    return results;
  }
}; // namespace molresponse_lib
