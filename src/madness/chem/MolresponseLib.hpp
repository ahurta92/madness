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
  // -----------------------------------------------------------------------------
  // Container for structured JSON fragments produced by the workflow
  // -----------------------------------------------------------------------------
  struct Results {
    nlohmann::json metadata;             // convergence metadata per state
    nlohmann::json properties;           // computed α, β, Raman property tables
    nlohmann::json vibrational_analysis; // vibrational analysis results
    nlohmann::json raman_spectra;        // Raman spectra results
    nlohmann::json debug_log;            // debug log of response calculations
  };
  static constexpr const char *label() { return "molresponse"; }

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
    // --- configure the ground-state archive location ---
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
    auto scf_json = nlohmann::json{};

    // moldft calc_info.json
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

    const auto molecule = read_molecule(moldft_checkpt);
    print("Read molecule with ", molecule.natom(), " atoms.");
    GroundStateData ground(world, relative_archive.string(), molecule);
    ResponseManager response_manager(world, calc_params);

    if (world.rank() == 0) {
      print("dipole.frequencies: ", response_params.dipole_frequencies());
      print("dipole.directions: ", response_params.dipole_directions());
      print("nuclear.frequencies: ", response_params.nuclear_frequencies());
      print("nuclear.directions: ", response_params.nuclear_directions());
      print("requested_properties: ", response_params.requested_properties());
    }
    // generate response states
    StateGenerator state_generator(molecule, calc_params.protocol(),
                                   ground.isSpinRestricted(), response_params);
    auto generated_states = state_generator.generateStates();
    if (world.rank() == 0) {
      GeneratedStateData::print_generated_state_map(generated_states.state_map);
    }
    world.gop.fence();

    // initialize metadata
    std::string meta_file = "response_metadata.json";

    ResponseRecord2 response_record(world, meta_file);
    response_record.initialize_states(generated_states.states);
    if (world.rank() == 0) {

      response_record.print_summary();
    }
    world.gop.fence();

    // debug logger
    ResponseDebugLogger debug_logger("response_log.json", true);

    // if any state and at any frequency needs solving at this protocol return
    // true
    auto needs_solving_at_protocol = [&](double protocol_thresh) {
      bool at_final_protocol =
          protocol_thresh ==
          calc_params.protocol().back(); // are we at the final protocol?
      for (const auto &state : generated_states.states) {
        for (size_t freq_idx = 0; freq_idx < state.frequencies.size();
             ++freq_idx) {
          bool is_saved = response_record.is_saved(state);
          bool should_solve =
              !is_saved ||
              (at_final_protocol && !response_record.is_converged(state));

          if (world.rank() == 0) {
            print("Checking state ", state.description(), " at thresh ",
                  protocol_thresh, " freq ", state.frequencies[freq_idx],
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

    // loop over thresholds
    for (double thresh : calc_params.protocol()) {

      if (!needs_solving_at_protocol(thresh)) {
        if (world.rank() == 0) {
          madness::print("✓ All states converged at thresh", thresh,
                         "skipping to next protocol.");
        }
        // advanced threshold for all states
        for (auto &state : generated_states.states) {
          state.advance_threshold();
        }
        continue;
      }

      response_manager.setProtocol(world, ground.getL(), thresh);
      ground.prepareOrbitals(world, FunctionDefaults<3>::get_k(), thresh);
      ground.computePreliminaries(world, *response_manager.getCoulombOp(),
                                  response_manager.getVtol(), fock_json_file);
      // if (world.rank() == 0)
      //   madness::print("hamiltonian:\n", ground.Hamiltonian);

      for (auto &state : generated_states.states) {
        //     if (state.is_converged || state.current_threshold() != thresh)
        //     continue;

        computeFrequencyLoop(world, response_manager, state, ground,
                             response_record, debug_logger);

        if (debug_logger.enabled() && world.rank() == 0) {
          debug_logger.write_to_disk();
        }
        state.advance_threshold();
      }
    }

    bool all_are_converged = std::all_of(
        generated_states.states.begin(), generated_states.states.end(),
        [response_record](const LinearResponseDescriptor &s) {
          return response_record.is_converged(s);
        });
    double thresh = calc_params.protocol().back();
    response_manager.setProtocol(world, ground.getL(), thresh);
    ground.prepareOrbitals(world, FunctionDefaults<3>::get_k(), thresh);
    ground.computePreliminaries(world, *response_manager.getCoulombOp(),
                                response_manager.getVtol(), fock_json_file);

    MADNESS_ASSERT(all_are_converged);
    VibrationalResults vib;
    RamanResults raman;

    // compute requested properties
    PropertyManager properties(world, "properties.json");
    std::string dip_dirs = response_params.dipole_directions();
    enum class PropertyType { Alpha, Beta, Raman };

    std::map<std::string, PropertyType> prop_map = {
        {"polarizability", PropertyType::Alpha},
        {"hyperpolarizability", PropertyType::Beta},
        {"raman", PropertyType::Raman}};

    for (const std::string &prop : response_params.requested_properties()) {
      // get rid of first and last characters
      auto prop_type = prop_map[prop.substr(1, prop.size() - 2)];
      if (prop_type == PropertyType::Alpha) {
        if (world.rank() == 0) {
          madness::print("▶️ Computing polarizability α...");
        }
        compute_alpha(world, generated_states.state_map, ground,
                      response_params.dipole_frequencies(),
                      response_params.dipole_directions(), properties);
        properties.save();
      } else if (prop_type == PropertyType::Beta) {
        if (world.rank() == 0) {

          madness::print("▶️ Computing hyperpolarizability β...");
        }

        compute_hyperpolarizability(world, ground,
                                    response_params.dipole_frequencies(),
                                    dip_dirs, properties);
        properties.save();
      } else if (prop_type == PropertyType::Raman) {
        if (world.rank() == 0) {
          madness::print("------------------------- Raman Computation "
                         "-------------------------");
        }

        vib = compute_hessian(world, generated_states.state_map, ground,
                              response_params.dipole_directions(), scf_calc);

        auto unit_amu_toau = std::sqrt(constants::atomic_mass_in_au);
        const double csg_factor = 142.9435756;
        Tensor<double> normal_modes = *vib.normalmodes_atomic;
        if(world.rank()==0){
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

        if (world.rank() == 0) {
          print(mode);
          print(mode[0], mode[mode.size() - 1]);
        }
        auto alpha_derivatives =
            compute_Raman(world, ground, response_params.dipole_frequencies(),
                          response_params.dipole_directions(),
                          response_params.nuclear_directions(), properties);
        raman.polarization_frequencies = response_params.dipole_frequencies();
        properties.save();
        world.gop.fence();
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
              beta2 += 0.5 * (3 * alpha(i, j) * alpha(i, j) -
                              alpha(i, i) * alpha(j, j));
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
        bool debug = world.rank() == 0 && rp_copy.print_level() > 1;
        if (debug) {
          for (int freq_idx = 0;
               freq_idx < raman.polarization_frequencies.size(); ++freq_idx) {
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

        auto print_row = [](std::ostream &os,
                            const RamanResults::RamanModeRow &r) {
          // Compute missing fields if not provided
          double pol =
              r.pol_int ? *r.pol_int : (45.0 * r.alpha2 + 4.0 * r.beta2);
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
             << std::setw(COL_FLOAT.width)
             << std::setprecision(COL_FLOAT.precision)
             << r.alpha2
             // Beta(a)**2
             << std::setw(COL_FLOAT.width)
             << std::setprecision(COL_FLOAT.precision)
             << r.beta2
             // Pol.Int.
             << std::setw(COL_FLOAT.width)
             << std::setprecision(COL_FLOAT.precision)
             << pol
             // Depol.Int.
             << std::setw(COL_FLOAT.width)
             << std::setprecision(COL_FLOAT.precision)
             << depi
             // Dep. Ratio
             << std::setw(COL_FLOAT.width)
             << std::setprecision(COL_FLOAT.precision) << rho << "\n";
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
        auto row_to_mode_data = [](const RamanResults::RamanModeRow &r) {
          return std::unordered_map<std::string, double>{
              {"mode", static_cast<double>(r.mode)},
              {"frequency_cm-1", r.freq_cm1},
              {"alpha2", r.alpha2},
              {"beta2", r.beta2},
              {"pol_int", r.pol_int ? *r.pol_int : 0.0},
              {"depol_int", r.depol_int ? *r.depol_int : 0.0},
              {"depol_ratio", r.dep_ratio ? *r.dep_ratio : 0.0},
          };
        };
        using RamanModeRow = RamanResults::RamanModeRow;

        for (int freq_idx = 0; freq_idx < raman.polarization_frequencies.size();
             ++freq_idx) {
          auto alpha_qi =
              raman.polarizability_derivatives_normal_modes[freq_idx];
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
            row.alpha2 = alpha2_au * csg_factor;
            auto depol_int_au = RamanIntensityL(alpha2_au, beta2_au);
            row.pol_int = depol_int_au * csg_factor;
            row.dep_ratio = DepolarizationRatio(alpha2_au, beta2_au);
            row.depol_int = (*row.dep_ratio) * (*row.depol_int);
            world.gop.fence();
            raman_rows.push_back(row);
          }
          raman.raman_spectra[raman.polarization_frequencies[freq_idx]] =
              raman_rows;
          world.gop.fence();
          if (world.rank() == 0) {
            print_block(std::cout, raman.polarization_frequencies[freq_idx],
                        raman_rows);
          }

          world.gop.fence();
        }
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
    results.metadata = response_record.to_json();
    results.properties = properties.to_json();
    results.debug_log = debug_logger.to_json();
    results.vibrational_analysis = vib.to_json();
    results.raman_spectra = raman.to_json();
    return results;
  }
}; // namespace molresponse_lib
