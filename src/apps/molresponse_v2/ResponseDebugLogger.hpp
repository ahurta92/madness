#pragma once
#include <fstream>
#include <iomanip>
#include <madness/external/nlohmann_json/json.hpp>

#include "ResponseState.hpp"

namespace fs = std::filesystem;
using json = nlohmann::json;

//==============================================================================
// ResponseDebugLogger
//   - accumulates debug timing/value data across runs
//   - preserves earlier entries for the same state.description()
//==============================================================================
class ResponseDebugLogger {
 public:
  // constructor: pass the filename you want to write to (e.g.
  // "responses/response_log.json")
  explicit ResponseDebugLogger(const std::string &filename,
                               bool enabled = false)
      : filename_(filename), enabled_(enabled) {
    // if an existing log file is there, load it
    if (fs::exists(filename_)) {
      std::ifstream in(filename_);
      in >> log_data_;
    }
  }

  [[nodiscard]] json to_json() const { return log_data_; }

  [[nodiscard]] bool enabled() const { return enabled_; }
  void set_enabled(bool on) { enabled_ = on; }

  // must call at the start of each new state
  void start_state(const LinearResponseDescriptor &state) {
    if (!enabled_) return;

    current_key_ = state.description();
    // if we've previously logged this state, start from its existing entry,
    // otherwise create a fresh template
    if (log_data_.contains(current_key_)) {
      current_entry_ = log_data_[current_key_];
    } else {
      current_entry_ = {{"iterations", json::array()}};
    }
  }

  // call at the start of each new iterate() call
  void begin_iteration(World &world, size_t iter_index) {
    if (!enabled_) return;
    if (world.rank() != 0) return;

    json values = {{"steps", json::object()}};
    json timings = {{"steps", json::object()}};

    current_iter_ = {
        {"iter", iter_index}, {"values", values}, {"timings", timings}};
  }

  // log a numeric value under a step name
  template <typename T>
  void log_value(World &world, const std::string &step_name, const T &value) {
    if (!enabled_) return;
    if (world.rank() != 0) return;
    current_iter_["values"]["steps"][step_name]["value"] = value;
  }

  // log wall + cpu times under a step name
  void log_timing(World &world, const std::string &step_name, double wall_time,
                  double cpu_time) {
    if (!enabled_) return;

    if (world.rank() == 0) {
      auto &S = current_iter_["timings"]["steps"][step_name];
      S["wall_time"] = wall_time;
      S["cpu_time"] = cpu_time;
    }
  }

  // convenience for both in one shot
  template <typename T>
  void log_value_and_time(World &world, const std::string &step_name,
                          const T &value, double wall_time, double cpu_time) {
    if (!enabled_) return;
    if (world.rank() == 0) {
      log_value(world, step_name, value);
      log_timing(world, step_name, wall_time, cpu_time);
    }
  }

  // call at the end of each iteration
  void end_iteration(World &world) {
    if (!enabled_) return;
    if (world.rank() == 0) {
      std::cout << "📊 Iteration summary for " << current_key_ << "\n";
      current_entry_["iterations"].push_back(current_iter_);
      std::cout << current_entry_.dump(2) << "\n";
    }
  }

  // call once per state, after the final iteration
  void finalize_state(World &world) {
    if (!enabled_) return;
    if (world.rank() == 0) {
      log_data_[current_key_] = current_entry_;
      print("📂 Wrote response log for", current_key_);
      std::cout << log_data_.dump(2) << "\n";

      write_to_disk(world);
    }
  }

  // write the entire accumulated JSON to disk
  void write_to_disk(World &world) const {
    if (!enabled_) return;
    if (world.rank() != 0) return;
    std::ofstream out(filename_);
    out << std::setw(2) << log_data_.dump(2) << "\n";
  }

  // ——— pretty‐printers ——————————————————————————————————————————————

  void print_timing_table(const std::string &description) const {
    if (!enabled_ || !log_data_.contains(description)) return;
    const auto &iter_data = log_data_.at(description)["iteration"]["timings"];

    // gather all step names
    std::set<std::string> step_names;
    for (auto &iter : iter_data) {
      for (auto &p : iter["steps"].items()) step_names.insert(p.key());
    }

    // map each to a 5-char header
    std::map<std::string, std::string> short_key;
    for (auto &name : step_names) {
      short_key[name] = name.substr(0, 5);
    }

    constexpr int W = 10;
    std::cout << "\n⏱️ Timing for " << description << "\n";
    std::cout << std::setw(6) << "Iter";
    for (auto &n : step_names) std::cout << std::setw(W) << short_key[n];
    std::cout << "\n" << std::string(6 + W * step_names.size(), '-') << "\n";

    for (auto &iter : iter_data) {
      std::cout << std::setw(6) << iter["iter"].get<int>();
      for (auto &n : step_names) {
        if (iter["steps"].contains(n)) {
          double w = iter["steps"][n]["wall_time"].get<double>();
          std::cout << std::setw(W) << std::fixed << std::setprecision(4) << w;
        } else {
          std::cout << std::setw(W) << "N/A";
        }
      }
      std::cout << "\n";
    }
  }

  void print_values_table(const std::string &description) const {
    if (!enabled_ || !log_data_.contains(description)) return;
    const auto &iter_data = log_data_.at(description)["iteration"]["values"];

    // gather step names
    std::set<std::string> step_names;
    for (auto &iter : iter_data) {
      for (auto &p : iter["steps"].items()) step_names.insert(p.key());
    }

    constexpr int W = 16;
    std::cout << "\n📋 Values for " << description << "\n";
    std::cout << std::setw(6) << "Iter";
    for (auto &n : step_names) std::cout << std::setw(W) << n;
    std::cout << "\n" << std::string(6 + W * step_names.size(), '-') << "\n";

    for (auto &iter : iter_data) {
      std::cout << std::setw(6) << iter["iter"].get<int>();
      for (auto &n : step_names) {
        if (iter["steps"].contains(n)) {
          double v = iter["steps"][n]["value"].get<double>();
          std::cout << std::setw(W) << std::scientific << std::setprecision(5)
                    << v;
        } else {
          std::cout << std::setw(W) << "N/A";
        }
      }
      std::cout << "\n";
    }
  }

 private:
  bool enabled_;
  std::string filename_;
  json log_data_;

  // per‐state temporary
  std::string current_key_;
  json current_entry_;
  json current_iter_;
};

class TimedValueLogger {
 public:
  TimedValueLogger(madness::World &world, const std::string &key,
                   ResponseDebugLogger *logger = nullptr)
      : world_(world), key_(key), logger_(logger) {
    world_.gop.fence();
    start_wall_ = madness::wall_time();
    start_cpu_ = madness::cpu_time();
  }

  void log() {
    double wall = madness::wall_time() - start_wall_;
    double cpu = madness::cpu_time() - start_cpu_;

    logger_->log_timing(world_, key_, wall, cpu);

    if (world_.rank() == 0) {
      std::cout << std::left << std::setw(30) << "[" + key_ + "]" << std::right
                << " | Wall: " << std::setw(7) << std::setprecision(3) << wall
                << "s | CPU: " << std::setw(7) << std::setprecision(3) << cpu
                << "s |" << std::endl;
    }
  }

  template <typename T>
  void log(const T &value) {
    double wall = madness::wall_time() - start_wall_;
    double cpu = madness::cpu_time() - start_cpu_;

    logger_->log_value_and_time(world_, key_, value, wall, cpu);

    if (world_.rank() == 0) {
      std::cout << std::left << std::setw(30) << "[" + key_ + "]" << std::right
                << " | Wall: " << std::setw(7) << std::setprecision(3) << wall
                << "s | CPU: " << std::setw(7) << std::setprecision(3) << cpu
                << "s | Value: " << std::setw(10) << std::setprecision(7)
                << std::fixed << value << " |" << std::endl;
    }
  }
  template <typename T>
  void log(const std::vector<T> &values) {
    double wall = madness::wall_time() - start_wall_;
    double cpu = madness::cpu_time() - start_cpu_;

    if (logger_ && world_.rank() == 0) {
      for (size_t i = 0; i < values.size(); ++i) {
        std::string key_i = key_ + "[" + std::to_string(i) + "]";
        logger_->log_value_and_time(key_i, values[i], wall, cpu);
      }
    }

    if (world_.rank() == 0) {
      std::cout << std::left << std::setw(20) << ("⏱️ [" + key_ + "]")
                << " | Values: ";
      for (const auto &v : values) {
        std::cout << std::right << std::setw(10) << std::setprecision(6)
                  << std::fixed << v << " ";
      }
      std::cout << "| Wall: " << std::setw(7) << std::setprecision(3) << wall
                << "s | CPU: " << std::setw(7) << std::setprecision(3) << cpu
                << "s |\n";
    }
  }

  template <typename T>
  void log_value(const T &value) {
    logger_->log_value(world_, key_, value);
  }

 private:
  madness::World &world_;
  std::string key_;
  ResponseDebugLogger *logger_;
  double start_wall_, start_cpu_;
};
