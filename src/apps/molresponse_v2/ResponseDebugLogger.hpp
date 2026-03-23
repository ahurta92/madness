#pragma once
#include <fstream>
#include <iomanip>
#include <optional>
#include <sstream>
#include <type_traits>
#include <madness/external/nlohmann_json/json.hpp>

#include "ResponseState.hpp"

namespace fs = std::filesystem;
using json = nlohmann::json;

class ResponseDebugLogger {
public:
  // Pass filename and whether logging is enabled.
  // If you're in MPI, pass rank0_only=true and only construct/use this on rank
  // 0, or set rank0_only=true and call write_to_disk() from all ranks safely.
  explicit ResponseDebugLogger(std::string filename, bool enabled = false)
      : enabled_(enabled), filename_(std::move(filename)) {
    if (!enabled_)
      return;
    if (fs::exists(filename_)) {
      std::ifstream inf(filename_);
      if (inf)
        inf >> log_data_;
    }
    ensure_root_();
  }

  [[nodiscard]] bool enabled() const { return enabled_; }
  void set_enabled(bool on) { enabled_ = on; }

  // Call at the start of each (state, protocol, frequency) block
  // before iteration begins.
  void start_state(const LinearResponsePoint &pt) {
    start_named_state(pt.perturbationDescription(), pt.threshold(),
                      pt.frequency());
  }

  // Start a logging context using an explicit state label and numeric protocol /
  // frequency tags. This is useful for excited-state bundle workflows that do not
  // naturally map to a LinearResponsePoint.
  void start_named_state(const std::string &state_key, double protocol,
                         double frequency) {
    if (!enabled_)
      return;
    state_key_ = state_key;
    proto_key_ = protocol_key_(protocol);
    freq_key_ = freq_key(frequency);

    auto &node = node_ref_();
    if (!node.contains("iteration_values"))
      node["iteration_values"] = json::array();
    if (!node.contains("iteration_timings"))
      node["iteration_timings"] = json::array();
  }

  // Start of each solver iterate() call
  void begin_iteration(std::size_t iter_index) {
    if (!enabled_)
      return;
    current_iter_values_ = {{"iter", iter_index}, {"steps", json::object()}};
    current_iter_timing_ = {{"iter", iter_index}, {"steps", json::object()}};
  }

  // Numeric value for a named step (keep T JSON-serializable; prefer
  // double/int)
  template <typename T>
  void log_value(const std::string &step_name, const T &value) {
    if (!enabled_)
      return;
    current_iter_values_["steps"][step_name]["value"] = value;
  }

  // Timing for a named step
  void log_timing(const std::string &step_name, double wall_time,
                  double cpu_time) {
    if (!enabled_)
      return;
    auto &s = current_iter_timing_["steps"][step_name];
    s["wall_time"] = wall_time;
    s["cpu_time"] = cpu_time;
  }

  // Convenience
  template <typename T>
  void log_value_and_time(const std::string &step_name, const T &value,
                          double wall_time, double cpu_time) {
    log_value(step_name, value);
    log_timing(step_name, wall_time, cpu_time);
  }

  // End of one iterate() call and append to arrays on the active node.
  void end_iteration() {
    if (!enabled_)
      return;
    auto &node = node_ref_();
    node["iteration_values"].push_back(current_iter_values_);
    node["iteration_timings"].push_back(current_iter_timing_);
  }

  // Optional: nothing to do right now, but good for symmetry
  void finalize_state() const {
    if (!enabled_)
      return;
    // no-op; everything is already written into log_data_
  }

  // Persist the entire log JSON
  // If rank0_only_ == true, call this only on rank 0 (or call from all ranks
  // with identical state).
  void write_to_disk() const {
    if (!enabled_) {
      return;
    }
    std::ofstream ofs(filename_);
    ofs << std::setw(2) << log_data_;
    if (ofs) {
      std::cout << "DEBUG_LOG_WRITTEN file=" << filename_ << "\n";
    } else {
      std::cerr << "WARN DEBUG_LOG_WRITE_FAILED file=" << filename_ << "\n";
    }
    // Make sure the directory exists
  }

  json to_json() const { return log_data_; }

  // Pretty printers

  void print_timing_table(const LinearResponsePoint &pt) const {
    if (!enabled_)
      return;
    const auto state_key = pt.perturbationDescription();
    const auto proto_key = protocol_key_(pt.threshold());
    const auto fkey = freq_key(pt.frequency());

    if (!exists_(state_key, proto_key, fkey))
      return;
    const auto &iter_data = log_data_.at(state_key)
                                .at("protocols")
                                .at(proto_key)
                                .at("freqs")
                                .at(fkey)
                                .at("iteration_timings");

    // gather all step names
    std::set<std::string> step_names;
    for (const auto &iter : iter_data) {
      for (auto &p : iter["steps"].items())
        step_names.insert(p.key());
    }

    constexpr int iter_width = 6;
    constexpr int col_width = 14;
    auto print_border = [&]() {
      std::cout << "+" << std::string(iter_width + 2, '-');
      for (size_t i = 0; i < step_names.size(); ++i) {
        std::cout << "+" << std::string(col_width + 2, '-');
      }
      std::cout << "+\n";
    };
    auto trunc_label = [&](const std::string &name) -> std::string {
      if (name.size() <= static_cast<size_t>(col_width)) {
        return name;
      }
      return name.substr(0, static_cast<size_t>(col_width - 1));
    };

    std::cout << "TIMING_TABLE state=" << state_key
              << " protocol=" << proto_key << " frequency=" << fkey << "\n";
    print_border();
    std::cout << "| " << std::setw(iter_width) << "iter" << " ";
    for (const auto &name : step_names) {
      std::cout << "| " << std::setw(col_width) << trunc_label(name) << " ";
    }
    std::cout << "|\n";
    print_border();

    for (const auto &iter : iter_data) {
      std::cout << "| " << std::setw(iter_width) << iter["iter"].get<int>()
                << " ";
      for (const auto &name : step_names) {
        if (iter["steps"].contains(name)) {
          const double wall = iter["steps"][name]["wall_time"].get<double>();
          std::ostringstream os;
          os << std::fixed << std::setprecision(4) << wall;
          std::cout << "| " << std::setw(col_width) << os.str() << " ";
        } else {
          std::cout << "| " << std::setw(col_width) << "N/A" << " ";
        }
      }
      std::cout << "|\n";
    }
    print_border();
  }

  void print_values_table(const LinearResponsePoint &pt) const {
    if (!enabled_)
      return;
    const auto state_key = pt.perturbationDescription();
    const auto proto_key = protocol_key_(pt.threshold());
    const auto fkey = freq_key(pt.frequency());

    if (!exists_(state_key, proto_key, fkey))
      return;
    const auto &iter_data = log_data_.at(state_key)
                                .at("protocols")
                                .at(proto_key)
                                .at("freqs")
                                .at(fkey)
                                .at("iteration_values");

    std::set<std::string> step_names;
    for (const auto &iter : iter_data) {
      for (auto &s : iter["steps"].items())
        step_names.insert(s.key());
    }

    constexpr int iter_width = 6;
    constexpr int col_width = 16;
    auto print_border = [&]() {
      std::cout << "+" << std::string(iter_width + 2, '-');
      for (size_t i = 0; i < step_names.size(); ++i) {
        std::cout << "+" << std::string(col_width + 2, '-');
      }
      std::cout << "+\n";
    };
    auto trunc_label = [&](const std::string &name) -> std::string {
      if (name.size() <= static_cast<size_t>(col_width)) {
        return name;
      }
      return name.substr(0, static_cast<size_t>(col_width - 1));
    };

    std::cout << "VALUE_TABLE state=" << state_key
              << " protocol=" << proto_key << " frequency=" << fkey << "\n";
    print_border();
    std::cout << "| " << std::setw(iter_width) << "iter" << " ";
    for (const auto &name : step_names) {
      std::cout << "| " << std::setw(col_width) << trunc_label(name) << " ";
    }
    std::cout << "|\n";
    print_border();

    for (const auto &iter : iter_data) {
      std::cout << "| " << std::setw(iter_width) << iter["iter"].get<int>()
                << " ";
      for (const auto &name : step_names) {
        if (iter["steps"].contains(name)) {
          const double value = iter["steps"][name]["value"].get<double>();
          std::ostringstream os;
          os << std::scientific << std::setprecision(5) << value;
          std::cout << "| " << std::setw(col_width) << os.str() << " ";
        } else {
          std::cout << "| " << std::setw(col_width) << "N/A" << " ";
        }
      }
      std::cout << "|\n";
    }
    print_border();
  }

private:
  bool enabled_ = false;
  std::string filename_;
  json log_data_;

  // Current tuple
  std::string state_key_;
  std::string proto_key_;
  std::string freq_key_;

  json current_iter_values_;
  json current_iter_timing_;

  // ---------- Canonical keys ----------
  static std::string freq_key(double f) {
    std::ostringstream os;
    os << std::fixed << std::setprecision(3) << f;
    return os.str(); // "0.500"
  }
  static std::string protocol_key_(double thr) {
    // Compact scientific like "1e-06"
    std::ostringstream os;
    os.setf(std::ios::scientific);
    os << std::setprecision(0) << thr;
    std::string s = os.str();
    // normalize "1.e-06" -> "1e-06"
    auto dot = s.find('.');
    if (dot != std::string::npos)
      s.erase(dot, 1);
    return s;
  }

  // ---------- Structure helpers ----------
  void ensure_root_() {
    if (!log_data_.is_object())
      log_data_ = json::object();
  }
  void ensure_state_(const std::string &s) {
    if (!log_data_.contains(s))
      log_data_[s] = json::object();
    if (!log_data_[s].contains("protocols"))
      log_data_[s]["protocols"] = json::object();
  }
  void ensure_proto_(const std::string &s, const std::string &p) {
    ensure_state_(s);
    if (!log_data_[s]["protocols"].contains(p))
      log_data_[s]["protocols"][p] = json::object();
    if (!log_data_[s]["protocols"][p].contains("freqs"))
      log_data_[s]["protocols"][p]["freqs"] = json::object();
  }
  void ensure_freq_(const std::string &s, const std::string &p,
                    const std::string &f) {
    ensure_proto_(s, p);
    auto &freqnode = log_data_[s]["protocols"][p]["freqs"];
    if (!freqnode.contains(f))
      freqnode[f] = json::object();
  }

  // Returns a **reference** to the active (state,proto,freq) node
  json &node_ref_() {
    ensure_freq_(state_key_, proto_key_, freq_key_);
    return log_data_[state_key_]["protocols"][proto_key_]["freqs"][freq_key_];
  }

  bool exists_(const std::string &s, const std::string &p,
               const std::string &f) const {
    if (!log_data_.contains(s))
      return false;
    if (!log_data_.at(s).contains("protocols"))
      return false;
    if (!log_data_.at(s).at("protocols").contains(p))
      return false;
    if (!log_data_.at(s).at("protocols").at(p).contains("freqs"))
      return false;
    if (!log_data_.at(s).at("protocols").at(p).at("freqs").contains(f))
      return false;
    return true;
  }
};
class TimedValueLogger {
public:
  static void set_iteration_context(size_t iter) {
    iteration_index_context_ = iter;
  }

  static void clear_iteration_context() { iteration_index_context_.reset(); }

  static void set_console_enabled(bool enabled) { console_enabled_ = enabled; }

  [[nodiscard]] static bool console_enabled() { return console_enabled_; }

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

    if (logger_ && world_.rank() == 0) {
      logger_->log_timing(key_, wall, cpu);
    }

    if (world_.rank() == 0 && console_enabled_) {
      print_console_row(key_, wall, cpu, "-");
    }
  }

  template <typename T> void log(const T &value) {
    double wall = madness::wall_time() - start_wall_;
    double cpu = madness::cpu_time() - start_cpu_;

    if (logger_ && world_.rank() == 0) {
      logger_->log_value_and_time(key_, value, wall, cpu);
    }

    if (world_.rank() == 0 && console_enabled_) {
      print_console_row(key_, wall, cpu, value_to_string(value));
    }
  }
  template <typename T> void log(const std::vector<T> &values) {
    double wall = madness::wall_time() - start_wall_;
    double cpu = madness::cpu_time() - start_cpu_;

    if (logger_ && world_.rank() == 0) {
      for (size_t i = 0; i < values.size(); ++i) {
        std::string key_i = key_ + "[" + std::to_string(i) + "]";
        logger_->log_value_and_time(key_i, values[i], wall, cpu);
      }
    }

    if (world_.rank() == 0 && console_enabled_) {
      if (values.empty()) {
        print_console_row(key_, wall, cpu, "-");
      } else {
        for (size_t i = 0; i < values.size(); ++i) {
          const std::string key_i = key_ + "[" + std::to_string(i) + "]";
          print_console_row(key_i, wall, cpu, value_to_string(values[i]));
        }
      }
    }
  }

  template <typename T> void log_value(const T &value) {
    if (logger_ && world_.rank() == 0) {
      logger_->log_value(key_, value);
    }
  }

private:
  static constexpr int k_iter_width = 6;
  static constexpr int k_key_width = 30;
  static constexpr int k_time_width = 10;
  static constexpr int k_value_width = 22;

  static std::string clip_cell(const std::string &text, size_t width) {
    if (text.size() <= width) {
      return text;
    }
    return text.substr(0, width);
  }

  template <typename T> static std::string value_to_string(const T &value) {
    std::ostringstream os;
    if constexpr (std::is_floating_point_v<T>) {
      os << std::scientific << std::setprecision(6) << value;
    } else if constexpr (std::is_same_v<T, bool>) {
      os << std::boolalpha << value;
    } else {
      os << value;
    }
    return os.str();
  }

  static void print_console_border() {
    std::cout << "+" << std::string(k_iter_width + 2, '-')
              << "+" << std::string(k_key_width + 2, '-')
              << "+" << std::string(k_time_width + 2, '-')
              << "+" << std::string(k_time_width + 2, '-')
              << "+" << std::string(k_value_width + 2, '-') << "+\n";
  }

  static void print_console_header_if_needed() {
    if (console_header_printed_) {
      return;
    }
    std::cout << "STEP_TIMING_TABLE_COLUMNS iter,key,wall_s,cpu_s,value\n";
    print_console_border();
    std::cout << "| " << std::setw(k_iter_width) << "iter"
              << " | " << std::left << std::setw(k_key_width) << "key"
              << " | " << std::right << std::setw(k_time_width) << "wall_s"
              << " | " << std::right << std::setw(k_time_width) << "cpu_s"
              << " | " << std::left << std::setw(k_value_width) << "value"
              << " |\n";
    print_console_border();
    console_header_printed_ = true;
  }

  static void print_console_row(const std::string &key, double wall, double cpu,
                                const std::string &value_text) {
    print_console_header_if_needed();

    std::ostringstream wall_stream;
    std::ostringstream cpu_stream;
    wall_stream << std::fixed << std::setprecision(4) << wall;
    cpu_stream << std::fixed << std::setprecision(4) << cpu;
    const std::string iter_cell = iteration_index_context_.has_value()
                                      ? std::to_string(*iteration_index_context_)
                                      : "-";

    std::cout << "| " << std::right << std::setw(k_iter_width)
              << clip_cell(iter_cell, k_iter_width)
              << " | " << std::left << std::setw(k_key_width)
              << clip_cell(key, k_key_width)
              << " | " << std::right << std::setw(k_time_width)
              << clip_cell(wall_stream.str(), k_time_width)
              << " | " << std::right << std::setw(k_time_width)
              << clip_cell(cpu_stream.str(), k_time_width)
              << " | " << std::left << std::setw(k_value_width)
              << clip_cell(value_text, k_value_width) << " |\n";

    ++console_rows_printed_;
    if (console_rows_printed_ % 20 == 0) {
      print_console_border();
    }
  }

  inline static bool console_header_printed_ = false;
  inline static size_t console_rows_printed_ = 0;
  inline static bool console_enabled_ = true;
  inline static std::optional<size_t> iteration_index_context_ =
      std::nullopt;

  madness::World &world_;
  std::string key_;
  ResponseDebugLogger *logger_;
  double start_wall_, start_cpu_;
};
