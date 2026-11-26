#pragma once

// The text summary of a run folder (summary.txt): results.json formatted for
// people. It reads only the parsed results.json it is given and computes no
// value; it only chooses which values to show and how to print them. Nothing
// ever reads the summary back.
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace madness::iteration {

namespace summary_text {

using json = nlohmann::json;

inline constexpr std::size_t max_width = 100;
inline constexpr double noise = 1e-10;
inline constexpr const char* above_target_note =
    "  * above target (driven states can pass under the scaled rule)";

inline std::string printf_text(const char* format, double x) {
    char text[64];
    std::snprintf(text, sizeof text, format, x);
    return text;
}

/// A convergence value or target: scientific, three significant figures.
inline std::string metric_text(double x) { return printf_text("%.2e", x); }

/// An energy in hartree: eight decimals. Numerical noise prints as 0.
inline std::string energy_text(double x) {
    if (std::abs(x) < noise) return "0";
    return printf_text("%.8f", x);
}

/// Any other value: six significant figures, fixed form unless very large or
/// small. Numerical noise prints as 0.
inline std::string value_text(double x) {
    if (std::abs(x) < noise) return "0";
    if (!std::isfinite(x)) return printf_text("%g", x);
    const int exponent = static_cast<int>(std::floor(std::log10(std::abs(std::stod(printf_text("%.5e", x))))));
    if (exponent < -4 || exponent > 5) return printf_text("%.5e", x);
    return printf_text(("%." + std::to_string(5 - exponent) + "f").c_str(), x);
}

/// A response property value: six decimals. Numerical noise prints as 0.
inline std::string decimals_text(double x) {
    if (std::abs(x) < noise) return "0";
    return printf_text("%.6f", x);
}

/// A response frequency, as the state folder names print it.
inline std::string frequency_text(double x) { return printf_text("%.5f", x); }

inline bool is_energy_unit(const json& unit) { return unit.is_string() && unit.get<std::string>() == "hartree"; }

/// A number as its kind prints: an energy, a response value (`decimals`), or
/// any other value. Integers print as they are; a missing value prints as -.
inline std::string number_text(const json& x, bool energy, bool decimals = false) {
    if (!x.is_number()) return x.is_null() ? "-" : x.is_string() ? x.get<std::string>() : x.dump();
    if (x.is_number_integer()) return x.dump();
    if (energy) return energy_text(x.get<double>());
    return decimals ? decimals_text(x.get<double>()) : value_text(x.get<double>());
}

/// Response property values print to six decimals.
inline bool decimal_values(const json& entry) { return entry.value("application", std::string()) == "molresponse"; }

/// The unit of `key` when `units` maps keys to units, or `units` itself.
inline json unit_of(const json& units, const std::string& key) {
    if (units.is_object()) return units.contains(key) ? units.at(key) : json();
    return units;
}

inline std::string unit_label(const json& units) {
    if (units.is_string()) return units.get<std::string>();
    return {};
}

/// Lines of `text` at most max_width long, broken at spaces; the first keeps
/// the leading spaces of `text`, the others start with `indent`.
inline std::vector<std::string> wrap(const std::string& text, const std::string& indent) {
    std::vector<std::string> lines;
    const auto lead = std::min(text.find_first_not_of(' '), text.size());
    std::string line = text.substr(0, lead);
    std::size_t prefix = line.size(); // the line holds no word yet while line.size() == prefix
    std::size_t start = lead;
    while (start < text.size()) {
        auto end = text.find(' ', start);
        if (end == std::string::npos) end = text.size();
        const auto word = text.substr(start, end - start);
        start = end + 1;
        if (word.empty()) continue;
        if (line.size() > prefix && line.size() + 1 + word.size() > max_width) {
            lines.push_back(line);
            line = indent;
            prefix = line.size();
        }
        line += (line.size() > prefix ? " " : "") + word;
    }
    if (line.size() > prefix) lines.push_back(line);
    return lines;
}

/// A left-aligned table: each column as wide as its widest cell plus two.
class Table {
    std::vector<std::vector<std::string>> rows_;
public:
    explicit Table(std::vector<std::string> header) { rows_.push_back(std::move(header)); }
    void add(std::vector<std::string> row) { row.resize(rows_.front().size()); rows_.push_back(std::move(row)); }
    [[nodiscard]] std::size_t width() const {
        std::size_t total = 0;
        // The two-space indent stands in for the padding the last column lacks.
        for (std::size_t c = 0; c < rows_.front().size(); ++c) total += column_width(c);
        return total;
    }
    [[nodiscard]] std::size_t column_width(std::size_t c) const {
        std::size_t w = 0;
        for (const auto& row : rows_) w = std::max(w, row[c].size());
        return w + 2;
    }
    void write(std::string& out) const {
        for (const auto& row : rows_) {
            std::string line = "  ";
            for (std::size_t c = 0; c < row.size(); ++c) {
                line += row[c];
                if (c + 1 < row.size()) line += std::string(column_width(c) - row[c].size(), ' ');
            }
            while (!line.empty() && line.back() == ' ') line.pop_back();
            out += line + "\n";
        }
    }
};

/// The targets that apply to a value: flat, or by kind when the entry names
/// kinds (a value keyed by root is an excited-state value).
inline json targets_for(const json& entry, bool excited) {
    auto targets = entry.value("targets", json::object());
    if (targets.contains("driven_states") || targets.contains("excited_states"))
        return targets.value(excited ? "excited_states" : "driven_states", json::object());
    return targets;
}

/// The metrics of the `convergence` objects of `items`, sorted as the file sorts them.
inline std::vector<std::string> convergence_metrics(const std::vector<json>& items) {
    std::vector<std::string> names;
    for (const auto& item : items) {
        const auto convergence = item.value("convergence", json::object());
        for (const auto& [name, value] : convergence.items())
            if (std::find(names.begin(), names.end(), name) == names.end()) names.push_back(name);
    }
    std::sort(names.begin(), names.end());
    return names;
}

/// One convergence cell, marked when above its target.
inline std::string convergence_cell(const json& convergence, const json& targets, const std::string& metric,
                                    bool& marked) {
    if (!convergence.is_object() || !convergence.contains(metric)) return "-";
    const double value = convergence.at(metric).get<double>();
    auto text = metric_text(value);
    if (targets.contains(metric) && targets.at(metric).is_number() && value > targets.at(metric).get<double>()) {
        text += "*";
        marked = true;
    }
    return text;
}

/// "metric value   metric value ..." for a convergence or targets object.
inline std::string metric_list(const json& values, const json& targets, bool& marked) {
    std::string text;
    for (const auto& [metric, value] : values.items()) {
        if (!value.is_number()) continue;
        if (!text.empty()) text += "   ";
        text += metric + " " + (targets.is_null() ? metric_text(value.get<double>())
                                                  : convergence_cell(values, targets, metric, marked));
    }
    return text;
}

/// The reasons of an availability object as (component, reason); a nested
/// object (reasons by coordinate, then component) joins its keys with a dot.
/// Depth first, in key order.
inline void availability_lines(const json& availability, std::vector<std::pair<std::string, std::string>>& lines) {
    if (!availability.is_object()) return;
    std::vector<std::pair<std::string, const json*>> pending;
    const auto push_children = [&](const std::string& prefix, const json& object) {
        std::vector<std::pair<std::string, const json*>> children;
        for (const auto& [key, reason] : object.items())
            children.emplace_back(prefix.empty() ? key : std::string(prefix).append(".").append(key), &reason);
        pending.insert(pending.end(), children.rbegin(), children.rend());
    };
    push_children({}, availability);
    while (!pending.empty()) {
        const auto [name, reason] = pending.back();
        pending.pop_back();
        if (reason->is_string()) lines.emplace_back(name, reason->get<std::string>());
        else if (reason->is_object() && reason->contains("reason") && reason->at("reason").is_string())
            lines.emplace_back(name, reason->at("reason").get<std::string>());
        else if (reason->is_object()) push_children(name, *reason);
        else lines.emplace_back(name, reason->dump());
    }
}

/// The distinct reasons of an availability object, in first-seen order,
/// joined with "; ".
inline std::string availability_reason(const json& availability) {
    std::vector<std::pair<std::string, std::string>> lines;
    availability_lines(availability, lines);
    std::vector<std::string> reasons;
    for (const auto& [name, reason] : lines)
        if (std::find(reasons.begin(), reasons.end(), reason) == reasons.end()) reasons.push_back(reason);
    if (reasons.empty()) return "no reason recorded";
    std::string text;
    for (const auto& reason : reasons) text += (text.empty() ? "" : "; ") + reason;
    return text;
}

/// The keys that say what an item belongs to, in the order they print.
inline const std::vector<std::string>& identifying_keys() {
    static const std::vector<std::string> keys{"root", "frequency", "frequency_b", "frequency_c"};
    return keys;
}

inline std::string item_label(const json& item) {
    std::string label;
    for (const auto& key : identifying_keys()) {
        if (!item.contains(key)) continue;
        if (!label.empty()) label += " ";
        const auto& v = item.at(key);
        label += key + " " + (v.is_number() ? frequency_text(v.get<double>()) : number_text(v, false));
    }
    return label;
}

/// The NOT COMPLETE lines, one per item with its reason: every property or
/// item that is not complete, and the run's stop reason, dropped work and
/// stalled states.
inline std::vector<std::string> not_complete(const json& results) {
    std::vector<std::string> lines;
    const auto add = [&](const std::string& head, const json& availability) {
        lines.push_back("  " + head + ": " + availability_reason(availability));
    };
    for (const auto& entry : results.value("calculations", json::array())) {
        const auto id = entry.value("calculation_id", std::string("?"));
        if (entry.contains("run")) {
            const auto& run = entry.at("run");
            const auto stop = run.value("stop_reason", std::string("unknown"));
            if (stop != "complete") lines.push_back(std::string("  ").append(id).append(": run stopped: ").append(stop));
            for (const auto& state : run.value("stalled", json::array()))
                lines.push_back("  " + id + ": stalled: " + (state.is_string() ? state.get<std::string>() : state.dump()));
            for (const auto& work : run.value("dropped_work", json::array()))
                lines.push_back("  " + id + ": dropped: " + work.value("id", std::string("?")) + ": " +
                                work.value("reason", std::string("no reason recorded")));
        }
        const auto properties = entry.value("properties", json::object());
        for (const auto& [name, property] : properties.items()) {
            const auto head = std::string(id).append(" ").append(name);
            if (property.contains("status") && property.at("status") != "complete")
                add(head, property.value("availability", json::object()));
            else if (!property.contains("status") && property.contains("availability") &&
                     !property.at("availability").empty())
                add(head, property.at("availability"));
            for (const auto& item : property.value("values", json::array()))
                if (item.value("status", std::string("complete")) != "complete")
                    add(head + " " + item_label(item), item.value("availability", json::object()));
        }
    }
    return lines;
}

inline std::string entry_title(const json& entry) {
    const auto id = entry.value("calculation_id", std::string("?"));
    std::string title = id == "ground_state" ? "Ground state" : id == "response" ? "Response" : id;
    std::string inner = entry.value("application", std::string());
    if (entry.contains("model") && entry.at("model").is_string())
        inner += (inner.empty() ? "" : ", ") + entry.at("model").get<std::string>();
    if (!inner.empty()) title += " (" + inner + ")";
    if (entry.contains("protocol")) {
        const auto& p = entry.at("protocol");
        if (p.contains("threshold")) title += ", level " + printf_text("%.1e", p.at("threshold").get<double>());
        if (p.contains("k")) title += " k " + std::to_string(p.at("k").get<int>());
    }
    return title;
}

/// Write `table` with its convergence columns; when that is wider than the
/// page, the convergence moves to a second table keyed by the item columns.
struct ConvergenceTable {
    std::vector<std::string> header;       // item columns, value columns
    std::size_t key_columns = 0;           // how many leading columns identify the item
    std::vector<std::string> metrics;
    std::vector<std::vector<std::string>> rows;
    std::vector<std::vector<std::string>> convergence; // per row; empty for a continuation row

    void write(std::string& out) const {
        Table joined([&] { auto h = header; h.insert(h.end(), metrics.begin(), metrics.end()); return h; }());
        for (std::size_t r = 0; r < rows.size(); ++r) {
            auto row = rows[r];
            const auto& c = convergence[r];
            row.insert(row.end(), c.begin(), c.end());
            joined.add(row);
        }
        if (metrics.empty() || joined.width() <= max_width) {
            joined.write(out);
            return;
        }
        Table values(header);
        for (const auto& row : rows) values.add(row);
        values.write(out);
        out += "\n";
        std::vector<std::string> h(header.begin(), header.begin() + static_cast<std::ptrdiff_t>(key_columns));
        h.insert(h.end(), metrics.begin(), metrics.end());
        Table split(h);
        for (std::size_t r = 0; r < rows.size(); ++r) {
            if (convergence[r].empty()) continue;
            std::vector<std::string> row(rows[r].begin(), rows[r].begin() + static_cast<std::ptrdiff_t>(key_columns));
            row.insert(row.end(), convergence[r].begin(), convergence[r].end());
            split.add(row);
        }
        split.write(out);
    }
};

/// The leading cells of an item's first row and of its continuation rows.
inline std::vector<std::string> key_cells(const json& item, const std::vector<std::string>& keys) {
    std::vector<std::string> cells;
    for (const auto& key : keys) {
        const auto& v = item.contains(key) ? item.at(key) : json();
        cells.push_back(v.is_number() ? frequency_text(v.get<double>()) : number_text(v, false));
    }
    return cells;
}

inline std::vector<std::string> present_keys(const json& items) {
    std::vector<std::string> keys;
    for (const auto& key : identifying_keys())
        for (const auto& item : items)
            if (item.contains(key)) {
                keys.push_back(key);
                break;
            }
    return keys;
}

/// Rows of one item: one per component (or one for a single value), with the
/// convergence on the first.
template <class Rows>
inline void item_rows(ConvergenceTable& table, const json& item, const std::vector<std::string>& keys,
                      const Rows& component_rows, const json& targets, bool& marked) {
    const auto keys_text = key_cells(item, keys);
    auto components = component_rows(item);
    if (components.empty()) components.push_back(std::vector<std::string>(table.header.size() - keys.size(), "-"));
    for (std::size_t i = 0; i < components.size(); ++i) {
        auto row = i == 0 ? keys_text : std::vector<std::string>(keys.size());
        row.insert(row.end(), components[i].begin(), components[i].end());
        std::vector<std::string> convergence;
        if (i == 0)
            for (const auto& metric : table.metrics)
                convergence.push_back(convergence_cell(item.value("convergence", json()), targets, metric, marked));
        table.rows.push_back(std::move(row));
        table.convergence.push_back(std::move(convergence));
    }
}

/// A components object printed as one row per component, if every component
/// is a number and there are at most nine (a 3 by 3 tensor).
inline bool small_components(const json& components) {
    if (!components.is_object() || components.size() > 9) return false;
    return std::all_of(components.begin(), components.end(), [](const json& v) { return v.is_number(); });
}

inline void write_title(std::string& out, const std::string& name, const std::string& unit,
                        const std::string& note = {}) {
    out += "\n";
    for (const auto& line :
         wrap("  " + name + (unit.empty() ? "" : " (" + unit + ")") + (note.empty() ? "" : ", " + note), "    "))
        out += line + "\n";
}

/// A property with a `values` list whose items hold `value` or a small
/// `components` (polarizability, hyperpolarizability, and unknown ones).
inline void write_values_property(std::string& out, const json& entry, const std::string& name,
                                  const json& property, std::vector<std::string>& not_printed, bool& marked) {
    const auto& items = property.at("values");
    const auto keys = present_keys(items);
    const auto unit = property.value("units", json());
    const bool energy = is_energy_unit(unit);
    const bool decimals = decimal_values(entry);
    bool components = false, large = false;
    for (const auto& item : items) {
        if (item.contains("components")) components = true;
        if (item.contains("components") && !small_components(item.at("components")) &&
            !item.at("components").empty())
            large = true;
        if (item.contains("value") && !item.at("value").is_number()) large = true;
    }
    ConvergenceTable table;
    table.header = keys;
    table.key_columns = keys.size();
    if (!large) {
        if (components) table.header.emplace_back("component");
        table.header.emplace_back("value");
    }
    table.metrics = convergence_metrics(std::vector<json>(items.begin(), items.end()));
    const auto rows = [&](const json& item) {
        std::vector<std::vector<std::string>> r;
        if (large) return r;
        if (item.contains("components")) {
            // Every requested component: those with a value, and those whose
            // availability says why they have none (printed as -).
            auto present = item.at("components");
            std::vector<std::pair<std::string, std::string>> missing;
            availability_lines(item.value("availability", json::object()), missing);
            for (const auto& [component, reason] : missing)
                if (!present.contains(component)) present[component] = nullptr;
            for (const auto& [component, value] : present.items())
                r.push_back({component, number_text(value, energy, decimals)});
        } else if (item.contains("value"))
            r.push_back(components ? std::vector<std::string>{"", number_text(item.at("value"), energy, decimals)}
                                   : std::vector<std::string>{number_text(item.at("value"), energy, decimals)});
        return r;
    };
    for (const auto& item : items)
        item_rows(table, item, keys, rows, targets_for(entry, item.contains("root")), marked);
    write_title(out, name, unit_label(unit));
    table.write(out);
    if (large) not_printed.push_back(name + " components");
}

/// excited_states: one row per root with its energy, oscillator strength and
/// convergence, then the transition dipoles. The quadrupoles are left out.
inline void write_excited_states(std::string& out, const json& entry, const json& property,
                                 std::vector<std::string>& not_printed, bool& marked) {
    const auto& items = property.at("values");
    const auto units = property.value("units", json::object());
    const bool decimals = decimal_values(entry);
    ConvergenceTable table;
    table.header = {"root", "energy", "osc. strength"};
    table.key_columns = 1;
    table.metrics = convergence_metrics(std::vector<json>(items.begin(), items.end()));
    const auto rows = [&](const json& item) {
        return std::vector<std::vector<std::string>>{
            {number_text(item.value("excitation_energy", json()), is_energy_unit(unit_of(units, "excitation_energy")),
                         decimals),
             number_text(item.value("oscillator_strength", json()), false, decimals)}};
    };
    for (const auto& item : items) item_rows(table, item, {"root"}, rows, targets_for(entry, true), marked);
    const auto energy_unit = unit_label(unit_of(units, "excitation_energy"));
    write_title(out, "excited states", energy_unit.empty() ? "" : "energy in " + energy_unit);
    table.write(out);

    Table dipoles({"root", "x", "y", "z"});
    bool any = false;
    for (const auto& item : items) {
        const auto d = item.value("transition_dipole", json::object());
        if (d.empty()) continue;
        any = true;
        dipoles.add({number_text(item.value("root", json()), false), number_text(d.value("x", json()), false, decimals),
                     number_text(d.value("y", json()), false, decimals),
                     number_text(d.value("z", json()), false, decimals)});
    }
    if (any) {
        write_title(out, "transition dipole", unit_label(unit_of(units, "transition_dipole")));
        dipoles.write(out);
    }
    bool quadrupole = false;
    for (const auto& item : items) quadrupole = quadrupole || item.contains("transition_quadrupole");
    if (quadrupole) not_printed.emplace_back("transition_quadrupole");
}

/// two_photon_absorption: per root, the transition tensor components and the
/// observables.
inline void write_two_photon(std::string& out, const json& entry, const json& property, bool& marked) {
    const auto& items = property.at("values");
    const auto units = property.value("units", json::object());
    const bool decimals = decimal_values(entry);
    ConvergenceTable table;
    table.header = {"root", "frequency", "component", "value"};
    table.key_columns = 2;
    table.metrics = convergence_metrics(std::vector<json>(items.begin(), items.end()));
    const auto rows = [&](const json& item) {
        std::vector<std::vector<std::string>> r;
        if (item.contains("energy"))
            r.push_back({"energy", number_text(item.at("energy"), is_energy_unit(unit_of(units, "energy")), decimals)});
        const auto components = item.value("components", json::object());
        const auto observables = item.value("observables", json::object());
        for (const auto& [component, value] : components.items())
            r.push_back({component, number_text(value, false, decimals)});
        for (const auto& [name, value] : observables.items()) r.push_back({name, number_text(value, false, decimals)});
        return r;
    };
    for (const auto& item : items) item_rows(table, item, {"root", "frequency"}, rows, targets_for(entry, true), marked);
    std::string unit;
    for (const auto& [key, u] : units.items())
        if (u.is_string()) unit += (unit.empty() ? "" : ", ") + key + " " + u.get<std::string>();
    const auto status = property.value("scientific_status", std::string());
    write_title(out, "two_photon_absorption", unit, status);
    table.write(out);
}

/// The unit text of a units entry, up to any explanation in parentheses.
inline std::string short_unit(const json& unit) {
    if (!unit.is_string()) return {};
    auto text = unit.get<std::string>();
    const auto open = text.find(" (");
    return open == std::string::npos ? text : text.substr(0, open);
}

/// raman_spectrum: per laser frequency, one row per mode.
inline void write_raman(std::string& out, const json& entry, const json& property, bool& marked) {
    const auto& items = property.at("values");
    const auto units = property.value("units", json::object());
    const bool decimals = decimal_values(entry);
    const std::vector<std::string> columns{"freq_cm1", "pol_int", "dep_ratio"};
    ConvergenceTable table;
    table.header = {"frequency", "mode"};
    table.header.insert(table.header.end(), columns.begin(), columns.end());
    table.key_columns = 1;
    table.metrics = convergence_metrics(std::vector<json>(items.begin(), items.end()));
    const auto rows = [&](const json& item) {
        std::vector<std::vector<std::string>> r;
        for (const auto& mode : item.value("value", json::array())) {
            std::vector<std::string> row{number_text(mode.value("mode", json()), false)};
            for (const auto& c : columns) row.push_back(number_text(mode.value(c, json()), false, decimals));
            r.push_back(std::move(row));
        }
        return r;
    };
    for (const auto& item : items) item_rows(table, item, {"frequency"}, rows, targets_for(entry, false), marked);
    std::string unit;
    for (const auto& c : columns)
        if (!short_unit(unit_of(units, c)).empty()) unit += (unit.empty() ? "" : ", ") + c + " " + short_unit(units.at(c));
    write_title(out, "raman_spectrum", unit);
    table.write(out);
}

/// polarizability_gradient: per laser frequency, convergence only; its
/// components (3N coordinates by 9) are left to results.json.
inline void write_values_convergence_only(std::string& out, const json& entry, const std::string& name,
                                          const json& property, std::vector<std::string>& not_printed, bool& marked) {
    const auto& items = property.at("values");
    const auto keys = present_keys(items);
    ConvergenceTable table;
    table.header = keys;
    table.key_columns = keys.size();
    table.metrics = convergence_metrics(std::vector<json>(items.begin(), items.end()));
    const auto rows = [](const json&) { return std::vector<std::vector<std::string>>{std::vector<std::string>{}}; };
    for (const auto& item : items)
        item_rows(table, item, keys, rows, targets_for(entry, item.contains("root")), marked);
    write_title(out, name, unit_label(property.value("units", json())));
    table.write(out);
    not_printed.push_back(name + " components");
}

/// A single-answer property whose value is not a number: a table of its
/// parts when they are small, otherwise named as not printed.
inline void write_structured_value(std::string& out, const json& entry, const std::string& name,
                                   const json& property, std::vector<std::string>& not_printed) {
    const auto unit = property.value("units", json());
    const bool energy = is_energy_unit(unit);
    const bool decimals = decimal_values(entry);
    const auto title = [&] { write_title(out, name, unit_label(unit)); };
    const auto& value = property.at("value");
    if (name == "vibrational_analysis" && value.is_object()) {
        const auto fu = unit_label(unit_of(unit, "frequencies"));
        write_title(out, name, fu.empty() ? "" : "frequency in " + fu);
        Table table({"mode", "frequency", "reduced_mass"});
        const auto f = value.value("frequencies", json::array());
        const auto m = value.value("reduced_mass", json::array());
        const bool energy_f = is_energy_unit(unit_of(unit, "frequencies"));
        for (std::size_t i = 0; i < f.size(); ++i)
            table.add({std::to_string(i), number_text(f.at(i), energy_f, decimals),
                       i < m.size() ? number_text(m.at(i), false, decimals) : "-"});
        table.write(out);
        for (const auto& [key, part] : value.items())
            if (key != "frequencies" && key != "reduced_mass") not_printed.push_back(std::string(name).append(" ").append(key));
        return;
    }
    if (value.is_object() && std::all_of(value.begin(), value.end(), [](const json& v) { return v.is_number(); })) {
        title();
        Table table({"term", "value"});
        // energy_components in the order the terms are usually read; any
        // other term follows in key order.
        std::vector<std::string> order;
        if (name == "energy_components")
            for (const auto* key : {"kinetic", "nuclear_attraction", "coulomb", "exchange_correlation", "nuclear_repulsion"})
                if (value.contains(key)) order.emplace_back(key);
        for (const auto& [key, v] : value.items())
            if (std::find(order.begin(), order.end(), key) == order.end()) order.push_back(key);
        for (const auto& key : order) table.add({key, number_text(value.at(key), energy, decimals)});
        table.write(out);
        return;
    }
    if (value.is_object() && std::all_of(value.begin(), value.end(), [](const json& v) {
            return v.is_array() && std::all_of(v.begin(), v.end(), [](const json& x) { return x.is_number(); });
        })) {
        title();
        Table table({"spin", "orbital", "value"});
        for (const auto& [key, list] : value.items())
            for (std::size_t i = 0; i < list.size(); ++i)
                table.add({key, std::to_string(i), number_text(list.at(i), energy, decimals)});
        table.write(out);
        return;
    }
    if (value.is_array() && !value.empty() && value.size() % 3 == 0 &&
        std::all_of(value.begin(), value.end(), [](const json& x) { return x.is_number(); })) {
        title();
        const bool per_atom = value.size() > 3;
        Table table(per_atom ? std::vector<std::string>{"atom", "x", "y", "z"} : std::vector<std::string>{"x", "y", "z"});
        for (std::size_t i = 0; i < value.size(); i += 3) {
            std::vector<std::string> row;
            if (per_atom) row.push_back(std::to_string(i / 3));
            for (std::size_t j = 0; j < 3; ++j) row.push_back(number_text(value.at(i + j), energy, decimals));
            table.add(row);
        }
        table.write(out);
        return;
    }
    not_printed.push_back(name);
}

/// A single-answer property shown in the shared table: its value is a number,
/// or it has no value (one line with -).
inline bool single_number(const json& property) {
    return !property.contains("value") || property.at("value").is_number();
}

inline void write_entry(std::string& out, const json& entry) {
    out += "\n" + entry_title(entry) + "\n";
    const auto targets = entry.value("targets", json::object());
    bool marked = false;
    if (entry.contains("convergence") && entry.at("convergence").is_object())
        out += "  reached   " + metric_list(entry.at("convergence"), targets_for(entry, false), marked) + "\n";
    if (targets.contains("driven_states") || targets.contains("excited_states")) {
        bool first = true;
        for (const auto* kind : {"driven_states", "excited_states"}) {
            if (!targets.contains(kind)) continue;
            const std::string label = std::string(kind) == "driven_states" ? "driven:  " : "excited: ";
            out += std::string(first ? "  targets   " : "            ") + label +
                   metric_list(targets.at(kind), json(), marked) + "\n";
            first = false;
        }
    } else if (!targets.empty()) {
        out += "  targets   " + metric_list(targets, json(), marked) + "\n";
    }

    std::vector<std::string> not_printed;
    const bool decimals = decimal_values(entry);
    const auto properties = entry.value("properties", json::object());
    // Single numbers share one table.
    Table scalars({"property", "value", "units"});
    bool any_scalar = false;
    for (const auto& [name, property] : properties.items()) {
        if (property.contains("values") || !single_number(property)) continue;
        any_scalar = true;
        const auto unit = property.value("units", json());
        scalars.add({name,
                     property.contains("value") ? number_text(property.at("value"), is_energy_unit(unit), decimals)
                                                : "-",
                     unit_label(unit)});
    }
    if (any_scalar) {
        out += "\n";
        scalars.write(out);
    }
    for (const auto& [name, property] : properties.items()) {
        if (property.contains("values")) {
            if (name == "excited_states") write_excited_states(out, entry, property, not_printed, marked);
            else if (name == "two_photon_absorption") write_two_photon(out, entry, property, marked);
            else if (name == "raman_spectrum") write_raman(out, entry, property, marked);
            else if (name == "polarizability_gradient")
                write_values_convergence_only(out, entry, name, property, not_printed, marked);
            else write_values_property(out, entry, name, property, not_printed, marked);
            continue;
        }
        if (!single_number(property)) write_structured_value(out, entry, name, property, not_printed);
    }
    if (!not_printed.empty()) {
        std::string list;
        for (const auto& n : not_printed) list += (list.empty() ? "" : ", ") + n;
        out += "\n";
        for (const auto& line : wrap("  not printed (see results.json): " + list, "    ")) out += line + "\n";
    }
    if (entry.contains("run")) {
        const auto& run = entry.at("run");
        const auto dropped = run.value("dropped_work", json::array()).size();
        const auto stalled = run.value("stalled", json::array()).size();
        const auto unconverged = run.value("unconverged", json::array()).size();
        const auto states = [](std::size_t n) { return std::to_string(n) + (n == 1 ? " state" : " states"); };
        const auto stop = run.value("stop_reason", std::string("unknown"));
        // The stop reason, with the count of states that ended without
        // meeting the request when there are any.
        const auto reason = unconverged == 0 ? stop
                            : stop == "unconverged" ? stop + " (" + states(unconverged) + ")"
                                                    : stop + ", unconverged (" + states(unconverged) + ")";
        out += "\n  run: " + reason + ". " +
               (dropped == 0 ? std::string("Nothing dropped.") : std::to_string(dropped) + " dropped.") + " " +
               (stalled == 0 ? std::string("Nothing stalled.") : std::to_string(stalled) + " stalled.") + "\n";
    }
    // One note for every * in the section.
    if (marked) {
        out += "\n";
        for (const auto& line : wrap(above_target_note, "    ")) out += line + "\n";
    }
}

/// The code version as people quote it: a full commit hash is shortened to
/// nine characters, keeping any suffix such as -dirty.
inline std::string short_code_version(const std::string& version) {
    const auto hex = version.find_first_not_of("0123456789abcdef");
    const auto length = hex == std::string::npos ? version.size() : hex;
    if (length != 40) return version;
    return version.substr(0, 9) + version.substr(40);
}

/// "2026-11-26T01:51:12.345Z" as "2026-11-26 01:51 UTC".
inline std::string launch_text(const std::string& at) {
    if (at.size() >= 16 && at[10] == 'T' && at.back() == 'Z') return at.substr(0, 10) + " " + at.substr(11, 5) + " UTC";
    return at;
}

} // namespace summary_text

/// The text of summary.txt for a parsed results.json. It reads nothing else.
inline std::string results_summary(const nlohmann::json& results) {
    using namespace summary_text;
    std::string out;
    const auto incomplete = not_complete(results);
    if (!incomplete.empty()) {
        out += "NOT COMPLETE\n";
        for (const auto& line : incomplete)
            for (const auto& wrapped : wrap(line, "        ")) out += wrapped + "\n";
        out += "\n";
    }
    out += "madqc summary\n";
    const auto provenance = results.value("provenance", json::object());
    if (provenance.contains("code_version"))
        out += "  code      " + short_code_version(provenance.at("code_version").get<std::string>()) + "\n";
    if (provenance.contains("launched_at"))
        out += "  launched  " + launch_text(provenance.at("launched_at").get<std::string>()) + "\n";
    if (provenance.contains("mpi_processes") && provenance.contains("threads_per_process")) {
        const auto np = provenance.at("mpi_processes").get<long>();
        out += "  run       " + std::to_string(np) + (np == 1 ? " process x " : " processes x ") +
               std::to_string(provenance.at("threads_per_process").get<long>()) + " threads\n";
    }
    for (const auto& entry : results.value("calculations", json::array())) write_entry(out, entry);
    return out;
}

} // namespace madness::iteration
