/*
 * QCCalculationParametersBase.cpp
 *
 *  Created on: 27 Jun 2019
 *      Author: fbischoff
 */

#include"QCCalculationParametersBase.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include"worldgop.h"



namespace madness {

namespace {
bool is_missing_datagroup_exception(const std::exception& ex) {
    constexpr const char* kPrefix = "position_stream: failed to locate ";
    const char* what = ex.what();
    if (what == nullptr) return false;
    const std::string message(what);
    return message.rfind(kPrefix, 0) == 0;
}
} // namespace




std::string QCParameter::print_line(const std::string& key, int key_width, int value_width) const {
    const int comment_column = key_width + 2 + value_width + 1;
    const std::string continuation = std::string(comment_column, ' ') + "# " + std::string(10, ' ');
    std::ostringstream output;
    output << std::right << std::setw(key_width) << key << "  "
           << std::left << std::setw(value_width) << get_value();
    // Values remain on one line so a printed deck is valid input. Long values
    // put their explanation on the next line instead of shifting the table.
    if (get_value().size() > static_cast<std::size_t>(value_width))
        output << '\n' << std::string(comment_column, ' ');
    else output << ' ';
    output << "# " << std::setw(10) << print_precedence();
    std::string description = get_comment();
    if (!allowed_values.empty()) {
        using madness::operators::operator<<;
        std::ostringstream allowed;
        allowed << allowed_values;
        description += "  Allowed: " + allowed.str();
    }
    const auto width = std::max(20, 120 - comment_column - 12);
    std::istringstream words(description);
    std::string word;
    int used = 0;
    while (words >> word) {
        if (used && used + 1 + static_cast<int>(word.size()) > width) {
            output << '\n' << continuation;
            used = 0;
        }
        if (used) { output << ' '; ++used; }
        output << word;
        used += static_cast<int>(word.size());
    }
    auto result = output.str();
    const auto last = result.find_last_not_of(' ');
    return result.substr(0, last + 1);
}

/// print all parameters
void QCCalculationParametersBase::print(const std::string header,
		const std::string footer) const {

	std::string body=print_to_string({"all"});
	if (header.size()>0) madness::print(header);
	::madness::print(body);
	if (footer.size()>0) madness::print(footer);
}

std::string QCCalculationParametersBase::print_to_string(const std::list<std::string> precedence) const {

	// sort parameters according to increasing print_order
	typedef std::tuple<int,std::string,QCParameter> keyvalT;
	std::list<keyvalT> list;
	for (auto& p : parameters)
        if (visible_input(p.first)) list.push_back(keyvalT(p.second.get_print_order(),p.first,p.second));
	list.sort([](const keyvalT& first, const keyvalT& second) {return std::get<0>(first) < std::get<0>(second);});

	// check if we have to print only non-default parameters
	bool print_all = std::find(precedence.begin(), precedence.end(), "all") != precedence.end();
	bool print_default = std::find(precedence.begin(), precedence.end(), "default") != precedence.end();
	bool print_derived = std::find(precedence.begin(), precedence.end(), "derived") != precedence.end();
	bool print_defined = std::find(precedence.begin(), precedence.end(), "defined") != precedence.end();

    int key_width = 20;
    for (const auto& p : list) key_width = std::max(key_width, static_cast<int>(std::get<1>(p).size()));
    std::stringstream ss;
    int counter=0;
	for (auto& p : list) {
		const QCParameter& param=std::get<2>(p);
		bool skip=
			(not print_default and param.precedence == QCParameter::def) or
			(not print_derived and param.precedence == QCParameter::derived) or
			(not print_defined and param.precedence == QCParameter::defined);
		if (print_all) skip=false;
		if (skip) continue;
		if ((counter++)>0) ss << std::endl; // no newline at the very end
        const auto row = param.print_line(std::get<1>(p), key_width);
        // A full listing documents defaults/derived settings without turning
        // them into explicit requests when copied into an input deck. Callers
        // explicitly selecting a precedence retain their serialization mode.
        if (print_all) {
            std::istringstream lines(row);
            std::string line;
            bool first = true;
            while (std::getline(lines, line)) {
                if (!first) ss << '\n';
                ss << (param.precedence == QCParameter::defined ? "  " : "# ") << line;
                first = false;
            }
        } else ss << row;
	}
	return ss.str();
}


bool QCCalculationParametersBase::file_exists(World& world, std::string filename) const {
    bool file_exists = true;
    if (world.rank() == 0) {
        std::ifstream ifs(filename);
        if (not ifs.is_open()) file_exists=false;
        ifs.close();
    }
    world.gop.broadcast_serializable(file_exists, 0);
    return file_exists;
}

void QCCalculationParametersBase::read_collectively(World& world, const std::function<void()>& reader) {
    std::string message;
    int error = 0; // 1: invalid_argument, 2: runtime_error
    if (world.rank() == 0) {
        const auto previous = parameters;
        try {
            reader();
        } catch (const std::invalid_argument& e) {
            message = e.what();
            error = 1;
        } catch (const std::exception& e) {
            message = e.what();
            error = 2;
        }
        if (error) parameters = previous;
    }
    world.gop.broadcast(error, 0);
    world.gop.broadcast_serializable(message, 0);
    if (error == 1) throw std::invalid_argument(message);
    if (error == 2) throw std::runtime_error(message);
    world.gop.broadcast_serializable(*this, 0);
}

void QCCalculationParametersBase::read_input_and_commandline_options(
        World& world, const commandlineparser& parser,
        const std::string tag) { // NOLINT(performance-unnecessary-value-param): preserve the existing signature.
    read_collectively(world, [&] {
        const auto filename = parser.value_raw("input");
        if (parser.key_exists("user_defined_input_file") && !std::ifstream(filename).is_open())
            throw std::invalid_argument("could not find user-defined input file: " + filename + "\n");
        read_input_on_root(world, filename, tag);
        read_commandline_on_root(world, parser, commandlineparser::tolower(tag));
    });
}

void QCCalculationParametersBase::read_input_on_root(
        World& world, const std::string& filename, const std::string& tag) {
    try {
        std::string contents, line;
        std::ifstream file(filename);
        while (std::getline(file, line)) contents += line + "\n";
        read_internal(world, contents, tag);
    } catch (const std::invalid_argument&) {
        throw;
    } catch (const std::exception& e) {
        if (is_missing_datagroup_exception(e)) {
            if (!throw_if_datagroup_not_found) return;
            throw std::runtime_error("could not read data group >>" + tag + "<< in file " + filename);
        }
        throw std::runtime_error("error while reading data group >>" + tag + "<< in file "
                                 + filename + ": " + e.what());
    }
}

void QCCalculationParametersBase::validate_commandline_occurrence(
        World& world, std::string contents, const std::string& tag) const {
    // Do not temporarily mutate the owning application: it may have state in
    // addition to the parameter map. A probe copies registrations and policies.
    struct Probe : QCCalculationParametersBase {
        const QCCalculationParametersBase& owner;
        explicit Probe(const QCCalculationParametersBase& source)
            : QCCalculationParametersBase(source), owner(source) {}
        [[nodiscard]] std::string get_tag() const override { return owner.get_tag(); }
        [[nodiscard]] std::set<std::string> case_sensitive_keys() const override { return owner.case_sensitive_keys(); }
        [[nodiscard]] std::set<std::string> full_precision_keys() const override { return owner.full_precision_keys(); }
        void validate_input_value(const std::string& key, const std::string& value) const override {
            owner.validate_input_value(key, value);
        }
    } probe(*this);
    probe.read_internal(world, contents, tag);
}

/// Apply only the last whole CLI group, after validating all original groups.
void QCCalculationParametersBase::read_commandline_options(World& world, const commandlineparser& parser,
                                                         const std::string tag) { // NOLINT(performance-unnecessary-value-param): preserve the public signature.
    const auto group = commandlineparser::tolower(tag);
    bool present = false;
    if (world.rank() == 0) {
        present = parser.key_exists(group);
        for (const auto& argument : parser.original_arguments())
            present = present || argument.key == group;
    }
    world.gop.broadcast(present, 0);
    if (!present) return;
    read_collectively(world, [&] { read_commandline_on_root(world, parser, group); });
}

void QCCalculationParametersBase::read_commandline_on_root(
        World& world, const commandlineparser& parser, const std::string& group) {
    const auto as_input = [&](std::string value) {
        std::replace(value.begin(), value.end(), '\'', ' ');
        std::replace(value.begin(), value.end(), ';', '\n');
        return group + "\n" + value + "\nend";
    };
    for (const auto& argument : parser.original_arguments()) {
        if (argument.key == group)
            validate_commandline_occurrence(world, as_input(argument.value), group);
    }
    if (parser.key_exists(group)) {
        auto contents = as_input(parser.value_raw(group));
        read_internal(world, contents, group);
    }
}


/// read the stream, starting from tag

/// only parameters that are defined in the constructor will be processed,
/// all others will be discarded.
void QCCalculationParametersBase::read_internal(World& world, std::string& filecontents, std::string tag) {
	std::stringstream f(filecontents);
	// the group starts at its header, the line whose first word is the tag. Taking the
	// tag as any word stopped at a value equal to it in an earlier group (`guess lcao`
	// in the dft group before the lcao group), and read that group's remaining lines
	position_stream_to_word(f, tag, '#', true, true, true);
	std::string line, key,value;

	// read input lines
	while (std::getline(f,line)) {

		// keep an original-case copy: values of case_sensitive_keys() (paths and
		// other free-form identifiers) must not be folded, or a deck cannot
		// round-trip them -- see case_sensitive_keys()
		std::string line_raw = line;

		// all in lower case
		std::transform(line.begin(), line.end(), line.begin(), ::tolower);

		// remove comments from line
		std::size_t last = line.find_first_of('#');
		line=line.substr(0,last);
        std::replace_copy(line.begin(), line.end(), line.begin(),'=', ' ');

		line_raw=line_raw.substr(0,line_raw.find_first_of('#'));
        std::replace_copy(line_raw.begin(), line_raw.end(), line_raw.begin(),'=', ' ');

		std::stringstream sline(line);

		// sline might be empty by now
		if (not (sline >> key)) continue;
		if (print_debug) ::madness::print("reading key ",key);

		// skip comment line
		if (key[0]=='#') continue;
		if (key=="end") break;

		// check if key exists in the initialized parameter list
		if (not (parameter_exists(key))) {
            if (not ignore_unknown_keys) {
                if (world.rank()==0) {
                    ::madness::print("found unknown key: ",key);
                    ::madness::print("in datagroup:      ",tag);
                }
                throw std::runtime_error("input error");
            }
            if ((not ignore_unknown_keys_silently)
                and (world.rank()==0)) madness::print("ignoring unknown parameter in input file: ",key);
			continue;
		}

		std::string word,line1;
		// for case-sensitive keys take the value from the original-case line
		std::stringstream svalues(keeps_case(key) ? line_raw : line);
		svalues >> word;                       // drop the key
		while (svalues >> word) {line1+=word+" ";}
		// trim result
		last = line1.find_last_not_of(' ');
		line1=line1.substr(0, last+1);


		// check which type to expect from the given key
		bool success=false;
		try {
			validate_input_value(key, line1);
			success = try_setting_any<all_parameter_types>(key, line1);
        } catch (std::invalid_argument& e) {
            throw;
		} catch (std::exception& e) {
			std::string errmsg="found an error for key >> "+key+" << \n" +e.what();
            throw std::runtime_error(errmsg);
		}
		if (not success) {
			madness::print("\n\ncould not assign the input parameter for key ",key);
			std::string requested_type=get_parameter(key).get_type();
			madness::print("\trequested type: ",requested_type,"\n");
			madness::print("add the corresponding type to QCCalculationsParametersBase.h\n\n");
			throw std::invalid_argument("add the corresponding type to QCCalculationsParametersBase.h");
		}
	}
};


bool operator==(const QCCalculationParametersBase& p1,
			   const QCCalculationParametersBase& p2) {
	if (p1.parameters.size() != p2.parameters.size()) return false;
	for (const auto& [key, value] : p1.parameters) {
		if (not p2.parameter_exists(key)) return false;
		if (not (value.get_value() == p2.get_parameter(key).get_value())) return false;
	}
	return true;
}


bool operator!=(const QCCalculationParametersBase& p1,
                const QCCalculationParametersBase& p2) {
    return !(p1 == p2);
}
} /* namespace madness */
