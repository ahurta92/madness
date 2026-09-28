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
*/

#ifndef MADNESS_CHEM_CONVERGENCELOG_H_INCLUDED
#define MADNESS_CHEM_CONVERGENCELOG_H_INCLUDED

/// \file ConvergenceLog.h
/// \brief Per-iteration convergence log of a solver: one CSV row per iteration.

#include <madness/world/MADworld.h>

#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace madness {

/// Value of a convergence target that the solve does not test.
inline constexpr double convergence_not_tested = std::numeric_limits<double>::quiet_NaN();

/// Append one iteration to `convergence/<name>.csv` in the current directory
/// (the task directory under madqc). `columns` are (name, value) pairs in a
/// fixed order; the header is written only into a new or empty file, so a
/// restarted solve appends to its earlier history. A target the solve does not
/// test is logged as nan (convergence_not_tested). Rank 0 writes; not collective.
inline void append_convergence_row(World& world, const std::string& name,
                                   const std::vector<std::pair<std::string, double>>& columns) {
    if (world.rank() != 0) return;
    std::error_code ec;
    std::filesystem::create_directories("convergence", ec);
    const std::string path = "convergence/" + name + ".csv";
    const bool fresh = !std::filesystem::exists(path, ec) || std::filesystem::file_size(path, ec) == 0;
    std::ofstream out(path, std::ios::app);
    if (!out) return;
    if (fresh) {
        for (std::size_t i = 0; i < columns.size(); ++i) out << (i ? "," : "") << columns[i].first;
        out << '\n';
    }
    out.precision(12);
    for (std::size_t i = 0; i < columns.size(); ++i) out << (i ? "," : "") << columns[i].second;
    out << '\n';
}

} // namespace madness

#endif // MADNESS_CHEM_CONVERGENCELOG_H_INCLUDED
