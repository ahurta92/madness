// ResponseSolver.cpp
//
// All per-type solver logic has been moved to ops/ sub-headers (see
// ResponseSolver.hpp).  This file is kept for the CMake target but contains
// no definitions — everything is inline in the ops headers.
//
// Per-type ops locations:
//   ops/StaticRestrictedOps.hpp   — static (ω=0) restricted
//   ops/TDARestrictedOps.hpp      — TDA excited-state restricted
//   ops/DynamicRestrictedOps.hpp  — full TDDFT restricted
//   ops/StaticUnrestrictedOps.hpp — static unrestricted (stubs)
//   ops/DynamicUnrestrictedOps.hpp — dynamic unrestricted (stubs)

#include "ResponseSolver.hpp"
