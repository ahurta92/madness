// LEGACY_PATCH: Force-included preamble for molresponse_legacy.
// The bsundahl codebase assumes 'using namespace madness' everywhere.
// Using -include ensures this fires before any source code, and including
// mra.h first guarantees namespace madness is defined before the directive.
#pragma once
#include <madness/mra/mra.h>
using namespace madness;
