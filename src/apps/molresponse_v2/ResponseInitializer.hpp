
#pragma once
#include "GroundStateData.hpp"
#include "ResponseState.hpp"
#include "ResponseVector.hpp"

ResponseVector initialize_guess_vector(World &world, const GroundStateData &gs,
                                       const LinearResponsePoint &pt);
