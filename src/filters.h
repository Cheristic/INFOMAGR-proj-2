#pragma once

#include <cmath>
#include "vec3.h"

inline double gaussianFilter(const point3& pos, const point3& photonPos, double dist, double alpha = 0.918, double beta = 1.953)
{
	double photonDist2 = (pos - photonPos).length_squared();
	double exponent = -beta * photonDist2 / (2 * dist);
	return alpha * (1 - ((1 - std::exp(exponent)) / (1 - std::exp(-beta))));
}