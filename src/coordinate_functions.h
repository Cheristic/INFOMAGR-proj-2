#pragma once

#include "vec3.h"
#include <cstdlib>

inline vec3& worldToLocal(const vec3& v, const vec3& localX, const vec3& localY, const vec3& localZ)
{
	return vec3(dot(v, localX), dot(v, localY), dot(v, localZ));
}

inline vec3& localToWorld(const vec3& v, const vec3& localX, const vec3& localY, const vec3& localZ)
{
	vec3& vWorld = vec3(0, 0, 0);
	for (int i = 0; i < 3; ++i)
	{
		vWorld[i] = v[0] * localX[i] + v[1] * localY[i] + v[2] * localZ[i];
	}

	return vWorld;
}

inline vec3& sphericalToCartesian(float theta, float phi)
{
	return vec3(std::cos(phi) * std::sin(theta), std::cos(theta), std::sin(phi) * std::cos(theta));
}

inline void orthonormalBasis(const vec3& n, vec3& t, vec3& b) {
	if (std::abs(n[1]) < 0.9f) {
		t = unit_vector(cross(n, vec3(0, 1, 0)));
	}
	else {
		t = unit_vector(cross(n, vec3(0, 0, -1)));
	}
	b = unit_vector(cross(t, n));
}

inline DirectionPair getTangentVectors(const vec3& normal)
{
	vec3 du, dv;
	orthonormalBasis(normal, du, dv);
	return DirectionPair(du, dv);
}

inline vec3 sampleCosineHemisphere(const float u, const float v, float& pdf) {
	const float max = (1.0f - 2.0f * u > -1) ? 1.0f - 2.0f * u : -1;
	const float clamped = (max < 1) ? max : 1;
	const float theta = 0.5f * std::acos(clamped);
	const float phi = pi * 2 * v;
	const float cosTheta = std::cos(theta);
	pdf = (1 / pi) * cosTheta;
	return sphericalToCartesian(theta, phi);
}