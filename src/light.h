#ifndef LIGHT_H
#define LIGHT_H

#include "vec3.h"
#include "ray.h"
#include "quad.h"
#include "hittable.h"
#include "hittable_list.h"
#include "rtweekend.h"


class light : public quad {
private:
    const vec3 le;
public:

    light(point3& Q, const vec3& u, const vec3& v, shared_ptr<material> mat, const vec3& le) : quad(Q, u, v, mat), le(le) {}

    vec3 Le() {
        return le;
    }

    // primarily used to sample points on light for photon map construction
    hit_record samplePoint(float& pdf) const {
        hit_record rec;
        rec.u = random_double(0, 1);
        rec.v = random_double(0, 1);
        rec.mat = mat;
        rec.p = Q + u * rec.u + v * rec.v;
        rec.normal = normal;

        // pdf = 1 / surface area
        pdf = 1.0 / cross(u, v).length();

        return rec;
    }

    // primarily used to sample points on light for photon map construction
    color sampleDirection(const hit_record& rec, float& pdf) const {

        const vec3 dir = sampleCosineHemisphere(random_double(0, 1), random_double(0, 1), pdf);

        // orthonormal basis
        vec3 b1, b2;
        orthonormalBasis(normal, b1, b2);

        return localToWorld(dir, b1, normal, b2);
    }
};

#endif