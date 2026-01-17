#ifndef MATERIAL_H
#define MATERIAL_H
//==============================================================================================
// Originally written in 2016 by Peter Shirley <ptrshrl@gmail.com>
//
// To the extent possible under law, the author(s) have dedicated all copyright and related and
// neighboring rights to this software to the public domain worldwide. This software is
// distributed without any warranty.
//
// You should have received a copy (see file COPYING.txt) of the CC0 Public Domain Dedication
// along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//==============================================================================================

#include "hittable.h"
#include "texture.h"
#include "vec3.h"

class material {
  public:
    virtual ~material() = default;

    virtual color emitted(double u, double v, const point3& p) const {
        return color(0,0,0);
    }
    virtual color evaluate(const vec3& r_in, const vec3& r_out) const = 0 {}
    virtual color sampleDirection(
        const vec3& r_in, const hit_record& rec, vec3& r_out, float& pdf
    ) const = 0 { }
    // schlick approximation of fresnel reflectance
    static float fresnel(float cosThetaI, float iorI, float iorT) {
        const float f0 =
            (iorI - iorT) * (iorI - iorT) / ((iorI + iorT) * (iorI + iorT));
        const auto pow5 = [](float x) { return x * x * x * x * x; };
        return f0 + (1.0f - f0) * pow5(std::max(1.0f - std::abs(cosThetaI), 0.0f));
    }
    static float cosTheta(const vec3& v) { return v[1]; }
    static float absCosTheta(const vec3& v) { return std::abs(cosTheta(v)); }

    virtual std::vector<DirectionPair> sampleAllDirections(const vec3& w0) const = 0 {}
};



class lambertian : public material {
  public:
    lambertian(const color& albedo) : tex(make_shared<solid_color>(albedo)) {}
    lambertian(shared_ptr<texture> tex) : tex(tex) {}

    color evaluate(const vec3& r_in, const vec3& r_out) const override {
        // when wo, wi is under the surface, return 0
        const float cosThetaO = cosTheta(r_in);
        const float cosThetaI = cosTheta(r_out);
        if (cosThetaO < 0 || cosThetaI < 0) return vec3(0, 0, 0);

        return tex->value(0, 0, vec3(0, 0, 0)) / pi; // return color / pi
    }

    color sampleDirection(const vec3& r_in, const hit_record& rec, vec3& r_out, float& pdf)
    const override {

        // orthonormal basis
        vec3 b1, b2;
        orthonormalBasis(rec.normal, b1, b2);

        const vec3 r_in_local = worldToLocal(r_in, b1, rec.normal, b2);

        vec3 r_out_local = sampleCosineHemisphere(random_double(0, 1), random_double(0, 1), pdf);

        r_out = localToWorld(r_out_local, b1, rec.normal, b2);

        return evaluate(r_in_local, r_out_local);
    }

    std::vector<DirectionPair> sampleAllDirections(const vec3& w0) const override
    {
        std::vector<DirectionPair> ret;
        return ret;
    }

  private:
    shared_ptr<texture> tex;
};


class metal : public material {
  public:
    metal(const color& albedo, double fuzz) : albedo(albedo), fuzz(fuzz < 1 ? fuzz : 1) {}

    color evaluate(const vec3& r_in, const vec3& r_out) const override {
        return vec3(0, 0, 0);
    }

    color sampleDirection(const vec3& r_in, const hit_record& rec, vec3& r_out, float& pdf)
    const override {
        // orthonormal basis
        vec3 b1, b2;
        orthonormalBasis(rec.normal, b1, b2);

        const vec3 r_in_local = worldToLocal(r_in, b1, rec.normal, b2);

        vec3 r_out_local = reflect(r_in, vec3(0, 1, 0));
        pdf = 1.0f;

        r_out = localToWorld(r_out_local, b1, rec.normal, b2);

        return albedo / absCosTheta(r_in_local);
    }

    std::vector<DirectionPair> sampleAllDirections(const vec3& w0) const override
    {
        std::vector<DirectionPair> ret;
        const vec3 wi = reflect(w0, vec3(0, 1, 0));
        ret.emplace_back(wi, albedo / absCosTheta(wi));
        return ret;
    }

  private:
    color albedo;
    double fuzz;
};


class diffuse_light : public material {
  public:
    diffuse_light(shared_ptr<texture> tex) : tex(tex) {}
    diffuse_light(const color& emit) : tex(make_shared<solid_color>(emit)) {}

    color emitted(double u, double v, const point3& p) const override {
        return tex->value(u, v, p);
    }

    vec3 evaluate(const vec3& r_in, const vec3& r_out) const override { return vec3(0, 0, 0); }

    color sampleDirection(const vec3& r_in, const hit_record& rec, vec3& r_out, float& pdf) const override {
        return color(0,0,0);
    }

    std::vector<DirectionPair> sampleAllDirections(const vec3& w0) const override { return std::vector<DirectionPair>(); }

  private:
    shared_ptr<texture> tex;
};



#endif
