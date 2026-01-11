#ifndef KDTREE_H
#define KDTREE_H

#include "vec3.h"
#include "hittable.h"
#include "ray.h"
#include "quad.h"
#include <vector>
#include "rtweekend.h"
#include "light.h"
#include "material.h"
#include <cstdlib>

class Photon {
public:
	vec3 pos;			// position
	color throughput;	// power packaed as 4 chars
	vec3 incident;		// compresed incident direction
	short flag;			// flag used in kdtree sorting + look-up
};

class KDTree {
private:
	std::vector<Photon> photons;
	struct KDNode {
		char axis; // x=0, y=1, z=2
		int idx; // index of median point
		int leftChildIdx;
		int rightChildIdx;

		KDNode() : axis(-1), idx(-1), leftChildIdx(-1), rightChildIdx(-1) {}
	};
public:
	void buildTree(const std::vector<Photon> photons, int nPoints) {
		this->photons = photons;
		std::vector<int> indices(nPoints);
		// replacing std::iota
		int value = 0;
		for (auto first = indices.begin(); first != indices.end(); ++first, ++value)
			*first = value;
		// build tree recursively

	}

	std::vector<int> queryKNearestPhotons(const vec3& p, int k,
		float& max_dist2) const {
		std::vector<int> a;
		return a;
	}

	const Photon& getIthPhoton(int i) const {
		return photons[i];
	}
};

class PhotonMap {
private:
	std::vector<Photon> photons;
	KDTree globalMap;
	KDTree causticsMap;
	shared_ptr<light> mainLight;

	ray sampleRayFromLight(const shared_ptr<light> light, vec3 throughput)
	{
		// sample point on light
		float light_pos_pdf;
		const hit_record light_surf = light->samplePoint(light_pos_pdf);
	
		// sample direction on light
		float light_dir_pdf;
		const vec3 dir = light->sampleDirection(light_surf, light_dir_pdf);
	
		// spawn ray
		ray r = ray(light_surf.p, dir);
		throughput = light->Le() / (light_pos_pdf * light_dir_pdf) * 
			std::abs(dot(dir, light_surf.normal));
		return r;
	}
public:
	int nPhotonsGlobal;
	int nPhotonsCaustic;
	int maxDepth;
	int finalGatheringDepth;
	int nEstimationGlobal;

	PhotonMap() {}

	static float cosTerm(const vec3& wo, const vec3& wi, const hit_record& surfaceInfo) {
		// function from yumcyaWiz repo has geometric normal and shading normal, don't know if that's needed
		const float wi_ns = dot(wi, surfaceInfo.normal);
		const float wi_ng = dot(wi, surfaceInfo.normal);
		const float wo_ns = dot(wo, surfaceInfo.normal);
		const float wo_ng = dot(wo, surfaceInfo.normal);
			// prevent light leaks
		if (wi_ng * wi_ns <= 0 || wo_ng * wo_ns <= 0) {
			return 0;
		}
		return std::abs(wo_ns) * std::abs(wi_ng) / std::abs(wo_ng);
	}

	void build(const hittable& world, const shared_ptr<light> light) {
		mainLight = light;
		for (int i = 0; i < nPhotonsGlobal; i++) {
			vec3 throughput;
			ray r = sampleRayFromLight(light, throughput);

			for (int k = 0; k < maxDepth; k++) {

				hit_record hit;
				if (world.hit(r, interval(0.001, infinity), hit, nullptr)) {
					// if is diffuse
					if (dynamic_cast<lambertian*>(hit.mat.get()) != nullptr) {
						photons.emplace_back(hit.p, throughput, -r.direction());
					}
				}

				if (k > 0) {
					const float russian_roulette_prob = std::min(
						std::max(throughput[0], std::max(throughput[1], throughput[2])),
							1.0);
					if (random_double(0, 1) >= russian_roulette_prob) { break; }
					throughput /= russian_roulette_prob;
				}

				// sample direction by BxDF
				ray scattered;
				color attenuation;
				hit.mat->scatter(r, hit, attenuation, scattered);

				throughput *= attenuation * cosTerm(-r.direction(), scattered.direction(), hit);

				r = ray(hit.p, scattered.direction());

			}
		}

		globalMap.buildTree(photons, photons.size());

		photons.clear();

		for (int i = 0; i < nPhotonsCaustic; i++) {
			vec3 throughput;
			ray r = sampleRayFromLight(light, throughput);

			bool prev_specular = false;
			for (int k = 0; k < maxDepth; k++) {
				hit_record hit;
				if (world.hit(r, interval(0.001, infinity), hit, nullptr)) {
					// if is diffuse and previous was not specular
					if (!prev_specular && dynamic_cast<lambertian*>(hit.mat.get()) != nullptr) {
						break;
					}

					// add photon when hitting diffuse surface after a specular
					if (prev_specular && dynamic_cast<lambertian*>(hit.mat.get()) != nullptr) {
						photons.emplace_back(hit.p, throughput, -r.direction());
						break;
					}

					prev_specular = dynamic_cast<metal*>(hit.mat.get()) != nullptr;

					if (k > 0) {
						const float russian_roulette_prob = std::min(
							std::max(throughput[0], std::max(throughput[1], throughput[2])),
							1.0);
						if (random_double(0, 1) >= russian_roulette_prob) { break; }
						throughput /= russian_roulette_prob;
					}

					// sample direction by BxDF
					ray scattered;
					color attenuation;
					hit.mat->scatter(r, hit, attenuation, scattered);

					throughput *= attenuation * cosTerm(-r.direction(), scattered.direction(), hit);

					r = ray(hit.p, scattered.direction());
				}
			}
		}

		causticsMap.buildTree(photons, photons.size());
	}

	color integrate(ray r, const hittable& world, int depth) const {
		if (depth >= maxDepth) return color(0, 0, 0);

		hit_record hit;
		shared_ptr<hittable> obj;
		// If the ray hits nothing, return the background color.
		if (!world.hit(r, interval(0.001, infinity), hit, obj))
			return color(0,0,0);

		if (dynamic_cast<diffuse_light*>(hit.mat.get()) != nullptr) {
			hittable* l = obj.get();
			return ((light*)l)->Le();
		}

		if (dynamic_cast<lambertian*>(hit.mat.get()) != nullptr) {
			if (depth >= finalGatheringDepth) {
				return computeRadianceWithPhotonMap(ray(r.origin(), - r.direction()), hit);
			}
			else {
				const vec3 Ld = computeDirectIllumination(world, ray(r.origin(), -r.direction()), hit);
				const vec3 Lc = computeCausticsWithPhotonMap(-r.direction(), hit);
				const vec3 Li = computeIndirectIllumination(world, -r.direction(), hit);
				return (Ld + Lc + Li);
			}
		}
		else if (dynamic_cast<metal*>(hit.mat.get()) != nullptr) {
			if (depth >= 3) {
				// sample direction by BxDF
				ray scattered;
				color attenuation;
				hit.mat->scatter(r, hit, attenuation, scattered);

				const vec3 throughput = attenuation * cosTerm(-r.direction(), scattered.direction(), hit);

				ray next_r = ray(hit.p, scattered.direction());
				return throughput * integrate(next_r, world, depth + 1);
			}
			else {
				// sample all directions
				const std::vector<DirectionPair> dir_pairs = hit.mat->scatter
			}
		}
		return color(0, 0, 0);
	}

	color computeRadianceWithPhotonMap(ray r, hit_record hit) const {
		float max_dist2;
		const std::vector<int> photon_indices = globalMap.queryKNearestPhotons(hit.p, nEstimationGlobal, max_dist2);

		vec3 Lo;
		for (const int photon_idx : photon_indices) {
			const Photon& photon = globalMap.getIthPhoton(photon_idx);
			ray scattered;
			color attenuation;
			hit.mat->scatter(r, hit, attenuation, scattered);

			Lo += attenuation * photon.throughput;
		}

		if (photon_indices.size() > 0) {
			Lo /= (nPhotonsGlobal * pi * max_dist2);
		}
		return Lo;
	}
	color computeDirectIllumination(const hittable& world, ray r_in, hit_record hit) const {
		vec3 Ld;

		// sample point on light
		float light_pos_pdf;
		const hit_record light_surf = mainLight->samplePoint(light_pos_pdf);

		const vec3 wi = unit_vector(light_surf.p - hit.p);
		const float r = (light_surf.p - hit.p).length();
		const float pdf_dir = light_pos_pdf * r * r / std::abs(dot(-wi, light_surf.normal));

		ray ray_shadow(hit.p, wi);
		ray_shadow.tmax = r - RAY_EPS;

		hit_record info_shadow;
		if (world.hit(ray_shadow, interval(0.001, infinity), info_shadow, nullptr)) {
			const vec3 Le = mainLight.get()->Le();
			// sample direction by BxDF
			ray scattered;
			color attenuation;
			hit.mat->scatter(r, hit, attenuation, scattered);
		}
	}
	color computeCausticsWithPhotonMap(vec3 dir, hit_record hit) const {

	}
	color computeIndirectIllumination(const hittable& world, vec3 dir, hit_record hit) const {

	}
};


#endif