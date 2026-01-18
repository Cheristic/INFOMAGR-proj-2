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
#include "filters.h"
#include <cstdlib>
#include <cmath>
#include <queue>
#include <fstream>

class Photon {
public:
	Photon(point3 pos, color throughput, vec3 incident) : pos(pos), throughput(throughput), incident(incident) {}

	point3 pos;			// position
	color throughput;	// power packaed as 4 chars
	vec3 incident;		// compresed incident direction
	short flag;			// flag used in kdtree sorting + look-up
};

class KDTree {
private:
	int numPhotons;

	struct KDNode {
		char axis; // x=0, y=1, z=2
		int idx; // index of median point
		int leftChildIdx;
		int rightChildIdx;

		KDNode() : axis(-1), idx(-1), leftChildIdx(-1), rightChildIdx(-1) {}
	};

	void BuildNode(int* idx, int nPoints, int depth)
	{
		if (nPoints <= 0) return;

		const int axis = depth % 3;

		std::sort(idx, idx + nPoints, [&](const int idx1, const int idx2) {return photons[idx1].pos[axis] < photons[idx2].pos[axis]; });

		const int mid = (nPoints - 1) / 2;

		const int parentIdx = nodes.size();
		KDNode node;
		node.axis = axis;
		node.idx = idx[mid];
		nodes.push_back(node);

		photons[idx[mid]].flag = axis;

		const int leftChildIdx = nodes.size();
		BuildNode(idx, mid, depth + 1);

		if (leftChildIdx == nodes.size())
		{
			nodes[parentIdx].leftChildIdx = -1;
		}
		else
		{
			nodes[parentIdx].leftChildIdx = leftChildIdx;
		}

		const int rightChildIdx = nodes.size();
		BuildNode(idx + mid + 1, nPoints - mid - 1, depth + 1);

		if (rightChildIdx == nodes.size())
		{
			nodes[parentIdx].leftChildIdx = -1;
		}
		else
		{
			nodes[parentIdx].rightChildIdx = rightChildIdx;
		}
	}

	using KNNQueue = std::priority_queue<std::pair<float, int>>;
	void SearchKNearestNode(int nodeIdx, const point3& queryPoint, int k, KNNQueue& queue) const
	{
		if (nodeIdx == -1 || nodeIdx >= nodes.size()) return;

		const KDNode& node = nodes[nodeIdx];

		const Photon& median = photons[node.idx];

		const float dist2 = (queryPoint - median.pos).length_squared();
		queue.emplace(dist2, node.idx);

		if (queue.size() > k)
		{
			queue.pop();
		}

		const bool isLower = queryPoint[node.axis] < median.pos[node.axis];
		if (isLower)
		{
			SearchKNearestNode(node.leftChildIdx, queryPoint, k, queue);
		}
		else
		{
			SearchKNearestNode(node.rightChildIdx, queryPoint, k, queue);
		}

		const float distToSiblings = median.pos[node.axis] - queryPoint[node.axis];
		if (queue.top().first > distToSiblings * distToSiblings)
		{
			if (isLower)
			{
				SearchKNearestNode(node.rightChildIdx, queryPoint, k, queue);
			}
			else
			{
				SearchKNearestNode(node.leftChildIdx, queryPoint, k, queue);
			}
		}
	}

public:
	KDTree() {}
	std::vector<KDNode> nodes;
	std::vector<Photon> photons;

	void SetPhotons(std::vector<Photon> photonArray, int nPhotons)
	{
		this->photons = photonArray;
		this->numPhotons = nPhotons;
	}

	void MakeTree()
	{
		std::vector<int> idx(numPhotons);
		// replacing std::iota
		int value = 0;
		for (auto first = idx.begin(); first != idx.end(); ++first, ++value)
			*first = value;

		BuildNode(idx.data(), numPhotons, 0);
	}

	std::vector<int> SearchKNearest(const point3 queryPoint, int k, float& maxDist2) const
	{
		KNNQueue queue;
		SearchKNearestNode(0, queryPoint, k, queue);

		std::vector<int> ret(queue.size());
		maxDist2 = 0;
		for (int i = 0; i < ret.size(); ++i)
		{
			const auto& p = queue.top();
			ret[i] = p.second;
			maxDist2 = std::max(maxDist2, p.first);
			queue.pop();
		}

		return ret;
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

	ray sampleRayFromLight(const shared_ptr<light> light, vec3& throughput)
	{
		// sample point on light
		float light_pos_pdf; // <-- this will always be the same value, same surface area
		const hit_record light_surf = light->samplePoint(light_pos_pdf);
	
		// sample direction on light
		float light_dir_pdf;
		const vec3 dir = light->sampleDirection(light_surf, light_dir_pdf);
		// spawn ray
		ray r = ray(light_surf.p, dir);

		auto v = (light->Le() / (light_pos_pdf * light_dir_pdf));
		auto v2 = std::abs(dot(dir, light_surf.normal));
		throughput = v * v2;

		return r;
	}
public:
	int nPhotonsGlobal;
	int nPhotonsCaustic;
	int maxDepth;
	int finalGatheringDepth;
	int nEstimationPhotons;

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
		int count = 0;
		for (int i = 0; i < nPhotonsGlobal; i++) {
			std::clog << "\rGlobal photons remaining: " << (nPhotonsGlobal - i) << ' ' << std::flush;
			vec3 throughput;
			//std::clog << "globalB_: " << (throughput) << "\n";
			ray r = sampleRayFromLight(light, throughput);
			//std::clog << "global: " << (throughput) << "\n";

			for (int k = 0; k < maxDepth; k++) {

				hit_record hit;
				//std::clog << "global: " << (throughput) << "\n";
				if (world.hit(r, interval(0.001, infinity), hit, nullptr)) {
					// if is diffuse
					if (dynamic_cast<lambertian*>(hit.mat.get()) != nullptr) {
						photons.emplace_back(hit.p, throughput, -r.direction());
					}
					if (k > 0) {
						const float russian_roulette_prob = std::min(
							std::max(throughput[0], std::max(throughput[1], throughput[2])),
							1.0);
						if (random_double(0, 1) >= russian_roulette_prob) { break; }
						throughput /= russian_roulette_prob;
					}

					// sample direction by BxDF
					vec3 r_out;
					float pdf_dir;
					vec3 col = hit.mat->sampleDirection(-r.direction(), hit, r_out, pdf_dir);

					throughput *= col * cosTerm(-r.direction(), r_out, hit) / pdf_dir;

					r = ray(hit.p, r_out);
					count++;
				}			
				else {
					break; // photon goes to sky
				}
			}
			std::clog << "\n";
		}
		globalMap.SetPhotons(photons, photons.size());
		globalMap.MakeTree();

		photons.clear();

		for (int i = 0; i < nPhotonsCaustic; i++) {
			vec3 throughput;
			ray r = sampleRayFromLight(light, throughput);

			bool prev_specular = false;
			for (int k = 0; k < maxDepth; k++) {
				hit_record hit;
				if (world.hit(r, interval(0.001, infinity), hit, nullptr)) {
					// if is diffuse and previous was not specular
					if (dynamic_cast<lambertian*>(hit.mat.get()) != nullptr) {
						if (prev_specular) photons.emplace_back(hit.p, throughput, -r.direction());
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
					vec3 r_out;
					float pdf_dir;
					vec3 col = hit.mat->sampleDirection(-r.direction(), hit, r_out, pdf_dir);

					throughput *= col * cosTerm(-r.direction(),r_out, hit) / pdf_dir;

					r = ray(hit.p, r_out);
				}
			}

			//std::clog << "caustics: " << (throughput) << "\n";
		}

		causticsMap.SetPhotons(photons, photons.size());
		causticsMap.MakeTree();
	}

	color integrate(ray r, const hittable& world, int depth) const {
		if (depth >= maxDepth) return color(0, 0, 0);

		hit_record hit;
		shared_ptr<hittable> obj;
		// If the ray hits nothing, return the background color.
		if (!world.hit(r, interval(0.001, infinity), hit, &obj))
			return color(0, 0, 0);

		if (dynamic_cast<diffuse_light*>(hit.mat.get()) != nullptr) {
			hittable* l = obj.get();
			return ((light*)l)->Le();
		}

		if (dynamic_cast<lambertian*>(hit.mat.get()) != nullptr) {
			if (depth >= finalGatheringDepth) {
				//std::clog << "1";
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
				vec3 r_out;
				float pdf_dir;
				vec3 col = hit.mat->sampleDirection(r.direction(), hit, r_out, pdf_dir);

				const vec3 throughput = col * cosTerm(r.direction(), r_out, hit) / pdf_dir;

				ray next_r = ray(hit.p, r_out);
				return throughput * integrate(next_r, world, depth + 1);
			}
			else {
				// sample all directions
				const std::vector<DirectionPair> dir_pairs = sampleAllBxDF(r.direction(), hit);

				// Recursively raytrace
				vec3 L0;
				for (const auto& dp : dir_pairs)
				{
					const vec3 dir = dp.first;
					const vec3 f = dp.second;

					const ray next_ray(hit.p, dir);
					const vec3 throughput = f * std::abs(dot(dir, hit.normal));

					L0 += throughput * integrate(next_ray, world, depth + 1);
				}
				return L0;
			}
		}
		else {
			std::clog << "\r[PhotonMap] invalid material type ";
			return color(0, 0, 0);
		}

		return color(0, 0, 0);
	}

	std::vector<DirectionPair> sampleAllBxDF(const vec3& w0, const hit_record& rec) const
	{
		DirectionPair tangentVecs = getTangentVectors(rec.normal);
		// Convert from world to local coords
		const vec3& w0Local = worldToLocal(w0, tangentVecs.first, rec.normal, tangentVecs.second);

		// Sample all directions in the tangent space
		std::vector<DirectionPair> dirPairs = rec.mat->sampleAllDirections(w0Local);

		// Convert from local to world coords
		for (auto& pair : dirPairs)
		{
			pair.first = localToWorld(pair.first, tangentVecs.first, rec.normal, tangentVecs.second);
		}

		return dirPairs;
	}

	color computeRadianceWithPhotonMap(ray r, hit_record hit) const {
		float max_dist2;
		const std::vector<int> photon_indices = globalMap.SearchKNearest(hit.p, nEstimationPhotons, max_dist2);

		vec3 Lo;
		for (const int photon_idx : photon_indices) {
			const Photon& photon = globalMap.getIthPhoton(photon_idx);

			// orthonormal basis
			vec3 b1, b2;
			orthonormalBasis(hit.normal, b1, b2);

			const vec3 r_in_local = worldToLocal(r.direction(), b1, hit.normal, b2);
			const vec3 r_out_local = worldToLocal(photon.incident, b1, hit.normal, b2);
			const color col = hit.mat->evaluate(r_in_local, r_out_local);

			Lo += col * photon.throughput;// *gaussianFilter(hit.p, photon.pos, max_dist2);
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

		// Convert positional pdf to directional pdf
		const vec3 wi = unit_vector(light_surf.p - hit.p);
		const float r = (light_surf.p - hit.p).length();
		const float pdf_dir = light_pos_pdf * r * r / std::abs(dot(-wi, light_surf.normal));

		// Make shadow ray
		ray ray_shadow(hit.p, wi);
		ray_shadow.tmax = r - RAY_EPS;

		// Trace ray to main light
		hit_record info_shadow;
		if (!world.hit(ray_shadow, interval(0.001, infinity), info_shadow, nullptr)) {
			const vec3 Le = mainLight.get()->Le();

			DirectionPair tangentVecs = getTangentVectors(hit.normal);
			const vec3 r_in_local = worldToLocal(r_in.direction(), tangentVecs.first, hit.normal, tangentVecs.second);
			const vec3 r_out_local = worldToLocal(wi, tangentVecs.first, hit.normal, tangentVecs.second);
			const vec3 col = hit.mat->evaluate(r_in_local, r_out_local);
			const float cos = std::abs(dot(wi, hit.normal));
			Ld = col * cos * Le / pdf_dir;
		}

		return Ld;
	}
	color computeCausticsWithPhotonMap(vec3 dir, hit_record hit) const 
	{
		// Get nearby photons
		float max_dist2;
		const std::vector<int> photonIdx = causticsMap.SearchKNearest(hit.p, nEstimationPhotons, max_dist2);

		vec3 Lo;
		for (const int idx : photonIdx)
		{
			const Photon& photon = causticsMap.getIthPhoton(idx);
			vec3 localDir, localwi;
			DirectionPair tangentVecs = getTangentVectors(hit.normal);
			localDir = worldToLocal(dir, tangentVecs.first, hit.normal, tangentVecs.second);
			localwi = worldToLocal(photon.incident, tangentVecs.first, hit.normal, tangentVecs.second);
			const vec3 f = hit.mat->evaluate(localDir, localwi);
			Lo += f * photon.throughput * gaussianFilter(hit.p, photon.pos, max_dist2);
		}
		if (photonIdx.size() > 0)
		{
			Lo /= (nPhotonsCaustic * pi * max_dist2);
		}

		return Lo;
	}

	color computeIndirectIlluminationRecursive(const hittable& world, const vec3& dir, hit_record hit, int depth) const
	{
		if (depth >= maxDepth) return vec3(0, 0, 0);
		
		vec3 Li;

		// sample direction by BxDF
		vec3 dir_out;
		float pdf_dir;
		color col = hit.mat->sampleDirection(dir, hit, dir_out, pdf_dir);
		const float cos = std::abs(dot(hit.normal, dir_out));

		// Trace final gathering ray
		ray rayFinalGathering(hit.p, dir_out);
		hit_record hitFinalGathering;
		if (world.hit(rayFinalGathering, interval(0.001, infinity), hitFinalGathering, nullptr))
		{
			// if hit lambertian, compute radiance with photon map
			if (dynamic_cast<lambertian*>(hitFinalGathering.mat.get()) != nullptr)
			{
				Li += col * cos * computeRadianceWithPhotonMap(
					ray(rayFinalGathering.origin(), -rayFinalGathering.direction()), 
					hitFinalGathering) / pdf_dir;
			}
			// If hit specular (metal), recursively call function
			else if (dynamic_cast<metal*>(hitFinalGathering.mat.get()) != nullptr)
			{
				Li += col * cos * computeIndirectIlluminationRecursive(
					world, -rayFinalGathering.direction(), hitFinalGathering, depth + 1);
			}
		}

		return Li;
	}

	color computeIndirectIllumination(const hittable& world, vec3 dir, hit_record hit) const {
		return computeIndirectIlluminationRecursive(world, dir, hit, 0);
	}

	const KDTree* getPhotonMapPtr() const { return &globalMap; }
};


#endif