#pragma once

#include "photon_antonio.h"
#include "vec3.h"
#include <algorithm>
#include <vector>
#include <iostream>
#include <queue>

class KDTree
{
private:
	struct Node
	{
		unsigned int idx;
		int leftChildIdx, rightChildIdx;
		int axis;

		Node() : axis(-1), idx(-1), leftChildIdx(-1), rightChildIdx(-1) {}
	};

	int numPhotons;

	void BuildNode(int* idx, int nPoints, int depth)
	{
		if (nPoints <= 0) return;

		const int axis = depth % 3;

		std::sort(idx, idx + nPoints, [&](const int idx1, const int idx2) {return photons[idx1].pos[axis] < photons[idx2].pos[axis]; });

		const int mid = (nPoints - 1) / 2;

		const int parentIdx = nodes.size();
		Node node;
		node.axis = axis;
		node.idx = idx[mid];
		nodes.push_back(node);

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

		const Node& node = nodes[nodeIdx];

		const photon& median = photons[node.idx];

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
	std::vector<Node> nodes;
	const photon* photons;

	void SetPhotons(const photon* photonArray, int nPhotons)
	{
		this->photons = photonArray;
		this->numPhotons = nPhotons;
	}

	void MakeTree()
	{
		std::vector<int> idx(numPhotons);
		for (int i = 0; i < numPhotons; i++)
		{
			idx[i] = i;
		}

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
};


