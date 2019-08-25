#include "bvh.h"

#include "CMU462/CMU462.h"
#include "static_scene/triangle.h"

#include <iostream>
#include <stack>

using namespace std;

namespace CMU462 {
namespace StaticScene {


BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {
  this->primitives = _primitives;

  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.

  BBox bb;
  for (size_t i = 0; i < primitives.size(); ++i) {
    bb.expand(primitives[i]->get_bbox());
  }

  root = new BVHNode(bb, 0, primitives.size());
  root->Partition(primitives, max_leaf_size);
}


BVHAccel::~BVHAccel() {
  // Implement a proper destructor for your BVH accelerator aggregate
  delete root;
}

BBox BVHAccel::get_bbox() const { return root->bb; }

bool BVHAccel::intersect(const Ray &ray) const {
  // TODO (PathTracer):
  // Implement ray - bvh aggregate intersection test. A ray intersects
  // with a BVH aggregate if and only if it intersects a primitive in
  // the BVH that is not an aggregate.

  bool hit = false;
  for (size_t p = 0; p < primitives.size(); ++p) {
    if (primitives[p]->intersect(ray)) hit = true;
  }

  return hit;
}

bool BVHAccel::intersect(const Ray &ray, Intersection *isect) const {
  // TODO (PathTracer):
  // Implement ray - bvh aggregate intersection test. A ray intersects
  // with a BVH aggregate if and only if it intersects a primitive in
  // the BVH that is not an aggregate. When an intersection does happen.
  // You should store the non-aggregate primitive in the intersection data
  // and not the BVH aggregate itself.

  bool hit = false;
  for (size_t p = 0; p < primitives.size(); ++p) {
    if (primitives[p]->intersect(ray, isect)) hit = true;
  }

  return hit;
}

BVHNode::~BVHNode()
{
	delete l;
	delete r;
}

void BVHNode::Partition(std::vector<Primitive*>& _primitives, size_t max_leaf_size)
{
	if (range <= max_leaf_size)
	{
		return;
	}

	const int BUCKET_SIZE = 8;

	vector<Primitive *> partitionA, partitionB;
	double sah = INF_D;

	for (size_t d = 0; d < 3; ++d)
	{
		vector<vector<Primitive *>> buckets(BUCKET_SIZE, vector<Primitive *>());
		vector<BBox> bucketBBoxes(BUCKET_SIZE);

		auto left = bb.min[d], right = bb.max[d];
		auto space = (right - left) / BUCKET_SIZE;

		for (size_t i = start; i < start + range; ++i)
		{
			auto p = _primitives[i];
			auto index = floorl((p->get_bbox().centroid()[d] - left) / space);
			buckets[index].push_back(p);
			bucketBBoxes[index].expand(p->get_bbox());
		}

		vector<BBox> leftBBoxes(BUCKET_SIZE);
		vector<BBox> rightBBoxes(BUCKET_SIZE);
		for (size_t i = 1; i <= BUCKET_SIZE - 1; ++i)
		{
			leftBBoxes[i] = leftBBoxes[i - 1];
			rightBBoxes[i] = rightBBoxes[i - 1];
			leftBBoxes[i].expand(bucketBBoxes[i - 1]);
			rightBBoxes[i].expand(bucketBBoxes[BUCKET_SIZE - i]);
		}

		size_t leftPartitionNum = 0, rightPartitionNum = 0;
		auto selectedPartition = -1;
		for (size_t i = 1; i <= BUCKET_SIZE - 1; ++i)
		{
			leftPartitionNum += buckets[i - 1].size();
			rightPartitionNum = range - leftPartitionNum;
			auto curSah = leftPartitionNum * leftBBoxes[i].surface_area() + rightPartitionNum * rightBBoxes[BUCKET_SIZE - 1].surface_area();
			if (curSah < sah)
			{
				sah = curSah;
				selectedPartition = i;
			}
		}

		if (selectedPartition != -1)
		{
			partitionA.clear();
			partitionB.clear();

			for (size_t i = 0; i < selectedPartition; ++i)
			{
				partitionA.insert(partitionA.end(), buckets[i].begin(), buckets[i].end());
			}
			for (size_t i = selectedPartition; i < BUCKET_SIZE; ++i)
			{
				partitionB.insert(partitionB.end(), buckets[i].begin(), buckets[i].end());
			}
		}
	}

	auto pIndex = start;
	for (auto p : partitionA)
	{
		_primitives[pIndex++] = p;
	}
	for (auto p : partitionB)
	{
		_primitives[pIndex++] = p;
	}

	if (partitionA.size() == range || partitionB.size() == range)
	{
		if (partitionB.size() == range)
		{
			swap(partitionA, partitionB);
		}
		auto half = partitionA.begin() + range / 2;
		partitionB.insert(partitionB.end(),
						  std::make_move_iterator(half),
						  std::make_move_iterator(partitionA.end()));
		partitionA.erase(half, partitionA.end());
	}

	{
		BBox bb;
		for (size_t i = 0; i < partitionA.size(); ++i)
		{
			bb.expand(partitionA[i]->get_bbox());
		}
		l = new BVHNode(bb, start, partitionA.size());
		l->Partition(_primitives, max_leaf_size);
	}

	{
		BBox bb;
		for (size_t i = 0; i < partitionB.size(); ++i)
		{
			bb.expand(partitionB[i]->get_bbox());
		}
		r = new BVHNode(bb, start + partitionA.size(), partitionB.size());
		r->Partition(_primitives, max_leaf_size);
	}
}

}  // namespace StaticScene
}  // namespace CMU462
