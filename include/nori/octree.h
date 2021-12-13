#pragma once
#include <nori/mesh.h>

NORI_NAMESPACE_BEGIN
struct treeNode;
class Octree{
    public:
    Octree(Mesh* mesh);
    BoundingBox3f WorldBound() const { return bounds;}
    bool IntersectkdTree(const Ray3f &ray_, Intersection &its, bool shadowRay) const;

    private:
    void buildTree(std::vector<int> meshIndex, BoundingBox3f bound, treeNode *parent, int nTri);
    Mesh* mesh;
    std::vector<int> meshIndices;
    BoundingBox3f bounds;
    std::vector<BoundingBox3f> tribounds;
    treeNode *root;
};
NORI_NAMESPACE_END
