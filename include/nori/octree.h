#pragma once
#include <nori/mesh.h>

NORI_NAMESPACE_BEGIN
struct treeNode;
class Octree{
    public:
    Octree(std::vector<Mesh *> &mesh, BoundingBox3f &bbox);
    BoundingBox3f WorldBound() const { return bounds;}
    bool IntersectOctree(Ray3f &ray_, Intersection &its, bool shadowRay,int &triIndex) const;

    private:
    int buildTree(std::vector<std::pair<int,int>> meshIndex, BoundingBox3f bound, treeNode *parent, int nTri);
    bool recursive(const treeNode &node, Ray3f &ray_, Intersection &its, BoundingBox3f bbox,int &triIndex) const;
    bool recursiveTest(const treeNode &node, Ray3f &ray_, BoundingBox3f bbox, int &triIndex) const;
    std::vector<Mesh *> mesh;
    std::vector<std::pair<int,int>> meshIndices;
    BoundingBox3f bounds;
    std::vector<std::vector<BoundingBox3f>> tribounds;
    treeNode *root;
};
NORI_NAMESPACE_END
