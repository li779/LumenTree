# pragma once
#include <embree3/rtcore.h>
#include <nori/mesh.h>

NORI_NAMESPACE_BEGIN

class accel_embree
{
private:
    RTCDevice embree_device = nullptr;
    RTCScene embree_scene = nullptr;
    std::vector<Mesh *> m_mesh;
    RTCGeometry embree_geometry(Mesh* mesh);
    
public:
    accel_embree(std::vector<Mesh *> &mesh);
    ~accel_embree();
    bool IntersectWithDetail(Ray3f &ray_, Intersection &its, int &triIndex) const;
    bool IntersectTest(Ray3f &ray_, Intersection &its,int &triIndex) const;
};

NORI_NAMESPACE_END