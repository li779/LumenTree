#include <nori/accel_embree.h>

NORI_NAMESPACE_BEGIN

accel_embree::accel_embree(std::vector<Mesh *> &mesh){
    this->m_mesh = mesh;
    embree_device = rtcNewDevice("");
    embree_scene = rtcNewScene(embree_device);
    rtcSetSceneFlags(embree_scene, RTC_SCENE_FLAG_DYNAMIC);
    for(uint32_t i = 0; i < m_mesh.size(); i++){
        rtcAttachGeometryByID(embree_scene,embree_geometry(m_mesh[i]),i);
    }
    rtcCommitScene(embree_scene);
}

accel_embree::~accel_embree(){rtcReleaseScene(embree_scene);}


RTCGeometry accel_embree::embree_geometry(Mesh* mesh){
    RTCGeometry geom = rtcNewGeometry(this->embree_device, RTC_GEOMETRY_TYPE_TRIANGLE);

    rtcSetSharedGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
                               mesh->getVertexPositions().data(), 0, 3 * sizeof(float),
                               mesh->getVertexCount());
    rtcSetSharedGeometryBuffer(geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
                               mesh->getIndices().data(), 0, 3 * sizeof(uint32_t),
                               mesh->getTriangleCount());

    rtcCommitGeometry(geom);
    return geom;
}

bool accel_embree::IntersectTest(Ray3f &ray_, Intersection &its,int &triIndex) const{
    RTCIntersectContext context;
    rtcInitIntersectContext(&context);
    RTCRay ray2;
    ray2.org_x = ray_.o.x();
    ray2.org_y = ray_.o.y();
    ray2.org_z = ray_.o.z();
    ray2.tnear = ray_.mint;
    ray2.dir_x = ray_.d.x();
    ray2.dir_y = ray_.d.y();
    ray2.dir_z = ray_.d.z();
    ray2.time = 0;
    ray2.tfar = ray_.maxt;
    ray2.mask = 0;
    ray2.id = 0;
    ray2.flags = 0;
    rtcOccluded1(embree_scene, &context, &ray2);
    return ray2.tfar != ray_.maxt;
}

bool accel_embree::IntersectWithDetail(Ray3f &ray_, Intersection &its, int &triIndex) const{
    RTCIntersectContext context;
    rtcInitIntersectContext(&context);
    RTCRayHit rh;
    rh.ray.org_x = ray_.o.x();
    rh.ray.org_y = ray_.o.y();
    rh.ray.org_z = ray_.o.z();
    rh.ray.tnear = ray_.mint;
    rh.ray.dir_x = ray_.d.x();
    rh.ray.dir_y = ray_.d.y();
    rh.ray.dir_z = ray_.d.z();
    rh.ray.time = 0;
    rh.ray.tfar = ray_.maxt;
    rh.ray.mask = 0;
    rh.ray.id = 0;
    rh.ray.flags = 0;
    rtcIntersect1(embree_scene, &context, &rh);
    if (rh.ray.tfar != ray_.maxt) {
        ray_.maxt = its.t = rh.ray.tfar;
        its.uv = Point2f(rh.hit.u, rh.hit.v);
        its.mesh = m_mesh[rh.hit.geomID];
        triIndex = rh.hit.primID;
        return true;
    }else{
        return false;
    }
}


NORI_NAMESPACE_END