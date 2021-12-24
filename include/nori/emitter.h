/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#pragma once

#include <nori/object.h>
#include <nori/ray.h>
#include <nori/mesh.h>

NORI_NAMESPACE_BEGIN

struct EmitterQueryRecord {
    Point3f intersection_point;
    Normal3f intersection_normal;
    Point3f emitter_samplingpoint;
    Normal3f emitter_samplingnormal;

    EmitterQueryRecord(Point3f intersection_point, Normal3f intersection_normal, Point3f emitter_samplingpoint, Normal3f emitter_samplingnormal)
        : intersection_point(std::move(intersection_point)),
        intersection_normal(std::move(intersection_normal)),
        emitter_samplingpoint(std::move(emitter_samplingpoint)),
        emitter_samplingnormal(std::move(emitter_samplingnormal))
    {}
};

/**
 * \brief Superclass of all emitters
 */
class Emitter : public NoriObject {
public:

    /**
     * \brief Return the type of object (i.e. Mesh/Emitter/etc.) 
     * provided by this instance
     * */
    EClassType getClassType() const { return EEmitter; }

    virtual Color3f eval(const EmitterQueryRecord& record) const = 0;
    virtual Color3f getRadiance() const = 0;
    virtual bool getIncomingRay(const EmitterQueryRecord& emitter_record, const Scene* scene, Ray3f &shadow_ray) =0;
    virtual float getLightPdf(int size, const Mesh* emitter_mesh, const EmitterQueryRecord& record)=0;
};

NORI_NAMESPACE_END
