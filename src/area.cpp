#include <nori/emitter.h>
#include <nori/color.h>
#include <nori/mesh.h>

NORI_NAMESPACE_BEGIN
class AreaLight : public Emitter
{
private:
    Color3f m_radiance;
public:
    AreaLight(const PropertyList &propList){
        m_radiance = propList.getColor("radiance");
    }

    Color3f getRadiance() const{
        return m_radiance;
    }

    Color3f eval(const EmitterQueryRecord &record) const override {
        Vector3f x_y = record.emitter_samplingpoint - record.intersection_point;
        Vector3f l_i = (x_y).normalized();
        float cos_theta_i = l_i.dot(record.intersection_normal);
        if (cos_theta_i <= 0)
            return {0.f};
        return m_radiance * cos_theta_i;
    }

    bool getIncomingRay(const EmitterQueryRecord& emitter_record, const Scene* scene, Ray3f shadow_ray) override{
        Vector3f wo = (emitter_record.emitter_samplingpoint - emitter_record.intersection_point);
        float dist = wo.norm();
        wo.normalize();
        if ((-wo).dot(emitter_record.emitter_samplingnormal) <= 0) {
            return false;
        }
        shadow_ray.d = emitter_record.intersection_point;
        shadow_ray.o = wo;
        shadow_ray.maxt = dist - Epsilon;
        return true;
    }

    float getLightPdf(int size, const Mesh* emitter_mesh, const EmitterQueryRecord& record) override{
        Vector3f light_dir = record.intersection_point - record.emitter_samplingpoint;
        float length_squared = light_dir.squaredNorm();
        light_dir.normalize();
        float cos_theta = abs(light_dir.dot(record.emitter_samplingnormal));
        return clamp(1.f / size * emitter_mesh->getPdf() * length_squared / cos_theta, 0.f, 100000000000.f);
    }

    std::string toString() const {
        return "AreaLight[ radiance: " + m_radiance.toString() + "]";
    }
};

NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END