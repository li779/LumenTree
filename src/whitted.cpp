#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN
class WhittedIntegrator : public Integrator{
    public:
    WhittedIntegrator(const PropertyList &props){}
    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const override{
        
        Intersection its;
        Color3f light_eval(0.f);
        Color3f emitted_light(0.f);

        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);
        if(its.mesh->getBSDF()->isDiffuse()){
            Mesh *emitter_mesh;
            const auto &emitter_meshes = scene->getEmitters();
            int emitter_sample_index = (int) (sampler->next1D() * emitter_meshes.size());
            float emitter_pdf = 1.f / (float) emitter_meshes.size();
            emitter_mesh = emitter_meshes[emitter_sample_index];

            Emitter *emitter = emitter_mesh->getEmitter();

            float emitter_samplingpdf;
            Normal3f emitter_samplingnormal;
            Point3f emitter_samplingpoint;
            emitter_mesh->samplingMeshSurface(sampler, emitter_samplingpoint, emitter_samplingnormal, emitter_samplingpdf);
            EmitterQueryRecord emitter_record(its.p, its.shFrame.n, emitter_samplingpoint, emitter_samplingnormal);

            if (its.mesh == emitter_mesh)
                    emitted_light += emitter_mesh->getEmitter()->getRadiance();

            Ray3f shadow_ray;    
            if(emitter->getIncomingRay(emitter_record, scene, shadow_ray)){
                if(!scene->rayIntersect(shadow_ray)){
                    light_eval += emitter->eval(emitter_record);
                }
            }

            Vector3f l_i = emitter_samplingpoint - its.p;
            float dist_2 = l_i.dot(l_i);
            l_i = l_i.normalized();
            float cos_theta_light = (-l_i).dot(emitter_samplingnormal);
            BSDFQueryRecord bsdf_record(its.toLocal(l_i), its.toLocal(-ray.d), ESolidAngle);
            Color3f bsdf_eval = its.mesh->getBSDF()->eval(bsdf_record);
            // convert probability measure
            emitter_samplingpdf *= dist_2 / cos_theta_light;
            return emitted_light + light_eval / emitter_samplingpdf / emitter_pdf * bsdf_eval;
        }else{
            BSDFQueryRecord bRec(its.toLocal(-ray.d));
            Color3f ref_color = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
            if (sampler->next1D() < 0.95 && ref_color.x() > 0.f) {
                return Li(scene, sampler, Ray3f(its.p, its.toWorld(bRec.wo))) / 0.95 * ref_color;
            } else {
                return Color3f(0.0f);
            }
        }

    }
    std::string toString() const {
        return "WhittedIntegrator[]";
    }
};
NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");
NORI_NAMESPACE_END
