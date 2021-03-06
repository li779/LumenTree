/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/bsdf.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF {
public:
    Dielectric(const PropertyList &propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);
    }

    Color3f eval(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    float pdf(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
        float cos_theta_i = Frame::cosTheta(bRec.wi);

        float kr = fresnel(cos_theta_i, m_extIOR, m_intIOR);
        bRec.measure = EDiscrete;

        if (sample.x() < kr) {
            // reflect
            bRec.wo = Vector3f(
                    -bRec.wi.x(),
                    -bRec.wi.y(),
                    bRec.wi.z()
            );
            bRec.eta = 1.f;
            // pdf weighted ray color, fresnel term reduces
            return {1.f};
        } else {
            //refract
            bRec.eta = cos_theta_i >= 0 ?  m_extIOR / m_intIOR : m_intIOR / m_extIOR;
            Normal3f n = cos_theta_i < 0 ? Normal3f(0.f, 0.f, -1.f) : Normal3f(0.f, 0.f, 1.f);
            cos_theta_i = abs(cos_theta_i);
            float cos_theta_o = sqrt(1 - bRec.eta * bRec.eta * fmax(0.f, 1.f - cos_theta_i * cos_theta_i));
            bRec.wo = bRec.eta * -bRec.wi + (bRec.eta * cos_theta_i - cos_theta_o) * n;
            bRec.wo.normalize();
            // pdf weighted ray color, fresnel term reduces
            return {bRec.eta * bRec.eta};
        }
    }

    std::string toString() const {
        return tfm::format(
            "Dielectric[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_intIOR, m_extIOR);
    }
private:
    float m_intIOR, m_extIOR;
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END
