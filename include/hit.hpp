#ifndef HIT_H
#define HIT_H

#include <vecmath.h>
#include "const.hpp"
#include "ray.hpp"

class Material;

class Hit {
public:

    // constructors
    Hit() {
        n = 0;
        t = INF;
        material = nullptr;
        radius = INIT_RADIUS;
        attenWeight = Vector3f(1);
        normal = flux = fluxLight = Vector3f::ZERO;
    }

    Hit(float _t, Material *m, const Vector3f &_normal) {
        n = 0;
        t = _t;
        material = m;
        normal = _normal;
        radius = INIT_RADIUS;
        attenWeight = Vector3f(1);
        flux = fluxLight = Vector3f::ZERO;
    }

    Hit(const Hit &h) {
        t = h.t;
        material = h.material;
        normal = h.normal;
    }

    // destructor
    ~Hit() = default;

    float getT() const {
        return t;
    }

    Material *getMaterial() const {
        return material;
    }

    const Vector3f &getNormal() const {
        return normal;
    }

    void set(float _t, Material *m, const Vector3f &n, const Vector3f _p) {
        t = _t;
        material = m;
        normal = n;
        p = _p;
    }

    int n;
    float t, radius; //长度，管辖半径
    Material *material; //材料
    Vector3f flux, fluxLight, attenWeight; //光通量、颜色
    Vector3f normal, p; //法向量，坐标
};

inline std::ostream &operator<<(std::ostream &os, const Hit &h) {
    os << "Hit <" << h.getT() << ", " << h.getNormal() << ">";
    return os;
}

#endif // HIT_H
