#ifndef MATERIAL_H
#define MATERIAL_H

#include <cassert>
#include <vecmath.h>

#include "ray.hpp"
#include "hit.hpp"
#include <iostream>

// TODO: Implement Shade function that computes Phong introduced in class.
class Material {
public:

    explicit Material(const Vector3f &d_color, const Vector3f &s_color = Vector3f::ZERO, float s = 0, 
                      const Vector3f &e_color = Vector3f::ZERO, float r = 1, Vector3f t = Vector3f(1, 0, 0)) :
            diffuseColor(d_color), specularColor(s_color), shininess(s),
            emission(e_color), reflect(r), type(t) {

            }

    virtual ~Material() = default;

    virtual Vector3f getDiffuseColor() const {
        return diffuseColor;
    }

    float clamp(const Vector3f& x, const Vector3f& y) {
        return (Vector3f::dot(x, y) >= 0) ? Vector3f::dot(x, y) : 0;
    }
    
    Vector3f diffuseColor;
    Vector3f specularColor;
    Vector3f emission;
    Vector3f type;
    float shininess;
    float reflect;
};


#endif // MATERIAL_H
