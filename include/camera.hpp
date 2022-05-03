#ifndef CAMERA_H
#define CAMERA_H

#include "ray.hpp"
#include <vecmath.h>
#include <float.h>
#include <cmath>


class Camera {
public:
    Camera(const Vector3f &center, const Vector3f &direction, const Vector3f &up, int imgW, int imgH) {
        this->center = center;
        this->direction = direction.normalized();
        this->horizontal = Vector3f::cross(this->direction, up);
		this->horizontal.normalize();
        this->up = Vector3f::cross(this->horizontal, this->direction);
        this->width = imgW;
        this->height = imgH;
    }

    // Generate rays for each screen-space coordinate
    virtual Ray generateRay(const Vector2f &point) = 0;
    virtual ~Camera() = default;

    int getWidth() const { return width; }
    int getHeight() const { return height; }

protected:
    // Extrinsic parameters
    Vector3f center;
    Vector3f direction;
    Vector3f up;
    Vector3f horizontal;
    // Intrinsic parameters
    int width;
    int height;
};

// TODO: Implement Perspective camera
// You can add new functions or variables whenever needed.
class PerspectiveCamera : public Camera {

public:
    PerspectiveCamera(const Vector3f &center, const Vector3f &direction,
            const Vector3f &up, int imgW, int imgH, float angle) : Camera(center, direction, up, imgW, imgH) {
        // angle is in radian.
        this->fx = (this->width/2)/tan(angle/2);
        this->fy = (this->height/2)/tan(angle/2);
    }

    virtual ~PerspectiveCamera() = default;

    Ray generateRay(const Vector2f &point) override {
        // generate ray emitted from point
        Vector3f ray_camera = Vector3f((point.x()-this->width/2)/fx, (this->height/2-point.y())/fy, 1).normalized();
        Matrix3f R(this->horizontal, -this->up, this->direction);
        Vector3f ray_world = R*ray_camera;
        return Ray(this->center, ray_world);
    }

protected:
    float fx;
    float fy;
};

#endif //CAMERA_H
