#ifndef BOUND_H
#define BOUND_H
#include <vecmath.h>
#include <vector>
#include "ray.hpp"
#include "const.hpp"
using std::vector;
class BoundBox {
   public:
    BoundBox() {
        bounds[0] = Vector3f(INF);
        bounds[1] = Vector3f(-INF);
    }

    BoundBox(const Vector3f &min, const Vector3f &max) {
        bounds[0] = min;
        bounds[1] = max;
    }

    void set(const Vector3f &lo, const Vector3f &hi) {
        bounds[0] = lo;
        bounds[1] = hi;
    }

    void updateBound(const Vector3f &vec) {
        for (int i = 0; i < 3; ++i) {
            bounds[0][i] = bounds[0][i] < vec[i] ? bounds[0][i] : vec[i];
            bounds[1][i] = bounds[1][i] < vec[i] ? vec[i] : bounds[1][i];
        }
    }

    bool intersect(const Ray &r, float &t_min) {
        Vector3f o(r.getOrigin()), invdir(1 / r.getDirection());
        vector<int> sgn = {invdir.x() < 0, invdir.y() < 0, invdir.z() < 0};
        t_min = INF;
        float tmin, tmax, tymin, tymax, tzmin, tzmax;
        tmin = (bounds[sgn[0]].x() - o.x()) * invdir.x();
        tmax = (bounds[1 - sgn[0]].x() - o.x()) * invdir.x();
        tymin = (bounds[sgn[1]].y() - o.y()) * invdir.y();
        tymax = (bounds[1 - sgn[1]].y() - o.y()) * invdir.y();
        if ((tmin > tymax) || (tymin > tmax)) return false;
        if (tymin > tmin) tmin = tymin;
        if (tymax < tmax) tmax = tymax;
        tzmin = (bounds[sgn[2]].z() - o.z()) * invdir.z();
        tzmax = (bounds[1 - sgn[2]].z() - o.z()) * invdir.z();
        if ((tmin > tzmax) || (tzmin > tmax)) return false;
        if (tzmin > tmin) tmin = tzmin;
        if (tzmax < tmax) tmax = tzmax;
        t_min = tmin;
        return true;
    }
    Vector3f bounds[2];
};

#endif  // !BOUND_H