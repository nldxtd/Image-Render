#ifndef OBJECT_H
#define OBJECT_H

#include "ray.hpp"
#include "hit.hpp"
#include "material.hpp"
#include "utils.hpp"
#include "const.hpp"
#include "bound.hpp"
#include <math.h>
#include <vector>
#include <tuple>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <utility>
#include <sstream>
#include <cmath>
#include <float.h>

class Object {
public:
    Material* material;
    Object(): material(nullptr) {}
    explicit Object(Material* material) {
        this->material = material;
    }
    virtual ~Object() = default;

    //ray intersect with object
    virtual bool intersect(const Ray &r, Hit &h, float tmin) = 0;
    //min, max, center
    virtual Vector3f vecmin() = 0;
    virtual Vector3f vecmax() = 0;
    virtual Vector3f veccenter() = 0;
    //shooting random ray
    virtual Ray shootRay(int axis=-1, long long int seed=0) {
        Vector3f zero = Vector3f::ZERO;
        return Ray(zero, zero);
    }
};

class Sphere : public Object {
public:
    Sphere() {
        // unit ball at the center
        this->Center = Vector3f(0, 0, 0);
        this->Radius = 1;
    }

    Sphere(const Vector3f &center, float radius, Material *material) : Object(material) {
        // giving center, radius, material
        this->Center = center;
        this->Radius = radius;
    }

    ~Sphere() override = default;

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        // Ray intersect with self at hit, with minimum t = tmin
        const float eps = 1e-10;

        // Computing params
        Vector3f dir = r.getDirection().normalized();
        Vector3f origin_minus_center = this->Center - r.getOrigin();
        // What about mistakes?
        float ans_t = 0;
        float radius_square = this->Radius*this->Radius;
        float midpoint = Vector3f::dot(dir, origin_minus_center);
        // Distance is r^2-d^2
        float distance = sqrt(radius_square - origin_minus_center.squaredLength() + midpoint * midpoint);
        if(distance >= 0) {
            // Inside sphere
            if(origin_minus_center.squaredLength() < radius_square) {
                ans_t = midpoint + distance;
            }
            // Outside sphere
            else if(origin_minus_center.squaredLength() > radius_square) {
                ans_t = midpoint - distance;
            }
            if(ans_t >= tmin && ans_t < h.getT()) {
                // Computing normal vector
                Vector3f hit_point = r.getOrigin() + dir*ans_t;
                Vector3f normal_vector = Vector3f::ZERO;
                if(origin_minus_center.squaredLength() < radius_square) 
                    normal_vector = (this->Center-hit_point).normalized();
                else 
                    normal_vector = (hit_point-this->Center).normalized();
                Vector3f point = r.pointAtParameter(ans_t);
                h.set(ans_t, this->material, normal_vector, point);
                return true;
            }
        }
        return false;
    }

    Ray shootRay(int axis=-1, long long int seed=0) override {
        float u = 2*random(axis, seed) - 1, v = 2*random(axis, seed) - 1;
        float r2 =u * u + v * v;
        while(r2>=1) {
            ++seed;
            u = 2*random(axis, seed) - 1;
            v = 2*random(axis, seed) - 1;
            r2 = u * u + v * v;
        }
        Vector3f dir(2*u*sqrtf(1-r2), 2*v*sqrt(1-r2),1-2*r2);
        dir.normalize();
        Vector3f ori = Center+Radius*dir;
        // if(ori.y() <= 8) {
        //     printf("generating ray shooting from %f, %f, %f\n", ori.x(), ori.y(), ori.z());
        // }
        return Ray(ori, dir);
    }

    Vector3f vecmin() override { 
        return Vector3f(Center.x() - Radius, Center.y() - Radius, Center.z() - Radius);
    }
    Vector3f vecmax() override { 
        return Vector3f(Center.x() + Radius, Center.y() + Radius, Center.z() + Radius);
    }
    Vector3f veccenter() override { return Center; }

protected:

    Vector3f Center;
    float Radius;
};

class Plane : public Object {
public:
    Plane() {
        // xy plane
        this->Normal = Vector3f(0, 0, 1);
        this->D = 0;
    }

    Plane(const Vector3f &normal, float d, Material *m) : Object(m) {
        this->Normal = normal;
        this->D = -d;
    }

    ~Plane() override = default;

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        Vector3f dir = r.getDirection().normalized();
        Vector3f ori = r.getOrigin();
        float ans_t = -((this->D + Vector3f::dot(ori, this->Normal))/Vector3f::dot(this->Normal, dir));
        if(ans_t > 0 && ans_t >= tmin && ans_t < h.getT()) {
            // Computing normal vector
            Vector3f normal_vector = (Vector3f::dot(dir, this->Normal) <= 0) ? this->Normal : -this->Normal;
            Vector3f point = r.pointAtParameter(ans_t);
            h.set(ans_t, this->material, this->Normal, point);
            return true;
        }
        return false;
    }

    Vector3f vecmin() override {
        return -INF * Vector3f(fabs(Normal.x()) < 1 - FLT_EPSILON,
                               fabs(Normal.y()) < 1 - FLT_EPSILON,
                               fabs(Normal.z()) < 1 - FLT_EPSILON) +
               Normal * D;
    }
    Vector3f vecmax() override {
        return INF * Vector3f(fabs(Normal.x()) < 1 - FLT_EPSILON,
                              fabs(Normal.y()) < 1 - FLT_EPSILON,
                              fabs(Normal.z()) < 1 - FLT_EPSILON) +
               Normal * D;
    }
    Vector3f veccenter() override { return Normal * D; }

protected:
    Vector3f Normal;
    float D;
};

class Triangle: public Object {

public:
	Triangle() = delete;

    // a b c are three vertex positions of the triangle
	Triangle( const Vector3f& a, const Vector3f& b, const Vector3f& c, Material* m) : Object(m) {
		this->vertices[0] = a;
		this->vertices[1] = b;
		this->vertices[2] = c;
		Vector3f param_e1 = this->vertices[0] - this->vertices[1];
		Vector3f param_e2 = this->vertices[0] - this->vertices[2];
		this->normal = Vector3f::cross(param_e1, param_e2).normalized();
        this->minvec = minE(minE(a, b), c);
        this->maxvec = maxE(maxE(a, b), c);
        this->center = (a + b + c) / 3;
	}

	bool intersect( const Ray& ray,  Hit& hit , float tmin) override {
		Vector3f param_e1 = this->vertices[0] - this->vertices[1];
		Vector3f param_e2 = this->vertices[0] - this->vertices[2];
		Vector3f param_s = this->vertices[0] - ray.getOrigin();
		Vector3f dir = ray.getDirection();
		float param_under = Matrix3f(dir, param_e1, param_e2).determinant();
		float param_anst = Matrix3f(param_s, param_e1, param_e2).determinant();
		float param_beta = Matrix3f(dir, param_s, param_e2).determinant();
		float param_gama = Matrix3f(dir, param_e1, param_s).determinant();
		float anst = param_anst/param_under;
		float beta = param_beta/param_under;
		float gama = param_gama/param_under;
		if(anst < hit.getT() && anst >= tmin) {
			if(gama >= 0 && beta >= 0 && gama+beta <= 1) {
				Vector3f normal_vector = (Vector3f::dot(dir, this->normal) <= 0) ? this->normal : -this->normal;
				Vector3f point = ray.pointAtParameter(anst);
				hit.set(anst, this->material, normal_vector, point);
				return true;
			}
		}
        return false;
	}

    Vector3f vecmin() override { return minvec; }
    Vector3f vecmax() override { return maxvec; }
    Vector3f veccenter() override { return center; }

	Vector3f normal;
	Vector3f vertices[3];
    Vector3f minvec, maxvec, center;
};

class Mesh : public Object {

public:
    struct TriangleIndex {
        TriangleIndex() {
            x[0] = 0; x[1] = 0; x[2] = 0;
        }
        int &operator[](const int i) { return x[i]; }
        // By Computer Graphics convention, counterclockwise winding is front face
        int x[3]{};
    };
    std::vector<Vector3f> v;
    std::vector<TriangleIndex> t;
    std::vector<Vector3f> n;
    BoundBox boundbox;
   
    Mesh(const char *filename, Material *material) : Object(material) {

        // Optional: Use tiny obj loader to replace this simple one.
        std::ifstream f;
        f.open(filename);
        if (!f.is_open()) {
            std::cout << "Cannot open " << filename << "\n";
            return;
        }
        std::string line;
        std::string vTok("v");
        std::string fTok("f");
        std::string texTok("vt");
        char bslash = '/', space = ' ';
        std::string tok;
        int texID;
        while (true) {
            std::getline(f, line);
            if (f.eof()) {
                break;
            }
            if (line.size() < 3) {
                continue;
            }
            if (line.at(0) == '#') {
                continue;
            }
            std::stringstream ss(line);
            ss >> tok;
            if (tok == vTok) {
                Vector3f vec;
                ss >> vec[0] >> vec[1] >> vec[2];
                v.push_back(vec);
                boundbox.updateBound(vec);
            } else if (tok == fTok) {
                if (line.find(bslash) != std::string::npos) {
                    std::replace(line.begin(), line.end(), bslash, space);
                    std::stringstream facess(line);
                    TriangleIndex trig;
                    facess >> tok;
                    for (int ii = 0; ii < 3; ii++) {
                        facess >> trig[ii] >> texID;
                        trig[ii]--;
                    }
                    t.push_back(trig);
                } else {
                    TriangleIndex trig;
                    for (int ii = 0; ii < 3; ii++) {
                        ss >> trig[ii];
                        trig[ii]--;
                    }
                    t.push_back(trig);
                }
            } else if (tok == texTok) {
                Vector2f texcoord;
                ss >> texcoord[0];
                ss >> texcoord[1];
            }
        }
        computeNormal();

        f.close();
    }

    void computeNormal() {
        n.resize(t.size());
        for (int triId = 0; triId < (int) t.size(); ++triId) {
            TriangleIndex& triIndex = t[triId];
            Vector3f a = v[triIndex[1]] - v[triIndex[0]];
            Vector3f b = v[triIndex[2]] - v[triIndex[0]];
            b = Vector3f::cross(a, b);
            n[triId] = b / b.length();
        }
    }

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        // Optional: Change this brute force method into a faster one.
        bool result = false;
        for (int triId = 0; triId < (int) t.size(); ++triId) {
            TriangleIndex& triIndex = t[triId];
            Triangle triangle(v[triIndex[0]],
                            v[triIndex[1]], v[triIndex[2]], material);
            triangle.normal = n[triId];
            result |= triangle.intersect(r, h, tmin);
        }
        return result;
    }

    Vector3f vecmin() override { return boundbox.bounds[0]; }
    Vector3f vecmax() override { return boundbox.bounds[1]; }
    Vector3f veccenter() override {
        return (boundbox.bounds[0] + boundbox.bounds[1]) / 2;
    }
};

struct CurvePoint {
    Vector3f V; // Vertex
    Vector3f T; // Tangent  (unit)
};

class Curve {
public:
    std::vector<Vector3f> controls;
    std::vector<CurvePoint> data;
    std::vector<double> knot;

    float mind, maxd, radius;
    float range[2];

    explicit Curve(std::vector<Vector3f> points) : controls(std::move(points)) {
        mind = INF, maxd = -INF, radius = 0;
        for(auto cp : controls) {
            mind = std::min(cp.y(), mind);
            maxd = std::max(cp.y(), maxd);
            radius = std::max(radius, std::max(fabs(cp.x()), fabs(cp.z())));
        }
    }

    bool intersect(const Ray &r, Hit &h, float tmin) {
        return false;
    }

    std::vector<Vector3f> &getControls() {
        return controls;
    }

    virtual void discretize(int resolution, std::vector<CurvePoint>& data) = 0;
};

class BezierCurve : public Curve {
public:
    explicit BezierCurve(const std::vector<Vector3f> &points) : Curve(points) {
        if (points.size() < 4 || points.size() % 3 != 1) {
            printf("Number of control points of BezierCurve must be 3n+1!\n");
            exit(0);
        }
        range[0] = 0;
        range[1] = 1;
        discretize(RESOLUTION, data);
    }

    int count(int a, int b) {
        if(a < b) return 0;
        int res = 1;
        for(int i = 1; i <= b; i++) {
            res *= (a - i + 1);
            res /= i;
        }
        return res;
    }

    double Bernstein(int i, int n, double t) {
        if(t == 0) return i == 0;
        if(t == 1) return i == n;
        return count(n, i) * pow(t, i) * pow(1-t, n-i);
    }

    void discretize(int resolution, std::vector<CurvePoint>& data) override {
        data.clear();
        int n = controls.size() - 1;
        for(int i = 0; i <= resolution; i++) {
            double t = (double)i / (resolution+1);
            knot.push_back(t);
            Vector3f point = Vector3f::ZERO;
            Vector3f tan = Vector3f::ZERO;
            for(int j = 0; j < n; j++) {
                point += controls[j] * Bernstein(j, n, t);
                tan += n * Bernstein(j, n - 1, t) * (controls[j+1]-controls[j]);
            }
            point += controls[n] * Bernstein(n, n, t);
            CurvePoint cp = {point, tan};
            data.push_back(cp);
        }
    }

protected:

};

class BsplineCurve : public Curve {

public:
    BsplineCurve(const std::vector<Vector3f> &points) : Curve(points) {
        if (points.size() < 4) {
            printf("Number of control points of BspineCurve must be more than 4!\n");
            exit(0);
        }
        int n = controls.size() - 1;
        int k = 3;
        range[0] = (float)k/(n+k);
        range[1] = (float)n/(n+k);
        discretize(RESOLUTION, data);
    }

    void discretize(int resolution, std::vector<CurvePoint>& data) override {
        data.clear();
        int n = controls.size() - 1;
        int k = 3;
        for(int i = 0; i < n+k+2; i++)
            knot.push_back((double)i / (n+k+1));
        for(int j = k; j <= n; j++) {
            for(int m = 1; m <= resolution; m++) {
                double t = m * (knot[j+1] - knot[j]) / (resolution+1) + knot[j];
                Vector3f point = Vector3f::ZERO;
                Vector3f tan = Vector3f::ZERO;
                vector<vector<double>>bpoint = vector<vector<double>>(n + k + 1, vector<double>(k + 1, 0.0));
                for(int bi = 0; bi < n + k + 1; bi++) {
                    if(t >= knot[bi] && t < knot[bi+1]) bpoint[bi][0] = 1;
                    else bpoint[bi][0] = 0;
                }
                for(int bj = 1; bj <= k; bj++) {
                    for(int bi = 0; bi <= n+k-bj; bi++) {
                        bpoint[bi][bj] = (t - knot[bi]) / (knot[bi + bj] - knot[bi]) * bpoint[bi][bj - 1] +
                    (knot[bi + bj + 1] - t) / (knot[bi + bj + 1] - knot[bi + 1]) * bpoint[bi + 1][bj - 1];
                    }
                }
                for (int i = 0; i <= n; i++) {
                    point += controls[i] * bpoint[i][k];
                    tan += controls[i] * k * (bpoint[i][k - 1] / (knot[i + k] - knot[i]) - bpoint[i + 1][k - 1] / (knot[i + k + 1] - knot[i + 1]));
                }
                CurvePoint cp = { point, tan };
                data.push_back(cp);
            }
        }
    }

protected:

};

class RevSurface : public Object {

    Curve *pCurve;
    BoundBox boundbox;

public:
    RevSurface(Curve *pCurve, Material* material) : pCurve(pCurve), Object(material) {
        // Check flat.
        for (const auto &cp : pCurve->getControls()) {
            if (cp.z() != 0.0) {
                printf("Profile of revSurface must be flat on xy plane.\n");
                exit(0);
            }
        }
        Vector3f lo = Vector3f(-pCurve->radius, pCurve->mind - 2, -pCurve->radius);
        Vector3f hi = Vector3f(pCurve->radius, pCurve->maxd + 2, pCurve->radius);
        boundbox.set(lo, hi);
    }

    ~RevSurface() override {
        delete pCurve;
    }

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        // (PA3 optional TODO): implement this for the ray-tracing routine using G-N iteration.
        float t;
        if(!boundbox.intersect(r, t) || t > h.getT()) return false;
        Vector3f normal;
        Vector3f point = r.pointAtParameter(t);
        float theta = M_PI + atan2(-point.z(), point.x());
        float mu = (pCurve->maxd-point.y())/(pCurve->maxd-pCurve->mind);
        bool isHit = false;
        //newton step
        Vector3f dtheta, dmu, hp;
        Quat4f rot;
        for(int i = 0; i < NEWTON_STEPS; i++) {
            theta = ((theta < 0 ? theta + 2 * M_PI : theta));
            if (theta >= 2 * M_PI) theta = fmod(theta, 2 * M_PI);
            if(mu >= 1) mu = 1.0 - FLT_EPSILON;
            if(mu <= 0) mu = FLT_EPSILON;
            rot.setAxisAngle(theta, Vector3f::UP);
            Matrix3f rotMat = Matrix3f::rotation(rot);
            //get point on spline
            int bpos = upper_bound(pCurve->knot.begin(), pCurve->knot.end(), mu) - pCurve->knot.begin() - 1;
            CurvePoint pt = pCurve->data[bpos];
            //rotate to scene
            hp = rotMat * pt.V;
            dmu = rotMat * pt.T;
            dtheta = Vector3f(-pt.V.x() * sin(theta), 0, -pt.V.x() * cos(theta));
            Vector3f hit2hit = point - hp;
            normal = Vector3f::cross(dmu, dtheta);
            if(hit2hit.squaredLength() < NEWTON_EPS) {
                return false;
            }
            float length  = Vector3f::dot(r.direction, normal);
            t -= Vector3f::dot(dmu, Vector3f::cross(dtheta, point)) / length;
            mu -= Vector3f::dot(r.direction, Vector3f::cross(dtheta, point)) / length;
            theta += Vector3f::dot(r.direction, Vector3f::cross(dmu, point)) / length;
            isHit = true;
        }
        if(!isnormal(mu) || !isnormal(theta) || !isnormal(t)) return false;
        if(t < 0 || mu < pCurve->range[0] || mu > pCurve->range[1] || t > h.getT()) return false;
        h.set(t, material, normal.normalized(), hp);
        return isHit;
    }
    Vector3f vecmin() override { return boundbox.bounds[0]; }
    Vector3f vecmax() override { return boundbox.bounds[1]; }
    Vector3f veccenter() override {
        return (boundbox.bounds[0] + boundbox.bounds[1]) / 2;
    }
};


#endif