#ifndef GROUP_H
#define GROUP_H

#include "bound.hpp"
#include "object.hpp"
#include "ray.hpp"
#include "hit.hpp"
#include "const.hpp"
#include <iostream>
#include <vector>
#include <map>
using std::vector;
using std::map;

class ObjectTreeNode {
public:
    Vector3f vecmin, vecmax;
    vector<Object*>* objs;
    ObjectTreeNode *lc, *rc;
    int l, r;
    bool isInside(Object* obj) {
        Vector3f objMin = obj->vecmin();
        Vector3f objMax = obj->vecmax();
        return (objMin.x() < objMax.x() || objMin.x() == objMax.x() && objMin.x() == objMax.x()) &&
               (objMax.x() > objMin.x() || objMax.x() == objMin.x() && objMin.x() == objMax.x()) &&
               (objMin.y() < objMax.y() || objMin.y() == objMax.y() && objMin.y() == objMax.y()) &&
               (objMax.y() > objMin.y() || objMax.y() == objMin.y() && objMin.y() == objMax.y()) &&
               (objMin.z() < objMax.z() || objMin.z() == objMax.z() && objMin.z() == objMax.z()) &&
               (objMax.z() > objMin.z() || objMax.z() == objMin.z() && objMin.z() == objMax.z());
    }
};

class ObjectKDTree {
public:
    int n;
    Vector3f** vertices;
    ObjectTreeNode* build(int depth, int dim, vector<Object*>* faces, const Vector3f& min, const Vector3f& max) {
        ObjectTreeNode* newNode = new ObjectTreeNode;
        newNode->vecmin = min;
        newNode->vecmax = max;
        Vector3f maxL, minR;
        if (dim == 0) {
            maxL = Vector3f((newNode->vecmin.x() + newNode->vecmax.x()) / 2, newNode->vecmax.y(), newNode->vecmax.z());
            minR = Vector3f((newNode->vecmin.x() + newNode->vecmax.x()) / 2, newNode->vecmin.y(), newNode->vecmin.z());
        } else if (dim == 1) {
            maxL = Vector3f(newNode->vecmax.x(), (newNode->vecmin.y() + newNode->vecmax.y()) / 2, newNode->vecmax.z());
            minR = Vector3f(newNode->vecmin.x(), (newNode->vecmin.y() + newNode->vecmax.y()) / 2, newNode->vecmin.z());
        } else {
            maxL = Vector3f(newNode->vecmax.x(), newNode->vecmax.y(), (newNode->vecmin.z() + newNode->vecmax.z()) / 2);
            minR = Vector3f(newNode->vecmin.x(), newNode->vecmin.y(), (newNode->vecmin.z() + newNode->vecmax.z()) / 2);
        }
        newNode->objs = new vector<Object*>;
        for (auto face : *faces)
            if (newNode->isInside(face)) newNode->objs->push_back(face);

        if (newNode->objs->size() > OBJ_FACES && depth < OBJ_DEPTH) {
            newNode->lc = build(depth + 1, (dim + 1) % 3, newNode->objs, min, maxL);
            newNode->rc = build(depth + 1, (dim + 1) % 3, newNode->objs, minR, max);

            vector<Object*>*faceL = newNode->lc->objs, *faceR = newNode->rc->objs;
            map<Object*, int> cnt;
            for (auto face : *faceL) cnt[face]++;
            for (auto face : *faceR) cnt[face]++;
            newNode->lc->objs = new vector<Object*>;
            newNode->rc->objs = new vector<Object*>;
            newNode->objs->clear();
            for (auto face : *faceL)
                if (cnt[face] == 1)
                    newNode->lc->objs->push_back(face);
                else
                    newNode->objs->push_back(face);
            for (auto face : *faceR)
                if (cnt[face] == 1) newNode->rc->objs->push_back(face);
        } else
            newNode->lc = newNode->rc = nullptr;
        return newNode;
    }

    void getFaces(ObjectTreeNode* p, vector<Object*>* faces) {
        p->l = faces->size();
        for (auto face : *(p->objs)) faces->push_back(face);
        p->r = faces->size();
        if (p->lc) getFaces(p->lc, faces);
        if (p->rc) getFaces(p->rc, faces);
    }

    ObjectTreeNode* root;
    vector<Object*>* faces;
    ObjectKDTree(vector<Object*>* objs) {
        Vector3f vecmin = Vector3f(INF, INF, INF);
        Vector3f vecmax = -vecmin;
        for (auto obj : *objs) {
            vecmin = minE(vecmin, obj->vecmin());
            vecmax = maxE(vecmax, obj->vecmax());
        }
        root = build(1, 0, faces, vecmin, vecmax);
        this->faces = new vector<Object*>;
        getFaces(root, this->faces);
    }

    float interWithBox(ObjectTreeNode* p, const Ray& ray) const {
        float t = INF;
        if (!p) return t;
        BoundBox(p->vecmin, p->vecmax).intersect(ray, t);
        return t;
    }

    bool intersect(ObjectTreeNode* p, const Ray& ray, Object*& nextFace, Hit& hit) const {
        bool flag = false;
        for (int i = 0; i < p->objs->size(); ++i)
            if ((*p->objs)[i]->intersect(ray, hit, 1e-6)) {
                nextFace = (*p->objs)[i];
                flag = true;
            }
        float tl = interWithBox(p->lc, ray), tr = interWithBox(p->rc, ray);
        if (tl < tr) {
            if (hit.t <= tl) return flag;
            if (p->lc) flag |= intersect(p->lc, ray, nextFace, hit);
            if (hit.t <= tr) return flag;
            if (p->rc) flag |= intersect(p->rc, ray, nextFace, hit);
        } else {
            if (hit.t <= tr) return flag;
            if (p->rc) flag |= intersect(p->rc, ray, nextFace, hit);
            if (hit.t <= tl) return flag;
            if (p->lc) flag |= intersect(p->lc, ray, nextFace, hit);
        }
        return flag;
    }
};

// TODO: Implement Group - add data structure to store a list of Object*
class Group{

public:

    Group() {
        this->groupSize = 0;
    }

    explicit Group (int num_objects) {
        this->groupSize = num_objects;
        this->vec_of_obj.resize(num_objects);
    }

    ~Group() {
        for(auto obj : this->vec_of_obj) {
            obj->~Object();
        }
    }

    bool intersect(const Ray &r, Hit &h, float tmin) {
        Object* nextface = nullptr;
        return tree->intersect(tree->root, r, nextface, h);
    }

    void addObject(int index, Object *obj) {
        vec_of_obj[index] = obj;
    }

    int getGroupSize() {
        return this->vec_of_obj.size();
    }

    std::vector<Object*> getShineObj() {
        std::vector<Object*> shinyObj;
        for(int i = 0; i < groupSize; i++) {
            if(vec_of_obj[i]->material->emission != Vector3f::ZERO)
                shinyObj.push_back(vec_of_obj[i]);
        }
        return shinyObj;
    }
    
    int groupSize;
    std::vector<Object*> vec_of_obj;
    ObjectKDTree* tree;
};

#endif
	
