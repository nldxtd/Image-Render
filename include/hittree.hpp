#ifndef HIT_TREE_H
#define HIT_TREE_H

//need rewrite
#define ALPHA 0.7
#include <algorithm>
#include <iostream>
#include <vector>

#include "hit.hpp"

#define minf(a, b) (((a) < (b)) ? (a) : (b))
#define maxf(a, b) (((a) > (b)) ? (a) : (b))

using std::vector;

class HitTreeNode {
public:
    Hit* hit;
    HitTreeNode* lc, *rc;
    Vector3f vecmax, vecmin;
    float maxradius;
};

class HitTree {
public:
    int numHits;
    Hit** hits;
    HitTreeNode* root;

    static bool compByX(Hit* a, Hit* b) { return a->p.x() < b->p.x(); }
    static bool compByY(Hit* a, Hit* b) { return a->p.y() < b->p.y(); }
    static bool compByZ(Hit* a, Hit* b) { return a->p.z() < b->p.z(); }

    HitTree(const vector<Hit *> &hitpoints) {
        numHits = hitpoints.size();
        hits = new Hit *[numHits]; 
        for(int i = 0; i < numHits; i++)
            hits[i] = hitpoints[i];
        root = build(0, numHits, 0);
        std::cout << "KDTree of hits point built succeed\n";
    }

    ~HitTree() {
        if(!root) return;
        del(root);
        delete [] hits;
    }

    void update(HitTreeNode* node, const Vector3f &photon, const Vector3f &attenWeight, const Vector3f &d) {
        if(!node) return;
        float mind = 0;
        if (photon.x() > node->vecmax.x()) mind += (photon.x() - node->vecmax.x()) * (photon.x() - node->vecmax.x());
        if (photon.x() < node->vecmin.x()) mind += (node->vecmin.x() - photon.x()) * (node->vecmin.x() - photon.x());
        if (photon.y() > node->vecmax.y()) mind += (photon.y() - node->vecmax.y()) * (photon.y() - node->vecmax.y());
        if (photon.y() < node->vecmin.y()) mind += (node->vecmin.y() - photon.y()) * (node->vecmin.y() - photon.y());
        if (photon.z() > node->vecmax.z()) mind += (photon.z() - node->vecmax.z()) * (photon.z() - node->vecmax.z());
        if (photon.z() < node->vecmin.z()) mind += (node->vecmin.z() - photon.z()) * (node->vecmin.z() - photon.z());
        if (mind > node->maxradius) return;
        if ((photon - node->hit->p).squaredLength() <= node->hit->radius) {
            Hit *hp = node->hit;
            float factor = (hp->n * ALPHA + ALPHA) / (hp->n * ALPHA + 1.);
            Vector3f dr = d - hp->normal * (2 * Vector3f::dot(d, hp->normal));
            hp->n++;
            hp->radius *= factor;
            hp->flux = (hp->flux + hp->attenWeight * attenWeight) * factor;
        }
        if (node->lc) update(node->lc, photon, attenWeight, d);
        if (node->rc) update(node->rc, photon, attenWeight, d);
        node->maxradius = node->hit->radius;
        if (node->lc && node->lc->hit->radius > node->maxradius) node->maxradius = node->lc->hit->radius;
        if (node->rc && node->rc->hit->radius > node->maxradius) node->maxradius = node->rc->hit->radius;
    }

    void del(HitTreeNode* ptr) {
        if(ptr->lc) del(ptr->lc);
        if(ptr->rc) del(ptr->rc);
        delete ptr;
    }

    HitTreeNode* build(int l, int r, int dim) {
        //left close right open, between[l, r)
        if(l >= r) { 
            return nullptr;
        }
        HitTreeNode* newNode = new HitTreeNode;
        newNode->vecmin = Vector3f(INF);
        newNode->vecmax = -newNode->vecmin;
        newNode->maxradius = 0;
        for(int i = l; i < r; i++) {
            newNode->vecmin = minE(newNode->vecmin, hits[i]->p);
            newNode->vecmax = maxE(newNode->vecmax, hits[i]->p);
            newNode->maxradius = std::max(newNode->maxradius, hits[i]->radius);
        }
        int m = (l+r)>>1;
        if(dim == 0) std::nth_element(hits+l, hits+m, hits+r, compByX);
        else if(dim == 1) std::nth_element(hits+l, hits+m, hits+r, compByY);
        else std::nth_element(hits+l, hits+m, hits+r, compByZ);
        newNode->hit = hits[m];
        newNode->lc = build(l, m, (dim+1)%3);
        newNode->rc = build(m+1, r, (dim+1)%3);

        return newNode;
    }

    Vector3f minE(const Vector3f& v1, const Vector3f& v2) {
        return Vector3f(minf(v1.x(), v2.x()), minf(v1.y(), v2.y()),
                    minf(v1.z(), v2.z()));
    }
    Vector3f maxE(const Vector3f& v1, const Vector3f& v2) {
        return Vector3f(maxf(v1.x(), v2.x()), maxf(v1.y(), v2.y()),
                        maxf(v1.z(), v2.z()));
    }
};

#endif
