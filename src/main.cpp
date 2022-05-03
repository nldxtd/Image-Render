#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <random>

#include "scene_parser.hpp"
#include "image.hpp"
#include "camera.hpp"
#include "group.hpp"
#include "light.hpp"
#include "hittree.hpp"
#include "utils.hpp"

using namespace std;

class SPPM {
public:
    int numRounds, numPhotons;    
    SceneParser parser;
    Camera* camera;
    Group* group;
    vector<Object*> shiny;
    vector<Hit*> hits;
    HitTree* hitTree;

    SPPM(const char* input_file, int _numRounds, int _numPhotons): parser(input_file) {
        camera = parser.getCamera();
        group = parser.getGroup();
        shiny = group->getShineObj();
        numRounds = _numRounds;
        numPhotons = _numPhotons;
        hitTree = nullptr;
        int w = camera->getWidth();
        int h = camera->getHeight();
        for(int u = 0; u < w; u++) {
            for(int v = 0; v < h; v++) {
                hits.push_back(new Hit());
            }
        }
        cout << "constructing an camera with width " << w << " and height " << h << endl;
    }

    ~SPPM() {
        int w = camera->getWidth();
        int h = camera->getHeight();
        for(int u = 0; u < w; u++) {
            for(int v = 0; v < h; v++) {
                delete hits[u*w+v];
            }
        }
        delete hitTree;
    }

    void render() {
        int w = camera->getWidth();
        int h = camera->getHeight();
        int cnt = 0;
        time_t start = time(NULL);
        for(int round = 1; round <= numRounds; round++) {
            float elapsed = (time(NULL) - start), progress = round + 0. / numRounds;
            fprintf(stderr, "\rRendering (%d/%d Rounds) %5.2f%% Time: %.2f/%.2f sec\n", round, numRounds, progress * 100., elapsed, elapsed / progress);
            //from each pixel shoot a ray
            for(int x = 0; x < w; x++) {
                for(int y = 0; y < h; y++) {
                    int index = x*h+y;
                    Ray camRay = camera->generateRay(Vector2f(x + RND, y + RND));
                    // cout << camRay.getDirection().x() << " " << camRay.getDirection().y() << " " << camRay.getDirection().z() << endl;
                    //ray started to shoot
                    int forward_depth = TRACE_DEPTH;
                    Vector3f attenWeight(1, 1, 1);
                    while(true) {
                        if(attenWeight.max() < MIN_ATTEN || --forward_depth <= 0) break;
                        hits[index]->t = INF;
                        if(!group->intersect(camRay, *hits[index], 1e-6)) {
                            hits[index]->fluxLight += hits[index]->attenWeight * parser.getBackgroundColor();
                            break;
                        }
                        cnt++;
                        camRay.origin = camRay.pointAtParameter(hits[index]->t);
                        Material* material = hits[index]->getMaterial();
                        Vector3f hitnorm = hits[index]->getNormal();
                        float type = RND2;
                        if (type <= material->type.x()) {  // Diffuse
                            hits[index]->attenWeight = attenWeight * material->diffuseColor;
                            hits[index]->fluxLight += hits[index]->attenWeight * material->emission;
                            if(x+y == 1) {
                                printf("hit material with color %f %f %f\n", material->diffuseColor.x(), material->diffuseColor.y(), material->diffuseColor.z());
                                printf("hit point with plot %f %f %f\n", hits[index]->p.x(), hits[index]->p.y(), hits[index]->p.z());
                            }
                            break;
                        } else if (type <= material->type.x() + material->type.y()) {
                            float cos = Vector3f::dot(camRay.direction, hitnorm);
                            camRay.direction = (camRay.direction - hitnorm * (cos * 2)).normalized();
                        } else {
                            float reflect = material->reflect;
                            float reflect0 = ((1.0 - reflect) * (1.0 - reflect)) / ((1.0 + reflect) * (1.0 + reflect));
                            if (Vector3f::dot(hitnorm, camRay.direction) > 0) {  // inside the medium
                                hitnorm.negate();
                                reflect = 1 / reflect;
                            }
                            reflect = 1 / reflect;
                            float cos1 = -Vector3f::dot(hitnorm, camRay.direction);
                            float cos2 = 1.0 - reflect * reflect * (1.0 - cos1 * cos1);
                            float Rprob = reflect0 + (1.0 - reflect0) * pow(1.0 - cos1, 5.0);   // Schlick-approximation
                            if (cos2 > 0 && RND2 > Rprob) {  // refraction direction
                                camRay.direction = ((camRay.direction * reflect) + (hitnorm * (reflect * cos1 - sqrt(cos2)))).normalized();
                            } else {  // reflection direction
                                camRay.direction = (camRay.direction + hitnorm * (cos1 * 2)).normalized();
                            }
                        }
                        attenWeight = attenWeight * material->diffuseColor;
                    }
                }
            }

            if(hitTree) delete hitTree;
            hitTree = new HitTree(hits);

            int shinySize = shiny.size();
            int photonMess = numPhotons / shinySize;
            for(int i = 0; i < shinySize; i++) {
                for(int j = 0; j < photonMess; j++) {
                    Ray randomRay = shiny[i]->shootRay(-1, (long long)round * numPhotons + (round + 1) * w * h + j);
                    long long seed = (long long)round * numPhotons + j;
                    int backwardDepth = TRACE_DEPTH;
                    Vector3f attenWeight = shiny[i]->material->diffuseColor * Vector3f(250, 250, 250);
                    // if(randomRay.origin.y() <= 8) {
                    //     printf("shoot ray from %f %f %f\n", randomRay.origin.x(), randomRay.origin.y(), randomRay.origin.z());
                    // }
                    while(true) {
                        if(attenWeight.max() < MIN_ATTEN || --backwardDepth <= 0) break;
                        Hit hit;
                        if(!group->intersect(randomRay, hit, 1e-6)) break;
                        randomRay.origin = randomRay.pointAtParameter(hit.t);
                        Material* material = hit.getMaterial();
                        Vector3f hitnorm = hit.getNormal();
                        float type = RND2;
                        if (type <= material->type.x()) {  // Diffuse
                            hitTree->update(hitTree->root, hit.p, attenWeight, randomRay.direction);
                            randomRay.direction = diffDir(hitnorm, -1, seed);
                        }
                        else if(type <= material->type.x() + material->type.y()) {
                            float cos = Vector3f::dot(randomRay.direction, hitnorm);
                            randomRay.direction = (randomRay.direction - hitnorm * (cos * 2)).normalized();
                        }
                        else {
                            float reflect = material->reflect;
                            float reflect0 = ((1.0 - reflect) * (1.0 - reflect)) / ((1.0 + reflect) * (1.0 + reflect));
                            if (Vector3f::dot(hitnorm, randomRay.direction) > 0) {  // inside the medium
                                hitnorm.negate();
                                reflect = 1 / reflect;
                            }
                            reflect = 1 / reflect;
                            float cos1 = -Vector3f::dot(hitnorm, randomRay.direction);
                            float cos2 = 1.0 - reflect * reflect * (1.0 - cos1 * cos1);
                            float Rprob = reflect0 + (1.0 - reflect0) * pow(1.0 - cos1, 5.0);   // Schlick-approximation
                            if (cos2 > 0 && RND2 > Rprob) {  // refraction direction
                                randomRay.direction = ((randomRay.direction * reflect) + (hitnorm * (reflect * cos1 - sqrt(cos2)))).normalized();
                            } else {  // reflection direction
                                randomRay.direction = (randomRay.direction + hitnorm * (cos1 * 2)).normalized();
                            }
                        }
                        attenWeight = attenWeight * hit.material->diffuseColor;
                    }
                }
            }
            float elapsed2 = (time(NULL) - start);
            fprintf(stderr, "\rRenderingTime: %.2f sec\n", elapsed2);
        }
    }
};

int main(int argc, char *argv[]) {
    for (int argNum = 1; argNum < argc; ++argNum) {
        std::cout << "Argument " << argNum << " is: " << argv[argNum] << std::endl;
    }

    if (argc != 5) {
        cout << "Usage: ./bin/PA1 <input scene file> <output bmp file> <num of rounds> <num of photons>" << endl;
        return 1;
    }
    string inputFile = argv[1];
    string outputFile = argv[2];
    int _numRounds = atoi(argv[3]);
    int _numPhotons = atoi(argv[4]);

    SPPM sppm(argv[1], _numRounds, _numPhotons);
    sppm.render();
    int w = sppm.camera->getWidth();
    int h = sppm.camera->getHeight();
    Image img(w, h);
    for(int i = 0; i < w; i++) {
        for(int j = 0; j < h; j++) {
            Hit* hit = sppm.hits[i*h+j];
            img.SetPixel(i, j, hit->flux/(M_PI * hit->radius * _numPhotons * _numRounds) + hit->fluxLight / _numRounds);
        }
    }
    img.SaveBMP(argv[2]);

    cout << "Hello! Computer Graphics!" << endl;
    return 0;
}

