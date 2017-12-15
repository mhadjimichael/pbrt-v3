
/*

 /$$$$$$$  /$$$$$$$$                   /$$$$$$ 
| $$__  $$|__  $$__/                  /$$__  $$
| $$  \ $$   | $$  /$$$$$$   /$$$$$$ | $$  \__/
| $$$$$$$/   | $$ |____  $$ |____  $$|  $$$$$$ 
| $$____/    | $$  /$$$$$$$  /$$$$$$$ \____  $$
| $$         | $$ /$$__  $$ /$$__  $$ /$$  \ $$
| $$         | $$|  $$$$$$$|  $$$$$$$|  $$$$$$/
|__/         |__/ \_______/ \_______/ \______/ 
                                               
                                               
*/

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_PTAAS_H
#define PBRT_INTEGRATORS_PTAAS_H

// integrators/ptaas.h*
#include "pbrt.h"
#include "integrator.h"
#include "integrators/mlt.h"
#include "integrators/path.h"
#include "scene.h"
#include "camera.h"

namespace pbrt {

// PtaasIntegrator Declarations
class PtaasIntegrator : public Integrator {
  public:
    // PtaasIntegrator Public Methods
    PtaasIntegrator(const ParamSet &params, 
            std::shared_ptr<Sampler> sampler,
            std::shared_ptr<const Camera> camera,
            std::shared_ptr<const Camera> camera2)
        {
            mlt = CreateMLTIntegrator(params, camera);
            path = CreatePathIntegrator(params, sampler, camera2);
            this->cam1 = camera;
            this->cam2 = camera2;
        }
    void Render(const Scene &scene) {

        mlt->Render(scene);
        path->Render(scene);

        LOG(ERROR) << "Bout to do some mergin'";

        cam1->film->MergeImage(0.5f, &(cam2->film->rgb[0]));
    }

  private:
    // PtaasIntegrator Private Data
    MLTIntegrator *mlt;
    PathIntegrator *path;
    std::shared_ptr<const Camera> cam1;
    std::shared_ptr<const Camera> cam2;

/*
    // MLTIntegrator Private Data
    std::shared_ptr<const Camera> camera;
    const int maxDepth;
    const int nBootstrap;
    const int nChains;
    const int mutationsPerPixel;
    const Float sigma, largeStepProbability;
    */
};

PtaasIntegrator *CreatePtaasIntegrator(
    const ParamSet &params, std::shared_ptr<Sampler> sampler,
    std::shared_ptr<const Camera> camera, std::shared_ptr<const Camera> cam2);

}  // namespace pbrt

#endif  // PBRT_INTEGRATORS_PTAAS_H
