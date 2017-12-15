
/*

 /$$$$$$$  /$$$$$$$$                   /$$$$$$ 
| $$__  $$|__  $$__/                  /$$__  $$
| $$  \ $$   | $$  /$$$$$$   /$$$$$$ | $$  \__/
| $$$$$$$/   | $$ |____  $$ |____  $$|  $$$$$$ 
| $$____/    | $$  /$$$$$$$  /$$$$$$$ \____  $$
| $$         | $$ /$$__  $$ /$$__  $$ /$$  \ $$
| $$         | $$|  $$$$$$$|  $$$$$$$|  $$$$$$/
|__/         |__/ \_______/ \_______/ \______/ 
                                               
Based on original MLT.

*/

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_ADAPTIVE_MLT_H
#define PBRT_INTEGRATORS_ADAPTIVE_MLT_H

// integrators/mlt.h*
#include "pbrt.h"
#include "integrator.h"
#include "sampler.h"
#include "spectrum.h"
#include "film.h"
#include "rng.h"
#include <unordered_map>

namespace pbrt {

// AdaptiveMLTSampler Declarations
class AdaptiveMLTSampler : public Sampler {
  public:
    // AdaptiveMLTSampler Public Methods
    AdaptiveMLTSampler(int mutationsPerPixel, int rngSequenceIndex, Float sigma,
               Float largeStepProbability, int streamCount)
        : Sampler(mutationsPerPixel),
          rng(rngSequenceIndex),
          sigma(sigma),
          largeStepProbability(largeStepProbability),
          streamCount(streamCount) {}
    Float Get1D();
    Point2f Get2D();
    std::unique_ptr<Sampler> Clone(int seed);
    void StartIteration();
    void Accept();
    void Reject();
    void StartStream(int index);
    int GetNextIndex() { return streamIndex + streamCount * sampleIndex++; }

    // PTaaS
    void logAcceptRatio(Float a, Spectrum b);

    Float n_s = 0; // the average probability that a small perturbation is accepted
    Float n_l = 0; // the average probability that a large perturbation is accepted
    Float n_0 = 0; // s the average probability that a large perturbation generates a non-zero contribution sample

    // num
    int num_s = 0;
    int num_l = 0;
    int num_0 = 0;
    int iters = 0;

    Float largeStepProbability;


  protected:
    // AdaptiveMLTSampler Private Declarations
    struct PrimarySample {
        Float value = 0;
        // PrimarySample Public Methods
        void Backup() {
            valueBackup = value;
            modifyBackup = lastModificationIteration;
        }
        void Restore() {
            value = valueBackup;
            lastModificationIteration = modifyBackup;
        }

        // PrimarySample Public Data
        int64_t lastModificationIteration = 0;
        Float valueBackup = 0;
        int64_t modifyBackup = 0;
    };

    // AdaptiveMLTSampler Private Methods
    void EnsureReady(int index);

    // AdaptiveMLTSampler Private Data
    RNG rng;
    const Float sigma;
    //Float largeStepProbability;
    const int streamCount;
    std::vector<PrimarySample> X;
    int64_t currentIteration = 0;
    bool largeStep = true;
    int64_t lastLargeStepIteration = 0;
    int streamIndex, sampleIndex;
};

// AdaptiveMLTIntegrator Declarations
class AdaptiveMLTIntegrator : public Integrator {
  public:
    // AdaptiveMLTIntegrator Public Methods
    AdaptiveMLTIntegrator(std::shared_ptr<const Camera> camera, int maxDepth,
                  int nBootstrap, int nChains, int mutationsPerPixel,
                  Float sigma, Float largeStepProbability)
        : camera(camera),
          maxDepth(maxDepth),
          nBootstrap(nBootstrap),
          nChains(nChains),
          mutationsPerPixel(mutationsPerPixel),
          sigma(sigma),
          largeStepProbability(largeStepProbability) {}
    void Render(const Scene &scene);
    Spectrum L(const Scene &scene, MemoryArena &arena,
               const std::unique_ptr<Distribution1D> &lightDistr,
               const std::unordered_map<const Light *, size_t> &lightToIndex,
               AdaptiveMLTSampler &sampler, int k, Point2f *pRaster);

  private:
    // AdaptiveMLTIntegrator Private Data
    std::shared_ptr<const Camera> camera;
    const int maxDepth;
    const int nBootstrap;
    const int nChains;
    const int mutationsPerPixel;
    const Float sigma;
    Float largeStepProbability;
};

AdaptiveMLTIntegrator *CreateAdaptiveMLTIntegrator(const ParamSet &params,
                                   std::shared_ptr<const Camera> camera);

}  // namespace pbrt

#endif  // PBRT_INTEGRATORS_ADAPTIVE_MLT_H
