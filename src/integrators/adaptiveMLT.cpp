
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


// integrators/adaptiveMLT.cpp*
#include "integrators/adaptiveMLT.h"
#include "integrators/bdpt.h"
#include "scene.h"
#include "film.h"
#include "sampler.h"
#include "integrator.h"
#include "camera.h"
#include "stats.h"
#include "filters/box.h"
#include "paramset.h"
#include "sampling.h"
#include "progressreporter.h"

namespace pbrt {

STAT_PERCENT("Integrator/Acceptance rate", acceptedMutations, totalMutations);

// AdaptiveMLTSampler Constants
static const int cameraStreamIndex = 0;
static const int lightStreamIndex = 1;
static const int connectionStreamIndex = 2;
static const int nSampleStreams = 3;

// AdaptiveMLTSampler Method Definitions
Float AdaptiveMLTSampler::Get1D() {
    ProfilePhase _(Prof::GetSample);
    int index = GetNextIndex();
    EnsureReady(index);
    return X[index].value;
}

Point2f AdaptiveMLTSampler::Get2D() { return {Get1D(), Get1D()}; }

std::unique_ptr<Sampler> AdaptiveMLTSampler::Clone(int seed) {
    LOG(FATAL) << "AdaptiveMLTSampler::Clone() is not implemented";
    return nullptr;
}

void AdaptiveMLTSampler::StartIteration() {
    currentIteration++;
    largeStep = rng.UniformFloat() < largeStepProbability;
}

void AdaptiveMLTSampler::Accept() {
    if (largeStep) lastLargeStepIteration = currentIteration;
}

void AdaptiveMLTSampler::logAcceptRatio(Float a, Spectrum contribution) {
    #define LIMITTTTT 1000
    if (iters > LIMITTTTT)  return;

    if (largeStep) {
        //num_l++;
        //n_l += a;
        if (a == 1.f) n_l++;
        else num_l++;

        num_0++;
        if (contribution == Spectrum(0.f)) {
            n_0 += 1;
        }
    }
    else {
        //num_s++;
        //n_s += a;
        if (a == 1.f) n_s++;
        else num_s++;
    }

    

    if (iters == LIMITTTTT && num_l > 0 && num_s > 0) {
        Float avg_n_l = n_l/(Float)num_l;
        Float avg_n_0 = n_0/(Float)num_0;
        Float avg_n_s = n_s/(Float)num_s;

        if (avg_n_l/avg_n_0 > 0.1f) {
        //if (avg_n_l*10.f > *avg_n_0) {
            largeStepProbability = std::min(avg_n_s/(2.f*(avg_n_s - avg_n_l)), 1.f);
            //std::cout << "Case 1\n";
        }
        else {
            largeStepProbability = 0.25f;
            //std::cout << "Case 2\n";
        }
        
        if (iters == LIMITTTTT) {
            std::cout << "Large Step Probability set to: " << largeStepProbability << "\n";
            std::cout << "avg_n_l: " << avg_n_l << "\n" <<
                        "avg_n_s: " << avg_n_s << "\n" <<
                        "avg_n_0 " << avg_n_0 << "\n" << 
                        "ratiooooo: " << avg_n_l/avg_n_0 << "\n";
        }
    }
    iters++;
}


void AdaptiveMLTSampler::EnsureReady(int index) {
    // Enlarge _AdaptiveMLTSampler::X_ if necessary and get current $\VEC{X}_i$
    if (index >= X.size()) X.resize(index + 1);
    PrimarySample &Xi = X[index];

    // Reset $\VEC{X}_i$ if a large step took place in the meantime
    if (Xi.lastModificationIteration < lastLargeStepIteration) {
        Xi.value = rng.UniformFloat();
        Xi.lastModificationIteration = lastLargeStepIteration;
    }

    // Apply remaining sequence of mutations to _sample_
    Xi.Backup();
    if (largeStep) {
        Xi.value = rng.UniformFloat();
    } else {
        int64_t nSmall = currentIteration - Xi.lastModificationIteration;
        // Apply _nSmall_ small step mutations

        // Sample the standard normal distribution $N(0, 1)$
        Float normalSample = Sqrt2 * ErfInv(2 * rng.UniformFloat() - 1);

        // Compute the effective standard deviation and apply perturbation to
        // $\VEC{X}_i$
        Float effSigma = sigma * std::sqrt((Float)nSmall);
        Xi.value += normalSample * effSigma;
        // Keep within [0,1)
        Xi.value -= std::floor(Xi.value);
    }
    Xi.lastModificationIteration = currentIteration;
}

void AdaptiveMLTSampler::Reject() {
    for (auto &Xi : X)
        if (Xi.lastModificationIteration == currentIteration) Xi.Restore();
    --currentIteration;
}

void AdaptiveMLTSampler::StartStream(int index) {
    CHECK_LT(index, streamCount);
    streamIndex = index;
    sampleIndex = 0;
}

// MLT Method Definitions
Spectrum AdaptiveMLTIntegrator::L(const Scene &scene, MemoryArena &arena,
                          const std::unique_ptr<Distribution1D> &lightDistr,
                          const std::unordered_map<const Light *, size_t> &lightToIndex,
                          AdaptiveMLTSampler &sampler, int depth, Point2f *pRaster) {
    sampler.StartStream(cameraStreamIndex);
    // Determine the number of available strategies and pick a specific one
    int s, t, nStrategies;
    if (depth == 0) {
        nStrategies = 1;
        s = 0;
        t = 2;
    } else {
        nStrategies = depth + 2;
        s = std::min((int)(sampler.Get1D() * nStrategies), nStrategies - 1);
        t = nStrategies - s;
    }

    // Generate a camera subpath with exactly _t_ vertices
    Vertex *cameraVertices = arena.Alloc<Vertex>(t);
    Bounds2f sampleBounds = (Bounds2f)camera->film->GetSampleBounds();
    *pRaster = sampleBounds.Lerp(sampler.Get2D());
    if (GenerateCameraSubpath(scene, sampler, arena, t, *camera, *pRaster,
                              cameraVertices) != t)
        return Spectrum(0.f);

    // Generate a light subpath with exactly _s_ vertices
    sampler.StartStream(lightStreamIndex);
    Vertex *lightVertices = arena.Alloc<Vertex>(s);
    if (GenerateLightSubpath(scene, sampler, arena, s, cameraVertices[0].time(),
                             *lightDistr, lightToIndex, lightVertices) != s)
        return Spectrum(0.f);

    // Execute connection strategy and return the radiance estimate
    sampler.StartStream(connectionStreamIndex);
    return ConnectBDPT(scene, lightVertices, cameraVertices, s, t, *lightDistr,
                       lightToIndex, *camera, sampler, pRaster) *
           nStrategies;
}

void AdaptiveMLTIntegrator::Render(const Scene &scene) {
    std::unique_ptr<Distribution1D> lightDistr =
        ComputeLightPowerDistribution(scene);

    // Compute a reverse mapping from light pointers to offsets into the
    // scene lights vector (and, equivalently, offsets into
    // lightDistr). Added after book text was finalized; this is critical
    // to reasonable performance with 100s+ of light sources.
    std::unordered_map<const Light *, size_t> lightToIndex;
    for (size_t i = 0; i < scene.lights.size(); ++i)
        lightToIndex[scene.lights[i].get()] = i;


    // calculate the prob
    {

        int i = 0;
        Film &film = *camera->film;
        int64_t nTotalMutations =
            (int64_t)mutationsPerPixel * (int64_t)film.GetSampleBounds().Area();

        // Generate bootstrap samples and compute normalization constant $b$
        int nBootstrapSamples = nBootstrap * (maxDepth + 1);
        std::vector<Float> bootstrapWeights(nBootstrapSamples, 0);
        if (scene.lights.size() > 0) {
            ProgressReporter progress(nBootstrap / 256,
                                    "Generating bootstrap paths");
            std::vector<MemoryArena> bootstrapThreadArenas(MaxThreadIndex());
            int chunkSize = Clamp(nBootstrap / 128, 1, 8192);
            ParallelFor([&](int i) {
                // Generate _i_th bootstrap sample
                MemoryArena &arena = bootstrapThreadArenas[ThreadIndex];
                for (int depth = 0; depth <= maxDepth; ++depth) {
                    int rngIndex = i * (maxDepth + 1) + depth;
                    AdaptiveMLTSampler sampler(mutationsPerPixel, rngIndex, sigma,
                                    largeStepProbability, nSampleStreams);
                    Point2f pRaster;
                    bootstrapWeights[rngIndex] =
                        L(scene, arena, lightDistr, lightToIndex, sampler, depth, &pRaster).y();
                    arena.Reset();
                }
                if ((i + 1) % 256 == 0) progress.Update();
            }, nBootstrap, chunkSize);
            progress.Done();
        }
        Distribution1D bootstrap(&bootstrapWeights[0], nBootstrapSamples);
        Float b = bootstrap.funcInt * (maxDepth + 1);

        int64_t nChainMutations =
            std::min((i + 1) * nTotalMutations / nChains, nTotalMutations) -
            i * nTotalMutations / nChains;
        // Follow {i}th Markov chain for _nChainMutations_
        MemoryArena arena;

        // Select initial state from the set of bootstrap samples
        RNG rng(i);
        int bootstrapIndex = bootstrap.SampleDiscrete(rng.UniformFloat());
        int depth = bootstrapIndex % (maxDepth + 1);

        // Initialize local variables for selected state
        AdaptiveMLTSampler sampler(mutationsPerPixel, bootstrapIndex, sigma,
                            largeStepProbability, nSampleStreams);
        Point2f pCurrent;
        Spectrum LCurrent =
            L(scene, arena, lightDistr, lightToIndex, sampler, depth, &pCurrent);

        // Run the Markov chain for _nChainMutations_ steps
        for (int64_t j = 0; j < nChainMutations; ++j) {
            sampler.StartIteration();
            Point2f pProposed;
            Spectrum LProposed =
                L(scene, arena, lightDistr, lightToIndex, sampler, depth, &pProposed);
            // Compute acceptance probability for proposed sample
            Float accept = std::min((Float)1, LProposed.y() / LCurrent.y());
            //sampler.logAcceptRatio(accept, LProposed.y()); //ptaas

            // Accept or reject the proposal
            if (rng.UniformFloat() < accept) {
                pCurrent = pProposed;
                LCurrent = LProposed;
                sampler.Accept();
                ++acceptedMutations;
                sampler.logAcceptRatio(1.0f, LProposed); //ptaas
            } else {
                sampler.Reject();
                sampler.logAcceptRatio(0.0f, LProposed); //ptaas
            }
            ++totalMutations;
            arena.Reset();
        }

        largeStepProbability = sampler.largeStepProbability;
    }

    // Generate bootstrap samples and compute normalization constant $b$
    int nBootstrapSamples = nBootstrap * (maxDepth + 1);
    std::vector<Float> bootstrapWeights(nBootstrapSamples, 0);
    if (scene.lights.size() > 0) {
        ProgressReporter progress(nBootstrap / 256,
                                  "Generating bootstrap paths");
        std::vector<MemoryArena> bootstrapThreadArenas(MaxThreadIndex());
        int chunkSize = Clamp(nBootstrap / 128, 1, 8192);
        ParallelFor([&](int i) {
            // Generate _i_th bootstrap sample
            MemoryArena &arena = bootstrapThreadArenas[ThreadIndex];
            for (int depth = 0; depth <= maxDepth; ++depth) {
                int rngIndex = i * (maxDepth + 1) + depth;
                AdaptiveMLTSampler sampler(mutationsPerPixel, rngIndex, sigma,
                                   largeStepProbability, nSampleStreams);
                Point2f pRaster;
                bootstrapWeights[rngIndex] =
                    L(scene, arena, lightDistr, lightToIndex, sampler, depth, &pRaster).y();
                arena.Reset();
            }
            if ((i + 1) % 256 == 0) progress.Update();
        }, nBootstrap, chunkSize);
        progress.Done();
    }
    Distribution1D bootstrap(&bootstrapWeights[0], nBootstrapSamples);
    Float b = bootstrap.funcInt * (maxDepth + 1);

    // Run _nChains_ Markov chains in parallel
    Film &film = *camera->film;
    int64_t nTotalMutations =
        (int64_t)mutationsPerPixel * (int64_t)film.GetSampleBounds().Area();
    if (scene.lights.size() > 0) {
        const int progressFrequency = 32768;
        ProgressReporter progress(nTotalMutations / progressFrequency,
                                  "Rendering");

        int64_t temp = std::min(nTotalMutations / nChains, nTotalMutations);
        std::cout << temp << "\n";
            
        ParallelFor([&](int i) {
            int64_t nChainMutations =
                std::min((i + 1) * nTotalMutations / nChains, nTotalMutations) -
                i * nTotalMutations / nChains;
            // Follow {i}th Markov chain for _nChainMutations_
            MemoryArena arena;

            // Select initial state from the set of bootstrap samples
            RNG rng(i);
            int bootstrapIndex = bootstrap.SampleDiscrete(rng.UniformFloat());
            int depth = bootstrapIndex % (maxDepth + 1);

            // Initialize local variables for selected state
            AdaptiveMLTSampler sampler(mutationsPerPixel, bootstrapIndex, sigma,
                               largeStepProbability, nSampleStreams);
            Point2f pCurrent;
            Spectrum LCurrent =
                L(scene, arena, lightDistr, lightToIndex, sampler, depth, &pCurrent);

            // Run the Markov chain for _nChainMutations_ steps
            for (int64_t j = 0; j < nChainMutations; ++j) {
                sampler.StartIteration();
                Point2f pProposed;
                Spectrum LProposed =
                    L(scene, arena, lightDistr, lightToIndex, sampler, depth, &pProposed);
                // Compute acceptance probability for proposed sample
                Float accept = std::min((Float)1, LProposed.y() / LCurrent.y());
                //sampler.logAcceptRatio(accept, LProposed.y()); //ptaas

                // Splat both current and proposed samples to _film_
                if (accept > 0)
                    film.AddSplat(pProposed,
                                  LProposed * accept / LProposed.y());
                film.AddSplat(pCurrent, LCurrent * (1 - accept) / LCurrent.y());

                // Accept or reject the proposal
                if (rng.UniformFloat() < accept) {
                    pCurrent = pProposed;
                    LCurrent = LProposed;
                    sampler.Accept();
                    ++acceptedMutations;
                    //sampler.logAcceptRatio(1.f, LProposed); //ptaas
                } else {
                    sampler.Reject();
                    //sampler.logAcceptRatio(0.0f, LProposed); //ptaas
                }
                ++totalMutations;
                if ((i * nTotalMutations / nChains + j) % progressFrequency ==
                    0)
                    progress.Update();
                arena.Reset();
            }
        }, nChains);
        progress.Done();
    }

    // Store final image computed with MLT
    camera->film->WriteImage(b / mutationsPerPixel);
}

AdaptiveMLTIntegrator *CreateAdaptiveMLTIntegrator(const ParamSet &params,
                                   std::shared_ptr<const Camera> camera) {
    int maxDepth = params.FindOneInt("maxdepth", 5);
    int nBootstrap = params.FindOneInt("bootstrapsamples", 100000);
    int64_t nChains = params.FindOneInt("chains", 1000);
    int mutationsPerPixel = params.FindOneInt("mutationsperpixel", 100);
    Float largeStepProbability =
        params.FindOneFloat("largestepprobability", 0.5f);
    Float sigma = params.FindOneFloat("sigma", .01f);
    if (PbrtOptions.quickRender) {
        mutationsPerPixel = std::max(1, mutationsPerPixel / 16);
        nBootstrap = std::max(1, nBootstrap / 16);
    }
    return new AdaptiveMLTIntegrator(camera, maxDepth, nBootstrap, nChains,
                             mutationsPerPixel, sigma, largeStepProbability);
}

}  // namespace pbrt
