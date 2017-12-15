
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

#include "integrators/ptaas.h"
#include "interaction.h"
#include "camera.h"
#include "film.h"
#include "paramset.h"

namespace pbrt {




PtaasIntegrator *CreatePtaasIntegrator(
    const ParamSet &params, std::shared_ptr<Sampler> sampler,
    std::shared_ptr<const Camera> camera, std::shared_ptr<const Camera> cam2) {

    return new PtaasIntegrator(params, sampler, camera, cam2);
}

}  // namespace pbrt
