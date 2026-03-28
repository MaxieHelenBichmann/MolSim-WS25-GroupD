#include "particles/generators/DiscGenerator.h"

#include "utils/MaxwellBoltzmannDistribution.h"

namespace mol_sim {

void DiscGenerator::generateParticles(ContainerRef particles, bool use_brownian_motion, bool use_init_temp) {
    if (radius == 0) {
        return;
    }
    particles.reserve(radius * radius * 4);
    const long loop_radius = static_cast<long>(radius);
    const double avg_v = use_init_temp ? sqrt(init_temp / mass) : avg_velo;

    for (long i = -loop_radius; i <= loop_radius; i++) {
        for (long j = -loop_radius; j <= loop_radius; j++) {
            if ((i * i) + (j * j) <= loop_radius * loop_radius) {
                R3 curr_pos = {position[0] + (static_cast<double>(i) * distance),
                               position[1] + (static_cast<double>(j) * distance), position[2]};
                R3 velo = {0., 0., 0.};
                if (use_brownian_motion) {
                    velo = maxwellBoltzmannDistributedVelocity(avg_v, 2);
                }
                particles.addParticle(curr_pos, velocity + velo, mass, epsilon, sigma);
            }
        }
    }
}
}  // namespace mol_sim
