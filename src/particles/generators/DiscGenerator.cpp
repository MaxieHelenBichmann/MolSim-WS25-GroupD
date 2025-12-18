#include "particles/generators/DiscGenerator.h"

#include "utils/MaxwellBoltzmannDistribution.h"

namespace mol_sim {

void DiscGenerator::generateParticles(ContainerRef particles) {
    if (radius == 0) {
        return;
    }
    // TODO Smarter reserve of particles (this is a bad upper bound)
    particles.reserve(radius * radius * 4);
    long loop_radius = static_cast<long>(radius);

    for (long i = -loop_radius; i <= loop_radius; i++) {
        for (long j = -loop_radius; j <= loop_radius; j++) {
            if ((i * i) + (j * j) <= loop_radius * loop_radius) {
                R3 curr_pos = {position[0] + (static_cast<double>(i) * distance),
                               position[1] + (static_cast<double>(j) * distance), position[2]};
                R3 velo = maxwellBoltzmannDistributedVelocity(avg_velo, 2);
                particles.addParticle(curr_pos, velocity + velo, mass, epsilon, sigma);
            }
        }
    }
}
}  // namespace mol_sim
