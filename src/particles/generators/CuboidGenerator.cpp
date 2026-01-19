#include "particles/generators/CuboidGenerator.h"

#include <algorithm>

#include "utils/MaxwellBoltzmannDistribution.h"

namespace mol_sim {

void CuboidGenerator::generateParticles(ContainerRef particles, bool use_brownian_motion, bool use_init_temp) {
    const double avg_v = use_init_temp ? sqrt(init_temp / mass) : avg_velo;
    particles.reserve(num_particles[0] * num_particles[1] * num_particles[2]);
    const size_t dimensions = num_particles[2] > 1 ? 3 : 2;
    for (size_t i = 0; i < num_particles[2]; i++) {
        for (size_t j = 0; j < num_particles[1]; j++) {
            for (size_t k = 0; k < num_particles[0]; k++) {
                R3 curr_pos = {position[0] + (static_cast<double>(k) * distance),
                               position[1] + (static_cast<double>(j) * distance),
                               position[2] + (static_cast<double>(i) * distance)};
                R3 velo = {0., 0., 0.};
                if (use_brownian_motion) {
                    velo = maxwellBoltzmannDistributedVelocity(avg_v, dimensions);
                }
                const N3 curr = {k, j, i};
                int type = (std::ranges::find(targets, curr) != targets.end()) ? 3 : 0;
                particles.addParticle(curr_pos, velocity + velo, mass, epsilon, sigma, type);
            }
        }
    }
}
}  // namespace mol_sim
