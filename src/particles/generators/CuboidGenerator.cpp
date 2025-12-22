#include "particles/generators/CuboidGenerator.h"

#include "utils/MaxwellBoltzmannDistribution.h"

namespace mol_sim {

void CuboidGenerator::generateParticles(ContainerRef particles, bool use_init_temp) {
    particles.reserve(num_particles[0] * num_particles[1] * num_particles[2]);
    for (size_t i = 0; i < num_particles[2]; i++) {
        for (size_t j = 0; j < num_particles[1]; j++) {
            for (size_t k = 0; k < num_particles[0]; k++) {
                R3 curr_pos = {position[0] + (static_cast<double>(k) * distance),
                               position[1] + (static_cast<double>(j) * distance),
                               position[2] + (static_cast<double>(i) * distance)};
                // TODO: Look at this more carefully, what do they mean with optional brownian.
                if (use_init_temp) {
                    avg_velo = sqrt(init_temp / mass);
                }

                R3 velo = maxwellBoltzmannDistributedVelocity(avg_velo, 2);
                particles.addParticle(curr_pos, velocity + velo, mass, epsilon, sigma);
            }
        }
    }
}
}  // namespace mol_sim
