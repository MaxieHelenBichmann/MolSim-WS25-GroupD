#include "particles/generators/CuboidGenerator.h"

#include "utils/MaxwellBoltzmannDistribution.h"

namespace mol_sim {

void CuboidGenerator::generateParticles(ContainerRef particles) {
    particles.reserve(num_particles[0] * num_particles[1] * num_particles[2]);
    for (size_t i = 0; i < num_particles[2]; i++) {
        for (size_t j = 0; j < num_particles[1]; j++) {
            for (size_t k = 0; k < num_particles[0]; k++) {
                R3 curr_pos = {position[0] + (k * distance), position[1] + (j * distance),
                               position[2] + (i * distance)};
                R3 velo = maxwellBoltzmannDistributedVelocity(avg_velo, 2);
                particles.addParticle(curr_pos, velocity + velo, mass);
            }
        }
    }
}
}  // namespace mol_sim
