#include "particles/generators/CuboidGenerator.h"

#include <cstddef>

#include "utils/MaxwellBoltzmannDistribution.h"

namespace mol_sim {

void CuboidGenerator::generateParticles(ContainerRef particles) {
    particles.reserve(num_particles[0] * num_particles[1] * num_particles[3]);
    R3 curr_pos = position;
    for (size_t i = 0; i < num_particles[2]; i++) {
        for (size_t j = 0; j < num_particles[1]; j++) {
            for (size_t k = 0; k < num_particles[0]; k++) {
                R3 velo = maxwellBoltzmannDistributedVelocity(avg_velo, 2);
                particles.addParticle(curr_pos, velo, mass);
                curr_pos[0] += distance;
            }
            curr_pos[1] += distance;
        }
        curr_pos[2] += distance;
    }
}
}  // namespace mol_sim
