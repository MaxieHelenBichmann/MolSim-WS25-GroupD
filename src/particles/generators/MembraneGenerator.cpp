#include "particles/generators/MembraneGenerator.h"

#include "utils/MaxwellBoltzmannDistribution.h"

namespace mol_sim {

void MembraneGenerator::generateParticles(ContainerRef particles, bool use_init_temp) {
    size_t base = particles.size();
    particles.reserve(num_particles[0] * num_particles[1]);
    for (size_t j = 0; j < num_particles[1]; j++) {
        for (size_t k = 0; k < num_particles[0]; k++) {
            R3 curr_pos = {position[0] + (static_cast<double>(k) * distance),
                           position[1] + (static_cast<double>(j) * distance), position[2]};

            if (use_init_temp) {
                avg_velo = sqrt(init_temp / mass);
            }
            N2 curr = {k, j};
            int type = (std::ranges::find(targets, curr) != targets.end()) ? 4 : 2;
            R3 velo = maxwellBoltzmannDistributedVelocity(avg_velo, 2);
            particles.addParticle(curr_pos, velocity + velo, mass, epsilon, sigma, type);
        }
    }
    const size_t row_size = num_particles[0];
    for (size_t index = 0; index < particles.size() - base; index++) {
        const size_t curr_elem = base + index;
        if (index % row_size != 0) {
            particles[curr_elem].getNeighbors()[0] = &particles[curr_elem - 1];
            particles[curr_elem - 1].getNeighbors()[1] = &particles[curr_elem];
        }
        // Check if bottom neighbor exists (previous row), if so set it and set current particle as its top neighbor
        if (index >= row_size) {
            particles[curr_elem].getNeighbors()[2] = &particles[curr_elem - row_size];
            particles[curr_elem - row_size].getNeighbors()[3] = &particles[curr_elem];
        }
        // Check if bottom-left diagonal neighbor exists, if so set it and set current particle as its top-right diagonal
        if (index >= row_size && index % row_size != 0) {
            particles[curr_elem].getNeighbors()[4] = &particles[curr_elem - row_size - 1];
            particles[curr_elem - row_size - 1].getNeighbors()[7] = &particles[curr_elem];
        }
        // Check if bottom-right diagonal neighbor exists, if so set it and set current particle as its top-left diagonal
        if (index >= num_particles[0] && (index + 1) % row_size != 0) {
            particles[curr_elem].getNeighbors()[5] = &particles[curr_elem - row_size + 1];
            particles[curr_elem - row_size + 1].getNeighbors()[6] = &particles[curr_elem];
        }
    }
}
}  // namespace mol_sim
