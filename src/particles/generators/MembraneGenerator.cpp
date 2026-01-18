#include "particles/generators/MembraneGenerator.h"

#include "utils/MaxwellBoltzmannDistribution.h"

namespace mol_sim {

void MembraneGenerator::generateParticles(ContainerRef particles, bool use_init_temp) {
    const size_t base = particles.size();
    particles.reserve(base + (num_particles[0] * num_particles[1]));

    const double avg_v = use_init_temp ? sqrt(init_temp / mass) : avg_velo;
    const size_t row_size = num_particles[0];

    for (size_t j = 0; j < num_particles[1]; j++) {
        for (size_t k = 0; k < num_particles[0]; k++) {
            R3 curr_pos = {position[0] + (static_cast<double>(k) * distance),
                           position[1] + (static_cast<double>(j) * distance), position[2]};

            const N2 curr = {k, j};
            const int type = (std::ranges::find(targets, curr) != targets.end()) ? 4 : 2;
            const R3 velo = maxwellBoltzmannDistributedVelocity(avg_v, 2);
            particles.addParticle(curr_pos, velocity + velo, mass, epsilon, sigma, type);
        }
    }

    for (size_t index = 0; index < num_particles[0] * num_particles[1]; index++) {
        const size_t curr_elem = base + index;
        auto& curr_neighbors = particles[curr_elem].getNeighbors();

        if (index % row_size != 0) {
            curr_neighbors[0] = curr_elem - 1;
            particles[curr_elem - 1].getNeighbors()[1] = curr_elem;
        }
        if (index >= row_size) {
            curr_neighbors[2] = curr_elem - row_size;
            particles[curr_elem - row_size].getNeighbors()[3] = curr_elem;
        }
        if (index >= row_size && index % row_size != 0) {
            curr_neighbors[4] = curr_elem - row_size - 1;
            particles[curr_elem - row_size - 1].getNeighbors()[7] = curr_elem;
        }
        if (index >= row_size && (index + 1) % row_size != 0) {
            curr_neighbors[5] = curr_elem - row_size + 1;
            particles[curr_elem - row_size + 1].getNeighbors()[6] = curr_elem;
        }
    }
}
}  // namespace mol_sim
