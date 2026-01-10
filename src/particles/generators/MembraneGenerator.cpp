#include <vector>

#include "particles/generators/CuboidGenerator.h"
#include "utils/MaxwellBoltzmannDistribution.h"

namespace mol_sim {

void CuboidGenerator::generateParticles(ContainerRef particles, bool use_init_temp) {
    std::vector<Particle> particles_vec;
    size_t index = 0;
    particles_vec.reserve(num_particles[0] * num_particles[1]);
    particles.reserve(num_particles[0] * num_particles[1]);
    for (size_t j = 0; j < num_particles[1]; j++) {
        for (size_t k = 0; k < num_particles[0]; k++) {
            R3 curr_pos = {position[0] + (static_cast<double>(k) * distance),
                           position[1] + (static_cast<double>(j) * distance), position[2]};

            if (use_init_temp) {
                avg_velo = sqrt(init_temp / mass);
            }

            R3 velo = maxwellBoltzmannDistributedVelocity(avg_velo, 2);
            Particle p = {curr_pos, velocity + velo, mass, epsilon, sigma, 2};
            // Check if left neighbor exists, if so set it and set current particle as its right neighbor
            if (index > 0) {
                p.getNeighbors()[0] = &particles_vec[index - 1];
                particles_vec[index - 1].getNeighbors()[1] = &p;
            }
            // Check if top neighbor exists, if so set it and set current particle as its bottom neighbor
            if (index > num_particles[0]) {
                p.getNeighbors()[3] = &particles_vec[index - num_particles[0]];
                particles_vec[index - num_particles[0]].getNeighbors()[2] = &p;
            }
            // Check if top left neighbor exist and we are not on border, if so set it and set current particle as its
            // bottom right neighbor
            if (index > num_particles[0] && index % num_particles[0] != 0) {
                p.getNeighbors()[6] = &particles_vec[index - num_particles[0] - 1];
                particles_vec[index - num_particles[0] - 1].getNeighbors()[5] = &p;
            }
            // Check if top right neighbor exist and we are not on border, if so set it and set current particle as its
            // bottom right neighbor
            if (index > num_particles[0] && (index + 1) % num_particles[0] != 0) {
                p.getNeighbors()[7] = &particles_vec[index - num_particles[0] + 1];
                particles_vec[index - num_particles[0] + 1].getNeighbors()[4] = &p;
            }

            particles_vec.push_back(p);
            index++;
        }
    }
    for (const auto& p : particles_vec) {
        particles.addParticle(p);
    }
}
}  // namespace mol_sim
