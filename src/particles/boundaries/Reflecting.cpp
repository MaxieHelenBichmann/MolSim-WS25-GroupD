#include "particles/boundaries/Reflecting.h"

namespace mol_sim {

Reflecting::Reflecting(LinkedCellContainer& particles)
: BoundaryCondition(particles) {
    max = particles.getDomainSize();
    use_particle_sigma_epsilon = true;
}

Reflecting::Reflecting(LinkedCellContainer& particles, double counter_sigma, double counter_epsilon) 
: BoundaryCondition(particles), counter_epsilon(counter_epsilon), counter_sigma(counter_sigma) {
    max = particles.getDomainSize();
    use_particle_sigma_epsilon = false;
}

void Reflecting::addCounterParticle(int sign, size_t coordinate, LinkedCellContainer& particles, Particle& p) {
    if (use_particle_sigma_epsilon) {
        counter_sigma = p.getSigma();
        counter_epsilon = p.getEpsilon();
    }
    R3 check = zero;
    check[coordinate] = sign * pow(2, 1.0/6.0) * counter_sigma;
    if (!particles.fitsContainer(p.getX() + check)) {
        R3 pos_counter_particle = sign < 0 ? zero : max;
        for (size_t i = 0; i < 3; i++) {
            if (i == coordinate) { continue; }
            pos_counter_particle[i] = p.getX()[i];
        }
        Particle ghost(pos_counter_particle, zero, .0, counter_epsilon, counter_sigma, -1);
        particles.addParticle(ghost);
        ghost_particles.push_back(&ghost);
    } 
}

void Reflecting::clean() {
    ghost_particles.clear();
}

void Reflecting::boundaryStrategy(Particle& p) {
    addCounterParticle( 1, 0, particles, p); 
    addCounterParticle(-1, 0, particles, p); 
    addCounterParticle( 1, 1, particles, p); 
    addCounterParticle(-1, 1, particles, p); 
    addCounterParticle( 1, 2, particles, p); 
    addCounterParticle(-1, 2, particles, p); 
}
} // namespace mol_sim
