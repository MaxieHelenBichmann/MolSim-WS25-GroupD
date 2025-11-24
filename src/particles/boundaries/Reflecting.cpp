#include "particles/boundaries/Reflecting.h"

namespace mol_sim {

template<ParticleContainer containerType>
Reflecting<containerType>::Reflecting(containerType& particles, std::optional<double> counter_sigma, std::optional<double> counter_epsilon) 
: BoundaryCondition(particles), counter_epsilon(counter_epsilon), counter_sigma(counter_sigma) {
    if (!dimension.has_value()) {
        SPDLOG_ERROR("Reflecting boundary condition specified but dimension for domain not provided! Defaulting dimension to {1.0, 1.0, 1.0}!");
    }
    max = dimension.value_or({1., 1., 1.});
}

template<ParticleContainer containerType>
void Reflecting<containerType>::addCounterParticle(int sign, size_t coordinate, containerType& particles, Particle& p) {
    double sigma = counter_sigma.value_or(p.getSigma());  
    double epsilon = counter_epsilon.value_or(p.getEpsilon());
    R3 check = zero;
    check[coordinate] = sign * pow(2, 1.0/6.0) * sigma;
    if (!particles.fitsContainer(p.getX() + check)) {
        R3 pos_counter_particle = sign < 0 ? zero : max;
        for (size_t i = 0; i < 3; i++) {
            if (i == coordinate) { continue; }
            pos_counter_particle[i] = p.getX()[i];
        }
        Particle ghost(pos_counter_particle, zero, .0, epsilon, sigma, -1);
        particles.addParticle(ghost);
        ghost_particles.push_back(&ghost);
    } 
}

template<ParticleContainer containerType>
void Reflecting<containerType>::boundaryStrategy(Particle& p) {
    addCounterParticle( 1, 0, particles, p); 
    addCounterParticle(-1, 0, particles, p); 
    addCounterParticle( 1, 1, particles, p); 
    addCounterParticle(-1, 1, particles, p); 
    addCounterParticle( 1, 2, particles, p); 
    addCounterParticle(-1, 2, particles, p); 
}
} // namespace mol_sim
