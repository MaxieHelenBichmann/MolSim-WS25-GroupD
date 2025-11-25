#include "particles/boundaries/Reflecting.h"

#include <spdlog/spdlog.h>

namespace mol_sim {

template <ParticleContainer containerType>
Reflecting<containerType>::Reflecting(containerType& particles, R3 dimension, BoundaryLocation location, std::optional<double> counter_sigma, std::optional<double> counter_epsilon) 
: BoundaryCondition<containerType>(particles, location), counter_epsilon(counter_epsilon), counter_sigma(counter_sigma), max(dimension) {}

template <ParticleContainer containerType>
Reflecting<containerType>::Reflecting(containerType* particles, R3 dimension, BoundaryLocation location, std::optional<double> counter_sigma, std::optional<double> counter_epsilon) 
: BoundaryCondition<containerType>(*particles, location), Reflecting<containerType>(*particles, dimension, location, counter_sigma, counter_epsilon) {}

template <ParticleContainer containerType>
void Reflecting<containerType>::addCounterParticle(int sign, size_t coordinate, containerType& particles, Particle& p) {
    double sigma = counter_sigma.value_or(p.getSigma());
    double epsilon = counter_epsilon.value_or(p.getEpsilon());
    R3 check = zero;
    check[coordinate] = sign * pow(2, 1.0 / 6.0) * sigma;
    check[coordinate] = sign * pow(2, 1.0 / 6.0) * sigma;
    if (!particles.fitsContainer(p.getX() + check)) {
        R3 pos_counter_particle = sign < 0 ? zero : max;
        for (size_t i = 0; i < 3; i++) {
            if (i == coordinate) { 
                continue; 
            }
            pos_counter_particle[i] = p.getX()[i];
        }
        Particle ghost(pos_counter_particle, zero, .0, epsilon, sigma, -1);
        particles.addParticle(ghost);
        ghost_particles.push_back(&ghost);
    }
}

template <ParticleContainer containerType>
void Reflecting<containerType>::boundaryStrategy(Particle& p) {
    addCounterParticle(1, 0, this->particles, p); 
    addCounterParticle(-1, 0, this->particles, p); 
    addCounterParticle(1, 1, this->particles, p); 
    addCounterParticle(-1, 1, this->particles, p); 
    addCounterParticle(1, 2, this->particles, p); 
    addCounterParticle(-1, 2, this->particles, p); 
}
}  // namespace mol_sim
