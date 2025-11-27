#include "particles/boundaries/Reflecting.h"

#include <spdlog/spdlog.h>

namespace mol_sim {

Reflecting::Reflecting(BoundaryLocation location, R3 dimension, std::optional<double> counter_sigma, std::optional<double> counter_epsilon) 
: BoundaryCondition(location, BoundaryType::REFLECTING), counter_epsilon(counter_epsilon), counter_sigma(counter_sigma), max(dimension) {}

Reflecting::Reflecting(BoundaryLocation location, R3 dimension, ContainerRef particles, std::optional<double> counter_sigma, std::optional<double> counter_epsilon) 
: BoundaryCondition(location, BoundaryType::REFLECTING), particles(particles), counter_epsilon(counter_epsilon), counter_sigma(counter_sigma), max(dimension) {}

Reflecting::~Reflecting() = default;

bool Reflecting::fitsDomain(R3 v) {
    return (v[0] > 0 || v[1] > 0 || v[2] > 0) || 
           (v[0] < max[0] || v[1] < max[1] || v[2] < max[2]);
}

void Reflecting::addCounterParticle(int sign, size_t coordinate, Particle& p) {
    double sigma = counter_sigma.value_or(p.getSigma());
    double epsilon = counter_epsilon.value_or(p.getEpsilon());
    R3 check = zero;
    check[coordinate] = sign * pow(2, 1.0 / 6.0) * sigma;
    check[coordinate] = sign * pow(2, 1.0 / 6.0) * sigma;
    if (!fitsDomain(p.getX() + check)) {
        R3 pos_counter_particle = sign < 0 ? zero : max;
        for (size_t i = 0; i < 3; i++) {
            if (i == coordinate) { 
                continue; 
            }
            pos_counter_particle[i] = p.getX()[i];
        }
        Particle ghost(pos_counter_particle, zero, .0, epsilon, sigma, -1);
        particles.addParticle(ghost);
    }
}

void Reflecting::setParticles(ContainerRef particles) { this->particles = particles; }

void Reflecting::boundaryStrategy(Particle& p) {
    switch (location) {
        case BoundaryLocation::LEFT:
            addCounterParticle(-1, 0, p);
            break;
        case BoundaryLocation::RIGHT: 
            addCounterParticle(1, 0, p); 
            break;
        case BoundaryLocation::FRONT:
            addCounterParticle(-1, 2, p); 
            break;
        case BoundaryLocation::BACK:
            addCounterParticle(1, 2, p); 
            break;
        case BoundaryLocation::UPPER:
            addCounterParticle(1, 1, p); 
            break;
        case BoundaryLocation::LOWER:
            addCounterParticle(-1, 1, p); 
            break;
        default:
            SPDLOG_ERROR("Unknown boundary type! Expected (LEFT, RIGHT, UPPER, LOWER, FRONT, BACK)!");
            exit(-1); 
    }
}
}  // namespace mol_sim
