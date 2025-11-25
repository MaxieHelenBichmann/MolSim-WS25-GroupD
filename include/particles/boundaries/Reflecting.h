#ifndef REFLECTING_H
#define REFLECTING_H

#include "particles/boundaries/BoundaryCondition.h"
#include "utils/Vector.h"

namespace mol_sim {

template <ParticleContainer containerType>
class Reflecting : public BoundaryCondition<containerType> {
    R3 zero = {.0, .0, .0};
    R3 max;
    std::optional<double> counter_epsilon = std::nullopt;
    std::optional<double> counter_sigma = std::nullopt;
    std::vector<Particle*> ghost_particles;

    void addCounterParticle(int sign, size_t coordinate, containerType& particles, Particle& p);
    bool boundaryConditionApplies(Particle& p) override;
    void boundaryStrategy(Particle& p) override;

   public: 
    Reflecting(containerType& particles, R3 dimension, BoundaryLocation location, std::optional<double> counter_sigma, std::optional<double> counter_epsilon);
    Reflecting(containerType* particles, R3 dimension, BoundaryLocation location, std::optional<double> counter_sigma, std::optional<double> counter_epsilon);
    ~Reflecting() override = default;
};

}  // namespace mol_sim

#endif