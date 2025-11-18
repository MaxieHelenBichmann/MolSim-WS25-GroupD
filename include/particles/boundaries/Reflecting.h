#ifndef REFLECTING_H
#define REFLECTING_H

#include <cstdint>
#include "particles/boundaries/BoundaryCondition.h"

namespace mol_sim {
using Reflecting = class Reflecting;

class Reflecting : public BoundaryCondition {
    R3 zero = {.0, .0, .0};
    R3 max;  
    // for now sigma + epsilon is that of the particle. but we could also do set up 
    // like add a setting in settings that specifies the sigma and epsilon 
    double counter_epsilon = 0;
    double counter_sigma = 0;
    bool use_particle_sigma_epsilon = false;
    std::vector<Particle*> ghost_particles;

    void addCounterParticle(int sign, size_t coordinate, LinkedCellContainer& particles, Particle& p);
    void boundaryStrategy(Particle& p) override;
    void clean() override;

    public: 
    Reflecting(LinkedCellContainer& particles);
    Reflecting(LinkedCellContainer& particles, double counter_sigma, double counter_epsilon);
    ~Reflecting() override = default;
};

} // namespace mol_sim

#endif