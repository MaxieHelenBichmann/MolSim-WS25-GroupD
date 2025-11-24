#include "particles/boundaries/Outflow.h"

namespace mol_sim {

template <ParticleContainer containerType>
Outflow<containerType>::Outflow(containerType* particles, BoundaryType location) {
    this->particles = particles;
    this->location = location;
}

template <ParticleContainer containerType>
bool Outflow<containerType>::boundaryConditionApplies(Particle& p) {
    (void)p;
    return false;
}

template <ParticleContainer containerType>
void Outflow<containerType>::boundaryStrategy(Particle& p) {
    (void)p;
}
}  // namespace mol_sim
