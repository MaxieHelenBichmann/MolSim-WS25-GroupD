#include "particles/boundaries/Outflow.h"

namespace mol_sim {

Outflow::Outflow(LinkedCellContainer& particles)
: BoundaryCondition(particles) {}

void Outflow::boundaryStrategy(Particle& p) {
    particles.eraseParticle(&p);
}

void Outflow::clean() {}

} // namespace mol_sim

