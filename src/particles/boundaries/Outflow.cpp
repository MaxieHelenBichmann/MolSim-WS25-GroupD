#include "particles/boundaries/Outflow.h"

namespace mol_sim {

Outflow::Outflow(BoundaryLocation location)
: BoundaryCondition(location, BoundaryType::OUTFLOW) {}

void Outflow::boundaryStrategy(Particle& p) {
    (void)p;
}
}  // namespace mol_sim

