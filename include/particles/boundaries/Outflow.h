#ifndef OUTFLOW_H
#define OUTFLOW_H

#include "particles/boundaries/BoundaryCondition.h"

namespace mol_sim {

template <ParticleContainer containerType>
class Outflow : public BoundaryCondition<containerType> {
    bool boundaryConditionApplies(Particle& p) override;
    void boundaryStrategy(Particle& p) override;

   public:
    Outflow(containerType* particles, BoundaryType location);
    ~Outflow() override = default;
};
}  // namespace mol_sim

#endif