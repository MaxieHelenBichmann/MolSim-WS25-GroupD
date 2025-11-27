#ifndef OUTFLOW_H
#define OUTFLOW_H

#include "particles/boundaries/BoundaryCondition.h"

namespace mol_sim {

class Outflow : public BoundaryCondition {
    public:
    Outflow(BoundaryLocation location);
    ~Outflow() override = default;
    void boundaryStrategy(Particle& p) override;
};
}  // namespace mol_sim

#endif