#ifndef OUTFLOW_H
#define OUTFLOW_H

#include "particles/boundaries/BoundaryCondition.h"

namespace mol_sim {
using Outflow = class Outflow;

class Outflow : public BoundaryCondition {
    void boundaryStrategy(Particle& p) override;
    void clean() override;

    public:
    Outflow(LinkedCellContainer& particles);
    ~Outflow() override = default;
};
} // namespace mol_sim

#endif