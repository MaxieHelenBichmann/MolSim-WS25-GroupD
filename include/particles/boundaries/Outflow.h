#ifndef OUTFLOW_H
#define OUTFLOW_H

#include "particles/boundaries/Boundary.h"
#include "physics/ForceSource.h"

namespace mol_sim {

class Outflow : public Boundary {
   public:
    explicit Outflow(BoundaryLocation location) : Boundary(location, BoundaryType::OUTFLOW) {}
    ~Outflow() override = default;

    void applyBoundary(Particle& p, const ForceSource& force) override {
        (void)p;
        (void)force;
    }
};

}  // namespace mol_sim

#endif
