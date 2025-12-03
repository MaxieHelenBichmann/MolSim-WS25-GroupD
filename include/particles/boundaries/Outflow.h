#ifndef OUTFLOW_H
#define OUTFLOW_H

#include "particles/boundaries/Boundary.h"

namespace mol_sim {

class Outflow : public Boundary {
   public:
    Outflow(BoundaryLocation location, R3 domain_size) noexcept
        : Boundary(location, BoundaryType::OUTFLOW, domain_size) {}
    ~Outflow() override = default;

    void applyBoundary([[maybe_unused]] Particle& p,
                       [[maybe_unused]] const ForceSource& force) const noexcept override {}
};

}  // namespace mol_sim

#endif
