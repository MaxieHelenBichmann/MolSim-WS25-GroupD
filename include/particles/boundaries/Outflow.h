#ifndef OUTFLOW_H
#define OUTFLOW_H

#include "particles/boundaries/Boundary.h"

namespace mol_sim {

class Outflow : public Boundary {
   public:
    Outflow(BoundaryLocation location, R3 domain_size) noexcept
        : Boundary(location, BoundaryType::OUTFLOW, domain_size) {}
    ~Outflow() override = default;

    /**
     * @brief The applyBoundary routine of Outflow. It doesn't do anything since any particle that lies on or beyond
     * the boundary of the simulation gets deleted automatically in the simulation pipeline.
     *
     * @param p The particle the boundary is to be applied on (in this case irrelevant).
     * @param force The force source that should be used in the boundary condition (in this case irrelevant).
     */
    void applyBoundary([[maybe_unused]] Particle& p,
                       [[maybe_unused]] const ForceSource& force) const noexcept override {}
};

}  // namespace mol_sim

#endif
