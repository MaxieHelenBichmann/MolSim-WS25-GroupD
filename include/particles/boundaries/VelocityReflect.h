#ifndef VELOCITYREFLECT_H
#define VELOCITYREFLECT_H

#include "particles/Particle.h"
#include "particles/boundaries/Boundary.h"
#include "utils/Vector.h"

namespace mol_sim {
class VelocityReflect : public Boundary {
   public:
    VelocityReflect(BoundaryLocation location, R3 domain_size) noexcept;
    ~VelocityReflect() override = default;

    void applyBoundary(Particle& p, [[maybe_unused]] const ForceSource& force) const noexcept override;
};
}  // namespace mol_sim
#endif