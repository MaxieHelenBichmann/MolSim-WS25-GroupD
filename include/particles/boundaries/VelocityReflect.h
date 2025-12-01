#ifndef VELOCITYREFLECT_H
#define VELOCITYREFLECT_H

#include "particles/boundaries/Boundary.h"
#include "particles/Particle.h"
#include "utils/Vector.h"

namespace mol_sim {
class VelocityReflect : public Boundary {
    R3 domain_size;
    
    public:
    VelocityReflect(BoundaryLocation location, R3 domain_size);
    ~VelocityReflect() override = default;

    [[nodiscard]] std::optional<Particle> applyBoundary(Particle& p) override; 
};
}  // namespace mol_sim
#endif