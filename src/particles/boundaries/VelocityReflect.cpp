#include "particles/boundaries/VelocityReflect.h"

namespace mol_sim {
VelocityReflect::VelocityReflect(BoundaryLocation location, R3 domain_size, bool angular_reflect)
    : Boundary(location, BoundaryType::VELOCITYREFLECT), domain_size(domain_size), angular_reflect(angular_reflect) {}

std::optional<Particle> VelocityReflect::applyBoundary(Particle& p) {
    if (angular_reflect) {
        size_t axis = getAxis();
        p.getV()[axis] = -1 * p.getV()[axis];
    } else {
        p.getV() = -1 * p.getV();
    }
    return std::nullopt;
}

}  // namespace mol_sim
