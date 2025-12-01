#include "particles/boundaries/VelocityReflect.h"

namespace mol_sim {
VelocityReflect::VelocityReflect(BoundaryLocation location, R3 domain_size)
    : Boundary(location, BoundaryType::VELOCITYREFLECT), domain_size(domain_size) {}

std::optional<Particle> VelocityReflect::applyBoundary(Particle& p) {
    size_t axis = getAxis();
    if (p.getX()[axis] < 0) {
            p.getX()[axis] *= -1;
        } else if (p.getX()[axis] > domain_size[axis]) {
            p.getX()[axis] -= 2 * p.getX()[axis];
        }
    return std::nullopt;
}

}  // namespace mol_sim
