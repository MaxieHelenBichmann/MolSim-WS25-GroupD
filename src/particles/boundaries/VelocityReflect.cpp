#include "particles/boundaries/VelocityReflect.h"

namespace mol_sim {
VelocityReflect::VelocityReflect(BoundaryLocation location, R3 domain_size) noexcept
    : Boundary(location, BoundaryType::VELOCITYREFLECT, domain_size) {}

std::optional<std::vector<Particle>> VelocityReflect::applyBoundary(
    Particle& p, [[maybe_unused]] const PairwiseForceSource& force) noexcept {
    size_t axis = getAxis();
    int sign = getSign();
    if (sign < 0 && p.getX()[axis] < 0) {
        p.getX()[axis] *= -1;
        p.getV()[axis] *= -1;
    } else if (sign > 0 && p.getX()[axis] > domain_size[axis]) {
        p.getX()[axis] -= 2 * (p.getX()[axis] - domain_size[axis]);
        p.getV()[axis] *= -1;
    }

    return std::nullopt;
}

}  // namespace mol_sim
