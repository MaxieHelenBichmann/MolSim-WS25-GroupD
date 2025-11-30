#include "particles/boundaries/VelocityReflect.h"

#include <spdlog/spdlog.h>

namespace mol_sim {
VelocityReflect::VelocityReflect(BoundaryLocation location, R3 domain_size, bool angular_reflect) 
    : Boundary(location, BoundaryType::VELOCITYREFLECT),
      domain_size(domain_size), 
      angular_reflect(angular_reflect) {}

std::optional<Particle> VelocityReflect::applyBoundary(Particle& p) {
    if (angular_reflect) {
        size_t axis = 0;
        switch (location) {
            case BoundaryLocation::LEFT:
            case BoundaryLocation::RIGHT:
                axis = 0;
                break;
            case BoundaryLocation::UPPER:
            case BoundaryLocation::LOWER:
                axis = 1;
                break;
            case BoundaryLocation::FRONT:
            case BoundaryLocation::BACK:
                axis = 2;
                break;
            default:
                SPDLOG_ERROR("Unrecognized boundary location!");
        }
        p.getV()[axis] = -1 * p.getV()[axis];
    } else {
        p.getV() = -1 * p.getV();
    }
    return std::nullopt;
}

}  // namespace mol_sim