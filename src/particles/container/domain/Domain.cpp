#include "particles/container/domain/Domain.h"

#include <spdlog/spdlog.h>

namespace mol_sim {
template <ParticleContainer containerType>
void Domain<containerType>::setBoundary(BoundaryConditionDeclaration boundary) {
    std::unique_ptr<BoundaryCondition<containerType>> ptr;
    switch (boundary.getLocation()) {
        case BoundaryType::LEFT:
            ptr = left_boundary;
            break;
        case BoundaryType::RIGHT:
            ptr = right_boundary;
            break;
        case BoundaryType::UPPER:
            ptr = upper_boundary;
            break;
        case BoundaryType::LOWER:
            ptr = lower_boundary;
            break;
        case BoundaryType::FRONT:
            ptr = front_boundary;
            break;
        case BoundaryType::BACK:
            ptr = back_boundary;
            break;
        default:
            SPDLOG_ERROR("Unknown boundary type! Expected (LEFT, RIGHT, UPPER, LOWER, FRONT, BACK)!");
    }

    switch (boundary.getType()) {
        case BoundaryConditionType::OUTFLOW:
            ptr.reset(new Outflow<containerType>(particles, boundary.getLocation()));
            break;
        case BoundaryConditionType::REFLECTING:
            ptr.reset(new Reflecting<containerType>(particles, boundary.getLocation(), boundary.getCounterSigma(),
                                                    boundary.getCounterEpsilon()));
            break;
        default:
            SPDLOG_ERROR("Unknown boundary condition!");
            break;
    }
}

template <ParticleContainer containerType>
std::unique_ptr<BoundaryCondition<containerType>> Domain<containerType>::getBoundary(BoundaryType boundary) {
    switch (boundary) {
        case BoundaryType::LEFT:
            return left_boundary;
        case BoundaryType::RIGHT:
            return right_boundary;
        case BoundaryType::UPPER:
            return upper_boundary;
        case BoundaryType::LOWER:
            return lower_boundary;
        case BoundaryType::FRONT:
            return front_boundary;
        case BoundaryType::BACK:
            return back_boundary;
        default:
            SPDLOG_ERROR("Unknown boundary type! Expected (LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK)!");
            return nullptr;
    }
}
}  // namespace mol_sim