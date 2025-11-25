#include "particles/container/domain/Domain.h"

#include <spdlog/spdlog.h>

namespace mol_sim {
    template <ParticleContainer containerType>
    void Domain<containerType>::setBoundary(BoundaryConditionDeclaration boundary) {
        std::unique_ptr<BoundaryCondition<containerType>> ptr;
        switch (boundary.getLocation()) {
            case BoundaryLocation::LEFT:
                ptr = left_boundary;
                break; 
            case BoundaryLocation::RIGHT:
                ptr = right_boundary;
                break; 
            case BoundaryLocation::UPPER:
                ptr = upper_boundary;
                break; 
            case BoundaryLocation::LOWER:
                ptr = lower_boundary;
                break; 
            case BoundaryLocation::FRONT:
                ptr = front_boundary;
                break; 
            case BoundaryLocation::BACK:
                ptr = back_boundary;
                break;
            default:
                SPDLOG_ERROR("Unknown boundary type! Expected (LEFT, RIGHT, UPPER, LOWER, FRONT, BACK)!"); 
        }
        
        switch (boundary.getType()) {
            case BoundaryType::OUTFLOW:
                ptr.reset(new Outflow<containerType>(particles, boundary.getLocation()));
                break;
            case BoundaryType::REFLECTING:
                ptr.reset(new Reflecting<containerType>(particles, dimension, boundary.getLocation(), boundary.getCounterSigma(), boundary.getCounterEpsilon()));
                break;
            default:
                SPDLOG_ERROR("Unknown boundary condition!");
                break;
        }
    }

    template <ParticleContainer containerType>
    std::unique_ptr<BoundaryCondition<containerType>> Domain<containerType>::getBoundary(BoundaryLocation location) {
        switch(location) {
            case BoundaryLocation::LEFT: return left_boundary;
            case BoundaryLocation::RIGHT: return right_boundary;
            case BoundaryLocation::UPPER: return upper_boundary; 
            case BoundaryLocation::LOWER: return lower_boundary;
            case BoundaryLocation::FRONT: return front_boundary;
            case BoundaryLocation::BACK: return back_boundary;
            default:
                SPDLOG_ERROR("Unknown boundary type! Expected (LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK)!");
                return nullptr; 
        }
    }
} //namespace mol_sim