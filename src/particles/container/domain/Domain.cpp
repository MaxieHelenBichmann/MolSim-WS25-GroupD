#include "particles/container/domain/Domain.h"

namespace mol_sim {
    template<ParticleContainer containerType>
    void Domain<containerType>::set_boundary(BoundaryConditionDeclaration boundary) {
        std::unique_ptr<BoundaryCondition<containertype>> ptr;
        switch(boundary.get_location()) {
            case LEFT:
                ptr = left_boundary;
                break; 
            case RIGHT:
                ptr = right_boundary;
                break; 
            case TOP:
                ptr = top_boundary;
                break; 
            case BOTTOM:
                ptr = bottom_boundary;
                break; 
            case FRONT:
                ptr = front_boundary;
                break; 
            case BACK:
                ptr = back_boundary;
                break;
            default:
                SPDLOG_ERROR("Unknown boundary type! Expected (LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK)!"); 
        }
        
        switch (boundary.get_type()) {
            case OUTFLOW:
                ptr.reset(new Outflow<containerType>(particles, boundary.get_location()));
                break;
            case REFLECTING:
                ptr.reset(new Reflecting<containerType>(particles, boundary.get_location(), boundary.get_counter_sigma(), boundary.get_counter_epsilon()));
                break;
            default:
                SPDLOG_ERROR("Unknown boundary condition!");
                break;
        }
    }

    template<ParticleContainer containerType>
    std::unique_ptr<BoundaryCondition<containerType>> Domain<containerType>::get_boundary(BoundaryType boundary) {
        switch(boundary) {
            case LEFT: return left_boundary;
            case RIGHT: return right_boundary;
            case UPPER: return top_boundary; 
            case LOWER: return bottom_boundary;
            case FRONT: return front_boundary;
            case BACK: return back_boundary;
            default:
                SPDLOG_ERROR("Unknown boundary type! Expected (LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK)!");
                return nullptr; 
        }
    }
} //namespace mol_sim