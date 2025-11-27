#include "particles/container/domain/Domain.h"

#include <spdlog/spdlog.h>

namespace mol_sim {
    Domain::Domain() = default;
    
    Domain::Domain(Domain& domain) noexcept {
        this->dimension = domain.dimension;
        left_boundary = domain.left_boundary;
        right_boundary = domain.right_boundary;
        upper_boundary = domain.upper_boundary;
        lower_boundary = domain.lower_boundary;
        front_boundary = domain.front_boundary;
        back_boundary = domain.back_boundary;
    }

    Domain::Domain(Domain&& domain) noexcept {
        this->dimension = domain.dimension;
        left_boundary = domain.left_boundary;
        right_boundary = domain.right_boundary;
        upper_boundary = domain.upper_boundary;
        lower_boundary = domain.lower_boundary;
        front_boundary = domain.front_boundary;
        back_boundary = domain.back_boundary;
    }

    Domain::Domain(R3 dimension, std::vector<std::optional<BoundaryCondition*>> boundaries) : dimension(dimension) {
        for (auto& boundary : boundaries) {
            if (boundary.has_value()) { 
                getBoundary(boundary.value()->getLocation()) = boundary.value();
            }
        }
    }

    Domain::Domain(R3 dimension, std::vector<std::optional<BoundaryCondition&>> boundaries) : dimension(dimension) {
        for (auto& boundary : boundaries) {
            if (boundary.has_value()) { 
                getBoundary(boundary.value().getLocation()) = &boundary.value();
            }
        }
    }

    R3 Domain::getDimension() { return dimension; }

    BoundaryCondition*& Domain::getBoundary(BoundaryLocation location) {
        switch(location) {
            case BoundaryLocation::LEFT: return left_boundary;
            case BoundaryLocation::RIGHT: return right_boundary;
            case BoundaryLocation::FRONT: return front_boundary;
            case BoundaryLocation::BACK: return back_boundary;
            case BoundaryLocation::UPPER: return upper_boundary; 
            case BoundaryLocation::LOWER: return lower_boundary;
            default:
                SPDLOG_ERROR("Unknown boundary type! Expected (LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK)!");
                exit(-1); 
        }
    }

    const BoundaryCondition* Domain::getBoundary(BoundaryLocation location) const {
        switch(location) {
            case BoundaryLocation::LEFT: return left_boundary;
            case BoundaryLocation::RIGHT: return right_boundary;
            case BoundaryLocation::FRONT: return front_boundary;
            case BoundaryLocation::BACK: return back_boundary;
            case BoundaryLocation::UPPER: return upper_boundary; 
            case BoundaryLocation::LOWER: return lower_boundary;
            default:
                SPDLOG_ERROR("Unknown boundary type! Expected (LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK)!");
                exit(-1); 
        }
    }

    void Domain::applyBoundary(Particle& p) { 
        left_boundary->boundaryStrategy(p);
        right_boundary->boundaryStrategy(p);
        front_boundary->boundaryStrategy(p);
        back_boundary->boundaryStrategy(p);
        upper_boundary->boundaryStrategy(p);
        lower_boundary->boundaryStrategy(p); 
    }
}  // namespace mol_sim