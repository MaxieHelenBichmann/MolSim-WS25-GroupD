#ifndef DOMAIN_H
#define DOMAIN_H

#include <memory>

#include "particles/boundaries/BoundaryCondition.h"
#include "particles/boundaries/Outflow.h"
#include "particles/boundaries/Reflecting.h"
#include "utils/Vector.h"

namespace mol_sim {
/**
 * @brief Provides a wrapper for domain configuration
 *
 */

template <ParticleContainer containerType>
class Domain {
    /**
     * @brief Dimensions/size of the domain
     *
     */
    R3 dimension = {1., 1., 1.};
    /**
     * @brief The container containing the actual particles
     *
     */
    containerType* particles = nullptr;
    /**
     * @brief Boundary Type of the left (x) boundary
     *
     */
    std::unique_ptr<BoundaryCondition<containerType>> left_boundary = 
        std::make_unique<Outflow<containerType>>(particles, BoundaryLocation::LEFT);
    /**
     * @brief Boundary Type of the right (x) boundary
     *
     */
    std::unique_ptr<BoundaryCondition<containerType>> right_boundary = 
        std::make_unique<Outflow<containerType>>(particles, BoundaryLocation::RIGHT);
    /**
     * @brief Boundary Type of the back (z) boundary
     *
     */
    std::unique_ptr<BoundaryCondition<containerType>> back_boundary = 
        std::make_unique<Outflow<containerType>>(particles, BoundaryLocation::BACK);
    /**
     * @brief Boundary Type of the front (z) boundary
     *
     */
    std::unique_ptr<BoundaryCondition<containerType>> front_boundary = 
        std::make_unique<Outflow<containerType>>(particles, BoundaryLocation::FRONT);
    /**
     * @brief Boundary Type of the lower (y) boundary
     *
     */
    std::unique_ptr<BoundaryCondition<containerType>> lower_boundary = 
        std::make_unique<Outflow<containerType>>(particles, BoundaryLocation::LOWER);
    /**
     * @brief Boundary Type of the upper (y) boundary
     *
     */
    std::unique_ptr<BoundaryCondition<containerType>> upper_boundary = 
        std::make_unique<Outflow<containerType>>(particles, BoundaryLocation::UPPER);

    /**
     * @brief Sets a specified boundary to a specified boundary condition.
     *
     * @param boundary The fully declared boundary condition.
     */
    void setBoundary(BoundaryConditionDeclaration boundary);
   
    public:
    Domain() = default; 

    Domain<containerType>& operator=(Domain<containerType>&& domain) noexcept {
        this->dimension = domain.dimension;
        this->particles = domain.particles;
        left_boundary = std::move(domain.left_boundary);
        right_boundary = std::move(domain.right_boundary);
        upper_boundary = std::move(domain.upper_boundary);
        lower_boundary = std::move(domain.lower_boundary);
        front_boundary = std::move(domain.front_boundary);
        back_boundary = std::move(domain.back_boundary);
        return *this;
    }

    Domain(Domain<containerType>& domain) noexcept {
        this->dimension = domain.dimension;
        this->particles = domain.particles;
        left_boundary = std::move(domain.left_boundary);
        right_boundary = std::move(domain.right_boundary);
        upper_boundary = std::move(domain.upper_boundary);
        lower_boundary = std::move(domain.lower_boundary);
        front_boundary = std::move(domain.front_boundary);
        back_boundary = std::move(domain.back_boundary);
    }

    Domain(Domain<containerType>&& domain) noexcept {
        this->dimension = domain.dimension;
        this->particles = domain.particles;
        left_boundary = std::move(domain.left_boundary);
        right_boundary = std::move(domain.right_boundary);
        upper_boundary = std::move(domain.upper_boundary);
        lower_boundary = std::move(domain.lower_boundary);
        front_boundary = std::move(domain.front_boundary);
        back_boundary = std::move(domain.back_boundary);
    }

    Domain(R3 dimension, std::vector<std::optional<BoundaryConditionDeclaration>> boundaries) : dimension(dimension) {
        for (auto& boundary : boundaries) {
            if (boundary.has_value()) { 
                setBoundary(boundary.value()); 
            }
        }
    }

    Domain(R3 dimension, containerType& particles, std::vector<std::optional<BoundaryConditionDeclaration>> boundaries) 
        : dimension(dimension), particles(particles) {
        for (auto& boundary : boundaries) {
            if (boundary.has_value()) { 
                setBoundary(boundary.value()); 
            }
        }
    }

    /**
     * @brief Get the boundary object
     *
     * @param location The location of the boundary (LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK) to be obtained.
     * @return std::unique_ptr<BoundaryCondition<containerType>> The unique pointer to the specified boundary.
     */
    std::unique_ptr<BoundaryCondition<containerType>> getBoundary(BoundaryLocation location);

    void setParticles(containerType& particles) { this->particles = &particles; }

    containerType& getParticles() { return particles; }

    R3 getDimension() { return dimension; }
};
}  // namespace mol_sim

#endif
