#ifndef DOMAIN_H
#define DOMAIN_H

#include "particles/boundaries/BoundaryCondition.h"

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
    containerType* particles;
    /**
     * @brief Boundary Type of the left (x) boundary
     *
     */
    std::unique_ptr<BoundaryCondition<containerType>> left_boundary = std::make_unique<Outflow<containerType>>(particles, BoundaryConditionDeclaration::BoundaryType::LEFT);
    /**
     * @brief Boundary Type of the right (x) boundary
     *
     */
    std::unique_ptr<BoundaryCondition<containerType>> right_boundary = std::make_unique<Outflow<containerType>>(particles, BoundaryConditionDeclaration::BoundaryType::RIGHT);
    /**
     * @brief Boundary Type of the back (z) boundary
     *
     */
    std::unique_ptr<BoundaryCondition<containerType>> back_boundary = std::make_unique<Outflow<containerType>>(particles, BoundaryConditionDeclaration::BoundaryType::BACK);
    /**
     * @brief Boundary Type of the front (z) boundary
     *
     */
    std::unique_ptr<BoundaryCondition<containerType>> front_boundary = std::make_unique<Outflow<containerType>>(particles, BoundaryConditionDeclaration::BoundaryType::FRONT);
    /**
     * @brief Boundary Type of the lower (y) boundary
     *
     */
    std::unique_ptr<BoundaryCondition<containerType>> lower_boundary = std::make_unique<Outflow<containerType>>(particles, BoundaryConditionDeclaration::BoundaryType::LOWER);
    /**
     * @brief Boundary Type of the upper (y) boundary
     *
     */
    std::unique_ptr<BoundaryCondition<containerType>> upper_boundary = std::make_unique<Outflow<containerType>>(particles, BoundaryConditionDeclaration::BoundaryType::UPPER);

    /**
     * @brief Sets a specified boundary to a specified boundary condition.
     *
     * @param boundary The fully declared boundary condition.
     */
    void set_boundary(BoundaryConditionDeclaration boundary);
   
    Domain<containerType>& operator=(const Domain<containerType>& domain) {
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

    public:
    Domain() = default;

    Domain(R3 dimension, std::vector<std::optional<BoundaryConditionDeclaration>> boundaries) 
    : dimension(dimension) {
        for (auto& boundary : boundaries) {
            if (boundary.has_value()) { set_boundary(boundary.value()); }
        }
    }

    Domain(R3 dimension, containerType& particles, std::vector<std::optional<BoundaryConditionDeclaration>> boundaries) 
    : dimension(dimension), particles(particles) {
        for (auto& boundary : boundaries) {
            if (boundary.has_value()) { set_boundary(boundary.value()); }
        }
    }

    /**
     * @brief Get the boundary object
     * 
     * @param boundary The boundary (LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK) to be obtained. 
     * @return std::unique_ptr<BoundaryCondition<containerType>> The unique pointer to the specified boundary.
     */
    std::unique_ptr<BoundaryCondition<containerType>> get_boundary(BoundaryType boundary);

    void set_particles(containerType& particles);

    containerType& get_particles();
};

} //namespace mol_sim

#endif
