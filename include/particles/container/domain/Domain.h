#ifndef DOMAIN_H
#define DOMAIN_H

#include "particles/boundaries/BoundaryCondition.h"
#include "particles/boundaries/Outflow.h"
#include "particles/boundaries/Reflecting.h"
#include "utils/Vector.h"

namespace mol_sim {
/**
 * @brief Provides a wrapper for domain configuration
 *
 */
class Domain {
    /**
     * @brief Dimensions/size of the domain
     *
     */
    R3 dimension = {1., 1., 1.};
    Outflow _left_default = Outflow(BoundaryLocation::LEFT);
    Outflow _right_default = Outflow(BoundaryLocation::RIGHT);
    Outflow _front_default = Outflow(BoundaryLocation::FRONT);
    Outflow _back_default = Outflow(BoundaryLocation::BACK);
    Outflow _upper_default = Outflow(BoundaryLocation::UPPER);
    Outflow _lower_default = Outflow(BoundaryLocation::LOWER);
    /**
     * @brief Boundary Type of the left (x) boundary
     *
     */
    BoundaryCondition* left_boundary = &_left_default;
    /**
     * @brief Boundary Type of the back (z) boundary
     *
     */
    BoundaryCondition* right_boundary = &_right_default;
    /**
     * @brief Boundary Type of the front (z) boundary
     *
     */
    BoundaryCondition* front_boundary = &_front_default;
    /**
     * @brief Boundary Type of the back (z) boundary
     *
     */
    BoundaryCondition* back_boundary = &_back_default;
    /**
     * @brief Boundary Type of the upper (y) boundary
     *
     */
    BoundaryCondition* upper_boundary = &_upper_default;
    /**
     * @brief Boundary Type of the lower (y) boundary
     *
     */
    BoundaryCondition* lower_boundary = &_lower_default;
   
    public:
    // constructors
    Domain() = default;
    Domain(Domain& domain) noexcept;
    Domain(Domain&& domain) noexcept;
    Domain(R3 dimension, std::vector<std::optional<BoundaryCondition*>> boundaries);
    Domain(R3 dimension, std::vector<std::optional<BoundaryCondition&>> boundaries);
    
    /**
     * @brief Get a non-const pointer to a specific BoundaryCondition object.
     *
     * @param location The location of the boundary (LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK) to be obtained.
     * @return BoundaryCondition* The pointer to the specified boundary.
     */
    BoundaryCondition*& getBoundary(BoundaryLocation location);
    /**
     * @brief Get a const pointer to a specific BoundaryCondition object.
     *
     * @param location The location of the boundary (LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK) to be obtained.
     * @return BoundaryCondition* The pointer to the specified boundary.
     */
    const BoundaryCondition* getBoundary(BoundaryLocation location) const;
    /**
     * @brief Get the domain's dimension vector.
     * 
     * @return R3 The vector specifiying the dimension of the domain. 
     */
    R3 getDimension(); 
    /**
     * @brief Applies the applicable boundary conditions to the given particle, i.e. exactly those, where the
     * particle is close enough. 
     * 
     * @param p The particle that the boundary condition should influence (if applicable). 
     * @param locations The locations where a boundary condition applies to the particle. 
     */
    void applyBoundary(Particle& p);
    
    Domain& operator=(const Domain& domain) noexcept {
        this->dimension = domain.dimension;
        left_boundary = domain.left_boundary;
        right_boundary = domain.right_boundary;
        upper_boundary = domain.upper_boundary;
        lower_boundary = domain.lower_boundary;
        front_boundary = domain.front_boundary;
        back_boundary = domain.back_boundary;
        return *this;
    }

    Domain& operator=(Domain&& domain) noexcept {
        this->dimension = domain.dimension;
        left_boundary = domain.left_boundary;
        right_boundary = domain.right_boundary;
        upper_boundary = domain.upper_boundary;
        lower_boundary = domain.lower_boundary;
        front_boundary = domain.front_boundary;
        back_boundary = domain.back_boundary;
        return *this;
    }
};
}  // namespace mol_sim

#endif
