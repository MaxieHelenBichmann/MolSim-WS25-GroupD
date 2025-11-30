#ifndef DOMAIN_H
#define DOMAIN_H

#include <array>
#include <memory>
#include <vector>

#include "particles/boundaries/Boundary.h"
#include "particles/boundaries/Outflow.h"
#include "particles/boundaries/Reflecting.h"
#include "utils/Vector.h"

namespace mol_sim {

/**
 * @brief Provides a wrapper for domain configuration
 */
class Domain {
    R3 dimension = {1., 1., 1.};

    /**
     * @brief Array of owned boundary conditions in order: LEFT, RIGHT, FRONT, BACK, UPPER, LOWER
     */
    std::array<std::unique_ptr<Boundary>, 6> boundaries;

    /**
     * @brief Maps BoundaryLocation enum to array index
     */
    static size_t locationToIndex(BoundaryLocation location);

   public:
    Domain();
    Domain(const Domain& other) = delete;
    Domain(Domain&& other) noexcept;

    /**
     * @brief Constructs a Domain with the given dimension and boundary conditions.
     *
     * @param dimension The dimensions of the domain (x, y, z).
     * @param boundaries Array of 6 boundary conditions in order: LEFT, RIGHT, FRONT, BACK, UPPER, LOWER.
     */
    Domain(R3 dimension, std::array<std::unique_ptr<Boundary>, 6> boundaries);

    /**
     * @brief Get a pointer to a specific Boundary object.
     *
     * @param location The location of the boundary.
     * @return Boundary* The pointer to the specified boundary.
     */
    Boundary* getBoundary(BoundaryLocation location);

    /**
     * @brief Get a const pointer to a specific Boundary object.
     *
     * @param location The location of the boundary.
     * @return const Boundary* The pointer to the specified boundary.
     */
    [[nodiscard]] const Boundary* getBoundary(BoundaryLocation location) const;

    /**
     * @brief Get the domain's dimension vector.
     *
     * @return R3 The dimensions of the domain.
     */
    [[nodiscard]] R3 getDimension() const;

    /**
     * @brief Applies boundary conditions to a particle. Also (if needed) computes ghost particles 
     * needed for a given particle.
     *
     * @param p The particle to apply boundaries to / compute ghost particles for.
     * @return Vector of ghost particles.
     */
    [[nodiscard]] std::vector<Particle> applyBoundary(Particle& p) const;

    Domain& operator=(const Domain& other) = delete;
    Domain& operator=(Domain&& other) noexcept;
};

}  // namespace mol_sim

#endif
