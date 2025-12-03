#ifndef DOMAIN_H
#define DOMAIN_H

#include <array>
#include <memory>
#include <vector>

#include "particles/boundaries/Boundary.h"
#include "physics/ForceSource.h"
#include "utils/Vector.h"

namespace mol_sim {

/**
 * @brief Provides a wrapper for domain configuration
 */
class Domain {
    R3 dimension = {100., 100., 100.};

    /**
     * @brief Array of owned boundary conditions in order: LEFT, RIGHT, FRONT, BACK, UPPER, LOWER
     */
    std::array<std::unique_ptr<Boundary>, 6> boundaries;

    /**
     * @brief Maps BoundaryLocation enum to array index
     * @throws BoundaryException if the given location is unknown.
     */
    static size_t locationToIndex(BoundaryLocation location);

   public:
    Domain();
    Domain(R3 dimension);
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
     * @brief Get a reference to a specific Boundary object.
     *
     * @param location The location of the boundary.
     * @return Boundary& The reference to the specified boundary.
     * @throws BoundaryException if the boundary at the given location is null.
     */
    Boundary& getBoundary(BoundaryLocation location);

    /**
     * @brief Get a const reference to a specific Boundary object.
     *
     * @param location The location of the boundary.
     * @return const Boundary& The reference to the specified boundary.
     * @throws BoundaryException if the boundary at the given location is null.
     */
    [[nodiscard]] const Boundary& getBoundary(BoundaryLocation location) const;

    /**
     * @brief Get the domain's dimension vector.
     *
     * @return R3 The dimensions of the domain.
     */
    [[nodiscard]] R3 getDimension() const noexcept;

    /**
     * @brief Applies boundary conditions to a particle. Also (if needed) computes ghost particles
     * needed for a given particle.
     *
     * @param p The particle to apply boundaries to / compute ghost particles for.
     * @param force The force source for ghost particle interactions.
     */
    void applyBoundary(Particle& p, const ForceSource& force) const noexcept;

    Domain& operator=(const Domain& other) = delete;
    Domain& operator=(Domain&& other) noexcept;
};

}  // namespace mol_sim

#endif
