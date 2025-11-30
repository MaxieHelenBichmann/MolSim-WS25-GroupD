#ifndef REFLECTING_H
#define REFLECTING_H

#include <cmath>
#include <optional>

#include "particles/boundaries/Boundary.h"
#include "utils/Vector.h"

namespace mol_sim {

class Reflecting : public Boundary {
    R3 domain_size;
    bool ghost_on_boundary = false;
    std::optional<double> boundary_epsilon;
    std::optional<double> boundary_sigma;

    /**
     * @brief Gets the axis index for the boundary location.
     * @return 0 for x-axis (LEFT/RIGHT), 1 for y-axis (UPPER/LOWER), 2 for z-axis (FRONT/BACK).
     */
    [[nodiscard]] size_t getAxis() const;

    /**
     * @brief Gets the sign/direction for the boundary.
     * @return -1 for boundaries at min (LEFT, LOWER, FRONT), +1 for boundaries at max (RIGHT, UPPER, BACK).
     */
    [[nodiscard]] int getSign() const;

    /**
     * @brief Computes the boundary position along its axis from domain size.
     * @return 0 for boundaries at min (LEFT, LOWER, FRONT), domain_size[axis] for boundaries at max.
     */
    [[nodiscard]] double getBoundaryPosition() const;

   public:
    /**
     * @brief Constructs a Reflecting boundary.
     *
     * @param location The boundary location (LEFT, RIGHT, etc.).
     * @param domain_size The dimensions of the entire domain (x, y, z).
     * @param ghost_on_boundary States that the ghost particles should be spawned right on the boundary if true 
     * instead of having them mirror.
     * @param sigma Optional sigma for ghost particle interactions.
     * @param epsilon Optional epsilon for ghost particle interactions.
     */
    Reflecting(BoundaryLocation location, R3 domain_size, bool ghost_on_boundary, std::optional<double> sigma = std::nullopt,
               std::optional<double> epsilon = std::nullopt);
    ~Reflecting() override = default;

    [[nodiscard]] std::optional<Particle> applyBoundary(Particle& p) override;
    [[nodiscard]] std::optional<double> getBoundarySigma() const { return boundary_sigma; }
    [[nodiscard]] std::optional<double> getBoundaryEpsilon() const { return boundary_epsilon; }
    [[nodiscard]] R3 getDomainSize() const { return domain_size; }
};

}  // namespace mol_sim

#endif
