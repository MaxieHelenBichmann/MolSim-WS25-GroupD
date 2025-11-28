#ifndef REFLECTING_H
#define REFLECTING_H

#include <cmath>
#include <optional>

#include "particles/boundaries/Boundary.h"
#include "utils/Vector.h"

namespace mol_sim {

class Reflecting : public Boundary {
    // position of this boundary along its axis
    double boundary_position;
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

   public:
    /**
     * @brief Constructs a Reflecting boundary.
     *
     * @param location The boundary location (LEFT, RIGHT, etc.).
     * @param boundary_position The coordinate of this boundary along its axis.
     * @param sigma Optional sigma for ghost particle interactions.
     * @param epsilon Optional epsilon for ghost particle interactions.
     */
    Reflecting(BoundaryLocation location, double boundary_position, std::optional<double> sigma = std::nullopt,
               std::optional<double> epsilon = std::nullopt);
    ~Reflecting() override = default;

    void applyBoundary(Particle& p) override;
    [[nodiscard]] std::optional<Particle> computeGhostParticle(const Particle& p) const override;

    [[nodiscard]] std::optional<double> getBoundarySigma() const { return boundary_sigma; }
    [[nodiscard]] std::optional<double> getBoundaryEpsilon() const { return boundary_epsilon; }
    [[nodiscard]] double getBoundaryPosition() const { return boundary_position; }
};

}  // namespace mol_sim

#endif
