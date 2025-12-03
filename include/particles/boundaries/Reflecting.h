#ifndef REFLECTING_H
#define REFLECTING_H

#include <optional>

#include "particles/boundaries/Boundary.h"

namespace mol_sim {

class Reflecting : public Boundary {
    bool ghost_on_boundary = false;
    std::optional<double> boundary_epsilon;
    std::optional<double> boundary_sigma;

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
    Reflecting(BoundaryLocation location, R3 domain_size, bool ghost_on_boundary,
               std::optional<double> sigma = std::nullopt, std::optional<double> epsilon = std::nullopt);
    ~Reflecting() override = default;

    void applyBoundary(Particle& p, const ForceSource& force) const override;
    [[nodiscard]] std::optional<double> getBoundarySigma() const { return boundary_sigma; }
    [[nodiscard]] std::optional<double> getBoundaryEpsilon() const { return boundary_epsilon; }
};

}  // namespace mol_sim

#endif
