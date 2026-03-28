#ifndef REFLECTING_H
#define REFLECTING_H

#include <optional>

#include "particles/boundaries/Boundary.h"

namespace mol_sim {

class Reflecting : public Boundary {
    /**
     * @brief Flag that specifies whether ghost particles should be spawned exactly on the
     * boundary or spawn like "mirror particles" like in the slides of Meeting/Lecture 3.
     *
     * If true then the ghost particles should be spawned on boundary, if false then we do the
     * "mirror particles".
     */
    bool ghost_on_boundary = false;
    /**
     * @brief The epsilon value the spawned ghost particles should have in this boundary.
     * If it is std::nullopt the epsilon of the particle which spawns the ghost will be used.
     */
    std::optional<double> boundary_epsilon;
    /**
     * @brief The sigma value the spawned ghost particles should have in this boundary.
     * If it is std::nullopt the sigma of the particle which spawns the ghost will be used.
     */
    std::optional<double> boundary_sigma;

   public:
    /**
     * @brief Constructs a Reflecting boundary.
     *
     * @param location The boundary location (LEFT, RIGHT, etc.).
     * @param domain_size The dimensions of the entire domain (x, y, z).
     * @param ghost_on_boundary States that the ghost particles should be spawned right on the boundary if true
     * instead of having them mirror.
     * @param sigma Optional sigma for ghost particle interactions. std::nullopt if you want to use the particle's
     * sigma.
     * @param epsilon Optional epsilon for ghost particle interactions. std::nullopt if you want to use the particle's
     * epsilon.
     */
    Reflecting(BoundaryLocation location, R3 domain_size, bool ghost_on_boundary,
               std::optional<double> sigma = std::nullopt, std::optional<double> epsilon = std::nullopt) noexcept;
    ~Reflecting() override = default;

    /**
     * @brief The applyBoundary routine of Reflecting. It checks if the given particle is close enough to the boundary
     * to where this particle would be in the repulsive range of the ghost particle. If that is the case, it spawns said
     * ghost particle to push the particle away from the boundary.
     *
     * @param p The particle which the boundary should be applied to.
     * @param force The force source to be used for the inter-force calculation between ghost and normal particle.
     * For now only supports Lennard-Jones-Force.
     * @return std::nullopt always, since ghost particles are created, applied (in terms of force calculation) and
     * destroyed entirely within this method.
     */
    void applyBoundary(Particle& p, const PairwiseForceSource& force) noexcept override;
    /**
     * @brief Returns the sigma of the ghost particles of this boundary.
     *
     * @return std::optional<double> The sigma of the ghost particles of this boundary. If
     * std::nullopt is returned that means the ghosts take the sigma of the normal particles.
     */
    [[nodiscard]] std::optional<double> getBoundarySigma() const noexcept { return boundary_sigma; }
    /**
     * @brief Returns the epsilon of the ghost particles of this boundary.
     *
     * @return std::optional<double> The epsilon of the ghost particles of this boundary. If
     * std::nullopt is returned that means the ghosts take the epsilon of the normal particles.
     */
    [[nodiscard]] std::optional<double> getBoundaryEpsilon() const noexcept { return boundary_epsilon; }
    /**
     * @brief Returns frag whether ghost particles are spawned on the boundary or mirrored.
     *
     * @return true if ghost particles are spawned on the boundary, false if they are mirrored.
     */
    [[nodiscard]] bool isGhostOnBoundary() const noexcept { return ghost_on_boundary; }
};

}  // namespace mol_sim

#endif
