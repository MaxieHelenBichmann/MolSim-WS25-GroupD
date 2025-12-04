#include "particles/boundaries/Reflecting.h"

#include <spdlog/spdlog.h>

#include <cmath>

namespace mol_sim {

Reflecting::Reflecting(BoundaryLocation location, R3 domain_size, bool ghost_on_boundary, std::optional<double> sigma,
                       std::optional<double> epsilon) noexcept
    : Boundary(location, BoundaryType::REFLECTING, domain_size),
      ghost_on_boundary(ghost_on_boundary),
      boundary_epsilon(epsilon),
      boundary_sigma(sigma) {}

void Reflecting::applyBoundary(Particle& p, const ForceSource& force) const noexcept {
    /**
     * Note that we spawn the ghost particle with either the specified sigma or the sigma of the passed particle,
     * BUT: we check whether or not the boundary condition applies with sigma_check. We do this because in the 
     * calculation of the Lennard-Jones-Force we use (particle.sigma + other.sigma)/2 as the sigma in the formula.
     * This means that if we want to be sure that the Lennard-Jones-Force actually pushes the particle away from the ghost
     * instead of having the ghost attract it, we need to make sure that the distance between the particle and the ghost is
     * < 2^(1/6)*((p.sigma + ghost.sigma)/2), instead of just 2^(1/6) * ghost.sigma(). This is because it could happen that:
     * 2^(1/6)*((p.sigma + ghost.sigma)/2) < 2^(1/6)*ghost.sigma. 
     * In other words this means if we used sigma_ghost instead of sigma_check for checking the applicability of the boundary 
     * condition, if p.sigma > ghost.sigma the ghost particle would actually attract the particle.
     * 
     * Again note that the whole reason we do this sigma = (p.sigma + ghost.sigma)/2 in LennardJonesForce.cpp is because
     * it's an intuitive way of combining two different sigmas in the LJforce calculation. This is essential if we want
     * to use the Newton's 3rd law optimization in 
     */
    
    double sigma_ghost = boundary_sigma.value_or(p.getSigma());
    double sigma_check = (boundary_sigma.value_or(p.getSigma()) + p.getSigma()) / 2;
    double epsilon = boundary_epsilon.value_or(p.getEpsilon());

    size_t axis = getAxis();
    int sign = getSign();
    double boundary_position = getBoundaryPosition();

    double threshold = ghost_on_boundary ? pow(2.0, 1.0 / 6.0) * sigma_check : pow(2.0, -5.0 / 6.0) * sigma_check;
    double particle_pos = p.getX()[axis];
    double distance_to_boundary = (sign < 0) ? (particle_pos - boundary_position) : (boundary_position - particle_pos);

    if (distance_to_boundary < threshold && distance_to_boundary >= 0) {
        R3 ghost_pos = p.getX();
        ghost_pos[axis] = boundary_position + (ghost_on_boundary ? 0 : sign * distance_to_boundary);

        SPDLOG_TRACE("Creating ghost particle at ({}, {}, {}) for particle at ({}, {}, {})", ghost_pos[0], ghost_pos[1],
                     ghost_pos[2], p.getX()[0], p.getX()[1], p.getX()[2]);
        p.getF() = p.getF() + force.applyForce(p, Particle(ghost_pos, {0.0, 0.0, 0.0}, p.getM(), epsilon, sigma_ghost, -1));
    }
}

}  // namespace mol_sim
