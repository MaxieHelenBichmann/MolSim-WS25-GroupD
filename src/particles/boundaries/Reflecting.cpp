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
    double sigma = boundary_sigma.value_or(p.getSigma());
    double epsilon = boundary_epsilon.value_or(p.getEpsilon());

    size_t axis = getAxis();
    int sign = getSign();
    double boundary_position = getBoundaryPosition();

    double threshold = ghost_on_boundary ? pow(2.0, 1.0 / 6.0) * sigma : pow(2.0, -5.0 / 6.0) * sigma;
    double particle_pos = p.getX()[axis];
    double distance_to_boundary = (sign < 0) ? (particle_pos - boundary_position) : (boundary_position - particle_pos);

    if (distance_to_boundary < threshold && distance_to_boundary >= 0) {
        R3 ghost_pos = p.getX();
        ghost_pos[axis] = boundary_position + (ghost_on_boundary ? 0 : sign * distance_to_boundary);

        SPDLOG_TRACE("Creating ghost particle at ({}, {}, {}) for particle at ({}, {}, {})", ghost_pos[0], ghost_pos[1],
                     ghost_pos[2], p.getX()[0], p.getX()[1], p.getX()[2]);
        p.getF() = p.getF() + force.applyForce(p, Particle(ghost_pos, {0.0, 0.0, 0.0}, p.getM(), epsilon, sigma, -1));
    }
}

}  // namespace mol_sim
