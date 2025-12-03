#include "particles/boundaries/Reflecting.h"

#include <spdlog/spdlog.h>

#include <cmath>

namespace mol_sim {

Reflecting::Reflecting(BoundaryLocation location, R3 domain_size, bool ghost_on_boundary, std::optional<double> sigma,
                       std::optional<double> epsilon)
    : Boundary(location, BoundaryType::REFLECTING),
      domain_size(domain_size),
      ghost_on_boundary(ghost_on_boundary),
      boundary_epsilon(epsilon),
      boundary_sigma(sigma) {}

int Reflecting::getSign() const {
    switch (location) {
        case BoundaryLocation::LEFT:
        case BoundaryLocation::LOWER:
        case BoundaryLocation::FRONT:
            return -1;
        case BoundaryLocation::RIGHT:
        case BoundaryLocation::UPPER:
        case BoundaryLocation::BACK:
            return 1;
        default:
            return 0;
    }
}

double Reflecting::getBoundaryPosition() const {
    size_t axis = getAxis();
    int sign = getSign();
    return (sign < 0) ? 0.0 : domain_size[axis];
}

std::optional<Particle> Reflecting::applyBoundary(Particle& p) {
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

        SPDLOG_DEBUG("New ghost particle on ( {} , {} , {} )", ghost_pos[0], ghost_pos[1], ghost_pos[2]);
        return Particle(ghost_pos, {0.0, 0.0, 0.0}, 0.0, epsilon, sigma, -1);
    }

    return std::nullopt;
}

}  // namespace mol_sim
