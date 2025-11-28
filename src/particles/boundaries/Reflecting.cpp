#include "particles/boundaries/Reflecting.h"

#include <cmath>

namespace mol_sim {

Reflecting::Reflecting(BoundaryLocation location, R3 domain_size, std::optional<double> sigma,
                       std::optional<double> epsilon)
    : Boundary(location, BoundaryType::REFLECTING),
      domain_size(domain_size),
      boundary_epsilon(epsilon),
      boundary_sigma(sigma) {}

size_t Reflecting::getAxis() const {
    switch (location) {
        case BoundaryLocation::LEFT:
        case BoundaryLocation::RIGHT:
            return 0;
        case BoundaryLocation::UPPER:
        case BoundaryLocation::LOWER:
            return 1;
        case BoundaryLocation::FRONT:
        case BoundaryLocation::BACK:
            return 2;
        default:
            return 0;
    }
}

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

void Reflecting::applyBoundary(Particle& p) { (void)p; }

std::optional<Particle> Reflecting::computeGhostParticle(const Particle& p) const {
    double sigma = boundary_sigma.value_or(p.getSigma());
    double epsilon = boundary_epsilon.value_or(p.getEpsilon());

    size_t axis = getAxis();
    int sign = getSign();
    double boundary_position = getBoundaryPosition();

    double threshold = std::pow(2.0, 1.0 / 6.0) * sigma;
    double particle_pos = p.getX()[axis];
    double distance_to_boundary = (sign < 0) ? (particle_pos - boundary_position) : (boundary_position - particle_pos);

    if (distance_to_boundary < threshold && distance_to_boundary >= 0) {
        R3 ghost_pos = p.getX();
        ghost_pos[axis] = boundary_position + (sign * (-distance_to_boundary));

        return Particle(ghost_pos, {0.0, 0.0, 0.0}, 0.0, epsilon, sigma, -1);
    }

    return std::nullopt;
}

}  // namespace mol_sim
