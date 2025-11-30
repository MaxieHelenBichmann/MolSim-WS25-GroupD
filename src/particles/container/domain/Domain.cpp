#include "particles/container/domain/Domain.h"

#include <spdlog/spdlog.h>

namespace mol_sim {

size_t Domain::locationToIndex(BoundaryLocation location) {
    switch (location) {
        case BoundaryLocation::LEFT:
            return 0;
        case BoundaryLocation::RIGHT:
            return 1;
        case BoundaryLocation::FRONT:
            return 2;
        case BoundaryLocation::BACK:
            return 3;
        case BoundaryLocation::UPPER:
            return 4;
        case BoundaryLocation::LOWER:
            return 5;
        default:
            SPDLOG_ERROR("Unknown boundary location!");
            return 0;
    }
}

Domain::Domain() {
    boundaries[0] = std::make_unique<Outflow>(BoundaryLocation::LEFT);
    boundaries[1] = std::make_unique<Outflow>(BoundaryLocation::RIGHT);
    boundaries[2] = std::make_unique<Outflow>(BoundaryLocation::FRONT);
    boundaries[3] = std::make_unique<Outflow>(BoundaryLocation::BACK);
    boundaries[4] = std::make_unique<Outflow>(BoundaryLocation::UPPER);
    boundaries[5] = std::make_unique<Outflow>(BoundaryLocation::LOWER);
}

Domain::Domain(Domain&& other) noexcept : dimension(other.dimension), boundaries(std::move(other.boundaries)) {}

Domain::Domain(R3 dimension, std::array<std::unique_ptr<Boundary>, 6> boundaries)
    : dimension(dimension), boundaries(std::move(boundaries)) {}

R3 Domain::getDimension() const { return dimension; }

Boundary* Domain::getBoundary(BoundaryLocation location) { 
    return boundaries[locationToIndex(location)].get(); 
}

const Boundary* Domain::getBoundary(BoundaryLocation location) const {
    return boundaries[locationToIndex(location)].get();
}

void Domain::applyBoundary(Particle& p) {
    for (auto& boundary : boundaries) {
        if (boundary) {
            boundary->applyBoundary(p);
        }
    }
}

std::vector<Particle> Domain::computeGhostParticles(const Particle& p) const {
    std::vector<Particle> ghosts;
    for (const auto& boundary : boundaries) {
        if (boundary) {
            if (auto ghost = boundary->computeGhostParticle(p)) {
                ghosts.push_back(std::move(ghost.value()));
            }
        }
    }
    return ghosts;
}

Domain& Domain::operator=(Domain&& other) noexcept {
    if (this != &other) {
        dimension = other.dimension;
        boundaries = std::move(other.boundaries);
    }
    return *this;
}

}  // namespace mol_sim
