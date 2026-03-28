#include "particles/container/domain/Domain.h"

#include <spdlog/spdlog.h>

#include "exceptions/BoundaryException.h"
#include "particles/boundaries/Outflow.h"
#include "particles/boundaries/Periodic.h"
#include "physics/pairwiseforces/LennardJonesForce.h"

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
            SPDLOG_ERROR("Unknown boundary location: {}", static_cast<int>(location));
            throw BoundaryException("Unknown boundary location");
    }
}

Domain::Domain() {
    boundaries[0] = std::make_unique<Outflow>(BoundaryLocation::LEFT, dimension);
    boundaries[1] = std::make_unique<Outflow>(BoundaryLocation::RIGHT, dimension);
    boundaries[2] = std::make_unique<Outflow>(BoundaryLocation::FRONT, dimension);
    boundaries[3] = std::make_unique<Outflow>(BoundaryLocation::BACK, dimension);
    boundaries[4] = std::make_unique<Outflow>(BoundaryLocation::UPPER, dimension);
    boundaries[5] = std::make_unique<Outflow>(BoundaryLocation::LOWER, dimension);
}

Domain::Domain(R3 dimension) : dimension(dimension) {
    boundaries[0] = std::make_unique<Outflow>(BoundaryLocation::LEFT, dimension);
    boundaries[1] = std::make_unique<Outflow>(BoundaryLocation::RIGHT, dimension);
    boundaries[2] = std::make_unique<Outflow>(BoundaryLocation::FRONT, dimension);
    boundaries[3] = std::make_unique<Outflow>(BoundaryLocation::BACK, dimension);
    boundaries[4] = std::make_unique<Outflow>(BoundaryLocation::UPPER, dimension);
    boundaries[5] = std::make_unique<Outflow>(BoundaryLocation::LOWER, dimension);
}

Domain::Domain(Domain&& other) noexcept : dimension(other.dimension), boundaries(std::move(other.boundaries)) {}

Domain::Domain(R3 dimension, std::array<std::unique_ptr<Boundary>, 6> boundaries)
    : dimension(dimension), boundaries(std::move(boundaries)) {}

Boundary& Domain::getBoundary(BoundaryLocation location) {
    auto& boundary = boundaries[locationToIndex(location)];
    if (!boundary) {
        SPDLOG_ERROR("Boundary at location {} is null", static_cast<int>(location));
        throw BoundaryException("Boundary at given location is null");
    }
    return *boundary;
}

const Boundary& Domain::getBoundary(BoundaryLocation location) const {
    const auto& boundary = boundaries[locationToIndex(location)];
    if (!boundary) {
        SPDLOG_ERROR("Boundary at location {} is null", static_cast<int>(location));
        throw BoundaryException("Boundary at given location is null");
    }
    return *boundary;
}

void Domain::applyBoundary(Particle& p) const noexcept {  // NOLINT
    const LennardJonesForce force;
    for (const auto& boundary : boundaries) {
        if (boundary) {
            auto new_particles = boundary->applyBoundary(p, force).value_or(std::vector<Particle>());
            particles.insert(particles.end(), new_particles.begin(), new_particles.end());
        }
    }

    // handle edge case (which can potentially cause instabilities relating to periodic boundaries)
    if (typeid(boundaries[0]) == typeid(std::unique_ptr<Periodic>) &&
        typeid(boundaries[1]) == typeid(std::unique_ptr<Periodic>)) {
        boundaries[0]->applyBoundary(p, force);
    }
    if (typeid(boundaries[2]) == typeid(std::unique_ptr<Periodic>) &&
        typeid(boundaries[3]) == typeid(std::unique_ptr<Periodic>)) {
        boundaries[2]->applyBoundary(p, force);
    }
    if (typeid(boundaries[4]) == typeid(std::unique_ptr<Periodic>) &&
        typeid(boundaries[5]) == typeid(std::unique_ptr<Periodic>)) {
        boundaries[4]->applyBoundary(p, force);
    }
}

Domain& Domain::operator=(Domain&& other) noexcept {
    if (this != &other) {
        dimension = other.dimension;
        boundaries = std::move(other.boundaries);
    }
    return *this;
}

}  // namespace mol_sim
