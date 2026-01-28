#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <spdlog/spdlog.h>

#include <cstdint>
#include <optional>
#include <string>

#include "particles/Particle.h"
#include "physics/pairwiseforces/PairwiseForceSource.h"
#include "utils/Vector.h"

namespace mol_sim {

/**
 * @brief Enum for boundary locations in the Linked-Cell Container.
 *
 * UPPER: +z direction (x-y plane at max z)
 * LOWER: -z direction (x-y plane at min z)
 * FRONT: -y direction (x-z plane at min y)
 * BACK: +y direction (x-z plane at max y)
 * LEFT: -x direction (y-z plane at min x)
 * RIGHT: +x direction (y-z plane at max x)
 */
enum class BoundaryLocation : std::uint8_t { UPPER, LOWER, FRONT, BACK, LEFT, RIGHT };

/**
 * @brief Enum of the supported boundary types.
 *
 * OUTFLOW: Outflow boundary condition: delete particles in halo cells
 * REFLECTING: Reflecting boundary condition: add ghost particles if particle gets too close to boundary
 * VELOCITYREFLECT: Reflecting boundary condition: reflect particle like a ball that flew against a flat surface
 * PERIODIC: Periodic boundary condition: imagine portals (like those in the Portal games) on each boundary.
 */
enum class BoundaryType : std::uint8_t { OUTFLOW, REFLECTING, VELOCITYREFLECT, PERIODIC };

class Boundary {
   protected:
    /**
     * @brief The location of the boundary (UPPER, LOWER, FRONT, BACK, LEFT, RIGHT)
     */
    BoundaryLocation location;
    /**
     * @brief The type of the boundary (OUTFLOW, REFLECTING, VELOCITYREFLECT, PERIODIC)
     */
    BoundaryType type;
    /**
     * @brief The size of the simulation domain (i.e the upper back right corner
     * assuming the lower front left corner is the origin)
     */
    R3 domain_size;

   public:
    Boundary(BoundaryLocation location, BoundaryType type, R3 domain_size) noexcept
        : location(location), type(type), domain_size(domain_size) {}
    virtual ~Boundary() = default;

    /**
     * @brief Applies the boundary condition to a particle (e.g., velocity reflection).
     * Computes a ghost particle if this boundary requires one for the given particle.
     * @param p The particle to check.
     * @param force The force source to use for ghost particle interactions.
     * @return Newly generated particles that need further processing beyond the scope of
     * this / an implementing class.
     */
    virtual void applyBoundary(Particle& p, const PairwiseForceSource& force) noexcept = 0;

    /**
     * @brief Returns the type of the boundary.
     *
     * @return BoundaryType The type of the boundary.
     */
    [[nodiscard]] BoundaryType getType() const noexcept { return type; }
    /**
     * @brief Returns the location of the boundary.
     *
     * @return BoundaryLocation The location of the boundary.
     */
    [[nodiscard]] BoundaryLocation getLocation() const noexcept { return location; }
    /**
     * @brief Returns the size of the domain the boundary belongs to.
     *
     * @return R3 The size of the domain the boundary belongs to.
     */
    [[nodiscard]] R3 getDomainSize() const noexcept { return domain_size; }

   protected:
    /**
     * @brief Gets the axis index for the boundary.
     * @return 0 for LEFT/RIGHT (x), 1 for FRONT/BACK (y), 2 for UPPER/LOWER (z).
     */
    [[nodiscard]] size_t getAxis() const noexcept {
        switch (location) {
            case BoundaryLocation::LEFT:
            case BoundaryLocation::RIGHT:
                return 0;
            case BoundaryLocation::FRONT:
            case BoundaryLocation::BACK:
                return 1;
            case BoundaryLocation::UPPER:
            case BoundaryLocation::LOWER:
                return 2;
            default:
                return 0;
        }
    }

    /**
     * @brief Gets the sign/direction for the boundary.
     * @return -1 for boundaries at min (LEFT, LOWER, FRONT), +1 for boundaries at max (RIGHT, UPPER, BACK).
     */
    [[nodiscard]] int getSign() const noexcept {
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

    /**
     * @brief Computes the boundary position along its axis from domain size.
     * @return 0 for boundaries at min (LEFT, LOWER, FRONT), domain_size[axis] for boundaries at max.
     */
    [[nodiscard]] double getBoundaryPosition() const noexcept {
        size_t axis = getAxis();
        int sign = getSign();
        return (sign < 0) ? 0.0 : domain_size[axis];
    }

    /**
     * @brief Returns the index of a BoundaryLocation.
     *
     * @return Index of the boundary location.
     * LEFT = 0, RIGHT = 1, FRONT = 2, BACK = 3, UPPER = 4, LOWER = 5
     *
     * @note This is needed in Periodic. Do not change this order.
     * We could also make BoundaryLocation into a non-class enum but that's bad practice (I think).
     */
    [[nodiscard]] inline size_t getBoundaryLocationIdx() const noexcept {
        switch (location) {  // could also make BoundaryLocation a non class enum but that's bad practice
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
                SPDLOG_ERROR("Unrecognized boundary location!");
        }
    }
};

/**
 * @brief Parses a string to a BoundaryType enum.
 *
 * @param type_str String representation (e.g., "OUTFLOW", "REFLECTING").
 * @return The corresponding BoundaryType, defaults to OUTFLOW if unknown.
 */
inline BoundaryType parseBoundaryType(const std::string& type_str) {
    if (type_str == "OUTFLOW" || type_str == "outflow" || type_str == "Outflow") {
        return BoundaryType::OUTFLOW;
    }
    if (type_str == "REFLECTING" || type_str == "reflecting" || type_str == "Reflecting") {
        return BoundaryType::REFLECTING;
    }
    if (type_str == "VELOCITYREFLECT" || type_str == "velocityreflect" || type_str == "VelocityReflect") {
        return BoundaryType::VELOCITYREFLECT;
    }
    if (type_str == "PERIODIC" || type_str == "periodic" || type_str == "Periodic") {
        return BoundaryType::PERIODIC;
    }
    return BoundaryType::OUTFLOW;
}

/**
 * @brief Parses a string to a BoundaryLocation enum.
 * @param location_str String representation (e.g., "left", "RIGHT").
 * @return The corresponding BoundaryLocation, defaults to LEFT if unknown.
 */
inline BoundaryLocation parseBoundaryLocation(const std::string& location_str) {
    if (location_str == "LEFT" || location_str == "left") {
        return BoundaryLocation::LEFT;
    }
    if (location_str == "RIGHT" || location_str == "right") {
        return BoundaryLocation::RIGHT;
    }
    if (location_str == "UPPER" || location_str == "upper") {
        return BoundaryLocation::UPPER;
    }
    if (location_str == "LOWER" || location_str == "lower") {
        return BoundaryLocation::LOWER;
    }
    if (location_str == "FRONT" || location_str == "front") {
        return BoundaryLocation::FRONT;
    }
    if (location_str == "BACK" || location_str == "back") {
        return BoundaryLocation::BACK;
    }
    return BoundaryLocation::LEFT;
}

}  // namespace mol_sim

#endif
