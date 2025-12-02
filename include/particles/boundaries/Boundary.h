#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <cstdint>
#include <optional>
#include <string>

#include "particles/Particle.h"

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
 */
enum class BoundaryType : std::uint8_t { OUTFLOW, REFLECTING, VELOCITYREFLECT };

class Boundary {
   protected:
    BoundaryLocation location;
    BoundaryType type;

   public:
    Boundary(BoundaryLocation location, BoundaryType type) : location(location), type(type) {}
    virtual ~Boundary() = default;

    /**
     * @brief Applies the boundary condition to a particle (e.g., velocity reflection).
     * Computes a ghost particle if this boundary requires one for the given particle.
     * @param p The particle to check.
     * @return The ghost particle if needed, nullopt otherwise.
     */
    [[nodiscard]] virtual std::optional<Particle> applyBoundary(Particle& p) = 0;
    [[nodiscard]] BoundaryType& getType() { return type; }
    [[nodiscard]] BoundaryLocation& getLocation() { return location; }
    [[nodiscard]] const BoundaryType& getType() const { return type; }
    [[nodiscard]] const BoundaryLocation& getLocation() const { return location; }

   protected:
    [[nodiscard]] size_t getAxis() const {
        switch (location) {
            case BoundaryLocation::LEFT:
            case BoundaryLocation::RIGHT:
                return 0;
            case BoundaryLocation::UPPER:
            case BoundaryLocation::LOWER:
                return 2;
            case BoundaryLocation::FRONT:
            case BoundaryLocation::BACK:
                return 1;
            default:
                return 0;
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
