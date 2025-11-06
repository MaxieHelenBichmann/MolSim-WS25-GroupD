#ifndef FORCE_SOURCE_H
#define FORCE_SOURCE_H

#include <concepts>
#include <cstdint>

#include "particles/Particle.h"
#include "utils/Vector.h"

namespace mol_sim {
/**
 * @brief Concept of a ForceSource.
 * fulfills the same function as the abstract Class approach.
 * @tparam C
 */
template <typename C>
concept ForceSource = requires(C c, const Particle& p1, const Particle& p2) {
    { c.applyForce(p1, p2) } -> std::same_as<Vector<double, 3>>;
};
/**
 * @brief Enum of all available Force Sources.
 *
 */
enum Force : std::uint8_t { GRAVITATIONAL, LENNARDJONES };

}  // namespace mol_sim

#endif
