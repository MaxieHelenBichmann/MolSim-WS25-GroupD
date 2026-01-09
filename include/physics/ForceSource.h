#ifndef FORCE_SOURCE_H
#define FORCE_SOURCE_H

#include <cstdint>

#include "particles/Particle.h"
#include "utils/Vector.h"

namespace mol_sim {
/**
 * @brief Abstract class of a ForceSource.
 * Each Force must implement this applyForce Method, which calculates the force between two particles.
 */
class ForceSource {
   public:
    [[nodiscard]] virtual Vector<double, 3> applyForce(const Particle& p1, const Particle& p2) const noexcept = 0;
    virtual ~ForceSource() = default;
};
/**
 * @brief Enum of all available Force Sources.
 *
 */
enum Force : std::uint8_t { GRAVITATIONAL, LENNARDJONES, S_LENNARDJONES };

}  // namespace mol_sim

#endif
