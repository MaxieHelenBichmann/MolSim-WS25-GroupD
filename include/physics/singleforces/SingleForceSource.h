#ifndef SINGLE_FORCE_SOURCE_H
#define SINGLE_FORCE_SOURCE_H

#include <cstdint>

#include "particles/Particle.h"
#include "utils/Vector.h"

namespace mol_sim {
/**
 * @brief Abstract class of a ForceSource.
 * Each Force must implement this applyForce Method, which calculates the force acting on a single particle.
 */
class SingleForceSource {
   public:
    [[nodiscard]] virtual Vector<double, 3> applyForce(const Particle& p1) const noexcept = 0;
    virtual ~SingleForceSource() = default;
};
/**
 * @brief Enum of all available Force Sources.
 *
 */
enum SingleForce : std::uint8_t { GRAV, HARMONIC };

}  // namespace mol_sim

#endif
