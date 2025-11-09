#ifndef FORCE_SOURCE_H
#define FORCE_SOURCE_H

#include <cstdint>

#include "particles/Particle.h"
#include "utils/Vector.h"

namespace mol_sim {

/**
 * @brief Interface for force source classes.
 *
 * Interface for force source classes.
 * These classes all implement the calculateForce(p1, p2) function.
 *
 */
class ForceSource {
   public:
    /**
     * @brief Virtual destructor for cleanup of derived classes.
     * Virtual destructor for cleanup of derived classes.
     */
    virtual ~ForceSource() = default;
    /**
     * @brief Calculates the force that a particle p1 exerts on a different particle p2.
     *
     * @param p1 Particle whos force is to be calculated.
     * @param p2 Particle which exerts force on p1.
     * @return Vector<double, 3> Force exerted on p1 by p2.
     *
     * Calculates the force that a particle p1 exerts on a different particle p2.
     */
    [[nodiscard]] virtual Vector<double, 3> calculateForce(const Particle& p1, const Particle& p2) const = 0;
};

/**
 * @brief Enum of all available Force Sources.
 *
 */
enum Force : std::uint8_t { GRAVITATIONAL };

}  // namespace mol_sim

#endif
