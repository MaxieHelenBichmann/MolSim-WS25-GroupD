#ifndef GRAV_FORCE_H
#define GRAV_FORCE_H

#include "physics/ForceSource.h"
#include "utils/Vector.h"

namespace mol_sim {

/**
 * @brief Class to calculate gravitational forces between two different objects. Implements the ForceSource interface.
 *
 * Class to calculate gravitational forces between two different objects.
 * Implements the ForceSource interface.
 */
class GravitationalForce : public ForceSource {
   public:
    /**
     * @brief Calculates the gravitational force a particle p2 exerts on a different particle p1.
     *
     * @param p1 Particle whos force is to be calculated.
     * @param p2 Particle which exerts force on p1.
     * @return Vector<double, 3> Force exerted on p1 by p2.
     *
     * Calculates the gravitational force a particle p2 exerts on a different particle p1.
     */
    [[nodiscard]] Vector<double, 3> calculateForce(const Particle& p1, const Particle& p2) const override;
};

}  // namespace mol_sim

#endif
