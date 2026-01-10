#ifndef HARMONIC_FORCE_H
#define HARMONIC_FORCE_H

#include "particles/Particle.h"
#include "physics/singleforces/SingleForceSource.h"
#include "utils/Vector.h"

namespace mol_sim {

/**
 * @brief Class to calculate gravitational forces between two different objects. Implements the ForceSource concept.
 *
 * Class to calculate gravitational forces between two different objects.
 * Implements the ForceSource concept.
 */
class HarmonicForce : public SingleForceSource {
   public:
    /**
     * @brief Calculates the gravitational force a particle p2 exerts on a different particle p1.
     *
     * @param p1 Particle whose force is to be calculated.
     * @param p2 Particle which exerts force on p1.
     * @return Vector<double, 3> Force exerted on p1 by p2.
     *
     * Calculates the gravitational force a particle p2 exerts on a different particle p1.
     */
    [[nodiscard]] Vector<double, 3> applyForce(const Particle& p1) const noexcept override { return {0., 0., 0.}; }
};
}  // namespace mol_sim

#endif
