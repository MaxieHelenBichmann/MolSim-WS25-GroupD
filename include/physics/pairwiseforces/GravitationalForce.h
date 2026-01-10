#ifndef PAIR_GRAV_FORCE_H
#define PAIR_GRAV_FORCE_H

#include "particles/Particle.h"
#include "physics/pairwiseforces/PairwiseForceSource.h"
#include "utils/Vector.h"

namespace mol_sim {

/**
 * @brief Class to calculate gravitational forces between two different objects. Implements the ForceSource concept.
 *
 * Class to calculate gravitational forces between two different objects.
 * Implements the ForceSource concept.
 */
class GravitationalForce : public PairwiseForceSource {
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
    [[nodiscard]] Vector<double, 3> applyForce(const Particle& p1, const Particle& p2) const noexcept override {
        const auto difference = p2.getX() - p1.getX();
        const double distance = difference.euclidNorm();
        if (distance == 0) {
            return {0.0, 0.0, 0.0};
        }
        const auto grav_force = p1.getM() * p2.getM() / (distance * distance * distance) * difference;
        return grav_force;
    }
};
}  // namespace mol_sim

#endif
