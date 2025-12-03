#ifndef LJ_FORCE_H
#define LJ_FORCE_H

#include "particles/Particle.h"
#include "physics/ForceSource.h"
#include "utils/Vector.h"

namespace mol_sim {

/**
 * @brief Class to calculate Lennard-Jones forces between two different objects. Implements the ForceSource Concept.
 *
 * Class to calculate Lennard-Jones forces between two different objects.
 * Implements the ForceSource concept.
 */
class LennardJonesForce : public ForceSource {
   public:
    /**
     * @brief Calculates the Lennard-Jones force a particle p2 exerts on a different particle p1.
     *
     * @param p1 Particle whos force is to be calculated.
     * @param p2 Particle which exerts force on p1.
     * @return Vector<double, 3> Force exerted on p1 by p2.
     *
     * Calculates the Lennard-Jones force a particle p2 exerts on a different particle p1.
     */
    [[nodiscard]] Vector<double, 3> applyForce(const Particle& p1, const Particle& p2) const noexcept override;
};

}  // namespace mol_sim

#endif
