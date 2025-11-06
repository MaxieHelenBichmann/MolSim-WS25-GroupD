#ifndef LJ_FORCE_H
#define LJ_FORCE_H

#include "physics/ForceSource.h"
#include "utils/Vector.h"

namespace mol_sim {

/**
 * @brief Class to calculate Lennard-Jones forces between two different objects. Implements the ForceSource interface.
 *
 * Class to calculate Lennard-Jones forces between two different objects.
 * Implements the ForceSource interface.
 */
class LennardJonesForce : public ForceSource {
   private:
    double epsilon;
    double sigma;

   public:
    /**
     * @brief Construct a new LennardJonesForce object, setting all necessary parameters for generation.
     *
     * @param epsilon Depth of the potential well
     * @param sigma Distance where Lennard-Jones force between two particles is 0
     */
    LennardJonesForce(double epsilon, double sigma) : epsilon(epsilon), sigma(sigma) {};
    /**
     * @brief Calculates the Lennard-Jones force a particle p2 exerts on a different particle p1.
     *
     * @param p1 Particle whos force is to be calculated.
     * @param p2 Particle which exerts force on p1.
     * @return Vector<double, 3> Force exerted on p1 by p2.
     *
     * Calculates the Lennard-Jones force a particle p2 exerts on a different particle p1.
     */
    [[nodiscard]] Vector<double, 3> applyForce(const Particle& p1, const Particle& p2) const override;
};

}  // namespace mol_sim

#endif
