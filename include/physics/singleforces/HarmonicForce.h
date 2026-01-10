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
   private:
    // stiffness constant
    const double K;
    // average bond length of molecule
    const double R_0;

   public:
    HarmonicForce(double k, double r_0) : K(k), R_0(r_0) {}
    /**
     * @brief Calculates the gravitational force a particle p2 exerts on a different particle p1.
     *
     * @param p1 Particle whose force is to be calculated.
     * @param p2 Particle which exerts force on p1.
     * @return Vector<double, 3> Force exerted on p1 by p2.
     *
     * Calculates the gravitational force a particle p2 exerts on a different particle p1.
     */
    [[nodiscard]] Vector<double, 3> applyForce(const Particle& p1) const noexcept override {
        Vector<double, 3> force = {0., 0., 0.};
        if (p1.getType() != 2 && p1.getType() != 4) {
            return force;
        }
        // compute direct neighbor influence
        for (size_t i = 0; i < 4; i++) {
            if ((p1.getNeighbors()[i]) == nullptr) {
                continue;
            }
            double dist = (p1.getX() - p1.getNeighbors()[i]->getX()).euclidNorm();
            double scalar = (K * 0.5 * (dist - R_0)) / dist;
            force += scalar * (p1.getNeighbors()[i]->getX() - p1.getX());
        }
        // compute diagonal neighbor influence
        for (size_t i = 4; i < 8; i++) {
            if ((p1.getNeighbors()[i]) == nullptr) {
                continue;
            }
            double dist = (p1.getX() - p1.getNeighbors()[i]->getX()).euclidNorm();
            double scalar = (K * 0.5 * (dist - std::numbers::sqrt2 * R_0)) / dist;
            force += scalar * (p1.getNeighbors()[i]->getX() - p1.getX());
        }
        return force;
    }
};
}  // namespace mol_sim

#endif
