#ifndef HARMONIC_FORCE_H
#define HARMONIC_FORCE_H

#include "particles/Particle.h"
#include "physics/singleforces/SingleForceSource.h"
#include "utils/Vector.h"

namespace mol_sim {

/**
 * @brief Class to calculate harmonic spring forces for membrane particles.
 *
 * Calculates forces based on harmonic potential between neighboring particles in a membrane structure.
 * Only applies to membrane particles (types 2 and 4).
 */
class HarmonicForce : public SingleForceSource {
   private:
    /**
     * @brief Stiffness constant of the harmonic spring.
     */
    const double K;
    /**
     * @brief Equilibrium bond length between direct neighbors.
     */
    const double R_0;

   public:
    /**
     * @brief Construct a new Harmonic Force calculator.
     *
     * @param k Stiffness constant.
     * @param r_0 Equilibrium bond length.
     */
    HarmonicForce(double k, double r_0) : K(k), R_0(r_0) {}
    /**
     * @brief Calculates the harmonic spring force acting on a membrane particle.
     *
     * @param p1 Particle whose force is to be calculated.
     * @return Vector<double, 3> Total harmonic force from all neighbors.
     *
     * Calculates forces from direct neighbors using distance R_0 and diagonal neighbors using sqrt(2)*R_0.
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
