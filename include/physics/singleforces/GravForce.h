#ifndef SINGLE_GRAV_FORCE_H
#define SINGLE_GRAV_FORCE_H

#include "particles/Particle.h"
#include "physics/singleforces/SingleForceSource.h"
#include "utils/Vector.h"

namespace mol_sim {

/**
 * @brief Class to calculate constant gravitational force on particles.
 *
 * Applies a uniform gravitational field to all particles based on mass.
 */
class GravForce : public SingleForceSource {
   private:
    /**
     * @brief Gravitational acceleration vector (e.g., {0, -9.81, 0} for downward gravity).
     */
    const R3 G_GRAV;

   public:
    /**
     * @brief Construct a new Gravitational Force calculator.
     *
     * @param g_grav Gravitational acceleration vector.
     */
    GravForce(R3 g_grav) : G_GRAV(g_grav) {}

    /**
     * @brief Calculates the gravitational force acting on a particle.
     *
     * @param p1 Particle whose force is to be calculated.
     * @return Vector<double, 3> Gravitational force (mass * gravitational acceleration).
     */
    [[nodiscard]] Vector<double, 3> applyForce(const Particle& p1) const noexcept override {
        return p1.getM() * G_GRAV;
    }
    [[nodiscard]] SingleForce getType() const override { return GRAV; }
};
}  // namespace mol_sim

#endif
