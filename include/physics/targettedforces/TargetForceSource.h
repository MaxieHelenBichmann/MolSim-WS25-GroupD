#ifndef TARGET_FORCE_SOURCE_H
#define TARGET_FORCE_SOURCE_H

#include "particles/Particle.h"
#include "utils/Vector.h"

namespace mol_sim {
/**
 * @brief Force source that applies a constant force to specific particle types for a limited number of iterations.
 * Used to apply external forces to target particles (types 3 and 4).
 */
class TargetForceSource {
   private:
    /**
     * @brief Direction in which the force is applied (unit vector recommended).
     */
    R3 direction;
    /**
     * @brief Magnitude of the applied force.
     */
    double magnitude;
    /**
     * @brief Maximum number of iterations to apply the force.
     */
    size_t max_iterations;
    /**
     * @brief Current iteration count.
     */
    size_t curr_iteration{};

   public:
    /**
     * @brief Construct a new Target Force Source.
     *
     * @param dir Direction in which the force is applied.
     * @param mag Magnitude of the applied force.
     * @param max_iter Maximum number of iterations to apply the force.
     */
    TargetForceSource(R3 dir, double mag, size_t max_iter) : direction(dir), magnitude(mag), max_iterations(max_iter) {}

    /**
     * @brief Apply force to a particle if it is a target particle (type 3 or 4) and within iteration limit.
     *
     * @param p Particle to apply force to.
     * @return Force vector, or zero if particle is not a target or iteration limit reached.
     */
    [[nodiscard]] Vector<double, 3> applyForce(const Particle& p) noexcept {
        if (p.getType() != 3 && p.getType() != 4) {
            return 0.0;
        }
        if (curr_iteration >= max_iterations) {
            return {0., 0., 0.};
        }
        curr_iteration++;

        return magnitude * direction;
    }

    /**
     * @brief Get the direction of the force.
     * @return Direction vector.
     */
    [[nodiscard]] R3 getDirection() const noexcept { return direction; }
    
    /**
     * @brief Get the magnitude of the force.
     * @return Force magnitude.
     */
    [[nodiscard]] double getMagnitude() const noexcept { return magnitude; }
    
    /**
     * @brief Get the maximum number of iterations.
     * @return Maximum iteration count.
     */
    [[nodiscard]] size_t getMaxIterations() const noexcept { return max_iterations; }
};

}  // namespace mol_sim

#endif
