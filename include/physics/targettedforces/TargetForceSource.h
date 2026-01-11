#ifndef TARGET_FORCE_SOURCE_H
#define TARGET_FORCE_SOURCE_H

#include "particles/Particle.h"
#include "utils/Vector.h"

namespace mol_sim {
/**
 * @brief Abstract class of a ForceSource.
 * Each Force must implement this applyForce Method, which calculates the force acting on a single particle.
 */
class TargetForceSource {
   private:
    R3 direction;           // Pull direction (e.g., {0, 1, 0} for upward)
    double magnitude;       // Force magnitude
    size_t max_iterations;  // Apply for first N iterations
    size_t curr_iteration{};

   public:
    TargetForceSource(R3 dir, double mag, size_t max_iter) : direction(dir), magnitude(mag), max_iterations(max_iter) {}

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
};

}  // namespace mol_sim

#endif
