#ifndef LJ_FORCE_H
#define LJ_FORCE_H

#include "particles/Particle.h"
#include "physics/pairwiseforces/PairwiseForceSource.h"
#include "utils/Vector.h"

namespace mol_sim {

/**
 * @brief Class to calculate Lennard-Jones forces between two different objects.
 *
 * A PairwiseForceSource that calculates the Lennard-Jones force between two particles.
 */
class LennardJonesForce : public PairwiseForceSource {
   public:
    /**
     * @brief Calculates the Lennard-Jones force a particle p2 exerts on a different particle p1.
     *
     * @param p1 Particle whose force is to be calculated.
     * @param p2 Particle which exerts force on p1.
     * @return Vector<double, 3> Force exerted on p1 by p2.
     *
     * Calculates the Lennard-Jones force a particle p2 exerts on a different particle p1.
     */
    [[nodiscard]] Vector<double, 3> applyForce(const Particle& p1, const Particle& p2) const noexcept override {
        const double sigma = (p1.getSigma() + p2.getSigma()) / 2;
        const double epsilon = std::sqrt(p1.getEpsilon() * p2.getEpsilon());
        const auto difference = p1.getX() - p2.getX();
        const double sqr_distance = difference.sqrEuclidNorm();

        // Guard against very small distances that would cause overflow
        constexpr double sqr_min_distance = 1e-20;
        if (sqr_distance < sqr_min_distance) {
            return {0.0, 0.0, 0.0};
        }

        const double inv_dist_sq = 1.0 / sqr_distance;
        const double sigma_div_distance_pow_2 = sigma * sigma * inv_dist_sq;
        const double sigma_div_distance_pow_6 =
            sigma_div_distance_pow_2 * sigma_div_distance_pow_2 * sigma_div_distance_pow_2;

        return ((-24.0 * epsilon * inv_dist_sq) *
                (sigma_div_distance_pow_6 - (2.0 * sigma_div_distance_pow_6 * sigma_div_distance_pow_6))) *
               difference;
    }
};

}  // namespace mol_sim

#endif
