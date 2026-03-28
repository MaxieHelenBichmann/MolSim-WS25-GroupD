#ifndef TRUNC_LJ_FORCE_H
#define TRUNC_LJ_FORCE_H

#include "particles/Particle.h"
#include "physics/pairwiseforces/PairwiseForceSource.h"
#include "utils/Vector.h"

namespace mol_sim {

/**
 * @brief Truncated Lennard-Jones force for membrane self-penetration prevention.
 *
 * This force is truncated at r_cutoff = 2^(1/6) * sigma, which is the minimum of the LJ potential.
 * Only the repulsive part of the Lennard-Jones potential is active.
 */
class TruncLennardJonesForce : public PairwiseForceSource {
   private:
    /**
     * @brief For squared distance comparison: r_cutoff^2 = (2^(1/6))^2 * sigma^2 = 2^(1/3) * sigma^2
     *
     */
    static constexpr double CUTOFF_FACTOR_SQUARED = 1.2599210498948731906665443602832965552806854248046875;  // 2^(1/3)

   public:
    /**
     * @brief Calculates the truncated Lennard-Jones force a particle p2 exerts on a different particle p1.
     *
     * @param p1 Particle whose force is to be calculated.
     * @param p2 Particle which exerts force on p1.
     * @return Vector<double, 3> Force exerted on p1 by p2 (zero if beyond cutoff or not membrane particles).
     */
    [[nodiscard]] Vector<double, 3> applyForce(const Particle& p1, const Particle& p2) const noexcept override {
        // Only apply to membrane particles (type 2)
        if ((p1.getType() != 2 && p1.getType() != 4) || (p2.getType() != 2 && p2.getType() != 4)) {
            return {0.0, 0.0, 0.0};
        }

        const double sigma = (p1.getSigma() + p2.getSigma()) / 2;
        const double epsilon = std::sqrt(p1.getEpsilon() * p2.getEpsilon());
        const auto difference = p1.getX() - p2.getX();
        const double sqr_distance = difference.sqrEuclidNorm();

        // Guard against very small distances
        constexpr double sqr_min_distance = 1e-20;
        if (sqr_distance < sqr_min_distance) {
            return {0.0, 0.0, 0.0};
        }

        // Cutoff check using squared distance: sqr_cutoff = 2^(1/3) * sigma^2
        const double sqr_cutoff = CUTOFF_FACTOR_SQUARED * sigma * sigma;
        if (sqr_distance >= sqr_cutoff) {
            return {0.0, 0.0, 0.0};
        }

        // Standard optimized LJ force calculation
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
