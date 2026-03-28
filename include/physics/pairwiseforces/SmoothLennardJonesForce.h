#ifndef SMOOTH_LJ_FORCE_H
#define SMOOTH_LJ_FORCE_H

#include <limits>

#include "particles/Particle.h"
#include "physics/pairwiseforces/PairwiseForceSource.h"
#include "utils/Vector.h"

namespace mol_sim {

/**
 * @brief Class to calculate smooth Lennard-Jones forces between two different objects.
 *
 * A PairwiseForceSource that calculates the smooth Lennard-Jones force between two particles.
 * It smoothens the force to zero between a smoothing radius rl and a cutoff radius rc.
 */
class SmoothLennardJonesForce : public PairwiseForceSource {
    /**
     * @brief Cutoff radius, beyond which no force is applied.
     */
    double rc = std::numeric_limits<double>::infinity();
    /**
     * @brief Smoothing radius, beyond which the force is smoothly reduced to zero.
     */
    double rl = std::numeric_limits<double>::infinity();

   public:
    /**
     * @brief Set the cutoff and smoothing radius for the force calculation.
     */
    void initForce(double rc, double rl) noexcept {
        this->rc = rc;
        this->rl = rl;
    };

    /**
     * @brief Calculates the smooth Lennard-Jones force a particle p2 exerts on a different particle p1.
     *
     * @param p1 Particle whose force is to be calculated.
     * @param p2 Particle which exerts force on p1.
     * @return Vector<double, 3> Force exerted on p1 by p2.
     *
     * Calculates the smooth Lennard-Jones force a particle p2 exerts on a different particle p1.
     */
    [[nodiscard]] Vector<double, 3> applyForce(const Particle& p1, const Particle& p2) const noexcept override {
        const auto difference = p1.getX() - p2.getX();
        const double sqr_distance = difference.sqrEuclidNorm();

        if (sqr_distance >= rc * rc || sqr_distance < 1e-20) {
            return {0.0, 0.0, 0.0};
        }

        const double sigma = (p1.getSigma() + p2.getSigma()) / 2;
        const double epsilon = std::sqrt(p1.getEpsilon() * p2.getEpsilon());

        if (sqr_distance <= rl * rl) {
            const double inv_dist_sq = 1.0 / sqr_distance;
            const double sigma_div_distance_pow_2 = sigma * sigma * inv_dist_sq;
            const double sigma_div_distance_pow_6 =
                sigma_div_distance_pow_2 * sigma_div_distance_pow_2 * sigma_div_distance_pow_2;

            return ((-24.0 * epsilon * inv_dist_sq) *
                    (sigma_div_distance_pow_6 - (2.0 * sigma_div_distance_pow_6 * sigma_div_distance_pow_6))) *
                   difference;
        }

        const double sigma_pow_2 = sigma * sigma;
        const double sigma_pow_6 = sigma_pow_2 * sigma_pow_2 * sigma_pow_2;

        const double distance = std::sqrt(sqr_distance);
        const double inv_dist_sq = 1.0 / sqr_distance;
        const double inv_dist_pow_6 = inv_dist_sq * inv_dist_sq * inv_dist_sq;
        const double inv_dist_pow_14 = inv_dist_pow_6 * inv_dist_pow_6 * inv_dist_sq;
        const double distance_pow_6 = sqr_distance * sqr_distance * sqr_distance;
        const double sigma_minus_distance = (2 * sigma_pow_6) - distance_pow_6;
        const double inv_rc_rl = 1 / (rc - rl);

        const double poly_term =
            (rc * (rc * sigma_minus_distance - ((3 * rl) - distance) * sigma_minus_distance)) +
            (distance * (rl * (5 * sigma_pow_6 - 2 * distance_pow_6) + distance * (-3 * sigma_pow_6 + distance_pow_6)));

        return difference * (24.0 * sigma_pow_6 * epsilon * inv_dist_pow_14 * inv_rc_rl * inv_rc_rl * inv_rc_rl *
                             (rc - distance) * poly_term);
    }
};

}  // namespace mol_sim

#endif
