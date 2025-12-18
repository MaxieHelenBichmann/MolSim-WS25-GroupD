#include "physics/LennardJonesForce.h"

#include "utils/Vector.h"

using namespace mol_sim;

Vector<double, 3> LennardJonesForce::applyForce(const Particle& p1, const Particle& p2) const noexcept {
    const double sigma = (p1.getSigma() + p2.getSigma()) / 2;
    const double epsilon = std::sqrt(p1.getEpsilon() * p2.getEpsilon());
    const auto difference = p1.getX() - p2.getX();
    const double sqr_distance = difference.sqrEuclidNorm();

    // Guard against very small distances that would cause overflow (TODO still correct with squared distance?)
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
