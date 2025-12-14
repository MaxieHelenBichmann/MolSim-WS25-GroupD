#include "physics/LennardJonesForce.h"

#include "utils/Vector.h"

using namespace mol_sim;

Vector<double, 3> LennardJonesForce::applyForce(const Particle& p1, const Particle& p2) const noexcept {
    const double sigma = (p1.getSigma() + p2.getSigma()) / 2;
    const double epsilon = std::sqrt(p1.getEpsilon() * p2.getEpsilon());
    const auto difference = p1.getX() - p2.getX();
    const double distance = difference.euclidNorm();

    // Guard against very small distances that would cause overflow
    constexpr double min_distance = 1e-10;
    if (distance < min_distance) {
        return {0.0, 0.0, 0.0};
    }

    const double sigma_div_distance = sigma / distance;
    const double sigma_div_distance_pow_6 = sigma_div_distance * sigma_div_distance * sigma_div_distance *
                                            sigma_div_distance * sigma_div_distance * sigma_div_distance;

    const auto llj_force = (((-24.0 * epsilon) / (distance * distance)) *
                            (sigma_div_distance_pow_6 - (2.0 * sigma_div_distance_pow_6 * sigma_div_distance_pow_6))) *
                           difference;
    return llj_force;
}
