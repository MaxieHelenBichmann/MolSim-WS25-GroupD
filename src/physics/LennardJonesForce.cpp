#include "physics/LennardJonesForce.h"

#include "utils/Vector.h"

using namespace mol_sim;

Vector<double, 3> LennardJonesForce::applyForce(const Particle& p1, const Particle& p2) const {
    const double sigma = (p1.getSigma() + p2.getSigma()) / 2;
    const double epsilon = std::sqrt(p1.getEpsilon() * p2.getEpsilon());
    const auto difference = p1.getX() - p2.getX();
    const double distance = difference.euclidNorm();
    if (distance == 0) {
        return {0.0, 0.0, 0.0};
    }
    const auto llj_force = (((-24.0 * epsilon) / (distance * distance)) *
                            (pow(sigma / distance, 6.0) - (2.0 * pow(sigma / distance, 12.0)))) *
                           -1.0 * difference;
    return llj_force;
}
