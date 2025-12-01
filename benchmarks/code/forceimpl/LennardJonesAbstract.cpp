#include "LennardJonesAbstract.h"
namespace mol_sim {
Vector<double, 3> lennardJonesAbstract::applyForce(const Particle& p1, const Particle& p2) const {
    const double sigma = (p1.getSigma() + p2.getSigma()) / 2;
    const double epsilon = std::sqrt(p1.getEpsilon() * p2.getEpsilon());
    const auto difference = p2.getX() - p1.getX();
    const double distance = difference.euclidNorm();
    if (distance == 0) {
        return {0.0, 0.0, 0.0};
    }
    const auto llj_force =
        (((-24 * epsilon) / (distance * distance)) * (pow(sigma / distance, 6) - (2 * pow(sigma / distance, 12)))) *
        -1 * difference;
    return llj_force;
}
}  // namespace mol_sim
