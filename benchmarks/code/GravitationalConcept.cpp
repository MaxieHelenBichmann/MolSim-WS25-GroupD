#include "GravitationalConcept.h"

#include "utils/Vector.h"

using namespace mol_sim;

Vector<double, 3> GravitationalConcept::applyForce(const Particle& p1, const Particle& p2) const {
    const auto difference = p2.getX() - p1.getX();
    const double distance = difference.euclidNorm();
    if (distance == 0) {
        return {0.0, 0.0, 0.0};
    }
    const auto grav_force = p1.getM() * p2.getM() / std::pow(distance, 3) * difference;
    return grav_force;
}
