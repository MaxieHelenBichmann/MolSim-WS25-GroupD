#include "physics/LennardJonesForce.h"
#include "utils/Vector.h"

using namespace mol_sim;

Vector<double, 3> LennardJonesForce::calculateForce(const Particle& p1, const Particle& p2) const {
    const auto difference = p2.getX() - p1.getX();
    const double distance = difference.euclidNorm();
    if (distance == 0) {
        return {0.0, 0.0, 0.0};
    }
    const auto llj_force = 
        (((-24 * epsilon) / (distance * distance)) * (pow(sigma/distance, 6) - (2 * pow(sigma/distance, 12)))) * -1 * difference;
    return llj_force;
}
