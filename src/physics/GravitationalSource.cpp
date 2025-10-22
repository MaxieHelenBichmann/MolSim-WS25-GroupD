#include "GravitationalForce.h"
#include "utils/Vector.h"

Vector<double, 3> GravitationalForce::calculateForce(const Particle& p1, const Particle& p2) const {
  const auto difference = p2.getX() - p1.getX();
  const double distance = difference.euclidNorm();
  const auto grav_force = p1.getM() * p2.getM() / std::pow(distance, 3) * difference;
  return grav_force;
}
