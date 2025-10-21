#include "GravitationalForce.h"

std::array<double, 3> GravitationalForce::calculateForce(const Particle& p1, const Particle& p2) const {
  const std::array<double, 3> diff = {p2.getX()[0] - p1.getX()[0], p2.getX()[1] - p1.getX()[1],
                                      p2.getX()[2] - p1.getX()[2]};
  double coeff = p1.getM() * p2.getM() / std::pow(ArrayUtils::L2Norm<std::array<double, 3>>(diff), 3);
  return {p1.getF()[0] - coeff * diff[0], p1.getF()[1] - coeff * diff[1], p1.getF()[2] - coeff * diff[2]};
}
