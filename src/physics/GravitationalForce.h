#ifndef GRAV_FORCE
#define GRAV_FORCE

#include "physics/ForceSource.h"
#include "utils/ArrayUtils.h"

/**
 * @brief Class to calculate gravitational forces between two different objects. Implements the ForceSource interface
 *
 * Class to calculate gravitational forces between two different objects.
 * Implements the ForceSource interface
 */
class GravitationalForce : public ForceSource {
 public:
  /**
   * @brief Calculates the gravitational force a particle p2 exerts on a different particle p1
   *
   * @param p1 Particle whos force is to be calculated
   * @param p2 Particle which exerts force on p1
   * @return std::array<double, 3> Force exerted on p1 by p2
   *
   * Calculates the gravitational force a particle p2 exerts on a different particle p1
   */
  [[nodiscard]] std::array<double, 3> calculateForce(const Particle& p1, const Particle& p2) const override;
};

#endif
