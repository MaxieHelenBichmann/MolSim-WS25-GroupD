#ifndef FORCE_SOURCE
#define FORCE_SOURCE

#include <array>

#include "particles/Particle.h"
/**
 * @brief Interface for force source classes
 *
 * Interface for force source classes
 * These classes all implement the calculateForce(p1, p2) function
 *
 */
class ForceSource {
 public:
  /**
   * @brief Virtual destructor for cleanup of derived classes
   * Virtual destructor for cleanup of derived classes
   */
  virtual ~ForceSource() = default;
  /**
   * @brief Calculates the force that a particle p1 exerts on a different particle p2
   *
   * @param p1 Particle whos force is to be calculated
   * @param p2 Particle which exerts force on p1
   * @return std::array<double, 3> Force exerted on p1 by p2
   *
   * Calculates the force that a particle p1 exerts on a different particle p2
   */
  [[nodiscard]] virtual std::array<double, 3> calculateForce(const Particle& p1, const Particle& p2) const = 0;
};

#endif
