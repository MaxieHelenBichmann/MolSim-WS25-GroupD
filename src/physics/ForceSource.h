#ifndef FORCE_SOURCE
#define FORCE_SOURCE

#include <array>

#include "particles/Particle.h"

class ForceSource {
 public:
  virtual ~ForceSource() = default;

  [[nodiscard]] virtual std::array<double, 3> calculateForce(const Particle& p1, const Particle& p2) const = 0;
};

#endif
