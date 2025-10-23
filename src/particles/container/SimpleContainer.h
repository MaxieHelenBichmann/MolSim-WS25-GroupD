#pragma once

#include <vector>

#include "../ParticleContainer.h"

namespace particle_containers {

class SimpleContainer : public std::vector<Particle> {
 public:
  // constructors
  using std::vector<Particle>::vector;

  // modify
  void addParticle(Particle&& value);
  void addParticle(const Particle& value);
  void addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg);
  void addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg, int type);
};
}  // namespace particle_containers