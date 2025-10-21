#pragma once

#include <vector>

#include "particles/Particle.h"

namespace particle_containers {

class ParticleContainer {
 public:
  virtual ~ParticleContainer() = default;  // ?

  // retrieve data
  virtual Particle& operator[](size_t idx) = 0;
  virtual const Particle& operator[](size_t idx) const = 0;

  virtual std::vector<Particle>& data() = 0;
  virtual const std::vector<Particle>& data() const = 0;

  virtual bool empty() const = 0;
  virtual size_t size() const = 0;

  // modify
  virtual void clear() = 0;
  virtual void reserve(size_t n) = 0;
  virtual void addParticle(Particle&& value) = 0;
  virtual void addParticle(const Particle& value) = 0;
  virtual void addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg) = 0;
  virtual void addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg, int type) = 0;

  // iterators
  virtual std::vector<Particle>::iterator begin() = 0;
  virtual std::vector<Particle>::const_iterator begin() const = 0;
  virtual std::vector<Particle>::iterator end() = 0;
  virtual std::vector<Particle>::const_iterator end() const = 0;
};

}  // namespace particle_containers