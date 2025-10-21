#pragma once

#include <vector>

#include "particles/Particle.h"

namespace particle_containers {

class ParticleContainer {
 public:
  virtual ~ParticleContainer() = default;

  // retrieve data
  virtual Particle& operator[](size_t idx);
  virtual const Particle& operator[](size_t idx) const;

  virtual std::vector<Particle>& data();
  virtual const std::vector<Particle>& data() const;

  virtual bool empty() const;
  virtual size_t size() const;

  virtual bool operator==(const ParticleContainer& other) const;

  // modify
  virtual void clear();
  virtual void reserve(size_t n);
  virtual void removeParticle(size_t idx);
  virtual void removeParticle(Particle& value);
  virtual void addParticle(Particle&& value);
  virtual void addParticle(const Particle& value);

  // iterators
  virtual std::vector<Particle>::iterator begin();
  virtual std::vector<Particle>::const_iterator begin() const;
  virtual std::vector<Particle>::iterator end();
  virtual std::vector<Particle>::const_iterator end() const;
};

}  // namespace particle_containers