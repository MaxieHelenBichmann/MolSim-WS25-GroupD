#pragma once

#include <vector>

#include "../ParticleContainer.h"

namespace particle_containers {

class SimpleContainer : ParticleContainer {
  std::vector<Particle> _data;

 public:
  // constructors
  SimpleContainer();
  SimpleContainer(size_t init_cap);
  SimpleContainer(std::initializer_list<Particle> init);

  // cpy constr
  SimpleContainer(const SimpleContainer& other);
  // cpy assignment
  SimpleContainer& operator=(const SimpleContainer& other);
  // mv constr
  SimpleContainer(SimpleContainer&& other);
  // mv assignment
  SimpleContainer& operator=(SimpleContainer&& other);
  // dstr
  ~SimpleContainer();

  // retrieve data
  Particle& operator[](size_t idx);
  const Particle& operator[](size_t idx) const;

  std::vector<Particle>& data();
  const std::vector<Particle>& data() const;

  bool empty() const;
  size_t size() const;

  bool operator==(const SimpleContainer& other) const;

  // modify
  void clear();
  void reserve(size_t n);
  void removeParticle(size_t idx);
  void removeParticle(Particle& value);
  void addParticle(Particle&& value);
  void addParticle(const Particle& value);

  // iterators
  std::vector<Particle>::iterator begin();
  std::vector<Particle>::const_iterator begin() const;
  std::vector<Particle>::iterator end();
  std::vector<Particle>::const_iterator end() const;
};
}  // namespace particle_containers