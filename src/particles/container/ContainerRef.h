#pragma once

#include <variant>

#include "SimpleContainer.h"

namespace particle_containers {

// only used for IO operations, not (!) in performance relevant code segments, then use template!!!
// has no ownership of container elements!!!
class ContainerRef {
  // more containers can be added in variant
  std::variant<SimpleContainer*> _instance;

 public:
  ContainerRef(SimpleContainer& c);
  operator SimpleContainer() const;

  // retrieve data
  Particle& operator[](size_t idx);
  const Particle& operator[](size_t idx) const;

  size_t size() const;
  bool empty() const;

  // modify
  void clear();
  void reserve(size_t n);

  void addParticle(Particle&& value);
  void addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg);
  void addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg, int type);

  // iterators
  std::vector<Particle>::iterator begin();
  std::vector<Particle>::const_iterator begin() const;
  std::vector<Particle>::const_iterator cbegin() const;
  std::vector<Particle>::iterator end();
  std::vector<Particle>::const_iterator end() const;
  std::vector<Particle>::const_iterator cend() const;
};
}  // namespace particle_containers