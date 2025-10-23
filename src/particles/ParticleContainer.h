#pragma once

#include <concepts>
#include <vector>

#include "particles/Particle.h"

namespace particle_containers {

using R3 = Vector<double, 3>;

template <typename C>
concept ParticleContainer = requires(C c) {
  // retrieve data
  { c[size_t()] } -> std::convertible_to<Particle>;
  { (*static_cast<const C*>(&c))[size_t()] } -> std::convertible_to<const Particle>;

  { c.size() } -> std::convertible_to<size_t>;
  { c.empty() } -> std::convertible_to<bool>;

  // modify
  { c.clear() };
  { c.reserve(size_t()) };

  { c.addParticle(Particle(0)) };
  { c.addParticle(R3(), R3(), 0.) };
  { c.addParticle(R3(), R3(), 0., 0) };

  // iterators
  { c.begin() } -> std::contiguous_iterator;
  { (*static_cast<const C*>(&c)).begin() } -> std::contiguous_iterator;
  { c.cbegin() } -> std::contiguous_iterator;
  { c.end() } -> std::contiguous_iterator;
  { (*static_cast<const C*>(&c)).end() } -> std::contiguous_iterator;
  { c.cend() } -> std::contiguous_iterator;
};

// template<ParticleContainer ty>
// imp(ty cont) {}

/*

class ParticleContainer {
 public:
  virtual ~ParticleContainer() = default;  // ?

  // retrieve data
  virtual Particle& operator[](size_t idx) = 0;
  virtual const Particle& operator[](size_t idx) const = 0;

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
};*/

}  // namespace particle_containers