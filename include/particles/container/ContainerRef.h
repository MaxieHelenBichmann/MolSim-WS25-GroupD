#ifndef CONTAINER_REF_H
#define CONTAINER_REF_H

#include <variant>

#include "particles/container/SimpleContainer.h"

namespace mol_sim {

/**
 * @brief Container Reference for Particles
 *
 * This container implements the concept ParticleContainer.
 * It essentially stores a pointer (not null) to an arbitrary ParticleContainer - currently only SimpleContainer
 * possible - but does not gain ownership over the elements in the Container. It calls the methods of the actual
 * SimpleContainer from which it was constructed. It is virtually a Generic Type for a ParticleContainer, which can be
 * instantiated.
 * Not being based on templates, it does not increase compile time as much, but checks for which function to call are
 * made at run time.
 *
 * Only use for UNCRITICAL functions/operations (e.g., IO). DO NOT USE IN PERFORMANCE-RELEVANT FUNCTIONS!
 *
 * Do not take References of a ContainerRef object, as it is a reference itself, and trivially copyable.
 * Alway pass-by-value (copy) - similar to std::span.
 *
 * (maybe will get erased later, if IO functions will be templated as well and compile time and code size acceptable)
 *
 */
class ContainerRef {
  // more containers can be added in variant
  std::variant<SimpleContainer*> _instance; //NOLINT

 public:
  // constructors
  ContainerRef(SimpleContainer& c);

  // retrieve data
  Particle& operator[](size_t idx);
  const Particle& operator[](size_t idx) const;

  [[nodiscard]] size_t size() const;
  [[nodiscard]] bool empty() const;

  // modify
  void clear();
  void reserve(size_t n);

  void addParticle(Particle&& value);
  void addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg);
  void addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg, int type);

  // iterators
  std::vector<Particle>::iterator begin();
  [[nodiscard]] std::vector<Particle>::const_iterator begin() const;
  [[nodiscard]] std::vector<Particle>::const_iterator cbegin() const;
  std::vector<Particle>::iterator end();
  [[nodiscard]] std::vector<Particle>::const_iterator end() const;
  [[nodiscard]] std::vector<Particle>::const_iterator cend() const;
};

}  // namespace mol_sim

#endif
