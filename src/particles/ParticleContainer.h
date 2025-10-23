#pragma once

#include <concepts>
#include <vector>

#include "particles/Particle.h"

namespace particle_containers {

using R3 = Vector<double, 3>;

/**
 * @brief Concept of a Particle Container
 *
 * This is a concept requiring all methods, that a arbitrary Particle Container should implement,
 * including access possibilities, modifiers, and iterators.
 * It functions like an interface or an abstract class, which specific particle containers should implement.
 * However, this implementation is with some certainty more performant than virtual function calls,
 * because their lookup happens at runtime, whereas this concepts makes lookup at compile time possible.
 *
 * To efficiently use this concept, the performance-relevant functions have to be templated:
 * ```
 * template<ParticleContainer containerType>
 * frequently_used_func(containerType& particles) { ... }
 * ```
 *
 * Then, this function can be called with any specific particle container, which fulfills this concept. For example:
 * ```
 * particle_containers::SimpleContainer cur_particles{p1, p2, p3};
 * frequently_used_func<particle_containers::SimpleContainer>(cur_particles);
 * ```
 *
 * For uncritical functions (e.g. Input/Output), a non-templated version can be used,
 * see particle_containers::ContainerRef.
 */
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

}  // namespace particle_containers