#ifndef PARTICLE_CONTAINER_H
#define PARTICLE_CONTAINER_H

#include <concepts>

#include "particles/Particle.h"

namespace mol_sim {

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
 * SimpleContainer cur_particles{p1, p2, p3};
 * frequently_used_func<SimpleContainer>(cur_particles);
 * ```
 *
 * For uncritical functions (e.g. Input/Output), a non-templated version can be used,
 * see ContainerRef.
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
    { c.addParticle(R3(), R3(), 0., 0., 0.) };
    { c.addParticle(R3(), R3(), 0., 0., 0., 0) };

    { c.updateParticlePosition(std::vector<Particle>::iterator(), R3()) };

    // iterators all particles
    { c.begin() } -> std::forward_iterator;
    { (*static_cast<const C*>(&c)).begin() } -> std::forward_iterator;
    { c.cbegin() } -> std::forward_iterator;
    { c.end() } -> std::forward_iterator;
    { (*static_cast<const C*>(&c)).end() } -> std::forward_iterator;
    { c.cend() } -> std::forward_iterator;

    // iterators particles in proximity
    { c.proximityBegin(R3(), double(), size_t()) } -> std::forward_iterator;
    { (*static_cast<const C*>(&c)).proximityBegin(R3(), double(), size_t()) } -> std::forward_iterator;
    { c.proximityEnd(R3(), double()) } -> std::forward_iterator;
    { (*static_cast<const C*>(&c)).proximityEnd(R3(), double()) } -> std::forward_iterator;
};

}  // namespace mol_sim

#endif
