#ifndef SIMPLE_CONTAINER_H
#define SIMPLE_CONTAINER_H

#include <vector>

#include "particles/ParticleContainer.h"

namespace mol_sim {

/**
 * @brief Simple Container for Particles
 *
 * This container implements the concept ParticleContainer.
 * It essentially uses a simple std::vector<Particle> to store the Particles.
 * Thus it is helpful for it to inherit all members (and constructors) directly from std::vector, immediately fulfilling
 * large part of the ParticleContainer concept. It implements the remaining functions for adding a Particle.
 *
 */
class SimpleContainer : public std::vector<Particle> {
   public:
    // constructors
    using std::vector<Particle>::vector;

    /**
     * @brief Adding an already existing Particle to the SimpleContainer by copying it.
     */
    void addParticle(Particle&& value);

    /**
     * @brief Adding an already existing const Particle to the SimpleContainer by copying it.
     */
    void addParticle(const Particle& value);

    /**
     * @brief Directly constructing a Particle with its required parameters in-place.
     */
    void addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg);

    /**
     * @brief Directly constructing a Particle with its required parameters, including type, in-place.
     */
    void addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg, int type);
};

}  // namespace mol_sim

#endif