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
     * @brief Adding an already existing Particle to the container by copying it.
     *
     * @param value Already constructed Particle.
     */
    void addParticle(Particle&& value);

    /**
     * @brief Adding an already existing const Particle to the container by copying it.
     *
     * @param value Already constructed Particle.
     */
    void addParticle(const Particle& value);

    /**
     * @brief Directly constructing a Particle with its required parameters, including type, in-place.
     *
     * @param x_arg Initial coordinates of the Particle.
     * @param v_arg Initial velocities of the Particle.
     * @param m_arg Mass of the Particle.
     * @param type Type of the Particle.
     */
    void addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg, double epsilon_arg,
                     double sigma_arg);
    void addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg, double epsilon_arg,
                     double sigma_arg, int type);
};

}  // namespace mol_sim

#endif
