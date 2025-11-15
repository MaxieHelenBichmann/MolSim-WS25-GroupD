#ifndef DISC_GENERATOR_H
#define DISC_GENERATOR_H

#include <cstddef>

#include "particles/Generator.h"
#include "utils/Vector.h"

namespace mol_sim {
/**
 * @brief Particle Generator, which generates a Cube of Particles with Brownian Motion.
 *
 */
class DiscGenerator : public Generator {
   private:
    /**
     * @brief Postion of the center of the disc of particles.
     *
     */
    R3 position;
    /**
     * @brief Initial velocity of the disc of particles.
     *
     */
    R3 velocity;
    /**
     * @brief number of particles along the radiues of the disc
     *
     */
    size_t radius;
    /**
     * @brief Mass of the particles.
     *
     */
    double mass;
    /**
     * @brief Mesh width of the grid of particles.
     *
     */
    double distance;
    /**
     * @brief Mean value of velocity of the Brownian Motion.
     *
     */
    double avg_velo;
    /**
     * @brief Epsilon value of the generated particles.
     *
     */
    double epsilon;
    /**
     * @brief Sigma value of the generated particles.
     *
     */
    double sigma;

   public:
    /**
     * @brief Function to generate Particles.
     *
     * @param particles Container to place Particles in.
     */
    void generateParticles(ContainerRef particles) override;
    /**
     * @brief Construct a new Cuboid Generator object, setting all necessary parameters for generation.
     *
     * @param position Postion of the lower left corner of the cube of particles.
     * @param velocity Initial velocity of the cube of particles.
     * @param num_particles Number of particles per dimension.
     * @param mass Mass of the particles.
     * @param distance Mesh width of the grid of particles.
     * @param avg_velo Mean value of velocity of the Brownian Motion.
     * @param epsilon Epsilon value of the generated particles.
     * @param sigma Sigma value of the generated particles.
     */
    DiscGenerator(R3 position, R3 velocity, size_t radius, double mass, double distance, double avg_velo,
                  double epsilon, double sigma)
        : position(position),
          velocity(velocity),
          radius(radius),
          mass(mass),
          distance(distance),
          avg_velo(avg_velo),
          epsilon(epsilon),
          sigma(sigma) {};
};
}  // namespace mol_sim

#endif
