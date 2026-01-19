#ifndef GENERATOR_H
#define GENERATOR_H

#include "particles/container/ContainerRef.h"
namespace mol_sim {
/**
 * @brief Abstract Class for particle Generators
 *
 */
class Generator {
   public:
    /**
     * @brief Generates particles and places them in the container.
     *
     * @param particles Container to place the generated particles in.
     * @param use_brownian_motion If true, apply Brownian motion to particle velocities.
     * @param use_init_temp If true, use init_temp to calculate avg_velo based on temperature.
     */
    virtual void generateParticles(ContainerRef particles, bool use_brownian_motion = true,
                                   bool use_init_temp = false) = 0;
};
}  // namespace mol_sim

#endif
