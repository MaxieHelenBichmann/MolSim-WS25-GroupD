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
     */
    virtual void generateParticles(ContainerRef particles) = 0;
};
}  // namespace mol_sim

#endif
