#ifndef GENERATOR_H
#define GENERATOR_H

#include "particles/container/ContainerRef.h"
namespace mol_sim {

class Generator {
    virtual void generateParticles(ContainerRef particles) = 0;
};
}  // namespace mol_sim

#endif
