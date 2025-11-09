#ifndef ABSTRACT_CONTAINER_H
#define ABSTRACT_CONTAINER_H

#include "particles/Particle.h"

namespace mol_sim {
/**
 * @deprecated ONLY USED FOR BENCHMARKING
 * @brief Abstract implementation of a Particle Container
 */
class AbstractContainer {
   public:
    virtual ~AbstractContainer() = default;

    virtual void addParticle(Particle&& value) = 0;
    virtual void addParticle(const Particle& value) = 0;
    virtual void addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg, double epsilon_arg,
                             double sigma_arg) = 0;

    virtual Particle* begin() = 0;
    virtual Particle* end() = 0;
};

}  // namespace mol_sim

#endif
