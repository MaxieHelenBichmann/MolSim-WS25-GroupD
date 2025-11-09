#ifndef CONTAINER_IMPL_H
#define CONTAINER_IMPL_H

#include <vector>

#include "AbstractContainer.h"

namespace mol_sim {
/**
 * @deprecated ONLY USED FOR BENCHMARKING
 * @brief Concrete Implementation of abstract container.
 */
class ContainerImpl : public AbstractContainer {
   private:
    std::vector<Particle> particles;

   public:
    void addParticle(Particle&& value) override;
    void addParticle(const Particle& value) override;
    void addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg, double epsilon_arg,
                     double sigma_arg) override;

    Particle* begin() override;
    Particle* end() override;
};

}  // namespace mol_sim

#endif
