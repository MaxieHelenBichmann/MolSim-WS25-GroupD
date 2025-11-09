#include "ContainerImpl.h"

namespace mol_sim {

void ContainerImpl::addParticle(Particle&& value) { particles.push_back(value); }

void ContainerImpl::addParticle(const Particle& value) { particles.push_back(value); }

void ContainerImpl::addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg) {
    particles.emplace_back(x_arg, v_arg, m_arg);
}

Particle* ContainerImpl::begin() { return particles.data(); }

Particle* ContainerImpl::end() { return particles.data() + particles.size(); }

}  // namespace mol_sim
