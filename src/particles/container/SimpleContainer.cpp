#include "particles/container/SimpleContainer.h"

using namespace mol_sim;

void SimpleContainer::addParticle(Particle&& value) { push_back(value); }
void SimpleContainer::addParticle(const Particle& value) { push_back(value); }
void SimpleContainer::addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg) {
    emplace_back(x_arg, v_arg, m_arg);
};
void SimpleContainer::addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg, int type) {
    emplace_back(x_arg, v_arg, m_arg, type);
};
void SimpleContainer::addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg, double epsilon_arg,
                                  double sigma_arg) {
    emplace_back(x_arg, v_arg, m_arg, epsilon_arg, sigma_arg);
};
void SimpleContainer::addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg, double epsilon_arg,
                                  double sigma_arg, int type) {
    emplace_back(x_arg, v_arg, m_arg, epsilon_arg, sigma_arg, type);
};

static_assert(ParticleContainer<SimpleContainer>);
