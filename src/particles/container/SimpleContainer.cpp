#include "particles/container/SimpleContainer.h"

using namespace mol_sim;

void SimpleContainer::addParticle(Particle&& value) { push_back(value); }
void SimpleContainer::addParticle(const Particle& value) { push_back(value); }

void SimpleContainer::addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg) {
    emplace_back(x_arg, v_arg, m_arg, epsilon_arg, sigma_arg);
};
void SimpleContainer::addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg, int type) {
    emplace_back(x_arg, v_arg, m_arg, epsilon_arg, sigma_arg, type);
};

void SimpleContainer::updateParticlePosition(std::vector<Particle>::iterator p, R3 new_x) { p->getX() = new_x; }

// proximity iterators
SimpleContainer::proximity_iterator SimpleContainer::proximityBegin(R3 center, double radius, size_t offset) {
    return proximity_iterator{center, radius, data(), data() + size(), offset};  // NOLINT
};

SimpleContainer::const_proximity_iterator SimpleContainer::proximityBegin(R3 center, double radius,
                                                                          size_t offset) const {
    return const_proximity_iterator{center, radius, data(), data() + size(), offset};  // NOLINT
};

SimpleContainer::proximity_iterator SimpleContainer::proximityEnd(R3 center, double radius) {
    return proximity_iterator{center, radius, data() + size(), data() + size(), 0};  // NOLINT
};

SimpleContainer::const_proximity_iterator SimpleContainer::proximityEnd(R3 center, double radius) const {
    return const_proximity_iterator{center, radius, data() + size(), data() + size(), 0};  // NOLINT
};

static_assert(ParticleContainer<SimpleContainer>);
