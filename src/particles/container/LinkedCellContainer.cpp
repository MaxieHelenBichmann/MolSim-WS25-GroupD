#include "particles/container/LinkedCellContainer.h"

using namespace mol_sim;

LinkedCellContainer::LinkedCellContainer() {}

Particle& LinkedCellContainer::operator[](size_t idx) { return data[idx]; }
const Particle& LinkedCellContainer::operator[](size_t idx) const { return data[idx]; }

size_t LinkedCellContainer::size() const { return data.size(); }
bool LinkedCellContainer::empty() const { return data.empty(); }

// modify
void LinkedCellContainer::clear() { data.clear(); }
void LinkedCellContainer::reserve(size_t n) { data.reserve(n); }

void LinkedCellContainer::addParticle(Particle&& value) { (void)value; }
void LinkedCellContainer::addParticle(const Particle& value) { (void)value; }

void LinkedCellContainer::addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg) {
    (void)x_arg;
    (void)v_arg;
    (void)m_arg;
    (void)epsilon_arg;
    (void)sigma_arg;
}
void LinkedCellContainer::addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg,
                                      int type) {
    (void)x_arg;
    (void)v_arg;
    (void)m_arg;
    (void)epsilon_arg;
    (void)sigma_arg;
    (void)type;
}

// normal iterators
std::vector<Particle>::iterator LinkedCellContainer::begin() { return data.begin(); }
std::vector<Particle>::const_iterator LinkedCellContainer::begin() const { return data.begin(); }
std::vector<Particle>::const_iterator LinkedCellContainer::cbegin() const { return data.cbegin(); }
std::vector<Particle>::iterator LinkedCellContainer::end() { return data.end(); }
std::vector<Particle>::const_iterator LinkedCellContainer::end() const { return data.end(); }
std::vector<Particle>::const_iterator LinkedCellContainer::cend() const { return data.cend(); }

// proximity iterators
LinkedCellContainer::proximity_iterator LinkedCellContainer::proximityBegin(R3 center, double radius, size_t offset) {
    (void)center;
    (void)radius;
    (void)offset;
    return {};
}

LinkedCellContainer::proximity_iterator LinkedCellContainer::proximityEnd(R3 center, double radius) {
    (void)center;
    (void)radius;
    return {};
}

LinkedCellContainer::const_proximity_iterator LinkedCellContainer::proximityBegin(R3 center, double radius,
                                                                                  size_t offset) const {
    (void)center;
    (void)radius;
    (void)offset;
    return {};
}
LinkedCellContainer::const_proximity_iterator LinkedCellContainer::proximityEnd(R3 center, double radius) const {
    (void)center;
    (void)radius;
    return {};
}

// static_assert(ParticleContainer<LinkedCellContainer>);