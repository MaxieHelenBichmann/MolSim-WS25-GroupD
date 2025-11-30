#include "particles/container/SimpleContainer.h"

#include <algorithm>

using namespace mol_sim;

SimpleContainer::SimpleContainer(R3 domain_size_arg, double cutoff_radius_arg)
    : domain_size(domain_size_arg), cutoff_radius(cutoff_radius_arg) {}

void SimpleContainer::addParticle(Particle&& value) { push_back(value); }
void SimpleContainer::addParticle(const Particle& value) { push_back(value); }

void SimpleContainer::addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg) {
    emplace_back(x_arg, v_arg, m_arg, epsilon_arg, sigma_arg);
};
void SimpleContainer::addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg, int type) {
    emplace_back(x_arg, v_arg, m_arg, epsilon_arg, sigma_arg, type);
};
std::vector<Particle>::iterator SimpleContainer::eraseParticle(std::vector<Particle>::iterator p) { return erase(p); };

SimpleContainer::proximity_iterator SimpleContainer::eraseParticle(SimpleContainer::proximity_iterator p) {
    Particle particle = *p;
    R3 center = p.getCenter();
    double radius = p.getRadius();
    bool prox = p.isProximity();
    std::set<BoundaryLocation> locations = p.getLocations();
    auto it = eraseParticle(std::vector<Particle>::iterator{&particle});
    return proximity_iterator{center, radius, data(), data() + it, 0, prox, locations};
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

// boundary and halo iterators
SimpleContainer::proximity_iterator SimpleContainer::haloBegin(const std::set<BoundaryLocation>& locations) {
    return proximity_iterator{domain_size, cutoff_radius, data(), data() + size(), 0, false, locations};  // NOLINT
};
SimpleContainer::proximity_iterator SimpleContainer::haloEnd(const std::set<BoundaryLocation>& locations) {
    return proximity_iterator{domain_size, cutoff_radius, data() + size(), data() + size(),  // NOLINT
                              0,           false,         locations};
};
SimpleContainer::const_proximity_iterator SimpleContainer::haloBegin(
    const std::set<BoundaryLocation>& locations) const {
    return const_proximity_iterator{domain_size, cutoff_radius, data(),   data() + size(),  // NOLINT
                                    0,           false,         locations};
};
SimpleContainer::const_proximity_iterator SimpleContainer::haloEnd(const std::set<BoundaryLocation>& locations) const {
    return const_proximity_iterator{domain_size, cutoff_radius, data() + size(), data() + size(),  // NOLINT
                                    0,           false,         locations};
};
SimpleContainer::proximity_iterator SimpleContainer::boundaryBegin(const std::set<BoundaryLocation>& locations) {
    return proximity_iterator{domain_size, -cutoff_radius, data(), data() + size(), 0, false, locations};  // NOLINT
};
SimpleContainer::proximity_iterator SimpleContainer::boundaryEnd(const std::set<BoundaryLocation>& locations) {
    return proximity_iterator{domain_size, -cutoff_radius, data() + size(), data() + size(),  // NOLINT
                              0,           false,          locations};
};
SimpleContainer::const_proximity_iterator SimpleContainer::boundaryBegin(
    const std::set<BoundaryLocation>& locations) const {
    return const_proximity_iterator{domain_size, -cutoff_radius, data(),   data() + size(),  // NOLINT
                                    0,           false,          locations};
};
SimpleContainer::const_proximity_iterator SimpleContainer::boundaryEnd(
    const std::set<BoundaryLocation>& locations) const {
    return const_proximity_iterator{domain_size, -cutoff_radius, data() + size(), data() + size(),  // NOLINT
                                    0,           false,          locations};
};

static_assert(ParticleContainer<SimpleContainer>);