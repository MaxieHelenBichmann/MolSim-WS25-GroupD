#include "SimpleContainerControl.h"

#include <cstddef>

using namespace mol_sim;

SimpleContainerControl::SimpleContainerControl(R3 domain_size_arg, double cutoff_radius_arg)
    : domain_size(domain_size_arg), cutoff_radius(cutoff_radius_arg) {}

void SimpleContainerControl::addParticle(Particle&& value) { push_back(std::move(value)); }
void SimpleContainerControl::addParticle(const Particle& value) { push_back(value); }

void SimpleContainerControl::addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg) {
    emplace_back(x_arg, v_arg, m_arg, epsilon_arg, sigma_arg);
};
void SimpleContainerControl::addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg,
                                         int type) {
    emplace_back(x_arg, v_arg, m_arg, epsilon_arg, sigma_arg, type);
};
std::vector<Particle>::iterator SimpleContainerControl::eraseParticle(std::vector<Particle>::iterator p) {
    return erase(p);
};

SimpleContainerControl::proximity_iterator SimpleContainerControl::eraseParticle(
    const SimpleContainerControl::proximity_iterator& p) {
    R3 center = p.getCenter();
    double radius = p.getRadius();
    bool prox = p.isProximity();
    std::set<BoundaryLocation> locations = p.getLocations();
    size_t idx = erase(begin() + (&(*p) - data())) - begin();
    return proximity_iterator{center, radius, data(), data() + idx, 0, prox, locations};
};

std::vector<Particle>::iterator SimpleContainerControl::updateParticlePosition(std::vector<Particle>::iterator p,
                                                                               R3 new_x) {
    p->getX() = new_x;
    return ++p;
}

// proximity iterators
SimpleContainerControl::proximity_iterator SimpleContainerControl::proximityBegin(R3 center,
                                                                                  [[maybe_unused]] size_t offset) {
    return proximity_iterator{center, cutoff_radius, data(), data() + size(), 0};  // NOLINT
};

SimpleContainerControl::const_proximity_iterator SimpleContainerControl::proximityBegin(
    R3 center, [[maybe_unused]] size_t offset) const {
    return const_proximity_iterator{center, cutoff_radius, data(), data() + size(), 0};  // NOLINT
};

SimpleContainerControl::proximity_iterator SimpleContainerControl::proximityEnd(R3 center) {
    return proximity_iterator{center, cutoff_radius, data() + size(), data() + size(), 0};  // NOLINT
};

SimpleContainerControl::const_proximity_iterator SimpleContainerControl::proximityEnd(R3 center) const {
    return const_proximity_iterator{center, cutoff_radius, data() + size(), data() + size(), 0};  // NOLINT
};

// boundary and halo iterators
SimpleContainerControl::proximity_iterator SimpleContainerControl::haloBegin(
    const std::set<BoundaryLocation>& locations) {
    return proximity_iterator{domain_size, cutoff_radius, data(), data() + size(), 0, false, locations};  // NOLINT
};
SimpleContainerControl::proximity_iterator SimpleContainerControl::haloEnd(
    const std::set<BoundaryLocation>& locations) {
    return proximity_iterator{domain_size, cutoff_radius, data() + size(), data() + size(),  // NOLINT
                              0,           false,         locations};
};
SimpleContainerControl::const_proximity_iterator SimpleContainerControl::haloBegin(
    const std::set<BoundaryLocation>& locations) const {
    return const_proximity_iterator{domain_size, cutoff_radius, data(),   data() + size(),  // NOLINT
                                    0,           false,         locations};
};
SimpleContainerControl::const_proximity_iterator SimpleContainerControl::haloEnd(
    const std::set<BoundaryLocation>& locations) const {
    return const_proximity_iterator{domain_size, cutoff_radius, data() + size(), data() + size(),  // NOLINT
                                    0,           false,         locations};
};
SimpleContainerControl::proximity_iterator SimpleContainerControl::boundaryBegin(
    const std::set<BoundaryLocation>& locations) {
    return proximity_iterator{domain_size, -cutoff_radius, data(), data() + size(), 0, false, locations};  // NOLINT
};
SimpleContainerControl::proximity_iterator SimpleContainerControl::boundaryEnd(
    const std::set<BoundaryLocation>& locations) {
    return proximity_iterator{domain_size, -cutoff_radius, data() + size(), data() + size(),  // NOLINT
                              0,           false,          locations};
};
SimpleContainerControl::const_proximity_iterator SimpleContainerControl::boundaryBegin(
    const std::set<BoundaryLocation>& locations) const {
    return const_proximity_iterator{domain_size, -cutoff_radius, data(),   data() + size(),  // NOLINT
                                    0,           false,          locations};
};
SimpleContainerControl::const_proximity_iterator SimpleContainerControl::boundaryEnd(
    const std::set<BoundaryLocation>& locations) const {
    return const_proximity_iterator{domain_size, -cutoff_radius, data() + size(), data() + size(),  // NOLINT
                                    0,           false,          locations};
};
