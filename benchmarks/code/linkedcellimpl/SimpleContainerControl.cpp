#include "SimpleContainerControl.h"

#include <particles/Particle.h>

#include <cstddef>

using namespace mol_sim;

SimpleContainerControl::SimpleContainerControl(R3 domain_size_arg, double cutoff_radius_arg)
    : domain_size(domain_size_arg), cutoff_radius(cutoff_radius_arg) {}

bool SimpleContainerControl::fitsDomain(R3 v) const noexcept {
    return (v[0] >= 0.0 && v[0] <= domain_size[0]) && (v[1] >= 0.0 && v[1] <= domain_size[1]) &&
           (v[2] >= 0.0 && v[2] <= domain_size[2]);
}

void SimpleContainerControl::addParticle(Particle&& value) { push_back(std::move(value)); }
void SimpleContainerControl::addParticle(const Particle& value) { push_back(value); }

void SimpleContainerControl::addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg) {
    emplace_back(x_arg, v_arg, m_arg, epsilon_arg, sigma_arg);
};
void SimpleContainerControl::addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg,
                                         int type) {
    emplace_back(x_arg, v_arg, m_arg, epsilon_arg, sigma_arg, type);
};
void SimpleContainerControl::addParticle(R3 x_arg, R3 old_x_arg, R3 v_arg, R3 f_arg, R3 old_f_arg, double m_arg,
                                         double epsilon_arg, double sigma_arg, int type) {
    emplace_back(x_arg, old_x_arg, v_arg, f_arg, old_f_arg, m_arg, epsilon_arg, sigma_arg, type);
};
std::vector<Particle>::iterator SimpleContainerControl::eraseParticle(std::vector<Particle>::iterator p) {
    return erase(p);
};

SimpleContainerControl::proximity_iterator<Particle> SimpleContainerControl::eraseParticle(
    const SimpleContainerControl::proximity_iterator<Particle>& p) {
    R3 center = p.getCenter();
    double radius = p.getRadius();
    bool prox = p.isProximity();
    std::set<BoundaryLocation> locations = p.getLocations();
    size_t idx = erase(begin() + (&(*p) - data())) - begin();
    return proximity_iterator<Particle>{center, radius, data(), data() + idx, 0, prox, locations};
};

std::vector<Particle>::iterator SimpleContainerControl::updateParticlePosition(std::vector<Particle>::iterator p,
                                                                               R3 new_x) {
    p->getX() = new_x;
    return ++p;
}

// proximity iterators
SimpleContainerControl::proximity_iterator<Particle> SimpleContainerControl::proximityBegin(
    R3 center, [[maybe_unused]] size_t offset) {
    return proximity_iterator<Particle>{center, cutoff_radius, data(), data() + size(), 0};  // NOLINT
};

SimpleContainerControl::proximity_iterator<const Particle> SimpleContainerControl::proximityBegin(
    R3 center, [[maybe_unused]] size_t offset) const {
    return proximity_iterator<const Particle>{center, cutoff_radius, data(), data() + size(), 0};  // NOLINT
};

SimpleContainerControl::proximity_iterator<Particle> SimpleContainerControl::proximityEnd(R3 center) {
    return proximity_iterator<Particle>{center, cutoff_radius, data() + size(), data() + size(), 0};  // NOLINT
};

SimpleContainerControl::proximity_iterator<const Particle> SimpleContainerControl::proximityEnd(R3 center) const {
    return proximity_iterator<const Particle>{center, cutoff_radius, data() + size(), data() + size(), 0};  // NOLINT
};

SimpleContainerControl::proximity_iterator<Particle> SimpleContainerControl::proximityBegin_no_N3L(R3 center) {
    return proximityBegin(center);
};

SimpleContainerControl::proximity_iterator<const Particle> SimpleContainerControl::proximityBegin_no_N3L(
    R3 center) const {
    return proximityBegin(center);
};

SimpleContainerControl::proximity_iterator<Particle> SimpleContainerControl::proximityEnd_no_N3L(R3 center) {
    return proximityEnd(center);
};

SimpleContainerControl::proximity_iterator<const Particle> SimpleContainerControl::proximityEnd_no_N3L(
    R3 center) const {
    return proximityEnd(center);
};

// boundary and halo iterators
SimpleContainerControl::proximity_iterator<Particle> SimpleContainerControl::haloBegin(
    const std::set<BoundaryLocation>& locations) {
    return proximity_iterator<Particle>{domain_size, cutoff_radius, data(),   data() + size(),  // NOLINT
                                        0,           false,         locations};
};
SimpleContainerControl::proximity_iterator<Particle> SimpleContainerControl::haloEnd(
    const std::set<BoundaryLocation>& locations) {
    return proximity_iterator<Particle>{domain_size, cutoff_radius, data() + size(), data() + size(),  // NOLINT
                                        0,           false,         locations};
};
SimpleContainerControl::proximity_iterator<const Particle> SimpleContainerControl::haloBegin(
    const std::set<BoundaryLocation>& locations) const {
    return proximity_iterator<const Particle>{domain_size, cutoff_radius, data(),   data() + size(),  // NOLINT
                                              0,           false,         locations};
};
SimpleContainerControl::proximity_iterator<const Particle> SimpleContainerControl::haloEnd(
    const std::set<BoundaryLocation>& locations) const {
    return proximity_iterator<const Particle>{domain_size, cutoff_radius, data() + size(), data() + size(),  // NOLINT
                                              0,           false,         locations};
};
SimpleContainerControl::proximity_iterator<Particle> SimpleContainerControl::boundaryBegin(
    const std::set<BoundaryLocation>& locations) {
    return proximity_iterator<Particle>{domain_size, -cutoff_radius, data(),   data() + size(),  // NOLINT
                                        0,           false,          locations};
};
SimpleContainerControl::proximity_iterator<Particle> SimpleContainerControl::boundaryEnd(
    const std::set<BoundaryLocation>& locations) {
    return proximity_iterator<Particle>{domain_size, -cutoff_radius, data() + size(), data() + size(),  // NOLINT
                                        0,           false,          locations};
};
SimpleContainerControl::proximity_iterator<const Particle> SimpleContainerControl::boundaryBegin(
    const std::set<BoundaryLocation>& locations) const {
    return proximity_iterator<const Particle>{domain_size, -cutoff_radius, data(),   data() + size(),  // NOLINT
                                              0,           false,          locations};
};
SimpleContainerControl::proximity_iterator<const Particle> SimpleContainerControl::boundaryEnd(
    const std::set<BoundaryLocation>& locations) const {
    return proximity_iterator<const Particle>{domain_size, -cutoff_radius, data() + size(), data() + size(),  // NOLINT
                                              0,           false,          locations};
};
