#include "particles/container/SimpleContainer.h"

#include <cstddef>

using namespace mol_sim;

SimpleContainer::SimpleContainer(R3 domain_size_arg, double cutoff_radius_arg)
    : domain_size(domain_size_arg), cutoff_radius(cutoff_radius_arg) {}

bool SimpleContainer::fitsContainer(R3 v) {
    return (v[0] >= -cutoff_radius && v[0] <= domain_size[0] + cutoff_radius) &&
           (v[1] >= -cutoff_radius && v[1] <= domain_size[1] + cutoff_radius) &&
           (v[2] >= -cutoff_radius && v[2] <= domain_size[2] + cutoff_radius);
}

void SimpleContainer::addParticle(Particle&& value) {
    if (!fitsContainer(value.getX())) {
        return;
    }
    push_back(std::move(value));
}
void SimpleContainer::addParticle(const Particle& value) {
    if (!fitsContainer(value.getX())) {
        return;
    }
    push_back(value);
}

void SimpleContainer::addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg) {
    if (!fitsContainer(x_arg)) {
        return;
    }
    emplace_back(x_arg, v_arg, m_arg, epsilon_arg, sigma_arg);
};
void SimpleContainer::addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg, int type) {
    if (!fitsContainer(x_arg)) {
        return;
    }
    emplace_back(x_arg, v_arg, m_arg, epsilon_arg, sigma_arg, type);
};
void SimpleContainer::addParticle(R3 x_arg, R3 old_x_arg, R3 v_arg, R3 f_arg, R3 old_f_arg, double m_arg,
                                  double epsilon_arg, double sigma_arg, int type) {
    if (!fitsContainer(x_arg)) {
        return;
    }
    emplace_back(x_arg, old_x_arg, v_arg, f_arg, old_f_arg, m_arg, epsilon_arg, sigma_arg, type);
};

std::vector<Particle>::iterator SimpleContainer::eraseParticle(std::vector<Particle>::iterator p) { return erase(p); };

std::vector<Particle>::iterator SimpleContainer::updateParticlePosition(std::vector<Particle>::iterator p, R3 new_x) {
    p->getX() = new_x;
    return ++p;
}

// proximity iterators
SimpleContainer::proximity_iterator<Particle> SimpleContainer::proximityBegin(R3 center, size_t offset) {
    return proximity_iterator<Particle>{center, cutoff_radius, data(), data() + size(),  // NOLINT
                                        offset == size() ? 0 : ++offset};
};

SimpleContainer::proximity_iterator<const Particle> SimpleContainer::proximityBegin(R3 center, size_t offset) const {
    return proximity_iterator<const Particle>{center, cutoff_radius, data(), data() + size(),  // NOLINT
                                              offset == size() ? 0 : ++offset};
};

SimpleContainer::proximity_iterator<Particle> SimpleContainer::proximityEnd(R3 center) {
    return proximity_iterator<Particle>{center, cutoff_radius, data() + size(), data() + size(), 0};  // NOLINT
};

SimpleContainer::proximity_iterator<const Particle> SimpleContainer::proximityEnd(R3 center) const {
    return proximity_iterator<const Particle>{center, cutoff_radius, data() + size(), data() + size(), 0};  // NOLINT
};

// boundary and halo iterators
SimpleContainer::proximity_iterator<Particle> SimpleContainer::haloBegin(const std::set<BoundaryLocation>& locations) {
    return proximity_iterator<Particle>{domain_size, cutoff_radius, data(),   data() + size(),  // NOLINT
                                        0,           false,         locations};
};
SimpleContainer::proximity_iterator<Particle> SimpleContainer::haloEnd(const std::set<BoundaryLocation>& locations) {
    return proximity_iterator<Particle>{domain_size, cutoff_radius, data() + size(), data() + size(),  // NOLINT
                                        0,           false,         locations};
};
SimpleContainer::proximity_iterator<const Particle> SimpleContainer::haloBegin(
    const std::set<BoundaryLocation>& locations) const {
    return proximity_iterator<const Particle>{domain_size, cutoff_radius, data(),   data() + size(),  // NOLINT
                                              0,           false,         locations};
};
SimpleContainer::proximity_iterator<const Particle> SimpleContainer::haloEnd(
    const std::set<BoundaryLocation>& locations) const {
    return proximity_iterator<const Particle>{domain_size, cutoff_radius, data() + size(), data() + size(),  // NOLINT
                                              0,           false,         locations};
};
SimpleContainer::proximity_iterator<Particle> SimpleContainer::boundaryBegin(
    const std::set<BoundaryLocation>& locations) {
    return proximity_iterator<Particle>{domain_size, -cutoff_radius, data(),   data() + size(),  // NOLINT
                                        0,           false,          locations};
};
SimpleContainer::proximity_iterator<Particle> SimpleContainer::boundaryEnd(
    const std::set<BoundaryLocation>& locations) {
    return proximity_iterator<Particle>{domain_size, -cutoff_radius, data() + size(), data() + size(),  // NOLINT
                                        0,           false,          locations};
};
SimpleContainer::proximity_iterator<const Particle> SimpleContainer::boundaryBegin(
    const std::set<BoundaryLocation>& locations) const {
    return proximity_iterator<const Particle>{domain_size, -cutoff_radius, data(),   data() + size(),  // NOLINT
                                              0,           false,          locations};
};
SimpleContainer::proximity_iterator<const Particle> SimpleContainer::boundaryEnd(
    const std::set<BoundaryLocation>& locations) const {
    return proximity_iterator<const Particle>{domain_size, -cutoff_radius, data() + size(), data() + size(),  // NOLINT
                                              0,           false,          locations};
};
