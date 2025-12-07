#include "particles/container/ContainerRef.h"

#include "particles/Particle.h"

using namespace mol_sim;

ContainerRef::ContainerRef() = default;
ContainerRef::ContainerRef(SimpleContainer& c) : instance(&c) {}
ContainerRef::ContainerRef(LinkedCellContainer& c) : instance(&c) {}
ContainerRef::ContainerRef(CONTAINER_REF& c) : instance(c) {}

Particle& ContainerRef::operator[](size_t idx) {
    return std::visit([idx](auto& c) -> Particle& { return (*c)[idx]; }, instance);
}
const Particle& ContainerRef::operator[](size_t idx) const {
    return std::visit([idx](const auto& c) -> const Particle& { return (*c)[idx]; }, instance);
}

size_t ContainerRef::size() const {
    return std::visit([](const auto& c) { return c->size(); }, instance);
}
bool ContainerRef::empty() const {
    return std::visit([](const auto& c) { return c->empty(); }, instance);
}

// modify
void ContainerRef::clear() {
    std::visit([](auto& c) { return c->clear(); }, instance);
}
void ContainerRef::reserve(size_t n) {
    std::visit([n](auto& c) { return c->reserve(n); }, instance);
}

void ContainerRef::addParticle(Particle&& value) {
    std::visit([&value](auto& c) { return c->addParticle(value); }, instance);
}
void ContainerRef::addParticle(const Particle& value) {
    std::visit([&value](auto& c) { return c->addParticle(value); }, instance);
}

void ContainerRef::addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg) {
    std::visit([&x_arg, &v_arg, m_arg, epsilon_arg,
                sigma_arg](auto& c) { return c->addParticle(x_arg, v_arg, m_arg, epsilon_arg, sigma_arg); },
               instance);
}
void ContainerRef::addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg, int type) {
    std::visit([&x_arg, &v_arg, m_arg, epsilon_arg, sigma_arg,
                type](auto& c) { return c->addParticle(x_arg, v_arg, m_arg, epsilon_arg, sigma_arg, type); },
               instance);
}

void ContainerRef::addParticle(R3 x_arg, R3 old_x_arg, R3 v_arg, R3 f_arg, R3 old_f_arg, double m_arg,
                               double epsilon_arg, double sigma_arg, int type) {
    std::visit(
        [&x_arg, &old_x_arg, &v_arg, &f_arg, &old_f_arg, m_arg, epsilon_arg, sigma_arg, type](auto& c) {
            return c->addParticle(x_arg, old_x_arg, v_arg, f_arg, old_f_arg, m_arg, epsilon_arg, sigma_arg, type);
        },
        instance);
}
std::vector<Particle>::iterator ContainerRef::eraseParticle(std::vector<Particle>::iterator p) {
    return std::visit([&p](auto& c) { return c->eraseParticle(p); }, instance);
}

ContainerRef::proximity_iterator<Particle, Cell> ContainerRef::eraseParticle(
    ContainerRef::proximity_iterator<Particle, Cell> p) {
    return std::visit([&p](auto& c) { return proximity_iterator<Particle, Cell>{c->eraseParticle(p)}; }, instance);
}

std::vector<Particle>::iterator ContainerRef::updateParticlePosition(std::vector<Particle>::iterator p, R3 new_x) {
    return std::visit([&p, &new_x](auto& c) { return c->updateParticlePosition(p, new_x); }, instance);
}

CONTAINER_REF ContainerRef::getInstance() { return instance; }

// normal iterators
std::vector<Particle>::iterator ContainerRef::begin() {
    return std::visit([](auto& c) { return c->begin(); }, instance);
}
std::vector<Particle>::const_iterator ContainerRef::begin() const {
    return std::visit([](const auto& c) { return c->begin(); }, instance);
}
std::vector<Particle>::const_iterator ContainerRef::cbegin() const {
    return std::visit([](const auto& c) { return c->cbegin(); }, instance);
}
std::vector<Particle>::iterator ContainerRef::end() {
    return std::visit([](auto& c) { return c->end(); }, instance);
}
std::vector<Particle>::const_iterator ContainerRef::end() const {
    return std::visit([](const auto& c) { return c->end(); }, instance);
}
std::vector<Particle>::const_iterator ContainerRef::cend() const {
    return std::visit([](const auto& c) { return c->cend(); }, instance);
}

// proximity iterators
ContainerRef::proximity_iterator<Particle, Cell> ContainerRef::proximityBegin(R3 center, size_t offset) {
    return std::visit(
        [center, offset](auto& c) { return proximity_iterator<Particle, Cell>{c->proximityBegin(center, offset)}; },
        instance);
}

ContainerRef::proximity_iterator<Particle, Cell> ContainerRef::proximityEnd(R3 center) {
    return std::visit([center](auto& c) { return proximity_iterator<Particle, Cell>{c->proximityEnd(center)}; },
                      instance);
}

ContainerRef::proximity_iterator<const Particle, const Cell> ContainerRef::proximityBegin(R3 center,
                                                                                          size_t offset) const {
    return std::visit(
        [center, offset](const auto& c) {
            return proximity_iterator<const Particle, const Cell>{std::as_const(*c).proximityBegin(center, offset)};
        },
        instance);
}
ContainerRef::proximity_iterator<const Particle, const Cell> ContainerRef::proximityEnd(R3 center) const {
    return std::visit(
        [center](const auto& c) {
            return proximity_iterator<const Particle, const Cell>{std::as_const(*c).proximityEnd(center)};
        },
        instance);
}

// boundary and halo iterators
ContainerRef::proximity_iterator<Particle, Cell> ContainerRef::haloBegin(const std::set<BoundaryLocation>& locations) {
    return std::visit([locations](auto& c) { return proximity_iterator<Particle, Cell>{c->haloBegin(locations)}; },
                      instance);
};
ContainerRef::proximity_iterator<Particle, Cell> ContainerRef::haloEnd(const std::set<BoundaryLocation>& locations) {
    return std::visit([locations](auto& c) { return proximity_iterator<Particle, Cell>{c->haloEnd(locations)}; },
                      instance);
};
ContainerRef::proximity_iterator<const Particle, const Cell> ContainerRef::haloBegin(
    const std::set<BoundaryLocation>& locations) const {
    return std::visit(
        [locations](const auto& c) {
            return proximity_iterator<const Particle, const Cell>{std::as_const(*c).haloBegin(locations)};
        },
        instance);
};
ContainerRef::proximity_iterator<const Particle, const Cell> ContainerRef::haloEnd(
    const std::set<BoundaryLocation>& locations) const {
    return std::visit(
        [locations](const auto& c) {
            return proximity_iterator<const Particle, const Cell>{std::as_const(*c).haloEnd(locations)};
        },
        instance);
};
ContainerRef::proximity_iterator<Particle, Cell> ContainerRef::boundaryBegin(
    const std::set<BoundaryLocation>& locations) {
    return std::visit([locations](auto& c) { return proximity_iterator<Particle, Cell>{c->boundaryBegin(locations)}; },
                      instance);
};
ContainerRef::proximity_iterator<Particle, Cell> ContainerRef::boundaryEnd(
    const std::set<BoundaryLocation>& locations) {
    return std::visit([locations](auto& c) { return proximity_iterator<Particle, Cell>{c->boundaryEnd(locations)}; },
                      instance);
};
ContainerRef::proximity_iterator<const Particle, const Cell> ContainerRef::boundaryBegin(
    const std::set<BoundaryLocation>& locations) const {
    return std::visit(
        [locations](const auto& c) {
            return proximity_iterator<const Particle, const Cell>{std::as_const(*c).boundaryBegin(locations)};
        },
        instance);
};
ContainerRef::proximity_iterator<const Particle, const Cell> ContainerRef::boundaryEnd(
    const std::set<BoundaryLocation>& locations) const {
    return std::visit(
        [locations](const auto& c) {
            return proximity_iterator<const Particle, const Cell>{std::as_const(*c).boundaryEnd(locations)};
        },
        instance);
};
