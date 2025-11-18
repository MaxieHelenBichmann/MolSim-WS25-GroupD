#include "particles/container/ContainerRef.h"

using namespace mol_sim;

ContainerRef::ContainerRef(SimpleContainer& c) : instance(&c) {}

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

void ContainerRef::updateParticlePosition(std::vector<Particle>::iterator p, R3 new_x) {
    std::visit([&p, &new_x](auto& c) { return c->updateParticlePosition(p, new_x); }, instance);
}

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
SimpleContainer::proximity_iterator ContainerRef::proximityBegin(R3 center, double radius, size_t offset) {
    return std::visit([center, radius, offset](auto& c) { return c->proximityBegin(center, radius, offset); },
                      instance);
}

SimpleContainer::proximity_iterator ContainerRef::proximityEnd(R3 center, double radius) {
    return std::visit([center, radius](auto& c) { return c->proximityEnd(center, radius); }, instance);
}

SimpleContainer::const_proximity_iterator ContainerRef::proximityBegin(R3 center, double radius, size_t offset) const {
    return std::visit(
        [center, radius, offset](const auto& c) { return std::as_const(*c).proximityBegin(center, radius, offset); },
        instance);
}
SimpleContainer::const_proximity_iterator ContainerRef::proximityEnd(R3 center, double radius) const {
    return std::visit([center, radius](const auto& c) { return std::as_const(*c).proximityEnd(center, radius); },
                      instance);
}

static_assert(ParticleContainer<ContainerRef>);
