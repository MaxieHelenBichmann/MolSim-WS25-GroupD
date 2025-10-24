#include "ContainerRef.h"

namespace particle_containers {

ContainerRef::ContainerRef(SimpleContainer& c) : _instance(&c) {}

Particle& ContainerRef::operator[](size_t idx) {
  return std::visit([idx](auto& c) -> Particle& { return (*c)[idx]; }, _instance);
}
const Particle& ContainerRef::operator[](size_t idx) const {
  return std::visit([idx](const auto& c) -> const Particle& { return (*c)[idx]; }, _instance);
}

size_t ContainerRef::size() const {
  return std::visit([](const auto& c) { return c->size(); }, _instance);
}
bool ContainerRef::empty() const {
  return std::visit([](const auto& c) { return c->empty(); }, _instance);
}

// modify
void ContainerRef::clear() {
  std::visit([](auto& c) { return c->clear(); }, _instance);
}
void ContainerRef::reserve(size_t n) {
  std::visit([n](auto& c) { return c->reserve(n); }, _instance);
}

void ContainerRef::addParticle(Particle&& value) {
  std::visit([&value](auto& c) { return c->addParticle(value); }, _instance);
}
void ContainerRef::addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg) {
  std::visit([&x_arg, &v_arg, m_arg](auto& c) { return c->addParticle(x_arg, v_arg, m_arg); }, _instance);
}
void ContainerRef::addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg, int type) {
  std::visit([&x_arg, &v_arg, m_arg, type](auto& c) { return c->addParticle(x_arg, v_arg, m_arg, type); }, _instance);
}

// iterators
std::vector<Particle>::iterator ContainerRef::begin() {
  return std::visit([](auto& c) { return c->begin(); }, _instance);
}
std::vector<Particle>::const_iterator ContainerRef::begin() const {
  return std::visit([](const auto& c) { return c->begin(); }, _instance);
}
std::vector<Particle>::const_iterator ContainerRef::cbegin() const {
  return std::visit([](const auto& c) { return c->cbegin(); }, _instance);
}
std::vector<Particle>::iterator ContainerRef::end() {
  return std::visit([](auto& c) { return c->end(); }, _instance);
}
std::vector<Particle>::const_iterator ContainerRef::end() const {
  return std::visit([](const auto& c) { return c->end(); }, _instance);
}
std::vector<Particle>::const_iterator ContainerRef::cend() const {
  return std::visit([](const auto& c) { return c->cend(); }, _instance);
}

static_assert(ParticleContainer<ContainerRef>);

}  // namespace particle_containers
