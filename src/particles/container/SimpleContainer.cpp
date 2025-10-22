#include "SimpleContainer.h"

namespace particle_containers {

SimpleContainer::SimpleContainer() = default;
SimpleContainer::SimpleContainer(size_t init_cap) { _data.reserve(init_cap); }
SimpleContainer::SimpleContainer(std::initializer_list<Particle> init) {
  _data.reserve(init.size());
  for (const Particle& elem : init) {
    _data.push_back(elem);
  }
}

// cpy constr
SimpleContainer::SimpleContainer(const SimpleContainer& other) = default;
// cpy assignment
SimpleContainer& SimpleContainer::operator=(const SimpleContainer& other) = default;
// mv constr
SimpleContainer::SimpleContainer(SimpleContainer&& other) = default;
// mv assignment
SimpleContainer& SimpleContainer::operator=(SimpleContainer&& other) = default;
// dstr
SimpleContainer::~SimpleContainer() = default;

// retrieve data
Particle& SimpleContainer::operator[](size_t idx) { return _data[idx]; }
const Particle& SimpleContainer::operator[](size_t idx) const { return _data[idx]; }

std::vector<Particle>& SimpleContainer::data() { return _data; }
const std::vector<Particle>& SimpleContainer::data() const { return _data; }

bool SimpleContainer::empty() const { return _data.empty(); }
size_t SimpleContainer::size() const { return _data.size(); }

bool SimpleContainer::operator==(const SimpleContainer& other) const { return _data == other._data; }

// modify
void SimpleContainer::clear() { _data.clear(); }
void SimpleContainer::reserve(size_t n) { _data.reserve(n); }
void SimpleContainer::addParticle(Particle&& value) { _data.push_back(value); }
void SimpleContainer::addParticle(const Particle& value) { _data.push_back(value); }
void SimpleContainer::addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg) {
  _data.emplace_back(x_arg, v_arg, m_arg);
};
void SimpleContainer::addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg, int type) {
  _data.emplace_back(x_arg, v_arg, m_arg, type);
};

// iterators
std::vector<Particle>::iterator SimpleContainer::begin() { return _data.begin(); }
std::vector<Particle>::const_iterator SimpleContainer::begin() const { return _data.begin(); }
std::vector<Particle>::iterator SimpleContainer::end() { return _data.end(); }
std::vector<Particle>::const_iterator SimpleContainer::end() const { return _data.end(); }

}  // namespace particle_containers
