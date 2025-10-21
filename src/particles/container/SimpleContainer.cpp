#include "SimpleContainer.h"

namespace particle_containers {

SimpleContainer::SimpleContainer() {}
SimpleContainer::SimpleContainer(size_t init_cap) { _data.reserve(init_cap); }
SimpleContainer::SimpleContainer(std::initializer_list<Particle> init) {}

// cpy constr
SimpleContainer::SimpleContainer(const SimpleContainer& other) {}
// cpy assignment
SimpleContainer& SimpleContainer::operator=(const SimpleContainer& other) {}
// mv constr
SimpleContainer::SimpleContainer(SimpleContainer&& other) : _data(std::move(other._data)) {}
// mv assignment
SimpleContainer& SimpleContainer::operator=(SimpleContainer&& other) { _data = std::move(other._data); }
// dstr
SimpleContainer::~SimpleContainer() {}
// -> nothing of this has to be custom

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
void SimpleContainer::removeParticle(size_t idx) {}
void SimpleContainer::removeParticle(Particle& value) {}
void SimpleContainer::addParticle(Particle&& value) {}
void SimpleContainer::addParticle(const Particle& value) {}

// iterators
std::vector<Particle>::iterator SimpleContainer::begin() { return _data.begin(); }
std::vector<Particle>::const_iterator SimpleContainer::begin() const { return _data.begin(); }
std::vector<Particle>::iterator SimpleContainer::end() { return _data.end(); }
std::vector<Particle>::const_iterator SimpleContainer::end() const { return _data.end(); }

}  // namespace particle_containers