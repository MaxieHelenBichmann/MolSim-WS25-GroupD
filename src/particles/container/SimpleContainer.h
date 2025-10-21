#pragma once

#include <vector>

#include "../ParticleContainer.h"

namespace particle_containers {

class SimpleContainer : public ParticleContainer {
  std::vector<Particle> _data;

 public:
  // constructors
  SimpleContainer();
  SimpleContainer(size_t init_cap);
  SimpleContainer(std::initializer_list<Particle> init);

  // cpy constr
  SimpleContainer(const SimpleContainer& other);
  // cpy assignment
  SimpleContainer& operator=(const SimpleContainer& other);
  // mv constr
  SimpleContainer(SimpleContainer&& other);
  // mv assignment
  SimpleContainer& operator=(SimpleContainer&& other);
  // dstr
  ~SimpleContainer();

  // retrieve data
  Particle& operator[](size_t idx) override;
  const Particle& operator[](size_t idx) const override;

  std::vector<Particle>& data() override;
  const std::vector<Particle>& data() const override;

  bool empty() const override;
  size_t size() const override;

  bool operator==(const SimpleContainer& other) const;

  // modify
  void clear() override;
  void reserve(size_t n) override;
  void addParticle(Particle&& value) override;
  void addParticle(const Particle& value) override;
  void addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg) override;
  void addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg, int type) override;

  // iterators
  std::vector<Particle>::iterator begin() override;
  std::vector<Particle>::const_iterator begin() const override;
  std::vector<Particle>::iterator end() override;
  std::vector<Particle>::const_iterator end() const override;
};
}  // namespace particle_containers