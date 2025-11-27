#ifndef REFLECTING_H
#define REFLECTING_H

#include "particles/boundaries/BoundaryCondition.h"
#include "utils/Vector.h"
#include "particles/container/ContainerRef.h"

namespace mol_sim {

class Reflecting : public BoundaryCondition {
    R3 zero = {.0, .0, .0};
    R3 max;
    std::optional<double> counter_epsilon = std::nullopt;
    std::optional<double> counter_sigma = std::nullopt;
    SimpleContainer _dummy_container;
    ContainerRef particles = ContainerRef(_dummy_container);

    void addCounterParticle(int sign, size_t coordinate, Particle& p);
    bool fitsDomain(R3 v);

   public:
    Reflecting(BoundaryLocation location, R3 dimension, std::optional<double> counter_sigma, std::optional<double> counter_epsilon);
    Reflecting(BoundaryLocation location, R3 dimension, ContainerRef particles, std::optional<double> counter_sigma, std::optional<double> counter_epsilon);
    ~Reflecting() override;
    void boundaryStrategy(Particle& p) override;
    void setParticles(ContainerRef particles);
};

}  // namespace mol_sim

#endif