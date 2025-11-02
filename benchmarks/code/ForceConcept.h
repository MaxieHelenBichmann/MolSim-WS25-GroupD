#ifndef FORCE_CONCEPT_H
#define FORCE_CONCEPT_H

#include <concepts>

#include "particles/Particle.h"

namespace mol_sim {

template <typename C>
concept ForceConcept = requires(C c, const Particle& p1, const Particle& p2) {
    { c.calculateForce(p1, p2) } -> std::same_as<Vector<double, 3>>;
};

}  // namespace mol_sim

#endif
