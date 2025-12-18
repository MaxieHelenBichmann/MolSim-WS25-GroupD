#ifndef LENNARDJONES_ABSTRACT_H
#define LENNARDJONES_ABSTRACT_H
#include "AbstractForce.h"
namespace mol_sim {

class lennardJonesAbstract : public ForceAbstract {
    [[nodiscard]] Vector<double, 3> applyForce(const Particle& p1, const Particle& p2) const override;
};
}  // namespace mol_sim
#endif
