#include "particles/boundaries/BoundaryCondition.h"

#include <set>
#include <typeinfo>

#include "particles/container/LinkedCellContainer.h"

namespace mol_sim {

//---------------------------------- BoundaryConditionDeclaration--------------------------------------------
BoundaryConditionDeclaration& BoundaryConditionDeclaration::setType(BoundaryConditionType type) {
    this->type = type;
    return *this;
}
BoundaryConditionDeclaration& BoundaryConditionDeclaration::setLocation(BoundaryType location) {
    this->location = location;
    return *this;
}
BoundaryConditionDeclaration& BoundaryConditionDeclaration::setCounterSigma(double counter_sigma) {
    this->counter_sigma.value() = counter_sigma;
    return *this;
}
BoundaryConditionDeclaration& BoundaryConditionDeclaration::setCounterEpsilon(double counter_epsilon) {
    this->counter_epsilon.value() = counter_epsilon;
    return *this;
}
BoundaryConditionType BoundaryConditionDeclaration::getType() { return type; }
BoundaryType BoundaryConditionDeclaration::getLocation() { return location; }
std::optional<double> BoundaryConditionDeclaration::getCounterSigma() { return counter_sigma; }
std::optional<double> BoundaryConditionDeclaration::getCounterEpsilon() { return counter_sigma; }
BoundaryConditionDeclaration::BoundaryConditionDeclaration(BoundaryType location, BoundaryConditionType type)
    : type(type), location(location) {}

//-----------------------------------------BoundaryCondition------------------------------------------------
template <ParticleContainer containerType>
BoundaryCondition<containerType>::BoundaryCondition(containerType& particles, BoundaryType location)
    : particles(particles), location(location) {}

template <ParticleContainer containerType>
void BoundaryCondition<containerType>::applyBoundary() {
    if (typeid(containerType) == typeid(LinkedCellContainer)) {
        const std::set<BoundaryType> boundary = {location};
        for (auto& it = particles->boundaryBegin(boundary); it != particles->boundaryEnd(boundary); it++) {
            if (boundaryConditionApplies(it)) {
                boundaryStrategy(it);
            }
        }
    }
}

}  // namespace mol_sim