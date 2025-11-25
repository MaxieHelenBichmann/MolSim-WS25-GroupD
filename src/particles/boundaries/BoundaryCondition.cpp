#include "particles/boundaries/BoundaryCondition.h"

#include <set>
#include <typeinfo>

#include "particles/container/LinkedCellContainer.h"

namespace mol_sim {

//---------------------------------- BoundaryConditionDeclaration--------------------------------------------    
BoundaryConditionDeclaration& BoundaryConditionDeclaration::setType(BoundaryType type) { 
    this->type = type; 
    return *this; 
}

BoundaryConditionDeclaration& BoundaryConditionDeclaration::setLocation(BoundaryLocation location) { 
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

BoundaryType BoundaryConditionDeclaration::getType() { return type; }
BoundaryLocation BoundaryConditionDeclaration::getLocation() { return location; }
std::optional<double> BoundaryConditionDeclaration::getCounterSigma() { return counter_sigma; }
std::optional<double> BoundaryConditionDeclaration::getCounterEpsilon() { return counter_sigma; }
BoundaryConditionDeclaration::BoundaryConditionDeclaration(BoundaryLocation location, BoundaryType type) 
    : location(location), type(type) {}
BoundaryConditionDeclaration::BoundaryConditionDeclaration(BoundaryLocation location, BoundaryType type, std::optional<double> counter_sigma, std::optional<double> counter_epsilon) 
    : location(location), type(type), counter_sigma(counter_sigma), counter_epsilon(counter_epsilon) {}

//-----------------------------------------BoundaryCondition------------------------------------------------   
template <ParticleContainer containerType>
BoundaryCondition<containerType>::BoundaryCondition(containerType& particles, BoundaryLocation location) 
    : particles(particles) {
        this->location = location;
    }

template <ParticleContainer containerType>
BoundaryCondition<containerType>::BoundaryCondition(containerType* particles, BoundaryLocation location) 
    : particles(*particles) {
        this->location = location;
    }

template <ParticleContainer containerType>
void BoundaryCondition<containerType>::applyBoundary() {
    const std::set<BoundaryLocation> boundary{location};
    for (auto& it = particles->boundaryBegin(boundary); it != particles->boundaryEnd(boundary); it++) {
        if (boundaryConditionApplies(it)) { 
            boundaryStrategy(it);
        } 
    }  
}

}  // namespace mol_sim