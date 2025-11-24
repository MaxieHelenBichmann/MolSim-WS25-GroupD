#include "particles/boundaries/BoundaryCondition.h"

namespace mol_sim {

//---------------------------------- BoundaryConditionDeclaration--------------------------------------------    
BoundaryConditionDeclaration& BoundaryConditionDeclaration::set_type(BoundaryConditionType type) { this->type = type; return *this; }
BoundaryConditionDeclaration& BoundaryConditionDeclaration::set_location(BoundaryType location) { this->location = location; return *this; }
BoundaryConditionDeclaration& BoundaryConditionDeclaration::set_counter_sigma(double counter_sigma) { this->counter_sigma.value() = counter_sigma; return *this; }
BoundaryConditionDeclaration& BoundaryConditionDeclaration::set_counter_epsilon(double counter_epsilon) { this->counter_epsilon.value() = counter_epsilon; return *this; }
BoundaryConditionType BoundaryConditionDeclaration::get_type() { return type; }
BoundaryType BoundaryConditionDeclaration::get_location() { return location; }
std::optional<double> BoundaryConditionDeclaration::get_counter_sigma() { return counter_sigma; }
std::optional<double> BoundaryConditionDeclaration::get_counter_epsilon() { return counter_sigma; }
BoundaryConditionDeclaration::BoundaryConditionDeclaration(BoundaryType location, BoundaryConditionType type) : location(location), type(type) {}

//-----------------------------------------BoundaryCondition------------------------------------------------   
template<ParticleContainer containerType>
BoundaryCondition<containerType>::BoundaryCondition(containerType& particles, BoundaryType location) : particles(particles), location(location) {}

template<ParticleContainer containerType>
void BoundaryCondition<containerType>::applyBoundary() {
    if (typeid(containerType) == typeid(LinkedCellContainer)) {
        const std::set<BoundaryType> boundary = { location };
        for (auto& it = particles.boundaryBegin(boundary); it != particles.boundaryEnd(boundary); it++) {
            if (boundaryConditionApplies(it)) { 
                boundaryStrategy(it);
            } 
        }  
    }
}

} //namespace mol_sim