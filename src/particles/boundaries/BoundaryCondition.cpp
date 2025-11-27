#include "particles/boundaries/BoundaryCondition.h"

#include <set>
#include <typeinfo>

#include "particles/container/LinkedCellContainer.h"

namespace mol_sim {

BoundaryCondition::BoundaryCondition(BoundaryLocation location, BoundaryType type) 
    : location(location), type(type) {}
BoundaryCondition::BoundaryCondition(BoundaryLocation location, BoundaryType type, std::optional<double> counter_sigma, std::optional<double> counter_epsilon) 
    : location(location), type(type), counter_sigma(counter_sigma), counter_epsilon(counter_epsilon) {}

BoundaryType& BoundaryCondition::getType() { return type; }
const BoundaryType& BoundaryCondition::getType() const { return type; }
BoundaryLocation& BoundaryCondition::getLocation() { return location; }
const BoundaryLocation& BoundaryCondition::getLocation() const { return location; }
std::optional<double>& BoundaryCondition::getCounterSigma() { return counter_sigma; }
const std::optional<double>& BoundaryCondition::getCounterSigma() const { return counter_sigma; }
std::optional<double>& BoundaryCondition::getCounterEpsilon() { return counter_sigma; }
const std::optional<double>& BoundaryCondition::getCounterEpsilon() const { return counter_sigma; }

}  // namespace mol_sim