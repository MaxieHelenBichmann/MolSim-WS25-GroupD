#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H

#include <vector>
#include "particles/container/LinkedCellContainer.h"

namespace mol_sim {

class BoundaryConditionDeclaration {
  public:
  /**
   * @brief Enum for boundary types in the Linked-Cell Container.
   *
   * UPPER: +z direction (x-y plane at max z)
   * LOWER: -z direction (x-y plane at min z)
   * FRONT: -y direction (x-z plane at min y)
   * BACK: +y direction (x-z plane at max y)
   * LEFT: -x direction (y-z plane at min x)
   * RIGHT: +x direction (y-z plane at max x)
  */
  enum class BoundaryType : std::uint8_t { UPPER, LOWER, FRONT, BACK, LEFT, RIGHT };
  enum class BoundaryConditionType : std::uint8_t { OUTFLOW, REFLECTING, PERIODIC }; 
  BoundaryConditionDeclaration& set_type(BoundaryConditionType type);
  BoundaryConditionDeclaration& set_location(BoundaryType location);
  BoundaryConditionDeclaration& set_counter_sigma(double counter_sigma);
  BoundaryConditionDeclaration& set_counter_epsilon(double counter_epsilon);

  BoundaryConditionType get_type();
  BoundaryType get_location();
  std::optional<double> get_counter_sigma();
  std::optional<double> get_counter_epsilon();

  BoundaryConditionDeclaration() = default;

  BoundaryConditionDeclaration(BoundaryType location, BoundaryConditionType type);

  ~BoundaryConditionDeclaration() = default;

  private:
  BoundaryConditionType type;
  BoundaryType location;
  std::optional<double> counter_sigma = std::nullopt;
  std::optional<double> counter_epsilon = std::nullopt;
};

using BoundaryType = BoundaryConditionDeclaration::BoundaryType;
using BoundaryConditionType = BoundaryConditionDeclaration::BoundaryConditionType;

template <ParticleContainer containerType>
class BoundaryCondition {
    virtual bool boundaryConditionApplies(Particle& p) = 0;
    virtual void boundaryStrategy(Particle& p) = 0;  
  protected:
    containerType& particles;
    BoundaryType location;

  public:
    BoundaryCondition(const BoundaryCondition& other) { 
      this->particles = other.particles;
      this->location = other.location;
      this->counter_sigma = other.counter_sigma;
      this->counter_epsilon = other.counter_epsilon;
    }
    BoundaryCondition(containerType& particles, BoundaryType location);
    ~BoundaryCondition() = default;
    /***
    * @brief Iterates over the boundary cells and applies the boundary condtion if necessary.
    */
    void applyBoundary();
};

} // namespace mol_sim

#endif