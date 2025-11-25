#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H

#include <optional>
#include <vector>
#include <cstdint>
#include <set>

#include "particles/Particle.h"
#include "particles/ParticleContainer.h"

namespace mol_sim {
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
enum class BoundaryLocation : std::uint8_t { UPPER, LOWER, FRONT, BACK, LEFT, RIGHT };
enum class BoundaryType : std::uint8_t { OUTFLOW, REFLECTING, PERIODIC }; 
class BoundaryConditionDeclaration {
  public:
  BoundaryConditionDeclaration& setType(BoundaryType type);
  BoundaryConditionDeclaration& setLocation(BoundaryLocation location);
  BoundaryConditionDeclaration& setCounterSigma(double counter_sigma);
  BoundaryConditionDeclaration& setCounterEpsilon(double counter_epsilon);

  BoundaryType getType();
  BoundaryLocation getLocation();
  std::optional<double> getCounterSigma();
  std::optional<double> getCounterEpsilon();

    BoundaryConditionDeclaration() = default;

  BoundaryConditionDeclaration(BoundaryLocation location, BoundaryType type);
  BoundaryConditionDeclaration(BoundaryLocation location, BoundaryType type, std::optional<double> counter_sigma, std::optional<double> counter_epsilon);

    ~BoundaryConditionDeclaration() = default;

  protected:
    BoundaryLocation location;
    BoundaryType type;
    std::optional<double> counter_sigma = std::nullopt;
    std::optional<double> counter_epsilon = std::nullopt;
};

template <ParticleContainer containerType>
class BoundaryCondition : public BoundaryConditionDeclaration {
    virtual bool boundaryConditionApplies(Particle& p) = 0;
    virtual void boundaryStrategy(Particle& p) = 0;  
  
  protected:
    containerType& particles;

  public:
    BoundaryCondition(containerType& particles, BoundaryLocation location);
    BoundaryCondition(containerType* particles, BoundaryLocation location);
    virtual ~BoundaryCondition() = 0;
    /***
     * @brief Iterates over the boundary cells and applies the boundary condtion if necessary.
     */
    void applyBoundary();
};

}  // namespace mol_sim

#endif