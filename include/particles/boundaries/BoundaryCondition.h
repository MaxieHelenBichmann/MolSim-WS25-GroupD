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
   * @brief Enum for boundary locations in the Linked-Cell Container.
   *
   * UPPER: +z direction (x-y plane at max z)
   * LOWER: -z direction (x-y plane at min z)
   * FRONT: -y direction (x-z plane at min y)
   * BACK: +y direction (x-z plane at max y)
   * LEFT: -x direction (y-z plane at min x)
   * RIGHT: +x direction (y-z plane at max x)
  */
enum class BoundaryLocation : std::uint8_t { UPPER, LOWER, FRONT, BACK, LEFT, RIGHT };

/**
 * @brief Enum of the supported boundary types.
 * 
 * OUTFLOW: Outflow boundary condition: delete particles in halo cells
 * REFLECTING: Reflecting boundary condition: add ghost particles if particle gets too close to boundary 
 */
enum class BoundaryType : std::uint8_t { OUTFLOW, REFLECTING }; 

class BoundaryCondition {
  protected:
    std::optional<double> counter_sigma = std::nullopt;
    std::optional<double> counter_epsilon = std::nullopt;
    BoundaryLocation location;
    BoundaryType type;
  
  public:
    BoundaryCondition(BoundaryLocation location, BoundaryType type);
    BoundaryCondition(BoundaryLocation location, BoundaryType type, std::optional<double> counter_sigma, std::optional<double> counter_epsilon);
    virtual ~BoundaryCondition() = 0;
    virtual void boundaryStrategy(Particle& p) = 0;  
    
    BoundaryType& getType();
    const BoundaryType& getType() const;
    BoundaryLocation& getLocation();
    const BoundaryLocation& getLocation() const;
    std::optional<double>& getCounterSigma();
    const std::optional<double>& getCounterSigma() const;
    std::optional<double>& getCounterEpsilon();
    const std::optional<double>& getCounterEpsilon() const;
};

}  // namespace mol_sim

#endif