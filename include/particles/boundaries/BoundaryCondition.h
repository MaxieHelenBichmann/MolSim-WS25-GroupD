#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H

#include <vector>
#include "particles/container/LinkedCellContainer.h"

namespace mol_sim {

enum BoundaryConditions : uint8_t { OUTFLOW, REFLECTING, PERIODIC };

class BoundaryCondition {
    protected:
    LinkedCellContainer& particles;
    virtual void boundaryStrategy(Particle& p) = 0;

    public:
    BoundaryCondition(LinkedCellContainer& particles) : particles(particles) {}
    virtual ~BoundaryCondition() = default;
    /***
    * @brief Iterates over the boundary cells and applies the boundary condtion if necessary.
    * 
    * (Comment below just for clarity. Will be removed in final product)
      we are only gonna iterate over the boundary cells. Yes i know that given
      a large enough velocity any particle in any inner cell could also "hit" or even
      cross the boundary within one timestep. However that can be fixed by making
      delta_t smaller. So yes that means more fiddling with delta_t for a given simulation
      so the boundary conditions hold but that to me seems like an ok tradeoff for the faster
      run-time. I think that on average this would actually by O(1) assuming constant density
      of the particles
    */
    void applyBoundary(Particle& p) {
        if (!particles.isOnBoundary(p)) { return; }
        boundaryStrategy(p); 
    }
    
    virtual void clean() = 0;
};

} // namespace mol_sim

#endif