#ifndef PERIODIC_H
#define PERIODIC_H

#include "particles/boundaries/Boundary.h"
#include "particles/container/LinkedCellContainer.h"

namespace mol_sim {

class Periodic : public Boundary {
    /**
     * @brief Array used to make sure we mirror the corners correctly (and not multiple times). 
     * This is relevant e.g. if only two boundaries are periodic and orthogonal to each other.
     * 0 means we mustn't copy this corner
     * 1 / -1 used to indicate whether this corner has already been copied in the current iteration. 
     * 
     * Mapping of array indicies to corners:
     * idx  (x,y,z) (where x = xmax, y = ymax, z = zmax)
     * -------------------------------------------------
     * 0    (0,0,0)
     * 1    (x,0,0)
     * 2    (0,y,0)
     * 3    (x,y,0)
     * 4    (0,0,z)
     * 5    (x,0,z)
     * 6    (0,y,z)
     * 7    (x,y,z)
     */
    static std::array<short, 6> corners;

    /**
     * @brief Marker to be used in the 'corners' array. Switches from -1
     * to 1 every iteration. Starts at -1. 
     */
    short mark = -1;
    
    /**
     * @brief The size of the corners in the domain. This is needed so the corners are copied correctly
     *
     * IMPORTANT: With this we assume that all boundaries on the same axis have the same width.
     * I.e. LEFT and RIGHT = width1, UPPER and LOWER = width2, FRONT and BACK = width3.
     */
    std::array<double, 3> corner_dimension;

    /**
     * @brief Returns true if the particle is in a corner of the simulation domain, false otherwise.
     * 
     * @param p The particle to be examined.
     * @return true if the particle is in a corner of the simulation domain.
     * @return false otherwise.
     */
    [[nodiscard]] bool isInCorner(Particle& p) const noexcept;

    /**
     * @brief Returns the vector indicating how a particle in the corner with index
     * 'corner_idx' must be shifted so it ends up at the mirrored corner.
     * 
     * @return R3 the vector indicating how a particle in the corner with index
     * 'corner_idx' must be shifted so it ends up at the mirrored corner.
     */
    [[nodiscard]] R3 getShift(size_t corner_idx) const noexcept;

   public:
    Periodic(BoundaryLocation location, R3 domain_size, double cutoff) noexcept;
    ~Periodic() override = default;

    /**
     * @brief The applyBoundary routine of Periodic. Puts particles in halo cells into respective
     * mirrored boundary cells and then copies particles in boundary cells to the respective halo cells.
     * 
     * @param p The particle the boundary is to be applied on.
     * @param force The force source that should be used in the boundary condition (in this case irrelevant).
     * @return The copied particles to be added to the halo cells of the container.
     */
    std::optional<std::vector<Particle>> applyBoundary(Particle& p,
                       [[maybe_unused]] const ForceSource& force) noexcept override;
};

}  // namespace mol_sim

#endif
