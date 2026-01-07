#ifndef PERIODIC_OLD_H
#define PERIODIC_OLD_H

#include <memory>

#include "particles/boundaries/Boundary.h"
#include "particles/container/LinkedCellContainer.h"

namespace mol_sim {

class PeriodicOld : public Boundary {
    /**
     * @brief Used to check if a particle lies within the halo of the domain.
     * 
     * @note If domain_size is not divisible by the cutoff radius in at least one dimension
     * doing it this way will flag some particles as "in the halo" of a SimpleContainer
     * even though they're not. This is additional overhead that can (and should) be eliminated
     * (especially if it creates a bunch of overhead). In the interest of readability and avoiding
     * more weird branches (if container == SimpleContainer...) and entanglement I decided to 
     * do it like this.
     */
    static std::unique_ptr<LinkedCellContainer> checker;
    
    /**
     * @brief The size of the corners in the domain. This is needed so the corners are copied correctly
     *
     * IMPORTANT: With this we assume that all boundaries on the same axis have the same width.
     * I.e. LEFT and RIGHT = width1, UPPER and LOWER = width2, FRONT and BACK = width3.
     */
    std::array<double, 3> halo_dimension;

    bool is2D = false; //NOLINT

    /**
     * @brief Returns true if the particle is in a corner of the simulation domain, false otherwise.
     * 
     * @param p The particle to be examined.
     * @return true if the particle is in a corner of the simulation domain.
     * @return false otherwise.
     */
    [[nodiscard]] bool isInHalo(R3 x) const noexcept;

    [[nodiscard]] bool isOnBoundary(R3 x, size_t axis, int sign) const noexcept;

    void teleportParticleIfOOB(Particle& p);

    std::vector<Particle> mirrorParticle(Particle& p);

    void updateOffset(R3& offset, size_t i);

   public:
    PeriodicOld(BoundaryLocation location, R3 domain_size, double cutoff, size_t dimensions) noexcept;
    ~PeriodicOld() override = default;

    /**
     * @brief The applyBoundary routine of PeriodicOld. Puts particles in halo cells into respective
     * mirrored boundary cells and then copies particles in boundary cells to the respective halo cells.
     * 
     * Particles that are:
     * - in the interior of the domain:        will be ignored
     * - in the boundary region of the domain: will be mirrored
     * - exactly on a boundary:                will be mirrored
     * - in the halo region of the domain:     will be teleported
     * - beyond the halo region of the domain: will be ignored
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
