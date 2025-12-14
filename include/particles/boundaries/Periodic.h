#ifndef PERIODIC_H
#define PERIODIC_H

#include "particles/boundaries/Boundary.h"
#include "particles/container/LinkedCellContainer.h"

namespace mol_sim {

class Periodic : public Boundary {
    double cutoff;
    /**
     * @brief The size of the corners in the domain. This is needed so the corners are copied correctly
     *
     * IMPORTANT: With this we assume that all boundaries on the same axis have the same width.
     * I.e. LEFT and RIGHT = width1, UPPER and LOWER = width2, FRONT and BACK = width3.
     */
    std::array<double, 3> halo_dimension;

    /**
     * @brief Returns true if the particle is in a corner of the simulation domain, false otherwise.
     * 
     * @param p The particle to be examined.
     * @return true if the particle is in a corner of the simulation domain.
     * @return false otherwise.
     */
    [[nodiscard]] bool isInHalo(R3 x) const noexcept;

    [[nodiscard]] bool isOnBoundary(R3 x) const noexcept;

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
