#ifndef PERIODIC_H
#define PERIODIC_H

#include <memory>

#include "particles/boundaries/Boundary.h"
#include "particles/container/LinkedCellContainer.h"
#include "utils/Vector.h"

namespace mol_sim {

class Periodic : public Boundary {
    /**
     * @brief The size of the corners in the domain. This is needed so the corners are copied correctly
     *
     * IMPORTANT: With this we assume that all boundaries on the same axis have the same width.
     * I.e. LEFT and RIGHT = width1, UPPER and LOWER = width2, FRONT and BACK = width3.
     */
    std::array<double, 3> halo_dimension;

    static constexpr const std::array<uint8_t, 62> MIRROR_IDX_LOOKUP = {
        //-----------------------------3D--------------------------------
        25, 17,  7, 15, 26,  8,  6, 24, //LEFT
        19, 11,  1,  9, 20,  2,  0, 18, //RIGHT
        25, 23, 19, 21, 26, 20, 18, 24, //FRONT
         7,  5,  1,  3,  8,  2,  0,  6, //BACK
        15, 21,  9,  3, 24, 18,  0,  6, //UPPER
        17, 23, 11,  5, 26, 20,  2,  8, //LOWER
        16, 10, 22, 4, 12, 15, //non-corner / edges mirror indicies in order: LEFT, RIGHT, FRONT, BACK, UPPER, LOWER
        
        //-----------------------------2D--------------------------------
        25, 19, 1, 7, //2D corners in order: lower-left, lower-right, upper-right, upper-left
        16, 10, 22, 4 //2D non-corner mirror indicies in order: LEFT, RIGHT, FRONT, BACK
    };

    /**
     * TODO: make this static const
     */
    //z inc -> x inc -> y inc
    std::array<R3, 27> shift_lookup;

    bool is2D = false; //NOLINT

    void teleportParticleIfOOB(Particle& p);
    
    [[nodiscard]] bool isOnBoundary(R3 x, size_t axis, int sign) const noexcept;
    
    [[nodiscard]] std::vector<Particle> mirrorParticle(Particle& p); 

    [[nodiscard]] size_t getIdx(Particle& p);
    
    static void addMirrorParticle(const R3& mirrorLocation, size_t mirrorIdx, Particle& p, std::vector<Particle>& mirrored_particles);

    void setShiftLookup(R3 domain_size);

   public:
    Periodic(BoundaryLocation location, R3 domain_size, double cutoff, bool is2D) noexcept;
    ~Periodic() override = default;

    /**
     * @brief The applyBoundary routine of Periodic. Puts particles in halo cells into respective
     * mirrored boundary cells and then copies particles in boundary cells to the respective halo cells.
     * 
     * Particles that are:
     * - in the interior of the domain:        will be ignored
     * - in the boundary region of the domain: will be mirrored
     * - exactly on a boundary:                will be mirrored (IDK!!!)
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
