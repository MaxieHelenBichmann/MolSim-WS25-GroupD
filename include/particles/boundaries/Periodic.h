#ifndef PERIODIC_H
#define PERIODIC_H

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

    /**
     * @brief Lookup array that indicates to which mirror locations the particle needs to
     * be mirrored. As to what each of the numbers represent, see the documentation in Particle.h/mirrorLocations
     * or see the slides of Assignment4/Periodic
     */
    static constexpr const std::array<uint8_t, 62> MIRROR_IDX_LOOKUP = {
        //-----------------------------3D--------------------------------
        25, 17, 7, 15, 26, 8, 6, 24,     // LEFT
        19, 11, 1, 9, 20, 2, 0, 18,      // RIGHT
        25, 23, 19, 21, 26, 20, 18, 24,  // FRONT
        7, 5, 1, 3, 8, 2, 0, 6,          // BACK
        15, 21, 9, 3, 24, 18, 0, 6,      // UPPER
        17, 23, 11, 5, 26, 20, 2, 8,     // LOWER
        16, 10, 22, 4, 12, 14,  // non-corner / edges mirror indicies in order: LEFT, RIGHT, FRONT, BACK, UPPER, LOWER

        //-----------------------------2D--------------------------------
        25, 19, 1, 7,  // 2D corners in order: lower-left, lower-right, upper-right, upper-left
        16, 10, 22, 4  // 2D non-corner mirror indicies in order: LEFT, RIGHT, FRONT, BACK
    };

    /**
     * @brief An array to look up the shifts one would need to add to p.getX()
     * to get to the coordinates of the i-th mirror location of p.
     *
     * For examples p.getX() + shift_lookup[i] = location of i-th mirror location of p.
     * TODO: make this static const
     */
    std::array<R3, 27> shift_lookup;

    /**
     * @brief Indicates the number of dimensions the simulation domain has.
     * Supported values are only 2 and 3.
     */
    size_t dimensions;

    /**
     * @brief Teleports the given particle to the opposite side of the domain if it's
     * inside the halo region, if need be.
     *
     * @param p The particle to be teleported, if need be.
     */
    void teleportParticleIfOOB(Particle& p);

    /**
     * @brief Determines whether a particle lies within a boundary region of the simulation domain.
     *
     * @param x The position of the particle.
     * @param axis The axis of the boundary {0 = x, 1 = y, 2 = z}.
     * @param sign The sign of the boundary {-1 = LEFT, LOWER, FRONT   1 = RIGHT, UPPER, BACK}
     * @return true If the particle is within the boundary region of the simulation domain.
     * @return false Otherwise.
     */
    [[nodiscard]] bool isOnBoundary(R3 x, size_t axis, int sign) const noexcept;

    /**
     * @brief Creates all necessary mirror particles for the given particle
     * if the particle lies within the simulation domain.
     *
     * CAREFUL: This method assumes that the passed particle lies within the
     * boundary region of the simulation domain
     *
     * @param p The particle to be mirrored. IMPORTANT: The particle must be within the boundary region of the domain.
     * @return std::vector<Particle> A vector containing the created mirror particles.
     */
    void mirrorParticle(Particle& p);

    /**
     * @brief Returns the index of the corner or edge if the passed particle is contained in either
     * of them, or 8, if it isn't.
     *
     * This is how the corners / edges are mapped to indicies for 2D / 3D:
     *
     *    2D (this is the entire simulation domain)
     *          3-------2
     *          |       |
     *          |       |
     *          |_______|
     *          0       1
     *
     *
     *    3D (this is just a single 3D boundary)
     *              3
     *          7-------6
     *          |       |
     *        0 |       | 2
     *          |_______|
     *          4   1   5
     *
     *    For 3D:
     *      If this boundary is LEFT:  horizontal axis = y, vertical axis = z
     *      If this boundary is RIGHT: horizontal axis = y, vertical axis = z
     *      If this boundary is FRONT: horizontal axis = x, vertical axis = z
     *      If this boundary is BACK:  horizontal axis = x, vertical axis = z
     *      If this boundary is UPPER: horizontal axis = x, vertical axis = z
     *      If this boundary is LOWER: horizontal axis = x, vertical axis = z
     *
     * @param p The particle to be mirrored
     * @return size_t The index of the corner / edge the particle is contained inside of.
     * 8 if its not in any corner / edge
     */
    [[nodiscard]] size_t getIdx(Particle& p);

    /**
     * @brief Adds a mirror particle at the given position to the
     * mirror_particles vector for the passed particle.
     *
     * @param mirrorLocation The location of the mirror particle
     * @param mirrorIdx The mirror location index of the mirror particle.
     * See Particle.h->mirrorLocations or Assignment4/Periodic slides for explanation.
     * @param p The particle the mirror particle belongs to.
     * @param mirrored_particles The vector of previously accumulated mirror particles.
     */
    void addMirrorParticle(const R3& mirrorLocation, size_t mirrorIdx, Particle& p);

    /**
     * @brief Initializes the shift lookup table.
     *
     * @param domain_size The domain size of the simulation domain.
     */
    void setShiftLookup(R3 domain_size);

   public:
    Periodic(BoundaryLocation location, R3 domain_size, double cutoff, size_t dimensions) noexcept;
    ~Periodic() override = default;

    /**
     * @brief The applyBoundary routine of Periodic. Puts particles in halo cells into respective
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
    void applyBoundary(Particle& p, [[maybe_unused]] const ForceSource& force) noexcept override;
};

}  // namespace mol_sim

#endif
