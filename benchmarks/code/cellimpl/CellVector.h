#ifndef CELL_VECTOR_H
#define CELL_VECTOR_H

#include <array>
#include <vector>

#include "utils/Vector.h"

namespace mol_sim {

/**
 * @brief Cell in the Linked-Cell container.
 *
 * Implemented with std::vector for storage of particle indices for benchmarking purposes.
 * Only implemented performance relevant methods, to compare to other Cell implementations.
 *
 * Expect better performance for iterators due to better cache locality, but worse performance for modifications.
 */
class CellVector {
    /**
     * std::vector storing the indices to all Particles in the cell.
     */
    std::vector<size_t> indices;

    /**
     * Boundaries of the cell. [xmin, xmax, ymin, ymax, zmin, zmax]
     */
    std::array<double, 6> bounds;

   public:
    /**
     * @brief Constructor, initializing a Cell with type and boundaries.
     */
    CellVector(std::array<double, 6> bounds);

    /**
     * @brief Adding an index to an existing Particle to the
     * cell.
     *
     * @param idx Index of already constructed Particle.
     */
    void addParticle(size_t idx);
    /**
     * @brief Remove a index to a Particle from the cell.
     *
     * @param idx Index of the pointer to remove.
     */
    void removeParticle(size_t idx);
    /**
     * @brief Update the index of a Particle in the cell.
     *
     * @param old_idx Old index of the Particle.
     * @param new_idx New index of the Particle.
     */
    void updateParticleIndex(size_t old_idx, size_t new_idx);
    /**
     * @brief Clear the entire cell, destructing no Particles.
     */
    void clear();
    /**
     * @brief Access the set of Particle pointers in the cell.
     *
     * @return Reference to the set of Particle pointers.
     */
    std::vector<size_t>& particles();

    /**
     * @brief Access the const set of Particle pointers in the cell.
     *
     * @return Reference to the set of Particle pointers.
     */
    [[nodiscard]] const std::vector<size_t>& particles() const;

    /**
     * @brief Check whether a Particle fits into the cell boundaries.
     *
     * @param x Coordinates to check.
     *
     * @return True if Particle is within cell boundaries.
     */
    [[nodiscard]] bool fits(R3 x) const;
    /**
     * @brief Returns the number of particles contained within the cell.
     *
     * @return size_t the number of particles contained within the cell.
     */
    size_t size();
};

}  // namespace mol_sim

#endif  // CELL_VECTOR_H