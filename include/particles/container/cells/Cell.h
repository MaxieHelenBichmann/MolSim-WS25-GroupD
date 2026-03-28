#ifndef CELL_H
#define CELL_H

#include <array>
#include <cstdint>
#include <unordered_set>
#include <vector>

#include "utils/Vector.h"

namespace mol_sim {

/**
 * @brief Enum for cell types in the Linked-Cell Container.
 */
enum class CellType : std::uint8_t { INNER, BOUNDARY, HALO };

/**
 * @brief Cell with specific boundaries and type in the Linked-Cell container.
 *
 * Unites methods for access and modification of Particles in the cell.
 */
class Cell {
    /**
     * std::unordered_set storing the indices (of the std::vector data from the LinkedCellContainer) to all Particles in
     * the cell.
     */
    std::unordered_set<size_t> indices;
    /**
     * Type of the cell. [INNER, BOUNDARY, HALO]
     */
    [[maybe_unused]] CellType type = CellType::INNER;
    /**
     * Boundaries of the cell. [xmin, xmax, ymin, ymax, zmin, zmax]
     */
    std::array<double, 6> bounds;

    /**
     * Cache for sorted particle indices to provide stable iteration.
     */
    mutable std::vector<size_t> sorted_cache;
    /**
     * Flag indicating whether the cache is dirty and needs to be updated.
     */
    mutable bool cache_dirty = true;

    /**
     * @brief Update the sorted cache if it is dirty.
     */
    void updateCache();

   public:
    /**
     * @brief Constructor, initializing a Cell with type and boundaries.
     */
    Cell(CellType cell_type, std::array<double, 6> bounds) noexcept;

    /**
     * @brief Adding an index (of the std::vector data from the LinkedCellContainer) to an existing Particle to the
     * cell.
     *
     * @param idx Index of already constructed Particle.
     */
    void addParticle(size_t idx);
    /**
     * @brief Remove a index (of the std::vector data from the LinkedCellContainer) to a Particle from the cell.
     *
     * @param idx Index of the Particle Index to remove.
     */
    void removeParticle(size_t idx) noexcept;
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
    void clear() noexcept;
    /**
     * @brief Access the set of Particle indices in the cell.
     *
     * @return Reference to the set of Particle indices.
     */
    std::unordered_set<size_t>& particles() noexcept { return indices; };

    /**
     * @brief Access the const set of Particle indices in the cell.
     *
     * @return Reference to the set of Particle indices.
     */
    [[nodiscard]] const std::unordered_set<size_t>& particles() const noexcept { return indices; };

    /**
     * @brief Stable iterator to the beginning of the sorted particle indices.
     *
     * @return Iterator to the beginning of the sorted particle indices.
     */
    std::vector<size_t>::iterator stableIteratorBegin() {
        if (cache_dirty) {
            updateCache();
        }
        return sorted_cache.begin();
    }

    /**
     * @brief Stable iterator to the end of the sorted particle indices.
     *
     * @return Iterator to the end of the sorted particle indices.
     */
    std::vector<size_t>::iterator stableIteratorEnd() {
        if (cache_dirty) {
            updateCache();
        }
        return sorted_cache.end();
    }

    /**
     * @brief Stable iterator to the beginning of the sorted particle indices.
     *
     * @return Iterator to the beginning of the sorted particle indices.
     */
    [[nodiscard]] std::vector<size_t>::const_iterator stableIteratorBegin() const {
        if (cache_dirty) {
            const_cast<Cell*>(this)->updateCache();
        }
        return sorted_cache.begin();
    }

    /**
     * @brief Stable iterator to the end of the sorted particle indices.
     *
     * @return Iterator to the end of the sorted particle indices.
     */
    [[nodiscard]] std::vector<size_t>::const_iterator stableIteratorEnd() const {
        if (cache_dirty) {
            const_cast<Cell*>(this)->updateCache();
        }
        return sorted_cache.end();
    }

    /**
     * @brief Check whether a Particle fits into the cell boundaries.
     *
     * @param x Coordinates to check.
     *
     * @return True if Particle is within cell boundaries.
     */
    [[nodiscard]] bool fits(R3 x) const noexcept;
    /**
     * @brief Returns the number of particles contained within the cell.
     *
     * @return size_t the number of particles contained within the cell.
     */
    [[nodiscard]] size_t size() const noexcept;
    /**
     * @brief Returns the type of the cell (INNER, BOUNDARY, HALO).
     *
     * @return CellType the type of the cell.
     */
    [[nodiscard]] CellType getType() const noexcept { return type; };
};

}  // namespace mol_sim

#endif  // CELL_H
