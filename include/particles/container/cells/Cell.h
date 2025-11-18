#ifndef CELL_H
#define CELL_H

#include <array>
#include <cstdint>
#include <vector>

#include "particles/Particle.h"
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
     * std::vector storing the indices (of the std::vector data from the LinkedCellContainer) to all Particles in
     * the cell.
     */
    std::vector<size_t> indices;
    /**
     * Type of the cell. [INNER, BOUNDARY, HALO]
     */
    [[maybe_unused]] CellType type = CellType::INNER;
    /**
     * Boundaries of the cell. [xmin, xmax, ymin, ymax, zmin, zmax]
     */
    std::array<double, 6> bounds;

   public:
    /**
     * @brief Constructor, initializing a Cell with type and boundaries.
     */
    Cell(CellType cell_type, std::array<double, 6> bounds);

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
     * @brief Access the vector of Particle pointers in the cell.
     *
     * @return Reference to the vector of Particle pointers.
     */
    std::vector<size_t>& particles();

    /**
     * @brief Access the const vector of Particle pointers in the cell.
     *
     * @return Reference to the vector of Particle pointers.
     */
    [[nodiscard]] const std::vector<size_t>& particles() const;

    size_t operator[](size_t idx);

    /**
     * @brief Check whether a Particle fits into the cell boundaries.
     *
     * @param x Coordinates to check.
     *
     * @return True if Particle is within cell boundaries.
     */
    [[nodiscard]] bool fits(R3 x) const;
};

}  // namespace mol_sim

#endif  // CELL_H