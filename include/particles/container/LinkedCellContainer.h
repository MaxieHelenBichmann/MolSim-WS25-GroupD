#ifndef LINKEDCELL_CONTAINER_H
#define LINKEDCELL_CONTAINER_H

#include <cstdint>
#include <vector>

#include "particles/ParticleContainer.h"

namespace mol_sim {

/**
 * @brief Linked-Cell Container for Particles
 *
 * This container implements the concept ParticleContainer.
 * It stored the Particles in linked cells to optimize proximity queries. However,
 * it is also possible to simply iterate over all Particles, because they are stored in a linearized vector.
 *
 */
class LinkedCellContainer {
    /**
     * @brief Cell with specific boundaries and type in the Linked-Cell container.
     *
     * Unites methods for access and modification of Particles in the cell.
     */
    struct Cell {
        enum class CellType : std::uint8_t { INNER, BOUNDARY, HALO };

        /**
         * @brief Constructor, initializing a Cell with type and boundaries.
         */
        Cell(CellType cell_type, std::array<double, 6> bounds);

        /**
         * @brief Adding a pointer to an existing Particle to the cell.
         *
         * @param value Pointer to already constructed Particle.
         */
        void addParticle(Particle* value);
        /**
         * @brief Remove a pointer at a specific index from the cell.
         *
         * @param idx Index of the pointer to remove.
         */
        void removeParticle(size_t idx);
        /**
         * @brief Remove a specific pointer from the cell.
         *
         * @param p Pointer to remove.
         */
        void removeParticle(Particle* p);
        /**
         * @brief Clear the entire cell, destructing no Particles.
         */
        void clear();
        /**
         * @brief Get the cell type (INNER, BOUNDARY, HALO)
         *  
         * @return CellType 
         */
        CellType getType();
        /**
         * @brief Returns the number of particles contained in the cell 
         * 
         * @return size_t 
         */
        size_t size();

        Particle* operator[](size_t idx);

        /**
         * @brief Check whether a Particle fits into the cell boundaries.
         *
         * @param p Pointer to Particle to check.
         *
         * @return True if Particle is within cell boundaries.
         */
        bool fits(Particle* p) const;

       private:
        /**
         * std::vector storing pointers to all Particles in the cell.
         */
        std::vector<Particle*> data;
        /**
         * Type of the cell. [INNER, BOUNDARY, HALO]
         */
        [[maybe_unused]] CellType type = CellType::INNER;
        /**
         * Boundaries of the cell. [xmin, xmax, ymin, ymax, zmin, zmax]
         */
        std::array<double, 6> bounds;
    };

    /**
     * @brief Switch a Particle from one cell to another. [HELPER FUNCTION]
     *
     * @param p Reference to Particle to switch.
     * @param old_cell_idx Index of the old cell.
     * @param new_cell_idx Index of the new cell.
     */
    void switchCell(Particle& p, size_t old_cell_idx, size_t new_cell_idx);

    /**
     * @brief Switch a Particle from one cell to another. [HELPER FUNCTION]
     *
     * @param p Pointer to Particle to switch.
     * @param old_cell_idx Index of the old cell.
     * @param new_cell_idx Index of the new cell.
     */
    void switchCell(Particle* p, size_t old_cell_idx, size_t new_cell_idx);

    /**
     * @brief Switch a Particle from one cell to another. [HELPER FUNCTION]
     *
     * @param p Index of the Particle in the old cell, which should be switched.
     * @param old_cell_idx Index of the old cell.
     * @param new_cell_idx Index of the new cell.
     */
    void switchCell(size_t p, size_t old_cell_idx, size_t new_cell_idx);

    /**
     * @brief Find the index of the cell, in which a Particle is located. [HELPER FUNCTION]
     *
     * @param p Pointer to the Particle.
     *
     * @return Index of the cell, in which the Particle is located.
     */
    size_t findCellIndex(Particle* p);

    /**
     * std::vector storing all Particles in the container.
     */
    std::vector<Particle> data;

    /**
     * std::vector storing all cells of the container.
     */
    std::vector<Cell> cells;

    /**
     * Vector storing the global bounds of the domain of this container.
     * The Linked-Cell container assumes a cuboidal domain from (0,0,0) to domain_size.
     */
    R3 domain_size;

    /**
     * Number of cells in each dimension.
     */
    std::array<size_t, 3> num_cells;

    /**
     * Cutoff radius used for cell size.
     */
    double cutoff_radius;

   public:
    // constructors

    /**
     * @brief Constructor, initializing a LinkedCellContainer.
     */
    LinkedCellContainer(R3 domain_size, double cutoff_radius);
    
    // retrieve data

    /**
     * @brief Check whether a (potential) Particle fits into the domain.
     *
     * @param v Coordinates of the Particle to check.
     *
     * @return True if position is within domain boundaries.
     */
    [[nodiscard]] bool fits(R3 v) const;

    Particle& operator[](size_t idx);
    const Particle& operator[](size_t idx) const;

    /**
     * @brief Returns size of container.
     *
     * @return Number of Particles in the container.
     */
    [[nodiscard]] size_t size() const;

    /**
     * @brief Queries whether the container is empty.
     *
     * @return True of container contains no Particles.
     */
    [[nodiscard]] bool empty() const;

    // modify

    /**
     * @brief Clear the entire container, destructing the Particles.
     */
    void clear();

    /**
     * @brief Reserve memory for n Particles, not yet constructing them.
     *
     * @param n Number of Particles, for which memory should be reserved.
     */
    void reserve(size_t n);

    /**
     * @brief Adding an already existing Particle to the container by copying it.
     *
     * @param value Already constructed Particle.
     */
    void addParticle(Particle&& value);

    /**
     * @brief Adding an already existing const Particle to the container by copying it.
     *
     * @param value Already constructed Particle.
     */
    void addParticle(const Particle& value);

    /**
     * @brief Directly constructing a Particle with its required parameters, including type, in-place.
     *
     * @param x_arg Initial coordinates of the Particle.
     * @param v_arg Initial velocities of the Particle.
     * @param m_arg Mass of the Particle.
     * @param epsilon_arg Epsilon of the Particle.
     * @param sigma_arg Sigma of the Particle.
     */
    void addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg);

    /**
     * @brief Directly constructing a Particle with its required parameters, including type, in-place.
     *
     * @param x_arg Initial coordinates of the Particle.
     * @param v_arg Initial velocities of the Particle.
     * @param m_arg Mass of the Particle.
     * @param epsilon_arg Epsilon of the Particle.
     * @param sigma_arg Sigma of the Particle.
     * @param type Type of the Particle.
     */
    void addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg, int type);

    // iterators

    /**
     * @brief Mutable iterator of the raw data (no logic).
     *
     * @return Mutable iterator to the first element of the container.
     */
    std::vector<Particle>::iterator begin();

    /**
     * @brief Const iterator of the raw data (no logic).
     *
     * @return Const iterator to the first element of the container.
     */
    [[nodiscard]] std::vector<Particle>::const_iterator begin() const;

    /**
     * @brief Const iterator of the raw data (no logic).
     *
     * @return Const iterator to the first element of the container.
     */
    [[nodiscard]] std::vector<Particle>::const_iterator cbegin() const;

    /**
     * @brief Mutable iterator of the raw data (no logic).
     *
     * @return Mutable past-the-end iterator of the container.
     */
    std::vector<Particle>::iterator end();

    /**
     * @brief Const iterator of the raw data (no logic).
     *
     * @return Const past-the-end iterator of the container.
     */
    [[nodiscard]] std::vector<Particle>::const_iterator end() const;

    /**
     * @brief Const iterator of the raw data (no logic).
     *
     * @return Const past-the-end iterator of the container
     */
    [[nodiscard]] std::vector<Particle>::const_iterator cend() const;

    class proximity_iterator {
       public:
        proximity_iterator() = default;
    };
    // static_assert(std::forward_iterator<proximity_iterator>);

    class const_proximity_iterator {
       public:
        const_proximity_iterator() = default;
    };
    // static_assert(std::forward_iterator<const_proximity_iterator>);

    /**
     * @brief Mutable Iterator over particles in proximity.
     *
     * @param center Center point to check proximity from (position of the particle).
     * @param radius Radius within which to consider particles in proximity.
     * @param offset Offset from the beginning of the container to start the iteration (used for N3L optimization).
     *
     * @return Mutable iterator to the first particle within the given radius of the center.
     */
    [[nodiscard]] proximity_iterator proximityBegin(R3 center, double radius, size_t offset = 0);

    /**
     * @brief Mutable Iterator over particles in proximity.
     *
     * @param center Center point to check proximity from (position of the particle).
     * @param radius Radius within which to consider particles in proximity.
     *
     * @return Mutable iterator after the last particle within the given radius of the center.
     */
    [[nodiscard]] proximity_iterator proximityEnd(R3 center, double radius);

    /**
     * @brief Const Iterator over particles in proximity.
     *
     * @param center Center point to check proximity from (position of the particle).
     * @param radius Radius within which to consider particles in proximity.
     * @param offset Offset from the beginning of the container to start the iteration (used for N3L optimization).
     *
     * @return Const iterator to the first particle within the given radius of the center.
     */
    [[nodiscard]] const_proximity_iterator proximityBegin(R3 center, double radius, size_t offset = 0) const;

    /**
     * @brief Const Iterator over particles in proximity.
     *
     * @param center Center point to check proximity from (position of the particle).
     * @param radius Radius within which to consider particles in proximity.
     *
     * @return Const iterator after the last particle within the given radius of the center.
     */
    [[nodiscard]] const_proximity_iterator proximityEnd(R3 center, double radius) const;

    /**
     * @brief Checks whether a given particle resides within a boundary cell.
     * 
     * @param p the Particle. 
     * @return true if the particle lies within a boundary cell. 
     * @return false if the particle does NOT lie within a boundary cell.
     */
    [[nodiscard]] bool isOnBoundary(Particle& p);

    /**
     * @brief Removes a given particle from the container. 
     * CAREFUL: This method assumes that no particle is contained more than once within the same container. 
     * 
     * @param p The particle to be removed. 
     */
    void removeParticle(Particle& p);

    /**
     * @brief Returns the domain size, i.e. the vector going from (0, 0, 0) to (x_max, y_max, z_max).
     * 
     * @return R3 The vector containing the domain size.
     */
    [[nodiscard]] R3 getDomainSize();

};

}  // namespace mol_sim

#endif
