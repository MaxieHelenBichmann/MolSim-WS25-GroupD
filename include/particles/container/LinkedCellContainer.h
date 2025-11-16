#ifndef LINKEDCELL_CONTAINER_H
#define LINKEDCELL_CONTAINER_H

#include <array>
#include <cstdint>
#include <set>
#include <vector>

#include "particles/ParticleContainer.h"

namespace mol_sim {

/**
 * @brief Enum for boundary types in the Linked-Cell Container.
 *
 * UPPER: +z direction (x-y plane at max z)
 * LOWER: -z direction (x-y plane at min z)
 * FRONT: -y direction (x-z plane at min y)
 * BACK: +y direction (x-z plane at max y)
 * LEFT: -x direction (y-z plane at min x)
 * RIGHT: +x direction (y-z plane at max x)
 */
enum class BoundaryType : std::uint8_t { UPPER, LOWER, FRONT, BACK, LEFT, RIGHT };

/**
 * @brief Enum for cell types in the Linked-Cell Container.
 */
enum class CellType : std::uint8_t { INNER, BOUNDARY, HALO };

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
         * @brief Access the vector of Particle pointers in the cell.
         *
         * @return Reference to the vector of Particle pointers.
         */
        std::vector<Particle*>& particles();

        /**
         * @brief Access the const vector of Particle pointers in the cell.
         *
         * @return Reference to the vector of Particle pointers.
         */
        [[nodiscard]] const std::vector<Particle*>& particles() const;

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
     * @param x Coordinates of the Particle.
     *
     * @return Index of the cell, in which the Particle is located.
     */
    [[nodiscard]] size_t findCellIndex(R3 vec) const;

    /**
     * @brief Find the adjacent cells of a cell. [HELPER FUNCTION]
     *
     * @param cell_idx Index of the cell.
     *
     * @return Vector of pointers to the adjacent cells, including the cell itself.
     */
    std::vector<Cell*> findAdjacentCells(size_t cell_idx);

    /**
     * @brief Find the adjacent cells of a cell. [HELPER FUNCTION]
     *
     * @param cell_idx Index of the cell.
     *
     * @return Vector of const pointers to the adjacent cells, including the cell itself.
     */
    [[nodiscard]] std::vector<const Cell*> findAdjacentCells(size_t cell_idx) const;

    /**
     * @brief Find the boundary or halo cells. [HELPER FUNCTION]
     *
     * @param type Type of the boundary.
     * @param boundary_cells Vector to store the found boundary cells.
     * @param offset Offset to consider for halo cells. (0 = halo cells, 1 = boundary cells)
     */
    void findBoundaryCells(BoundaryType type, std::vector<const Cell*>& boundary_cells, size_t offset = 0) const;

    /**
     * @brief Find the boundary or halo cells. [HELPER FUNCTION]
     *
     * @param type Type of the boundary.
     * @param boundary_cells Vector to store the found boundary cells.
     * @param offset Offset to consider for halo cells. (0 = halo cells, 1 = boundary cells)
     */
    void findBoundaryCells(BoundaryType type, std::vector<Cell*>& boundary_cells, size_t offset = 0);

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
     * Number of cells in each dimension, including halo cells.
     */
    std::array<size_t, 3> num_cells;

    /**
     * Effective edge length of a single cell in each spatial dimension.
     * At least as big as the given cutoff radius, to still only check the immediate neighbors.
     */
    std::array<double, 3> cell_length{};

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
    [[nodiscard]] bool fitsDomain(R3 v) const;

    /**
     * @brief Check whether a (potential) Particle fits into the container (so including halo cells).
     *
     * @param v Coordinates of the Particle to check.
     *
     * @return True if position is within container boundaries.
     */
    [[nodiscard]] bool fitsContainer(R3 v) const;

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

    /**
     * @brief Remove a Particle if it lies in a halo cell.
     *
     * @param p Pointer to the Particle to remove.
     */
    void removeParticle(Particle* p);

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

    /**
     * @brief Iterator that iterates over all particles that apply a force on a given particle.
     *
     * Satisfies the forward iterator concept.
     * Only iterates over particles within the given radius of the center, or all if the radius is infinite.
     *
     * Enables cutoff radius, and iterates only over the particles in the adjacent cells.
     */
    class proximity_iterator {
        Particle** cur;
        Particle** end;
        Particle** cell_end;
        std::vector<Cell*> cells;
        double radius;
        R3 center;

        void inc() {
            if (cur == cell_end) {  // reached end of current cell
                cells.erase(cells.begin());
                cur = cells.front()->particles().data();
                cell_end = cur + cells.front()->particles().size();  // NOLINT
            } else {                                                 // iterate within current cell
                ++cur;                                               // NOLINT
            }
        }

        void satisfy() {
            if (std::isinf(radius)) {
                return;
            }
            while (cur != end && !((center - (*cur)->getX()).euclidNorm() <= radius)) {
                inc();
            }
        }

       public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = Particle;
        using difference_type = std::ptrdiff_t;
        using pointer = Particle*;
        using reference = Particle&;

        proximity_iterator() noexcept : cur(nullptr), end(nullptr), cell_end(nullptr), radius(0.0) {}
        proximity_iterator(R3 center, double radius, Particle** cur, std::vector<Cell*> cells)
            : cur(cur),
              end(cells.back()->particles().data() + cells.back()->particles().size()),         // NOLINT
              cell_end(cells.front()->particles().data() + cells.front()->particles().size()),  // NOLINT
              cells(cells),
              radius(radius),
              center(center) {
            satisfy();
        }

        reference operator*() const { return **cur; }
        pointer operator->() const { return *cur; }

        proximity_iterator& operator++() {
            inc();
            satisfy();
            return *this;
        }

        proximity_iterator operator++(int) {
            proximity_iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        friend bool operator==(const proximity_iterator& a, const proximity_iterator& b) { return a.cur == b.cur; }
        friend bool operator!=(const proximity_iterator& a, const proximity_iterator& b) { return !(a == b); }
    };
    static_assert(std::forward_iterator<proximity_iterator>);

    /**
     * @brief Iterator that iterates over all particles that apply a force on a given particle.
     *
     * Satisfies the forward iterator concept.
     * Only iterates over particles within the given radius of the center, or all if the radius is infinite.
     *
     * Enables cutoff radius, and iterates only over the particles in the adjacent cells.
     */
    class const_proximity_iterator {
        const Particle* const* cur;
        const Particle* const* end;
        const Particle* const* cell_end;
        std::vector<const Cell*> cells;
        double radius;
        R3 center;

        void inc() {
            if (cur == cell_end) {  // reached end of current cell
                cells.erase(cells.begin());
                cur = cells.front()->particles().data();
                cell_end = cur + cells.front()->particles().size();  // NOLINT
            } else {                                                 // iterate within current cell
                ++cur;                                               // NOLINT
            }
        }

        void satisfy() {
            if (std::isinf(radius)) {
                return;
            }
            while (cur != end && !((center - (*cur)->getX()).euclidNorm() <= radius)) {
                inc();
            }
        }

       public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = const Particle;
        using difference_type = std::ptrdiff_t;
        using pointer = const Particle*;
        using reference = const Particle&;

        const_proximity_iterator() noexcept : cur(nullptr), end(nullptr), cell_end(nullptr), radius(0.0) {}
        const_proximity_iterator(R3 center, double radius, const Particle* const* cur, std::vector<const Cell*> cells)
            : cur(cur),
              end(cells.back()->particles().data() + cells.back()->particles().size()),         // NOLINT
              cell_end(cells.front()->particles().data() + cells.front()->particles().size()),  // NOLINT
              cells(cells),
              radius(radius),
              center(center) {
            satisfy();
        }

        reference operator*() const { return **cur; }
        pointer operator->() const { return *cur; }

        const_proximity_iterator& operator++() {
            inc();
            satisfy();
            return *this;
        }

        const_proximity_iterator operator++(int) {
            const_proximity_iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        friend bool operator==(const const_proximity_iterator& a, const const_proximity_iterator& b) {
            return a.cur == b.cur;
        }
        friend bool operator!=(const const_proximity_iterator& a, const const_proximity_iterator& b) {
            return !(a == b);
        }
    };
    static_assert(std::forward_iterator<const_proximity_iterator>);

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

    // boundary and halo iterators

    /**
     * @brief Iterator over particles in halo cells.
     *
     * @param boundary_types Boundary types to specify which halo cells to iterate over. Defaults to all sides.
     *
     * @return Iterator to the first particle within the given halo cells.
     */
    [[nodiscard]] proximity_iterator haloBegin(const std::set<BoundaryType>& boundary_types = {
                                                   BoundaryType::UPPER, BoundaryType::LOWER, BoundaryType::FRONT,
                                                   BoundaryType::BACK, BoundaryType::LEFT, BoundaryType::RIGHT});

    /**
     * @brief Const Iterator over particles in halo cells.
     *
     * @param boundary_types Boundary types to specify which halo cells to iterate over. Defaults to all sides.
     *
     * @return Const iterator to the first particle within the given halo cells.
     */
    [[nodiscard]] const_proximity_iterator haloBegin(const std::set<BoundaryType>& boundary_types = {
                                                         BoundaryType::UPPER, BoundaryType::LOWER, BoundaryType::FRONT,
                                                         BoundaryType::BACK, BoundaryType::LEFT,
                                                         BoundaryType::RIGHT}) const;

    /**
     * @brief Iterator over particles in halo cells.
     *
     * @param boundary_types Boundary types to specify which halo cells to iterate over. Defaults to all sides.
     *
     * @return Iterator after the last particle within the given halo cells.
     */
    [[nodiscard]] proximity_iterator haloEnd(const std::set<BoundaryType>& boundary_types = {
                                                 BoundaryType::UPPER, BoundaryType::LOWER, BoundaryType::FRONT,
                                                 BoundaryType::BACK, BoundaryType::LEFT, BoundaryType::RIGHT});

    /**
     * @brief Const Iterator over particles in halo cells.
     *
     * @param boundary_types Boundary types to specify which halo cells to iterate over. Defaults to all sides.
     *
     * @return Const iterator after the last particle within the given halo cells.
     */
    [[nodiscard]] const_proximity_iterator haloEnd(const std::set<BoundaryType>& boundary_types = {
                                                       BoundaryType::UPPER, BoundaryType::LOWER, BoundaryType::FRONT,
                                                       BoundaryType::BACK, BoundaryType::LEFT,
                                                       BoundaryType::RIGHT}) const;

    /**
     * @brief Iterator over particles in boundary cells.
     *
     * @param boundary_types Boundary types to specify which boundary cells to iterate over. Defaults to all sides.
     *
     * @return Iterator to the first particle within the given boundary cells.
     */
    [[nodiscard]] proximity_iterator boundaryBegin(const std::set<BoundaryType>& boundary_types = {
                                                       BoundaryType::UPPER, BoundaryType::LOWER, BoundaryType::FRONT,
                                                       BoundaryType::BACK, BoundaryType::LEFT, BoundaryType::RIGHT});

    /**
     * @brief Const Iterator over particles in boundary cells.
     *
     * @param boundary_types Boundary types to specify which boundary cells to iterate over. Defaults to all sides.
     *
     * @return Const iterator to the first particle within the given boundary cells.
     */
    [[nodiscard]] const_proximity_iterator boundaryBegin(const std::set<BoundaryType>& boundary_types = {
                                                             BoundaryType::UPPER, BoundaryType::LOWER,
                                                             BoundaryType::FRONT, BoundaryType::BACK,
                                                             BoundaryType::LEFT, BoundaryType::RIGHT}) const;

    /**
     * @brief Iterator over particles in boundary cells.
     *
     * @param boundary_types Boundary types to specify which boundary cells to iterate over. Defaults to all sides.
     *
     * @return Iterator after the last particle within the given boundary cells.
     */
    [[nodiscard]] proximity_iterator boundaryEnd(const std::set<BoundaryType>& boundary_types = {
                                                     BoundaryType::UPPER, BoundaryType::LOWER, BoundaryType::FRONT,
                                                     BoundaryType::BACK, BoundaryType::LEFT, BoundaryType::RIGHT});

    /**
     * @brief Const Iterator over particles in boundary cells.
     *
     * @param boundary_types Boundary types to specify which boundary cells to iterate over. Defaults to all sides.
     *
     * @return Const iterator after the last particle within the given boundary cells.
     */
    [[nodiscard]] const_proximity_iterator boundaryEnd(const std::set<BoundaryType>& boundary_types = {
                                                           BoundaryType::UPPER, BoundaryType::LOWER,
                                                           BoundaryType::FRONT, BoundaryType::BACK, BoundaryType::LEFT,
                                                           BoundaryType::RIGHT}) const;
};

}  // namespace mol_sim

#endif
