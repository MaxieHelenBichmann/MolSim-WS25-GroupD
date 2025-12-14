#ifndef LINKEDCELL_CONTAINER_H
#define LINKEDCELL_CONTAINER_H

#include <particles/Particle.h>
#include <spdlog/spdlog.h>

#include <array>
#include <cstddef>
#include <set>
#include <span>
#include <vector>

#include "particles/ParticleContainer.h"
#include "particles/boundaries/Boundary.h"
#include "particles/container/cells/Cell.h"

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
    std::vector<Cell*> findAdjacentCellsN3L(size_t cell_idx);

    /**
     * @brief Find the adjacent cells of a cell. [HELPER FUNCTION]
     *
     * @param cell_idx Index of the cell.
     *
     * @return Vector of const pointers to the adjacent cells, including the cell itself.
     */
    [[nodiscard]] std::vector<const Cell*> findAdjacentCellsN3L(size_t cell_idx) const;

    /**
     * @brief Find the boundary or halo cells. [HELPER FUNCTION]
     *
     * @param type Type of the boundary.
     * @param boundary_cells Vector to store the found boundary cells.
     * @param offset Offset to consider for halo cells. (0 = halo cells, 1 = boundary cells)
     */
    void findBoundaryCells(BoundaryLocation type, std::vector<const Cell*>& boundary_cells, size_t offset = 0) const;

    /**
     * @brief Find the boundary or halo cells. [HELPER FUNCTION]
     *
     * @param type Type of the boundary.
     * @param boundary_cells Vector to store the found boundary cells.
     * @param offset Offset to consider for halo cells. (0 = halo cells, 1 = boundary cells)
     */
    void findBoundaryCells(BoundaryLocation type, std::vector<Cell*>& boundary_cells, size_t offset = 0);

    /**
     * @brief Update the indices of Particles in cells after erasing a Particle. [HELPER FUNCTION]
     *
     * @param starting_idx Index from which to start updating.
     */
    void decreaseCellIndices(size_t starting_idx);

    /**
     * std::vector storing all Particles in the container.
     */
    std::vector<Particle> data;

    /**
     * std::vector storing all cells of the container.
     */
    std::vector<Cell> cells;

    /**
     * Cutoff radius with which the Container is initialized.
     * Used for proximity queries. Particles further apart than this radius are not considered.
     */
    double cutoff_radius;

    /**
     * Vector storing the global bounds of the domain of this container.
     * The Linked-Cell container assumes a cuboidal domain from (0,0,0) to domain_size.
     *
     * NOTE: All components specified here must be positive!
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
     * @brief Parses the number of cells (for LinkedCellContainer) and the cell_length / corner_length for
     * LinkedCellContainer and Periodic respectively based on the provided domain_size and cutoff_radius
     *
     * @param num_cells Where the computed number of cells in each dimension will be stored
     * @param cell_length Where the computed length of the cells / corners in each dimension will be stored
     * @param domain_size The size of the domain
     * @param cutoff_radius The cutoff radius
     * 
     * @note (remove this in final product) 
     * We need this code inside Periodic.cpp aswell because of the edge cases that arise when the domain_size
     * isn't divisible (in at least 1 dimension) by the cutoff radius for LCC.
     * Factoring out this code block to here avoids code dupcliation in LCC.cpp and Periodic.cpp.
     * Other ideas for this refactoring would've been:
     *  1) Putting this code into YAMLReader.cpp (and removing it from LCC.cpp) and then adjusting LCC.cpp constructor
     *  2) Just copy-pasting this code into Periodic.h (w/ minor changes)
     *  3) Making a new header file just for this (and other things?) seems overkill
     *  4) Make haloDimension a field in SimpleContainer aswell and then use a getter (plus pass a ContainerRef to 
     *     Periodic) -> entanglement
     *  5) Additional check in MolSim.cpp or Simulation.h (or YAMLReader.cpp) *after* parsing domain type about
     *     whether it is a LCC. If it is then do something like x_boundary.setCornerDimension(cell_length) 
     *     -> more random checks
     * This is the best refactoring I could come up with (except maybe option 1 but I didn't want to 
     * mess with LCC like that for now). If you can think of a nicer one please feel free to change it to that.
     */
    static void computeCellsOrCorners(std::array<size_t, 3>& num_cells, std::array<double, 3>& cell_length, 
        R3 domain_size, double cutoff_radius);

    /**
     * @brief Check whether a (potential) Particle fits into the domain.
     *
     * @param v Coordinates of the Particle to check.
     *
     * @return True if position is within domain boundaries.
     */
    [[nodiscard]] bool fitsDomain(R3 v) const noexcept;

    /**
     * @brief Check whether a (potential) Particle fits into the container (so including halo cells).
     *
     * @param v Coordinates of the Particle to check.
     *
     * @return True if position is within container boundaries.
     */
    [[nodiscard]] bool fitsContainer(R3 v) const noexcept;

    Particle& operator[](size_t idx) noexcept;
    const Particle& operator[](size_t idx) const noexcept;

    /**
     * @brief Returns size of container.
     *
     * @return Number of Particles in the container.
     */
    [[nodiscard]] size_t size() const noexcept;

    /**
     * @brief Queries whether the container is empty.
     *
     * @return True of container contains no Particles.
     */
    [[nodiscard]] bool empty() const noexcept;

    // modify

    /**
     * @brief Clear the entire container, destructing the Particles.
     */
    void clear() noexcept;

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
     * @brief Directly constructing a Particle with its required parameters, including type, in-place.
     * Only used for checkpointing.
     *
     * @param x_arg Initial coordinates of the Particle.
     * @param old_x_arg Old coordinates of the Particle (only used when restoring from checkpoint).
     * @param v_arg Initial velocities of the Particle.
     * @param f_arg Forces on the Particle (only used when restoring from checkpoint).
     * @param old_f_arg Old forces on the Particle (only used when restoring from checkpoint).
     * @param m_arg Mass of the Particle.
     * @param epsilon_arg Epsilon of the Particle.
     * @param sigma_arg Sigma of the Particle.
     * @param type Type of the Particle.
     */
    void addParticle(R3 x_arg, R3 old_x_arg, R3 v_arg, R3 f_arg, R3 old_f_arg, double m_arg, double epsilon_arg,
                     double sigma_arg, int type);

    /**
     * @brief Remove a Particle.
     *
     * @param p Iterator to to the Particle to remove.
     *
     * @return Iterator to the next Particle after the removed one.
     */
    std::vector<Particle>::iterator eraseParticle(std::vector<Particle>::iterator p);

    /**
     * @brief Update the position of a given Particle in the Container.
     *
     * @param p Iterator to the Particle to update.
     * @param new_x New position to set.
     *
     * @return Iterator to the next Particle.
     */
    std::vector<Particle>::iterator updateParticlePosition(std::vector<Particle>::iterator p, R3 new_x);

    // iterators

    /**
     * @brief Mutable iterator of the raw data (no logic).
     *
     * @return Mutable iterator to the first element of the container.
     */
    std::vector<Particle>::iterator begin() noexcept;

    /**
     * @brief Const iterator of the raw data (no logic).
     *
     * @return Const iterator to the first element of the container.
     */
    [[nodiscard]] std::vector<Particle>::const_iterator begin() const noexcept;

    /**
     * @brief Const iterator of the raw data (no logic).
     *
     * @return Const iterator to the first element of the container.
     */
    [[nodiscard]] std::vector<Particle>::const_iterator cbegin() const noexcept;

    /**
     * @brief Mutable iterator of the raw data (no logic).
     *
     * @return Mutable past-the-end iterator of the container.
     */
    std::vector<Particle>::iterator end() noexcept;

    /**
     * @brief Const iterator of the raw data (no logic).
     *
     * @return Const past-the-end iterator of the container.
     */
    [[nodiscard]] std::vector<Particle>::const_iterator end() const noexcept;

    /**
     * @brief Const iterator of the raw data (no logic).
     *
     * @return Const past-the-end iterator of the container
     */
    [[nodiscard]] std::vector<Particle>::const_iterator cend() const noexcept;

    /**
     * @brief Iterator that iterates over all particles that apply a force on a given particle.
     *
     * Satisfies the forward iterator concept.
     * Can be instantiated as mutable and const iterator with P = (const) Particle and C = (const) Cell.
     * @tparam P Either Particle or const Particle
     * @tparam C Either Cell or const Cell
     *
     * Only iterates over particles within the given radius of the center, or all if the radius is infinite.
     * Or is used for the iteration over halo and boundary particles.
     *
     * Enables cutoff radius, and iterates only over the particles in the adjacent cells.
     */
    template <typename P, typename C>
        requires((std::is_same_v<P, Particle> && std::is_same_v<C, Cell>) ||
                 (std::is_same_v<P, const Particle> && std::is_same_v<C, const Cell>))
    class proximity_iterator {
        std::set<size_t>::iterator cur;
        std::set<size_t>::iterator end;
        std::set<size_t>::iterator cell_end;
        std::vector<C*> cells;
        std::span<P> container_data;
        double radius;
        R3 center;
        size_t center_idx;

        void inc() {
            SPDLOG_DEBUG("Incrementing proximity iterator");
            if (cur != cell_end) {
                ++cur;
            }
            while (cur == cell_end && cells.size() > 1) {  // reached end of current cell
                cells.erase(cells.begin());
                cur = cells.front()->particles().begin();
                cell_end = cells.front()->particles().end();
            }
        }

        void satisfy() {
            while (cur != end && (cur == cell_end || !((center - container_data[*cur].getX()).euclidNorm() <= radius) ||
                                  (cells.size() == 1 && *cur <= center_idx && center_idx != container_data.size()))) {
                inc();
            }
        }

       public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = P;
        using difference_type = std::ptrdiff_t;
        using pointer = P*;
        using reference = P&;

        proximity_iterator() noexcept : radius(0.0) {}
        proximity_iterator(R3 center, double radius, std::set<size_t>::iterator cur, std::vector<C*> cells,
                           std::span<P> data, size_t center_idx)
            : cur(cur),
              end(cells.back()->particles().end()),
              cell_end(cells.front()->particles().end()),
              cells(cells),
              container_data(data),
              radius(radius),
              center(center),
              center_idx(center_idx) {
            satisfy();
        }

        reference operator*() const noexcept { return container_data[*cur]; }
        pointer operator->() const noexcept { return &container_data[*cur]; }

        proximity_iterator<P, C>& operator++() {
            inc();
            satisfy();
            return *this;
        }

        proximity_iterator<P, C> operator++(int) {
            proximity_iterator<P, C> tmp = *this;
            ++(*this);
            return tmp;
        }

        friend bool operator==(const proximity_iterator<P, C>& a, const proximity_iterator<P, C>& b) noexcept {
            return a.cur == b.cur;
        }
        friend bool operator!=(const proximity_iterator<P, C>& a, const proximity_iterator<P, C>& b) noexcept {
            return !(a == b);
        }

        [[nodiscard]] std::vector<C*> getCells() const { return cells; }
        [[nodiscard]] double getRadius() const noexcept { return radius; }
        [[nodiscard]] R3 getCenter() const { return center; }
        [[nodiscard]] size_t getCenterIdx() const noexcept { return center_idx; }
        [[nodiscard]] size_t getIdx() const noexcept { return *cur; }
    };
    static_assert(std::forward_iterator<proximity_iterator<Particle, Cell>>);
    static_assert(std::forward_iterator<proximity_iterator<const Particle, const Cell>>);

    /**
     * @brief Remove a Particle. Used when iterating with proximity_iterator, not used yet (thus not tested), but could
     * be useful.
     *
     * @param p Iterator to to the Particle to remove.
     *
     * @return Iterator to the next Particle after the removed one.
     */
    proximity_iterator<Particle, Cell> eraseParticle(proximity_iterator<Particle, Cell> p);

    /**
     * @brief Mutable Iterator over particles in proximity.
     *
     * @param center Center point to check proximity from (position of the particle).
     * @param offset Offset from the beginning of the container to start the iteration (used for N3L optimization).
     *
     * @return Mutable iterator to the first particle within the given radius of the center.
     */
    [[nodiscard]] proximity_iterator<Particle, Cell> proximityBegin(R3 center, size_t offset);

    /**
     * @brief Mutable Iterator over particles in proximity.
     *
     * @param center Center point to check proximity from (position of the particle).
     *
     * @return Mutable iterator after the last particle within the given radius of the center.
     */
    [[nodiscard]] proximity_iterator<Particle, Cell> proximityEnd(R3 center);

    /**
     * @brief Const Iterator over particles in proximity.
     *
     * @param center Center point to check proximity from (position of the particle).
     * @param offset Offset from the beginning of the container to start the iteration (used for N3L optimization).
     *
     * @return Const iterator to the first particle within the given radius of the center.
     */
    [[nodiscard]] proximity_iterator<const Particle, const Cell> proximityBegin(R3 center, size_t offset) const;

    /**
     * @brief Const Iterator over particles in proximity.
     *
     * @param center Center point to check proximity from (position of the particle).
     *
     * @return Const iterator after the last particle within the given radius of the center.
     */
    [[nodiscard]] proximity_iterator<const Particle, const Cell> proximityEnd(R3 center) const;

    // boundary and halo iterators

    /**
     * @brief Iterator over particles in halo cells.
     *
     * @param boundary_types Boundary types to specify which halo cells to iterate over. Defaults to all sides.
     *
     * @return Iterator to the first particle within the given halo cells.
     */
    [[nodiscard]] proximity_iterator<Particle, Cell> haloBegin(const std::set<BoundaryLocation>& boundary_types = {
                                                                   BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                                   BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                                   BoundaryLocation::LEFT, BoundaryLocation::RIGHT});

    /**
     * @brief Const Iterator over particles in halo cells.
     *
     * @param boundary_types Boundary types to specify which halo cells to iterate over. Defaults to all sides.
     *
     * @return Const iterator to the first particle within the given halo cells.
     */
    [[nodiscard]] proximity_iterator<const Particle, const Cell> haloBegin(
        const std::set<BoundaryLocation>& boundary_types = {BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                            BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                            BoundaryLocation::LEFT, BoundaryLocation::RIGHT}) const;

    /**
     * @brief Iterator over particles in halo cells.
     *
     * @param boundary_types Boundary types to specify which halo cells to iterate over. Defaults to all sides.
     *
     * @return Iterator after the last particle within the given halo cells.
     */
    [[nodiscard]] proximity_iterator<Particle, Cell> haloEnd(const std::set<BoundaryLocation>& boundary_types = {
                                                                 BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                                 BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                                 BoundaryLocation::LEFT, BoundaryLocation::RIGHT});

    /**
     * @brief Const Iterator over particles in halo cells.
     *
     * @param boundary_types Boundary types to specify which halo cells to iterate over. Defaults to all sides.
     *
     * @return Const iterator after the last particle within the given halo cells.
     */
    [[nodiscard]] proximity_iterator<const Particle, const Cell> haloEnd(
        const std::set<BoundaryLocation>& boundary_types = {BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                            BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                            BoundaryLocation::LEFT, BoundaryLocation::RIGHT}) const;

    /**
     * @brief Iterator over particles in boundary cells.
     *
     * @param boundary_types Boundary types to specify which boundary cells to iterate over. Defaults to all sides.
     *
     * @return Iterator to the first particle within the given boundary cells.
     */
    [[nodiscard]] proximity_iterator<Particle, Cell> boundaryBegin(const std::set<BoundaryLocation>& boundary_types = {
                                                                       BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                                       BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                                       BoundaryLocation::LEFT,
                                                                       BoundaryLocation::RIGHT});

    /**
     * @brief Const Iterator over particles in boundary cells.
     *
     * @param boundary_types Boundary types to specify which boundary cells to iterate over. Defaults to all sides.
     *
     * @return Const iterator to the first particle within the given boundary cells.
     */
    [[nodiscard]] proximity_iterator<const Particle, const Cell> boundaryBegin(
        const std::set<BoundaryLocation>& boundary_types = {BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                            BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                            BoundaryLocation::LEFT, BoundaryLocation::RIGHT}) const;

    /**
     * @brief Iterator over particles in boundary cells.
     *
     * @param boundary_types Boundary types to specify which boundary cells to iterate over. Defaults to all sides.
     *
     * @return Iterator after the last particle within the given boundary cells.
     */
    [[nodiscard]] proximity_iterator<Particle, Cell> boundaryEnd(const std::set<BoundaryLocation>& boundary_types = {
                                                                     BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                                     BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                                     BoundaryLocation::LEFT, BoundaryLocation::RIGHT});

    /**
     * @brief Const Iterator over particles in boundary cells.
     *
     * @param boundary_types Boundary types to specify which boundary cells to iterate over. Defaults to all sides.
     *
     * @return Const iterator after the last particle within the given boundary cells.
     */
    [[nodiscard]] proximity_iterator<const Particle, const Cell> boundaryEnd(
        const std::set<BoundaryLocation>& boundary_types = {BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                            BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                            BoundaryLocation::LEFT, BoundaryLocation::RIGHT}) const;
};
static_assert(ParticleContainer<LinkedCellContainer>);

}  // namespace mol_sim

#endif