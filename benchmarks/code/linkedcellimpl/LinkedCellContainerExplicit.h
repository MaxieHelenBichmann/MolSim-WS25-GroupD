#ifndef LINKEDCELL_CONTAINER_EXPLICIT_H
#define LINKEDCELL_CONTAINER_EXPLICIT_H

#include <spdlog/spdlog.h>

#include <array>
#include <cstdint>
#include <set>
#include <vector>

#include "particles/Particle.h"
#include "particles/container/LinkedCellContainer.h"
#include "particles/container/cells/Cell.h"

namespace mol_sim {

/**
 * @brief Cell as a management structure, keeping references to Particles in the form of indices.
 *
 * Implemented with std::vector to be comparable with CellDirect.
 * Implements all methods of the Cell (so documentation is analogous), but only used for benchmarking.
 */
class CellExplicit {
    /**
     * @brief std::vector storing the indices of Particles in the cell.
     */
    std::vector<size_t> indices;
    std::array<double, 6> bounds;
    [[maybe_unused]] CellType type = CellType::INNER;

   public:
    CellExplicit(CellType cell_type, std::array<double, 6> bounds);

    void addParticle(size_t idx);
    void removeParticle(size_t idx);
    void updateParticleIndex(size_t old_idx, size_t new_idx);
    void clear();
    std::vector<size_t>& particles();
    [[nodiscard]] const std::vector<size_t>& particles() const;
    [[nodiscard]] bool fits(R3 x) const;
    size_t size();
    CellType getType();
};

/**
 * @brief Linked-Cell Container for Particles
 *
 * This container DOES NOT implement the concept ParticleContainer.
 *
 * Implemented with an explicit storage for Particles and Cell management structure for benchmarking purposes.
 * Implements all methods of the LinkedCellContainer (so documentation is analogous), but only used for benchmarking.
 *
 * Expect better performance for modifications, but worse performance for iterators.
 */
class LinkedCellContainerExplicit {
    [[nodiscard]] size_t findCellIndex(R3 vec) const;
    std::vector<CellExplicit*> findAdjacentCells(size_t cell_idx);
    [[nodiscard]] std::vector<const CellExplicit*> findAdjacentCells(size_t cell_idx) const;
    void findBoundaryCells(BoundaryLocation type, std::vector<const CellExplicit*>& boundary_cells,
                           size_t offset = 0) const;
    void findBoundaryCells(BoundaryLocation type, std::vector<CellExplicit*>& boundary_cells, size_t offset = 0);
    void decreaseCellIndices(size_t starting_idx);

    /**
     * @brief Storage vector for actual Particle data.
     */
    std::vector<Particle> data;
    /**
     * @brief Management structure for Cells.
     */
    std::vector<CellExplicit> cells;
    R3 domain_size;
    std::array<size_t, 3> num_cells;
    std::array<double, 3> cell_length{};

   public:
    // constructors

    LinkedCellContainerExplicit(R3 domain_size, double cutoff_radius);

    // retrieve data

    [[nodiscard]] bool fitsDomain(R3 v) const;
    [[nodiscard]] bool fitsContainer(R3 v) const;
    Particle& operator[](size_t idx);
    const Particle& operator[](size_t idx) const;
    [[nodiscard]] size_t size() const;
    [[nodiscard]] bool empty() const;

    // modify

    void clear();
    void reserve(size_t n);
    void addParticle(Particle&& value);
    void addParticle(const Particle& value);
    void addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg);
    void addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg, int type);
    void eraseParticle(Particle* p);
    void updateParticlePosition(std::vector<Particle>::iterator p, R3 new_x);

    // iterators

    /**
     * @brief Iterator that iterates over all particles that apply a force on a given particle.
     *
     * Requires indirection over storage.
     */
    class proximity_iterator {
        std::vector<size_t>::iterator cur;
        std::vector<size_t>::iterator end;
        std::vector<size_t>::iterator cell_end;
        std::vector<CellExplicit*> cells;
        std::vector<Particle>* container_data;
        double radius;
        R3 center;

        void inc() {
            SPDLOG_DEBUG("Incrementing proximity iterator");
            if (cur != cell_end) {
                ++cur;
            }
            while (cur == cell_end && cur != end) {  // reached end of current cell
                cells.erase(cells.begin());
                cur = cells.front()->particles().begin();
                cell_end = cells.front()->particles().end();
            }
        }

        void satisfy() {
            while (cur != end &&
                   (cur == cell_end || !((center - (*container_data)[*cur].getX()).euclidNorm() <= radius))) {
                inc();
            }
        }

       public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = Particle;
        using difference_type = std::ptrdiff_t;
        using pointer = Particle*;
        using reference = Particle&;

        proximity_iterator() noexcept : container_data(nullptr), radius(0.0) {}
        proximity_iterator(R3 center, double radius, std::vector<size_t>::iterator cur,
                           std::vector<CellExplicit*> cells, std::vector<Particle>* data)
            : cur(cur),
              end(cells.back()->particles().end()),
              cell_end(cells.front()->particles().end()),
              cells(cells),
              container_data(data),
              radius(radius),
              center(center) {
            satisfy();
        }

        reference operator*() const { return (*container_data)[*cur]; }
        pointer operator->() const { return &(*container_data)[*cur]; }

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
     * @brief Cost Iterator that iterates over all particles that apply a force on a given particle.
     *
     * Requires indirection over storage.
     */
    class const_proximity_iterator {
        std::vector<size_t>::const_iterator cur;
        std::vector<size_t>::const_iterator end;
        std::vector<size_t>::const_iterator cell_end;
        std::vector<const CellExplicit*> cells;
        const std::vector<Particle>* container_data;
        double radius;
        R3 center;

        void inc() {
            SPDLOG_DEBUG("Incrementing const proximity iterator");
            if (cur != cell_end) {
                ++cur;
            }
            while (cur == cell_end && cur != end) {  // reached end of current cell
                cells.erase(cells.begin());
                cur = cells.front()->particles().begin();
                cell_end = cells.front()->particles().end();
            }
        }

        void satisfy() {
            if (std::isinf(radius)) {
                return;
            }
            while (cur != end &&
                   (cur == cell_end || !((center - (*container_data)[*cur].getX()).euclidNorm() <= radius))) {
                inc();
            }
        }

       public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = const Particle;
        using difference_type = std::ptrdiff_t;
        using pointer = const Particle*;
        using reference = const Particle&;

        const_proximity_iterator() noexcept : container_data(nullptr), radius(0.0) {}
        const_proximity_iterator(R3 center, double radius, std::vector<size_t>::const_iterator cur,
                                 std::vector<const CellExplicit*> cells, const std::vector<Particle>* container_data)
            : cur(cur),
              end(cells.back()->particles().end()),
              cell_end(cells.front()->particles().end()),
              cells(cells),
              container_data(container_data),
              radius(radius),
              center(center) {
            satisfy();
        }

        reference operator*() const { return (*container_data)[*cur]; }
        pointer operator->() const { return &(*container_data)[*cur]; }

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

    std::vector<Particle>::iterator begin();
    [[nodiscard]] std::vector<Particle>::const_iterator begin() const;
    [[nodiscard]] std::vector<Particle>::const_iterator cbegin() const;
    std::vector<Particle>::iterator end();
    [[nodiscard]] std::vector<Particle>::const_iterator end() const;
    [[nodiscard]] std::vector<Particle>::const_iterator cend() const;

    [[nodiscard]] proximity_iterator proximityBegin(R3 center, double radius, size_t offset = 0);
    [[nodiscard]] proximity_iterator proximityEnd(R3 center, double radius);
    [[nodiscard]] const_proximity_iterator proximityBegin(R3 center, double radius, size_t offset = 0) const;
    [[nodiscard]] const_proximity_iterator proximityEnd(R3 center, double radius) const;

    // boundary and halo iterators

    [[nodiscard]] proximity_iterator haloBegin(const std::set<BoundaryLocation>& boundary_types = {
                                                   BoundaryLocation::UPPER, BoundaryLocation::LOWER, BoundaryLocation::FRONT,
                                                   BoundaryLocation::BACK, BoundaryLocation::LEFT, BoundaryLocation::RIGHT});
    [[nodiscard]] const_proximity_iterator haloBegin(const std::set<BoundaryLocation>& boundary_types = {
                                                         BoundaryLocation::UPPER, BoundaryLocation::LOWER, BoundaryLocation::FRONT,
                                                         BoundaryLocation::BACK, BoundaryLocation::LEFT,
                                                         BoundaryLocation::RIGHT}) const;
    [[nodiscard]] proximity_iterator haloEnd(const std::set<BoundaryLocation>& boundary_types = {
                                                 BoundaryLocation::UPPER, BoundaryLocation::LOWER, BoundaryLocation::FRONT,
                                                 BoundaryLocation::BACK, BoundaryLocation::LEFT, BoundaryLocation::RIGHT});
    [[nodiscard]] const_proximity_iterator haloEnd(const std::set<BoundaryLocation>& boundary_types = {
                                                       BoundaryLocation::UPPER, BoundaryLocation::LOWER, BoundaryLocation::FRONT,
                                                       BoundaryLocation::BACK, BoundaryLocation::LEFT,
                                                       BoundaryLocation::RIGHT}) const;

    [[nodiscard]] proximity_iterator boundaryBegin(const std::set<BoundaryLocation>& boundary_types = {
                                                       BoundaryLocation::UPPER, BoundaryLocation::LOWER, BoundaryLocation::FRONT,
                                                       BoundaryLocation::BACK, BoundaryLocation::LEFT, BoundaryLocation::RIGHT});
    [[nodiscard]] const_proximity_iterator boundaryBegin(const std::set<BoundaryLocation>& boundary_types = {
                                                             BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                             BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                             BoundaryLocation::LEFT, BoundaryLocation::RIGHT}) const;
    [[nodiscard]] proximity_iterator boundaryEnd(const std::set<BoundaryLocation>& boundary_types = {
                                                     BoundaryLocation::UPPER, BoundaryLocation::LOWER, BoundaryLocation::FRONT,
                                                     BoundaryLocation::BACK, BoundaryLocation::LEFT, BoundaryLocation::RIGHT});
    [[nodiscard]] const_proximity_iterator boundaryEnd(const std::set<BoundaryLocation>& boundary_types = {
                                                           BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                           BoundaryLocation::FRONT, BoundaryLocation::BACK, BoundaryLocation::LEFT,
                                                           BoundaryLocation::RIGHT}) const;

    [[nodiscard]] bool isOnBoundary(Particle& p);
    [[nodiscard]] R3 getDomainSize();
};
// static_assert(ParticleContainer<LinkedCellContainerExplicit>); // TODO: Update to new concept

}  // namespace mol_sim

#endif
