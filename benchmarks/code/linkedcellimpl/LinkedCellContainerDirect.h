#ifndef LINKEDCELL_CONTAINER_DIRECT_H
#define LINKEDCELL_CONTAINER_DIRECT_H

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
 * @brief Cell with direct access to Particle data.
 *
 * Implemented with std::vector to be comparable with CellExplicit.
 * Implements all methods of the Cell (so documentation is analogous), but only used for benchmarking.
 */
class CellDirect {
    /**
     * @brief std::vector storing all Particles directly in the cell.
     */
    std::vector<Particle> data;
    std::array<double, 6> bounds;
    [[maybe_unused]] CellType type = CellType::INNER;

   public:
    CellDirect(CellType cell_type, std::array<double, 6> bounds);

    void addParticle(Particle&& value);
    void addParticle(const Particle& value);
    void addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg);
    void addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg, int type);
    Particle removeParticle(size_t idx);
    void clear();
    std::vector<Particle>& particles();
    [[nodiscard]] const std::vector<Particle>& particles() const;
    [[nodiscard]] bool fits(R3 x) const;
    size_t size();
    CellType getType();
};

/**
 * @brief Linked-Cell Container for Particles
 *
 * This container DOES NOT implement the concept ParticleContainer.
 *
 * Implemented with cells that have direct access to Particles for benchmarking purposes.
 * Implements all methods of the LinkedCellContainer (so documentation is analogous), but only used for benchmarking.
 *
 * Expect better performance for iterators (due to memory locality), but worse performance for modifications.
 */
class LinkedCellContainerDirect {
    [[nodiscard]] size_t findCellIndex(R3 vec) const;
    std::vector<CellDirect*> findAdjacentCells(size_t cell_idx);
    [[nodiscard]] std::vector<const CellDirect*> findAdjacentCells(size_t cell_idx) const;
    void findBoundaryCells(BoundaryType type, std::vector<const CellDirect*>& boundary_cells, size_t offset = 0) const;
    void findBoundaryCells(BoundaryType type, std::vector<CellDirect*>& boundary_cells, size_t offset = 0);

    /**
     * @brief Storage and Management as std::vector storing all cells of the container.
     */
    std::vector<CellDirect> cells;
    R3 domain_size;
    std::array<size_t, 3> num_cells;
    std::array<double, 3> cell_length{};

   public:
    // constructors

    LinkedCellContainerDirect(R3 domain_size, double cutoff_radius);

    // retrieve data

    [[nodiscard]] bool fitsDomain(R3 v) const;
    [[nodiscard]] bool fitsContainer(R3 v) const;
    [[nodiscard]] size_t size() const;
    [[nodiscard]] bool empty() const;

    // modify

    void clear();
    void addParticle(Particle&& value);
    void addParticle(const Particle& value);
    void addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg);
    void addParticle(R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg, int type);
    void eraseParticle(Particle* p);
    void updateParticlePosition(Particle& p, R3 new_x);

    // iterators

    /**
     * @brief Iterator that iterates over all particles that apply a force on a given particle.
     *
     * Does not require indirection over storage.
     */
    class proximity_iterator {
        std::vector<Particle>::iterator cur;
        std::vector<Particle>::iterator end;
        std::vector<Particle>::iterator cell_end;
        std::vector<CellDirect*> cells;
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
            while (cur != end && (cur == cell_end || !((center - (*cur).getX()).euclidNorm() <= radius))) {
                inc();
            }
        }

       public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = Particle;
        using difference_type = std::ptrdiff_t;
        using pointer = Particle*;
        using reference = Particle&;

        proximity_iterator() noexcept : radius(0.0) {}
        proximity_iterator(R3 center, double radius, std::vector<Particle>::iterator cur,
                           std::vector<CellDirect*> cells)
            : cur(cur),
              end(cells.back()->particles().end()),
              cell_end(cells.front()->particles().end()),
              cells(cells),
              radius(radius),
              center(center) {
            satisfy();
        }

        reference operator*() const { return *cur; }
        pointer operator->() const { return &(*cur); }

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
     * @brief Const Iterator that iterates over all particles that apply a force on a given particle.
     *
     * Does not require indirection over storage.
     */
    class const_proximity_iterator {
        std::vector<Particle>::const_iterator cur;
        std::vector<Particle>::const_iterator end;
        std::vector<Particle>::const_iterator cell_end;
        std::vector<const CellDirect*> cells;
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
            while (cur != end && (cur == cell_end || !((center - (*cur).getX()).euclidNorm() <= radius))) {
                inc();
            }
        }

       public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = const Particle;
        using difference_type = std::ptrdiff_t;
        using pointer = const Particle*;
        using reference = const Particle&;

        const_proximity_iterator() noexcept : radius(0.0) {}
        const_proximity_iterator(R3 center, double radius, std::vector<Particle>::const_iterator cur,
                                 std::vector<const CellDirect*> cells)
            : cur(cur),
              end(cells.back()->particles().end()),
              cell_end(cells.front()->particles().end()),
              cells(cells),
              radius(radius),
              center(center) {
            satisfy();
        }

        reference operator*() const { return *cur; }
        pointer operator->() const { return &(*cur); }

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

    proximity_iterator begin();
    [[nodiscard]] const_proximity_iterator begin() const;
    [[nodiscard]] const_proximity_iterator cbegin() const;
    proximity_iterator end();
    [[nodiscard]] const_proximity_iterator end() const;
    [[nodiscard]] const_proximity_iterator cend() const;

    [[nodiscard]] proximity_iterator proximityBegin(R3 center, double radius, size_t offset = 0);
    [[nodiscard]] proximity_iterator proximityEnd(R3 center, double radius);
    [[nodiscard]] const_proximity_iterator proximityBegin(R3 center, double radius, size_t offset = 0) const;
    [[nodiscard]] const_proximity_iterator proximityEnd(R3 center, double radius) const;

    // boundary and halo iterators

    [[nodiscard]] proximity_iterator haloBegin(const std::set<BoundaryType>& boundary_types = {
                                                   BoundaryType::UPPER, BoundaryType::LOWER, BoundaryType::FRONT,
                                                   BoundaryType::BACK, BoundaryType::LEFT, BoundaryType::RIGHT});
    [[nodiscard]] const_proximity_iterator haloBegin(const std::set<BoundaryType>& boundary_types = {
                                                         BoundaryType::UPPER, BoundaryType::LOWER, BoundaryType::FRONT,
                                                         BoundaryType::BACK, BoundaryType::LEFT,
                                                         BoundaryType::RIGHT}) const;
    [[nodiscard]] proximity_iterator haloEnd(const std::set<BoundaryType>& boundary_types = {
                                                 BoundaryType::UPPER, BoundaryType::LOWER, BoundaryType::FRONT,
                                                 BoundaryType::BACK, BoundaryType::LEFT, BoundaryType::RIGHT});
    [[nodiscard]] const_proximity_iterator haloEnd(const std::set<BoundaryType>& boundary_types = {
                                                       BoundaryType::UPPER, BoundaryType::LOWER, BoundaryType::FRONT,
                                                       BoundaryType::BACK, BoundaryType::LEFT,
                                                       BoundaryType::RIGHT}) const;

    [[nodiscard]] proximity_iterator boundaryBegin(const std::set<BoundaryType>& boundary_types = {
                                                       BoundaryType::UPPER, BoundaryType::LOWER, BoundaryType::FRONT,
                                                       BoundaryType::BACK, BoundaryType::LEFT, BoundaryType::RIGHT});
    [[nodiscard]] const_proximity_iterator boundaryBegin(const std::set<BoundaryType>& boundary_types = {
                                                             BoundaryType::UPPER, BoundaryType::LOWER,
                                                             BoundaryType::FRONT, BoundaryType::BACK,
                                                             BoundaryType::LEFT, BoundaryType::RIGHT}) const;
    [[nodiscard]] proximity_iterator boundaryEnd(const std::set<BoundaryType>& boundary_types = {
                                                     BoundaryType::UPPER, BoundaryType::LOWER, BoundaryType::FRONT,
                                                     BoundaryType::BACK, BoundaryType::LEFT, BoundaryType::RIGHT});
    [[nodiscard]] const_proximity_iterator boundaryEnd(const std::set<BoundaryType>& boundary_types = {
                                                           BoundaryType::UPPER, BoundaryType::LOWER,
                                                           BoundaryType::FRONT, BoundaryType::BACK, BoundaryType::LEFT,
                                                           BoundaryType::RIGHT}) const;

    [[nodiscard]] bool isOnBoundary(Particle& p);
    [[nodiscard]] R3 getDomainSize();
};

}  // namespace mol_sim

#endif
