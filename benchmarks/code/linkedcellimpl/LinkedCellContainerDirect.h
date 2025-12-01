#ifndef LINKEDCELL_CONTAINER_DIRECT_H
#define LINKEDCELL_CONTAINER_DIRECT_H

#include <spdlog/spdlog.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <set>
#include <vector>

#include "particles/Particle.h"
#include "particles/ParticleContainer.h"
#include "particles/container/LinkedCellContainer.h"
#include "particles/container/cells/Cell.h"
#include "utils/Vector.h"

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
 * This container implements the concept ParticleContainer.
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
    void findBoundaryCells(BoundaryLocation type, std::vector<const CellDirect*>& boundary_cells,
                           size_t offset = 0) const;
    void findBoundaryCells(BoundaryLocation type, std::vector<CellDirect*>& boundary_cells, size_t offset = 0);

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

    // iterators

    /**
     * @brief Iterator that iterates over all particles that apply a force on a given particle.
     *
     * Does not require indirection over storage.
     */
    class proximity_iterator {
        std::vector<Particle>::iterator cur;
        std::vector<Particle>::iterator end;
        std::vector<Particle>::iterator begin;
        std::vector<CellDirect*>::iterator cur_cell;
        std::vector<CellDirect*> cells;
        std::vector<Particle> skipped_particles;
        double radius;
        R3 center;

        void inc() {
            SPDLOG_DEBUG("Incrementing proximity iterator");
            if (cur != end) {
                ++cur;
            }
            while (cur_cell != cells.end() && cur == (*cur_cell)->particles().end()) {  // reached end of current cell
                cur_cell++;
                cur = (*cur_cell)->particles().begin();
            }
        }
        void dec() {
            SPDLOG_DEBUG("Decremeting proximity iterator");
            if (cur != begin) {
                --cur;
            }
            while (cur_cell != cells.begin() && cur == (*cur_cell)->particles().begin() &&
                   cur != cells.front()->particles().begin()) {  // reached start of current cell
                cur_cell--;
                cur = (*cur_cell)->particles().end();
            }
        }

        void satisfyInc() {
            while (cur != end &&
                   (cur == (*cur_cell)->particles().end() || !((center - (*cur).getX()).euclidNorm() <= radius) ||
                    std::ranges::find(skipped_particles.begin(), skipped_particles.end(), *cur) !=
                        skipped_particles.end())) {
                inc();
            }
        }
        void satisfyDec() {
            while (cur != begin &&
                   (cur == (*cur_cell)->particles().begin() || !((center - (*cur).getX()).euclidNorm() <= radius) ||
                    std::ranges::find(skipped_particles.begin(), skipped_particles.end(), *cur) !=
                        skipped_particles.end())) {
                dec();
            }
        }

       public:
        using iterator_category = std::random_access_iterator_tag;
        using value_type = Particle;
        using difference_type = std::ptrdiff_t;
        using pointer = Particle*;
        using reference = Particle&;

        proximity_iterator() noexcept : radius(0.0) {}
        proximity_iterator(R3 center, double radius, std::vector<Particle>::iterator cur,
                           std::vector<CellDirect*> cells, std::vector<Particle> skipped_particles = {})
            : cur(cur),
              end(cells.back()->particles().end()),
              begin(cells.front()->particles().begin()),
              cur_cell(cells.begin()),
              cells(cells),
              skipped_particles(std::move(skipped_particles)),
              radius(radius),
              center(center) {
            satisfyInc();
        }

        reference operator*() const { return *cur; }
        pointer operator->() const { return &(*cur); }

        proximity_iterator& operator++() {
            inc();
            satisfyInc();
            return *this;
        }

        proximity_iterator operator++(int) {
            proximity_iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        proximity_iterator& operator--() {
            dec();
            satisfyDec();
            return *this;
        }

        proximity_iterator operator--(int) {
            proximity_iterator tmp = *this;
            --(*this);
            return tmp;
        }

        proximity_iterator& operator+=(difference_type n) {
            if (n >= 0) {
                for (difference_type i = 0; i < n; ++i) {
                    inc();
                    satisfyInc();
                }
            } else {
                for (difference_type i = 0; i < -n; ++i) {
                    dec();
                    satisfyDec();
                }
            }
            return *this;
        }

        proximity_iterator& operator-=(difference_type n) { return *this += -n; }
        proximity_iterator operator+(difference_type n) const {
            proximity_iterator tmp = *this;
            tmp += n;
            return tmp;
        };
        friend proximity_iterator operator+(difference_type n, proximity_iterator it) {
            it += n;
            return it;
        }
        proximity_iterator operator-(difference_type n) const {
            proximity_iterator tmp = *this;
            tmp -= n;
            return tmp;
        };
        difference_type operator-([[maybe_unused]] const proximity_iterator& other) const {
            return 0;
        }  // dummy function
        value_type& operator[](difference_type idx) const {
            if (idx < 0) {
                idx = -idx;
            }
            proximity_iterator tmp = *this;
            tmp += idx;
            return *tmp;

        }  // dummy function
        bool operator<(const proximity_iterator& other) const { return cur_cell < other.cur_cell; }
        bool operator<=(const proximity_iterator& other) const { return cur_cell <= other.cur_cell; }
        bool operator>(const proximity_iterator& other) const { return cur_cell > other.cur_cell; }
        bool operator>=(const proximity_iterator& other) const { return cur_cell >= other.cur_cell; }

        friend bool operator==(const proximity_iterator& a, const proximity_iterator& b) { return a.cur == b.cur; }
        friend bool operator!=(const proximity_iterator& a, const proximity_iterator& b) { return !(a == b); }

        operator std::vector<Particle>::iterator() const { return cur; }
        [[nodiscard]] std::vector<CellDirect*> getCells() const { return cells; }
        [[nodiscard]] double getRadius() const { return radius; }
        [[nodiscard]] R3 getCenter() const { return center; }
        [[nodiscard]] std::vector<Particle> getSkipped() const { return skipped_particles; }
        void skipParticle(const Particle& p) { skipped_particles.push_back(p); }
    };
    static_assert(std::random_access_iterator<proximity_iterator>);

    /**
     * @brief Const Iterator that iterates over all particles that apply a force on a given particle.
     *
     * Does not require indirection over storage.
     */
    class const_proximity_iterator {
        std::vector<Particle>::const_iterator cur;
        std::vector<Particle>::const_iterator end;
        std::vector<Particle>::const_iterator begin;
        std::vector<const CellDirect*>::iterator cur_cell;
        std::vector<const CellDirect*> cells;
        std::vector<Particle> skipped_particles;
        double radius;
        R3 center;

        void inc() {
            SPDLOG_DEBUG("Incrementing proximity iterator");
            if (cur != end) {
                ++cur;
            }
            while (cur_cell != cells.end() && cur == (*cur_cell)->particles().end()) {  // reached end of current cell
                cur_cell++;
                cur = (*cur_cell)->particles().begin();
            }
        }
        void dec() {
            SPDLOG_DEBUG("Decremeting proximity iterator");
            if (cur != begin) {
                --cur;
            }
            while (cur_cell != cells.begin() && cur == (*cur_cell)->particles().begin() &&
                   cur != cells.front()->particles().begin()) {  // reached start of current cell
                cur_cell--;
                cur = (*cur_cell)->particles().end();
            }
        }

        void satisfyInc() {
            while (cur != end &&
                   (cur == (*cur_cell)->particles().end() || !((center - (*cur).getX()).euclidNorm() <= radius) ||
                    std::ranges::find(skipped_particles.begin(), skipped_particles.end(), *cur) !=
                        skipped_particles.end())) {
                inc();
            }
        }
        void satisfyDec() {
            while (cur != begin &&
                   (cur == (*cur_cell)->particles().begin() || !((center - (*cur).getX()).euclidNorm() <= radius) ||
                    std::ranges::find(skipped_particles.begin(), skipped_particles.end(), *cur) !=
                        skipped_particles.end())) {
                dec();
            }
        }

       public:
        using iterator_category = std::random_access_iterator_tag;
        using value_type = const Particle;
        using difference_type = std::ptrdiff_t;
        using pointer = const Particle*;
        using reference = const Particle&;

        const_proximity_iterator() noexcept : radius(0.0) {}
        const_proximity_iterator(R3 center, double radius, std::vector<Particle>::const_iterator cur,
                                 std::vector<const CellDirect*> cells, std::vector<Particle> skipped_particles = {})
            : cur(cur),
              end(cells.back()->particles().end()),
              begin(cells.front()->particles().begin()),
              cur_cell(cells.begin()),
              cells(cells),
              skipped_particles(std::move(skipped_particles)),
              radius(radius),
              center(center) {
            satisfyInc();
        }

        reference operator*() const { return *cur; }
        pointer operator->() const { return &(*cur); }

        const_proximity_iterator& operator++() {
            inc();
            satisfyInc();
            return *this;
        }

        const_proximity_iterator operator++(int) {
            const_proximity_iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        const_proximity_iterator& operator--() {
            dec();
            satisfyDec();
            return *this;
        }

        const_proximity_iterator operator--(int) {
            const_proximity_iterator tmp = *this;
            --(*this);
            return tmp;
        }

        const_proximity_iterator& operator+=(difference_type n) {
            if (n >= 0) {
                for (difference_type i = 0; i < n; ++i) {
                    inc();
                    satisfyInc();
                }
            } else {
                for (difference_type i = 0; i < -n; ++i) {
                    dec();
                    satisfyDec();
                }
            }
            return *this;
        }

        const_proximity_iterator& operator-=(difference_type n) { return *this += -n; }
        const_proximity_iterator operator+(difference_type n) const {
            const_proximity_iterator tmp = *this;
            tmp += n;
            return tmp;
        };
        friend const_proximity_iterator operator+(difference_type n, const_proximity_iterator it) {
            it += n;
            return it;
        }
        const_proximity_iterator operator-(difference_type n) const {
            const_proximity_iterator tmp = *this;
            tmp -= n;
            return tmp;
        };
        difference_type operator-([[maybe_unused]] const const_proximity_iterator& other) const {
            return 0;
        }  // dummy function
        value_type& operator[](difference_type idx) const {
            if (idx < 0) {
                idx = -idx;
            }
            const_proximity_iterator tmp = *this;
            tmp += idx;
            return *tmp;

        }  // dummy function
        bool operator<(const const_proximity_iterator& other) const { return cur_cell < other.cur_cell; }
        bool operator<=(const const_proximity_iterator& other) const { return cur_cell <= other.cur_cell; }
        bool operator>(const const_proximity_iterator& other) const { return cur_cell > other.cur_cell; }
        bool operator>=(const const_proximity_iterator& other) const { return cur_cell >= other.cur_cell; }

        friend bool operator==(const const_proximity_iterator& a, const const_proximity_iterator& b) {
            return a.cur == b.cur;
        }
        friend bool operator!=(const const_proximity_iterator& a, const const_proximity_iterator& b) {
            return !(a == b);
        }
        operator std::vector<Particle>::const_iterator() const { return cur; }
    };
    static_assert(std::random_access_iterator<const_proximity_iterator>);

    proximity_iterator eraseParticle(proximity_iterator p);
    proximity_iterator updateParticlePosition(proximity_iterator p, R3 new_x);

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

    [[nodiscard]] proximity_iterator haloBegin(const std::set<BoundaryLocation>& boundary_types = {
                                                   BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                   BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                   BoundaryLocation::LEFT, BoundaryLocation::RIGHT});
    [[nodiscard]] const_proximity_iterator haloBegin(const std::set<BoundaryLocation>& boundary_types = {
                                                         BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                         BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                         BoundaryLocation::LEFT, BoundaryLocation::RIGHT}) const;
    [[nodiscard]] proximity_iterator haloEnd(const std::set<BoundaryLocation>& boundary_types = {
                                                 BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                 BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                 BoundaryLocation::LEFT, BoundaryLocation::RIGHT});
    [[nodiscard]] const_proximity_iterator haloEnd(const std::set<BoundaryLocation>& boundary_types = {
                                                       BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                       BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                       BoundaryLocation::LEFT, BoundaryLocation::RIGHT}) const;

    [[nodiscard]] proximity_iterator boundaryBegin(const std::set<BoundaryLocation>& boundary_types = {
                                                       BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                       BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                       BoundaryLocation::LEFT, BoundaryLocation::RIGHT});
    [[nodiscard]] const_proximity_iterator boundaryBegin(const std::set<BoundaryLocation>& boundary_types = {
                                                             BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                             BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                             BoundaryLocation::LEFT, BoundaryLocation::RIGHT}) const;
    [[nodiscard]] proximity_iterator boundaryEnd(const std::set<BoundaryLocation>& boundary_types = {
                                                     BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                     BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                     BoundaryLocation::LEFT, BoundaryLocation::RIGHT});
    [[nodiscard]] const_proximity_iterator boundaryEnd(const std::set<BoundaryLocation>& boundary_types = {
                                                           BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                           BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                           BoundaryLocation::LEFT, BoundaryLocation::RIGHT}) const;

    [[nodiscard]] bool isOnBoundary(Particle& p);
    [[nodiscard]] R3 getDomainSize();
};
static_assert(ParticleContainer<LinkedCellContainerDirect>);

}  // namespace mol_sim

#endif
