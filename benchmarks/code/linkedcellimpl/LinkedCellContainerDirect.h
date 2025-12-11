#ifndef LINKEDCELL_CONTAINER_DIRECT_H
#define LINKEDCELL_CONTAINER_DIRECT_H

#include <spdlog/spdlog.h>

#include <array>
#include <cstddef>
#include <iterator>
#include <set>
#include <vector>

#include "particles/Particle.h"
#include "particles/ParticleContainer.h"
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
    std::array<double, 6>& getBounds() { return bounds; }
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
    double cutoff_radius;

   public:
    // constructors
    std::vector<CellDirect>& getCells() { return cells; };

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
     * Does not require indirection over storage, but is also used for iteration over all particles (without logic), so
     * has to implement the random-access iterator concept.
     */
    template <typename P, typename C>
        requires((std::is_same_v<P, Particle> && std::is_same_v<C, CellDirect>) ||
                 (std::is_same_v<P, const Particle> && std::is_same_v<C, const CellDirect>))
    class proximity_iterator {
        std::vector<C*> cells;
        size_t cur_cell;
        size_t cur;
        std::vector<Particle> skipped_particles;
        double radius;
        R3 center;

        void inc() {
            if (cells.empty() || cur_cell == cells.size() ||
                (cur_cell == cells.size() - 1 && cur == cells[cur_cell]->particles().size())) {
                SPDLOG_DEBUG("cur = {}, cur_cell = {}, cells size = {}", cur, cur_cell, cells.size());
                return;
            }
            SPDLOG_DEBUG("Incrementing proximity iterator");
            if (cur_cell != cells.size() && cells[cur_cell] != nullptr && cur != cells[cur_cell]->particles().size()) {
                ++cur;
                SPDLOG_DEBUG("cur = {}, cur_cell = {}, cells size = {}", cur, cur_cell, cells.size());
                return;
            }
            while (cur_cell != cells.size() && cells[cur_cell] != nullptr &&
                   cur == cells[cur_cell]->particles().size()) {
                ++cur_cell;
                if (cur_cell == cells.size()) {
                    cur = cells[cur_cell - 1]->particles().size();
                    break;
                }
                while (cur_cell != cells.size() && cells[cur_cell] == nullptr) {
                    ++cur_cell;
                }
                if (cur_cell != cells.size()) {
                    cur = 0;
                } else {
                    cur = cells[cur_cell - 1]->particles().size();
                }
            }
            SPDLOG_DEBUG("cur = {}, cur_cell = {}, cells size = {}", cur, cur_cell, cells.size());
        }
        void dec() {
            if (cells.empty() || (cur_cell == 0 && cur == 0)) {
                return;
            }
            SPDLOG_DEBUG("Decremeting proximity iterator");

            // Handle case when at end position (cur_cell points past last element)
            if (cur_cell == cells.size()) {
                --cur_cell;
                cur = cells[cur_cell]->particles().size();
            }

            if (cur != 0) {
                --cur;
                return;
            }
            while (cur_cell != 0 && cur == 0) {
                --cur_cell;
                cur = cells[cur_cell]->particles().size();
            }
            // After moving to end of previous cell, decrement to get last element
            if (cur != 0) {
                --cur;
            }
        }

        void satisfyInc() {
            while (!cells.empty() && !(cur_cell == cells.size() ||
                                       (cur_cell == cells.size() - 1 && cur == cells[cur_cell]->particles().size()))) {
                // Check if current cell is null or we're at end of current cell's particles
                if (cells[cur_cell] == nullptr || cur == cells[cur_cell]->particles().size()) {
                    inc();
                    continue;
                }
                if (!((center - cells[cur_cell]->particles()[cur].getX()).euclidNorm() <= radius) ||
                    std::ranges::find(skipped_particles.begin(), skipped_particles.end(),
                                      cells[cur_cell]->particles()[cur]) != skipped_particles.end()) {
                    inc();
                    continue;
                }
                break;
            }
        }
        void satisfyDec() {
            // Handle being at the end position first
            while (!cells.empty() && cur_cell == cells.size()) {
                dec();
            }
            while (!cells.empty() && cur_cell != 0) {
                if (cells[cur_cell] == nullptr || cur == 0) {
                    dec();
                    continue;
                }
                if (!((center - cells[cur_cell]->particles()[cur].getX()).euclidNorm() <= radius) ||
                    std::ranges::find(skipped_particles.begin(), skipped_particles.end(),
                                      cells[cur_cell]->particles()[cur]) != skipped_particles.end()) {
                    dec();
                    continue;
                }
                break;
            }
        }

       public:
        using iterator_category = std::random_access_iterator_tag;
        using value_type = P;
        using difference_type = std::ptrdiff_t;
        using pointer = P*;
        using reference = P&;

        proximity_iterator() noexcept : radius(0.0) {}
        proximity_iterator(R3 center, double radius, size_t cur, std::vector<C*> cells_input, size_t cur_cell_input,
                           std::vector<Particle> skipped_particles = {})
            : cells(std::move(cells_input)),
              cur_cell(cur_cell_input),
              cur(cur),
              skipped_particles(std::move(skipped_particles)),
              radius(radius),
              center(center) {
            satisfyInc();
        }

        reference operator*() const { return cells[cur_cell]->particles()[cur]; }
        pointer operator->() const { return &(cells[cur_cell]->particles()[cur]); }

        proximity_iterator<P, C>& operator++() {
            inc();
            satisfyInc();
            return *this;
        }

        proximity_iterator<P, C> operator++(int) {
            proximity_iterator<P, C> tmp = *this;
            ++(*this);
            return tmp;
        }

        proximity_iterator<P, C>& operator--() {
            dec();
            satisfyDec();
            return *this;
        }

        proximity_iterator<P, C> operator--(int) {
            proximity_iterator<P, C> tmp = *this;
            --(*this);
            return tmp;
        }

        proximity_iterator<P, C>& operator+=(difference_type n) {
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

        proximity_iterator<P, C>& operator-=(difference_type n) { return *this += -n; }
        proximity_iterator<P, C> operator+(difference_type n) const {
            proximity_iterator<P, C> tmp = *this;
            tmp += n;
            return tmp;
        };
        friend proximity_iterator<P, C> operator+(difference_type n, proximity_iterator<P, C> it) {
            it += n;
            return it;
        }
        proximity_iterator<P, C> operator-(difference_type n) const {
            proximity_iterator<P, C> tmp = *this;
            tmp -= n;
            return tmp;
        };
        difference_type operator-([[maybe_unused]] const proximity_iterator<P, C>& other) const {
            return 0;
        }  // dummy function
        value_type& operator[](difference_type idx) const {
            if (idx < 0) {
                idx = -idx;
            }
            proximity_iterator<P, C> tmp = *this;
            tmp += idx;
            return *tmp;

        }  // dummy function
        bool operator<(const proximity_iterator<P, C>& other) const { return cur_cell < other.cur_cell; }
        bool operator<=(const proximity_iterator<P, C>& other) const { return cur_cell <= other.cur_cell; }
        bool operator>(const proximity_iterator<P, C>& other) const { return cur_cell > other.cur_cell; }
        bool operator>=(const proximity_iterator<P, C>& other) const { return cur_cell >= other.cur_cell; }

        friend bool operator==(const proximity_iterator<P, C>& a, const proximity_iterator<P, C>& b) {
            SPDLOG_DEBUG("Comparing proximity iterators");
            SPDLOG_DEBUG("cur: {} vs {} --- cur_cell: {} vs {}", a.cur, b.cur, a.cur_cell, b.cur_cell);
            return (a.cur == b.cur && a.cur_cell == b.cur_cell) ||
                   (a.cur_cell == a.cells.size() && b.cur_cell == b.cells.size() - 1 &&
                    b.cur == b.cells[b.cur_cell]->particles().size()) ||
                   (a.cur_cell == a.cells.size() - 1 && a.cur == a.cells[a.cur_cell]->particles().size() &&
                    b.cur_cell == b.cells.size());
        }
        friend bool operator!=(const proximity_iterator<P, C>& a, const proximity_iterator<P, C>& b) {
            return !(a == b);
        }

        operator std::vector<Particle>::iterator() const { return cells[cur_cell]->particles().begin() + cur; }
        operator std::vector<Particle>::const_iterator() const { return cells[cur_cell]->particles().begin() + cur; }
        [[nodiscard]] std::vector<CellDirect*> getCells() const { return cells; }
        [[nodiscard]] size_t getCurCell() const { return cur_cell; }
        [[nodiscard]] size_t getCur() const { return cur; }
        [[nodiscard]] double getRadius() const { return radius; }
        [[nodiscard]] R3 getCenter() const { return center; }
        [[nodiscard]] std::vector<Particle> getSkipped() const { return skipped_particles; }
        void skipParticle(const Particle& p) {
            skipped_particles.push_back(p);
            satisfyInc();
        }
        [[nodiscard]] value_type getParticle() const { return cells[cur_cell]->particles()[cur]; }
    };
    static_assert(std::random_access_iterator<proximity_iterator<Particle, CellDirect>>);
    static_assert(std::random_access_iterator<proximity_iterator<const Particle, const CellDirect>>);

    proximity_iterator<Particle, CellDirect> eraseParticle(proximity_iterator<Particle, CellDirect> p);
    proximity_iterator<Particle, CellDirect> updateParticlePosition(proximity_iterator<Particle, CellDirect> p,
                                                                    R3 new_x);

    proximity_iterator<Particle, CellDirect> begin();
    [[nodiscard]] proximity_iterator<const Particle, const CellDirect> begin() const;
    [[nodiscard]] proximity_iterator<const Particle, const CellDirect> cbegin() const;
    proximity_iterator<Particle, CellDirect> end();
    [[nodiscard]] proximity_iterator<const Particle, const CellDirect> end() const;
    [[nodiscard]] proximity_iterator<const Particle, const CellDirect> cend() const;

    [[nodiscard]] proximity_iterator<Particle, CellDirect> proximityBegin(R3 center, size_t offset = 0);
    [[nodiscard]] proximity_iterator<Particle, CellDirect> proximityEnd(R3 center);
    [[nodiscard]] proximity_iterator<const Particle, const CellDirect> proximityBegin(R3 center,
                                                                                      size_t offset = 0) const;
    [[nodiscard]] proximity_iterator<const Particle, const CellDirect> proximityEnd(R3 center) const;

    // boundary and halo iterators

    [[nodiscard]] proximity_iterator<Particle, CellDirect> haloBegin(
        const std::set<BoundaryLocation>& boundary_types = {BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                            BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                            BoundaryLocation::LEFT, BoundaryLocation::RIGHT});
    [[nodiscard]] proximity_iterator<const Particle, const CellDirect> haloBegin(
        const std::set<BoundaryLocation>& boundary_types = {BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                            BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                            BoundaryLocation::LEFT, BoundaryLocation::RIGHT}) const;
    [[nodiscard]] proximity_iterator<Particle, CellDirect> haloEnd(const std::set<BoundaryLocation>& boundary_types = {
                                                                       BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                                       BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                                       BoundaryLocation::LEFT,
                                                                       BoundaryLocation::RIGHT});
    [[nodiscard]] proximity_iterator<const Particle, const CellDirect> haloEnd(
        const std::set<BoundaryLocation>& boundary_types = {BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                            BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                            BoundaryLocation::LEFT, BoundaryLocation::RIGHT}) const;

    [[nodiscard]] proximity_iterator<Particle, CellDirect> boundaryBegin(
        const std::set<BoundaryLocation>& boundary_types = {BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                            BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                            BoundaryLocation::LEFT, BoundaryLocation::RIGHT});
    [[nodiscard]] proximity_iterator<const Particle, const CellDirect> boundaryBegin(
        const std::set<BoundaryLocation>& boundary_types = {BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                            BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                            BoundaryLocation::LEFT, BoundaryLocation::RIGHT}) const;
    [[nodiscard]] proximity_iterator<Particle, CellDirect> boundaryEnd(
        const std::set<BoundaryLocation>& boundary_types = {BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                            BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                            BoundaryLocation::LEFT, BoundaryLocation::RIGHT});
    [[nodiscard]] proximity_iterator<const Particle, const CellDirect> boundaryEnd(
        const std::set<BoundaryLocation>& boundary_types = {BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                            BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                            BoundaryLocation::LEFT, BoundaryLocation::RIGHT}) const;

    [[nodiscard]] R3 getDomainSize();
};
static_assert(ParticleContainer<LinkedCellContainerDirect>);

}  // namespace mol_sim

#endif
