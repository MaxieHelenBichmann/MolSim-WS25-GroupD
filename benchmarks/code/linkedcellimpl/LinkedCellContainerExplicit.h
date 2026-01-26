#ifndef LINKEDCELL_CONTAINER_EXPLICIT_H
#define LINKEDCELL_CONTAINER_EXPLICIT_H

#include <spdlog/spdlog.h>

#include <array>
#include <set>
#include <span>
#include <vector>

#include "particles/Particle.h"
#include "particles/ParticleContainer.h"
#include "particles/container/cells/Cell.h"

namespace mol_sim {

/**
 * @brief Cell as a management structure, keeping references to Particles in the form of indices.
 *
 * Implemented with std::vector to be comparable with CellDirect.
 * Implements all methods of the Cell (so documentation is analogous), but only used for benchmarking.
 *
 * @note Certain optimizations and additions (dated after Worksheet 3) that were being made to the real
 * LinkedCellContainer NOT included.
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
 * This container implements the concept ParticleContainer.
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
    double cutoff_radius;

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
    void addParticle(R3 x_arg, R3 old_x_arg, R3 v_arg, R3 f_arg, R3 old_f_arg, double m_arg, double epsilon_arg,
                     double sigma_arg, int type);
    std::vector<Particle>::iterator eraseParticle(std::vector<Particle>::iterator p);
    std::vector<Particle>::iterator updateParticlePosition(std::vector<Particle>::iterator p, R3 new_x);

    // iterators

    /**
     * @brief Iterator that iterates over all particles that apply a force on a given particle.
     *
     * Requires indirection over storage.
     */
    template <typename P, typename C>
        requires((std::is_same_v<P, Particle> && std::is_same_v<C, CellExplicit>) ||
                 (std::is_same_v<P, const Particle> && std::is_same_v<C, const CellExplicit>))
    class proximity_iterator {
        std::vector<size_t>::const_iterator cur;
        std::vector<size_t>::const_iterator end;
        std::vector<size_t>::const_iterator cell_end;
        std::vector<C*> cells;
        std::span<P> container_data;
        double radius;
        R3 center;

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
            while (cur != end &&
                   (cur == cell_end || !((center - container_data[*cur].getX()).euclidNorm() <= radius))) {
                inc();
            }
        }

       public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = P;
        using difference_type = std::ptrdiff_t;
        using pointer = P*;
        using reference = P&;

        proximity_iterator() noexcept : container_data(nullptr), radius(0.0) {}
        proximity_iterator(R3 center, double radius, std::vector<size_t>::const_iterator cur, std::vector<C*> cells,
                           std::span<P> data)
            : cur(cur),
              end(cells.back()->particles().end()),
              cell_end(cells.front()->particles().end()),
              cells(cells),
              container_data(data),
              radius(radius),
              center(center) {
            satisfy();
        }

        reference operator*() const { return container_data[*cur]; }
        pointer operator->() const { return &container_data[*cur]; }

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

        friend bool operator==(const proximity_iterator<P, C>& a, const proximity_iterator<P, C>& b) {
            return a.cur == b.cur;
        }
        friend bool operator!=(const proximity_iterator<P, C>& a, const proximity_iterator<P, C>& b) {
            return !(a == b);
        }
    };
    static_assert(std::forward_iterator<proximity_iterator<Particle, CellExplicit>>);
    static_assert(std::forward_iterator<proximity_iterator<const Particle, const CellExplicit>>);

    std::vector<Particle>::iterator begin();
    [[nodiscard]] std::vector<Particle>::const_iterator begin() const;
    [[nodiscard]] std::vector<Particle>::const_iterator cbegin() const;
    std::vector<Particle>::iterator end();
    [[nodiscard]] std::vector<Particle>::const_iterator end() const;
    [[nodiscard]] std::vector<Particle>::const_iterator cend() const;

    [[nodiscard]] proximity_iterator<Particle, CellExplicit> proximityBegin(R3 center, size_t offset = 0);
    [[nodiscard]] proximity_iterator<Particle, CellExplicit> proximityEnd(R3 center);
    [[nodiscard]] proximity_iterator<const Particle, const CellExplicit> proximityBegin(R3 center,
                                                                                        size_t offset = 0) const;
    [[nodiscard]] proximity_iterator<const Particle, const CellExplicit> proximityEnd(R3 center) const;

    [[nodiscard]] proximity_iterator<Particle, CellExplicit> proximityBegin_no_N3L(R3 center);   // NOLINT
    [[nodiscard]] proximity_iterator<Particle, CellExplicit> proximityEnd_no_N3L(R3 center);     // NOLINT
    [[nodiscard]] proximity_iterator<const Particle, const CellExplicit> proximityBegin_no_N3L(  // NOLINT
        R3 center) const;
    [[nodiscard]] proximity_iterator<const Particle, const CellExplicit> proximityEnd_no_N3L(  // NOLINT
        R3 center) const;

    // boundary and halo iterators

    [[nodiscard]] proximity_iterator<Particle, CellExplicit> haloBegin(
        const std::set<BoundaryLocation>& boundary_types = {BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                            BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                            BoundaryLocation::LEFT, BoundaryLocation::RIGHT});
    [[nodiscard]] proximity_iterator<const Particle, const CellExplicit> haloBegin(
        const std::set<BoundaryLocation>& boundary_types = {BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                            BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                            BoundaryLocation::LEFT, BoundaryLocation::RIGHT}) const;
    [[nodiscard]] proximity_iterator<Particle, CellExplicit> haloEnd(
        const std::set<BoundaryLocation>& boundary_types = {BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                            BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                            BoundaryLocation::LEFT, BoundaryLocation::RIGHT});
    [[nodiscard]] proximity_iterator<const Particle, const CellExplicit> haloEnd(
        const std::set<BoundaryLocation>& boundary_types = {BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                            BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                            BoundaryLocation::LEFT, BoundaryLocation::RIGHT}) const;

    [[nodiscard]] proximity_iterator<Particle, CellExplicit> boundaryBegin(
        const std::set<BoundaryLocation>& boundary_types = {BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                            BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                            BoundaryLocation::LEFT, BoundaryLocation::RIGHT});
    [[nodiscard]] proximity_iterator<const Particle, const CellExplicit> boundaryBegin(
        const std::set<BoundaryLocation>& boundary_types = {BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                            BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                            BoundaryLocation::LEFT, BoundaryLocation::RIGHT}) const;
    [[nodiscard]] proximity_iterator<Particle, CellExplicit> boundaryEnd(
        const std::set<BoundaryLocation>& boundary_types = {BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                            BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                            BoundaryLocation::LEFT, BoundaryLocation::RIGHT});
    [[nodiscard]] proximity_iterator<const Particle, const CellExplicit> boundaryEnd(
        const std::set<BoundaryLocation>& boundary_types = {BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                            BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                            BoundaryLocation::LEFT, BoundaryLocation::RIGHT}) const;

    [[nodiscard]] R3 getDomainSize();

    void prepareForParallelIteration() const {}  // No-op for benchmark container
};
static_assert(ParticleContainer<LinkedCellContainerExplicit>);

}  // namespace mol_sim

#endif
