#ifndef SIMPLE_CONTAINER_H
#define SIMPLE_CONTAINER_H

#include <limits>
#include <set>
#include <vector>

#include "particles/ParticleContainer.h"
#include "particles/boundaries/Boundary.h"

namespace mol_sim {

/**
 * @brief Simple Container for Particles
 *
 * This container implements the concept ParticleContainer.
 * It essentially uses a simple std::vector<Particle> to store the Particles.
 * Thus it is helpful for it to inherit all members (and constructors) directly from std::vector, immediately fulfilling
 * large part of the ParticleContainer concept. It implements the remaining functions for adding a Particle.
 *
 */
class SimpleContainer : public std::vector<Particle> {
    R3 domain_size = R3{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity()};
    double cutoff_radius = std::numeric_limits<double>::infinity();

   public:
    // constructors
    using std::vector<Particle>::vector;

    SimpleContainer(R3 domain_size_arg, double cutoff_radius_arg);

    /**
     * @brief Checks that a given position fits within the domain + the halo region (1 cutoff radius)
     * around the domain.
     *
     * @param v The position to be checked.
     * @return true if the position lies within the domain + halo region
     * @return false else
     */
    inline bool fitsContainer(R3 v);

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
     * @brief Update the position of a given Particle in the Container.
     *
     * @param p Iterator to the Particle to update.
     * @param new_x New position to set.
     *
     * @return Iterator to the next Particle.
     */
    std::vector<Particle>::iterator updateParticlePosition(std::vector<Particle>::iterator p, R3 new_x);

    /**
     * @brief Removes a given particle from the container.
     *
     * @param p Iterator to the Particle to be removed.
     *
     * @return Iterator to the next Particle after the removed one.
     */
    std::vector<Particle>::iterator eraseParticle(std::vector<Particle>::iterator p);

    /**
     * @brief Iterator template that iterates over all particles that apply a force on a given particle.
     *
     * Satisfies the forward iterator concept.
     * Can be instantiated as mutable and const iterator with P = (const) Particle.
     * @tparam P Either Particle or const Particle
     *
     * Only iterates over particles within the given radius of the center, or all if the radius is infinite.
     * Or is used for the iteration over halo and boundary particles.
     *
     * Enables cutoff radius, BUT still iterates over all particles when checking proximity/boundary/halo.
     * To optimize that behavior and enable optimizations, use the LinkedCellContainer.
     * Can be used for the N3L optimization via the offset.
     *
     * This class is mainly used to fulfill the interface with the calculateF() function in Simulation.
     */
    template <typename P>
        requires(std::is_same_v<P, Particle> || std::is_same_v<P, const Particle>)
    class proximity_iterator {
        P* cur;
        P* end;
        double radius;
        R3 center_or_domain;
        bool prox;
        std::set<BoundaryLocation> locations;

        bool fitBoundary() {  // NOLINT
            double effective_radius = -radius;
            if (cur->getX()[0] > center_or_domain[0] || cur->getX()[1] > center_or_domain[1] ||
                cur->getX()[2] > center_or_domain[2] || cur->getX()[0] < 0.0 || cur->getX()[1] < 0.0 ||
                cur->getX()[2] < 0.0) {
                return false;
            }
            for (auto location : locations) {
                switch (location) {
                    case BoundaryLocation::UPPER: {
                        if (cur->getX()[2] >= center_or_domain[2] - effective_radius &&
                            cur->getX()[2] <= center_or_domain[2]) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::LOWER: {
                        if (cur->getX()[2] >= 0.0 && cur->getX()[2] <= effective_radius) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::FRONT: {
                        if (cur->getX()[1] >= 0.0 && cur->getX()[1] <= effective_radius) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::BACK: {
                        if (cur->getX()[1] >= center_or_domain[1] - effective_radius &&
                            cur->getX()[1] <= center_or_domain[1]) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::LEFT: {
                        if (cur->getX()[0] >= 0.0 && cur->getX()[0] <= effective_radius) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::RIGHT: {
                        if (cur->getX()[0] >= center_or_domain[0] - effective_radius &&
                            cur->getX()[0] <= center_or_domain[0]) {
                            return true;
                        }
                        break;
                    }
                    default:
                        break;
                }
            }
            return false;
        }
        bool fitHalo() {  // NOLINT
            if (cur->getX()[0] > center_or_domain[0] + radius || cur->getX()[1] > center_or_domain[1] + radius ||
                cur->getX()[2] > center_or_domain[2] + radius || cur->getX()[0] < -radius || cur->getX()[1] < -radius ||
                cur->getX()[2] < -radius) {
                return false;
            }
            for (auto location : locations) {
                switch (location) {
                    case BoundaryLocation::UPPER: {
                        if (cur->getX()[2] > center_or_domain[2] && cur->getX()[2] <= center_or_domain[2] + radius) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::LOWER: {
                        if (cur->getX()[2] >= -radius && cur->getX()[2] < 0.0) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::FRONT: {
                        if (cur->getX()[1] >= -radius && cur->getX()[1] < 0.0) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::BACK: {
                        if (cur->getX()[1] > center_or_domain[1] && cur->getX()[1] <= center_or_domain[1] + radius) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::LEFT: {
                        if (cur->getX()[0] >= -radius && cur->getX()[0] < 0.0) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::RIGHT: {
                        if (cur->getX()[0] > center_or_domain[0] && cur->getX()[0] <= center_or_domain[0] + radius) {
                            return true;
                        }
                        break;
                    }
                    default:
                        break;
                }
            }
            return false;
        }
        void satisfy() {
            if (prox) {  // proximity check
                if (std::isinf(radius)) {
                    return;
                }
                while (cur != end && !((center_or_domain - cur->getX()).euclidNorm() <= radius)) {
                    ++cur;
                }
            } else {
                if (radius < 0.0) {  // boundary check
                    while (cur != end && !(fitBoundary())) {
                        ++cur;
                    }
                } else {  // halo check
                    while (cur != end && !(fitHalo())) {
                        ++cur;
                    }
                }
            }
        }

       public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = P;
        using difference_type = std::ptrdiff_t;
        using pointer = P*;
        using reference = P&;

        proximity_iterator() noexcept : cur(nullptr), end(nullptr), radius(0.0), prox(true) {}
        proximity_iterator(R3 center_or_domain, double radius, P* cur, P* end, size_t offset, bool prox = true,
                           const std::set<BoundaryLocation>& locations = {})
            : cur(cur + offset),
              end(end),
              radius(radius),
              center_or_domain(center_or_domain),
              prox(prox),
              locations(locations) {
            satisfy();
        }

        reference operator*() const noexcept { return *cur; }
        pointer operator->() const noexcept { return cur; }

        proximity_iterator<P>& operator++() {
            ++cur;  // NOLINT
            satisfy();
            return *this;
        }

        proximity_iterator<P> operator++(int) {
            proximity_iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        friend bool operator==(const proximity_iterator<P>& a, const proximity_iterator<P>& b) noexcept {
            return a.cur == b.cur;
        }
        friend bool operator!=(const proximity_iterator<P>& a, const proximity_iterator<P>& b) noexcept {
            return !(a == b);
        }

        [[nodiscard]] R3 getCenter() const noexcept { return center_or_domain; }
        [[nodiscard]] double getRadius() const noexcept { return radius; }
        [[nodiscard]] bool isProximity() const noexcept { return prox; }
        [[nodiscard]] std::set<BoundaryLocation> getLocations() const noexcept { return locations; }
    };
    static_assert(std::forward_iterator<proximity_iterator<Particle>>);
    static_assert(std::forward_iterator<proximity_iterator<const Particle>>);

    /**
     * @brief Mutable Iterator over particles in proximity.
     *
     * @param center Center point to check proximity from (position of the particle).
     * @param offset Offset from the beginning of the container to start the iteration (used for N3L optimization).
     *
     * @return Mutable iterator to the first particle within the given radius of the center.
     */
    [[nodiscard]] proximity_iterator<Particle> proximityBegin(R3 center, size_t offset = 0);

    /**
     * @brief Const Iterator over particles in proximity.
     *
     * @param center Center point to check proximity from (position of the particle).
     * @param offset Offset from the beginning of the container to start the iteration (used for N3L optimization).
     *
     * @return Const iterator to the first particle within the given radius of the center.
     */
    [[nodiscard]] proximity_iterator<const Particle> proximityBegin(R3 center, size_t offset = 0) const;

    /**
     * @brief Mutable Iterator over particles in proximity.
     *
     * @param center Center point to check proximity from (position of the particle).
     *
     * @return Mutable iterator after the last particle within the given radius of the center.
     */
    [[nodiscard]] proximity_iterator<Particle> proximityEnd(R3 center);

    /**
     * @brief Const Iterator over particles in proximity.
     *
     * @param center Center point to check proximity from (position of the particle).
     *
     * @return Const iterator after the last particle within the given radius of the center.
     */
    [[nodiscard]] proximity_iterator<const Particle> proximityEnd(R3 center) const;

    /**
     * @brief Mutable Iterator over particles in proximity, but does not use the N3L optimization.
     *
     * @param center Center point to check proximity from (position of the particle).
     *
     * @return Mutable iterator to the first particle within the given radius of the center.
     */
    [[nodiscard]] proximity_iterator<Particle> proximityBegin_no_N3L(R3 center);  // NOLINT

    /**
     * @brief Mutable Iterator over particles in proximity, but does not use the N3L optimization.
     *
     * @param center Center point to check proximity from (position of the particle).
     *
     * @return Mutable iterator after the last particle within the given radius of the center.
     */
    [[nodiscard]] proximity_iterator<Particle> proximityEnd_no_N3L(R3 center);  // NOLINT

    /**
     * @brief Const Iterator over particles in proximity, but does not use the N3L optimization.
     *
     * @param center Center point to check proximity from (position of the particle).
     *
     * @return Const iterator to the first particle within the given radius of the center.
     */
    [[nodiscard]] proximity_iterator<const Particle> proximityBegin_no_N3L(R3 center) const;  // NOLINT

    /**
     * @brief Const Iterator over particles in proximity, but does not use the N3L optimization.
     *
     * @param center Center point to check proximity from (position of the particle).
     *
     * @return Const iterator after the last particle within the given radius of the center.
     */
    [[nodiscard]] proximity_iterator<const Particle> proximityEnd_no_N3L(R3 center) const;  // NOLINT

    // boundary and halo iterators

    /**
     * @brief Iterator over particles in halo, which is the same region around the domain as the halo cells in the
     * LinkedCellContainer.
     *
     * @param locations Boundary types to specify which halo boundary to iterate over. Defaults to all sides.
     *
     * @return Iterator to the first particle within the given halo.
     */
    [[nodiscard]] proximity_iterator<Particle> haloBegin(const std::set<BoundaryLocation>& locations = {
                                                             BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                             BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                             BoundaryLocation::LEFT, BoundaryLocation::RIGHT});

    /**
     * @brief Const Iterator over particles in halo, which is the same region around the domain as the halo cells in the
     * LinkedCellContainer.
     *
     * @param locations Boundary types to specify which halo boundary to iterate over. Defaults to all sides.
     *
     * @return Const iterator to the first particle within the given halo.
     */
    [[nodiscard]] proximity_iterator<const Particle> haloBegin(const std::set<BoundaryLocation>& locations = {
                                                                   BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                                   BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                                   BoundaryLocation::LEFT,
                                                                   BoundaryLocation::RIGHT}) const;

    /**
     * @brief Iterator over particles in halo, which is the same region around the domain as the halo cells in the
     * LinkedCellContainer.
     *
     * @param locations Boundary types to specify which halo boundary to iterate over. Defaults to all sides.
     *
     * @return Iterator after the last particle within the given halo.
     */
    [[nodiscard]] proximity_iterator<Particle> haloEnd(const std::set<BoundaryLocation>& locations = {
                                                           BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                           BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                           BoundaryLocation::LEFT, BoundaryLocation::RIGHT});

    /**
     * @brief Const Iterator over particles in halo, which is the same region around the domain as the halo cells in
     * the LinkedCellContainer.
     *
     * @param locations Boundary types to specify which halo boundary to iterate over. Defaults to all sides.
     *
     * @return Const iterator after the last particle within the given halo.
     */
    [[nodiscard]] proximity_iterator<const Particle> haloEnd(const std::set<BoundaryLocation>& locations = {
                                                                 BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                                 BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                                 BoundaryLocation::LEFT,
                                                                 BoundaryLocation::RIGHT}) const;

    /**
     * @brief Iterator over particles in boundary, which is the same region around the domain as the boundary cells in
     * the LinkedCellContainer.
     *
     * @param locations Boundary types to specify which boundary to iterate over. Defaults to all sides.
     *
     * @return Iterator to the first particle within the given boundary.
     */
    [[nodiscard]] proximity_iterator<Particle> boundaryBegin(const std::set<BoundaryLocation>& locations = {
                                                                 BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                                 BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                                 BoundaryLocation::LEFT, BoundaryLocation::RIGHT});

    /**
     * @brief Const Iterator over particles in boundary, which is the same region around the domain as the boundary
     * cells in the LinkedCellContainer.
     *
     * @param locations Boundary types to specify which boundary to iterate over. Defaults to all sides.
     *
     * @return Const iterator to the first particle within the given boundary.
     */
    [[nodiscard]] proximity_iterator<const Particle> boundaryBegin(const std::set<BoundaryLocation>& locations = {
                                                                       BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                                       BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                                       BoundaryLocation::LEFT,
                                                                       BoundaryLocation::RIGHT}) const;

    /**
     * @brief Iterator over particles in boundary, which is the same region around the domain as the boundary
     * cells in the LinkedCellContainer.
     *
     * @param locations Boundary types to specify which boundary to iterate over. Defaults to all sides.
     *
     * @return Iterator after the last particle within the given boundary.
     */
    [[nodiscard]] proximity_iterator<Particle> boundaryEnd(const std::set<BoundaryLocation>& locations = {
                                                               BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                               BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                               BoundaryLocation::LEFT, BoundaryLocation::RIGHT});

    /**
     * @brief Const Iterator over particles in boundary, which is the same region around the domain as the boundary
     * cells in the LinkedCellContainer.
     *
     * @param locations Boundary types to specify which boundary to iterate over. Defaults to all sides.
     *
     * @return Const iterator after the last particle within the given boundary.
     */
    [[nodiscard]] proximity_iterator<const Particle> boundaryEnd(const std::set<BoundaryLocation>& locations = {
                                                                     BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                                     BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                                     BoundaryLocation::LEFT,
                                                                     BoundaryLocation::RIGHT}) const;
};
static_assert(ParticleContainer<SimpleContainer>);

}  // namespace mol_sim

#endif