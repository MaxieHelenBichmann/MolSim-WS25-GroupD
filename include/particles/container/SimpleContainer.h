#ifndef SIMPLE_CONTAINER_H
#define SIMPLE_CONTAINER_H

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
    R3 domain_size;
    double cutoff_radius = 0.0;

   public:
    // constructors
    using std::vector<Particle>::vector;

    SimpleContainer(R3 domain_size_arg, double cutoff_radius_arg);

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
     * @brief Update the position of a given Particle in the Container.
     *
     * @param p Iterator to the Particle to update.
     * @param new_x New position to set.
     */
    void updateParticlePosition(std::vector<Particle>::iterator p, R3 new_x);

    /**
     * @brief Removes a given particle from the container.
     *
     * @param p Iterator to the Particle to be removed.
     *
     * @return Iterator to the next Particle after the removed one.
     */
    std::vector<Particle>::iterator eraseParticle(std::vector<Particle>::iterator p);

    /**
     * @brief Const iterator that iterates over all particles that apply a force on a given particle.
     *
     * Satisfies the forward iterator concept.
     * Only iterates over particles within the given radius of the center, or all if the radius is infinite.
     *
     * Enables cutoff radius, BUT still iterates over all particles when checking proximity.
     * To optimize that behavior and enable optimizations, use the LinkedCellContainer.
     *
     * This class is mainly used to fulfill the interface with the calculateF() function in Simulation.
     */
    class const_proximity_iterator {
        const Particle* cur;
        const Particle* end;
        double radius;
        R3 center_or_domain;
        bool prox;
        std::set<BoundaryLocation> locations;

        bool fitBoundary() {  // NOLINT
            if (cur->getX()[0] > center_or_domain[0] || cur->getX()[1] > center_or_domain[1] ||
                cur->getX()[2] > center_or_domain[2] || cur->getX()[0] < 0.0 || cur->getX()[1] < 0.0 ||
                cur->getX()[2] < 0.0) {
                return false;
            }
            for (auto location : locations) {
                switch (location) {
                    case BoundaryLocation::UPPER: {
                        if (cur->getX()[2] >= center_or_domain[2] + radius && cur->getX()[2] <= center_or_domain[2]) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::LOWER: {
                        if (cur->getX()[2] >= 0.0 && cur->getX()[2] <= -radius) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::FRONT: {
                        if (cur->getX()[1] >= 0.0 && cur->getX()[1] <= -radius) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::BACK: {
                        if (cur->getX()[1] >= center_or_domain[1] + radius && cur->getX()[1] <= center_or_domain[1]) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::LEFT: {
                        if (cur->getX()[0] >= 0.0 && cur->getX()[0] <= -radius) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::RIGHT: {
                        if (cur->getX()[0] >= center_or_domain[0] + radius && cur->getX()[0] <= center_or_domain[0]) {
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
                        if (cur->getX()[2] >= center_or_domain[2] && cur->getX()[2] <= center_or_domain[2] + radius) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::LOWER: {
                        if (cur->getX()[2] >= -radius && cur->getX()[2] <= 0.0) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::FRONT: {
                        if (cur->getX()[1] >= -radius && cur->getX()[1] <= 0.0) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::BACK: {
                        if (cur->getX()[1] >= center_or_domain[1] && cur->getX()[1] <= center_or_domain[1] + radius) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::LEFT: {
                        if (cur->getX()[0] >= -radius && cur->getX()[0] <= 0.0) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::RIGHT: {
                        if (cur->getX()[0] >= center_or_domain[0] && cur->getX()[0] <= center_or_domain[0] + radius) {
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
        using value_type = Particle;
        using difference_type = std::ptrdiff_t;
        using pointer = const Particle*;
        using reference = const Particle&;

        const_proximity_iterator() noexcept : cur(nullptr), end(nullptr), radius(0.0), prox(true) {}
        const_proximity_iterator(R3 center_or_domain, double radius, const Particle* cur, const Particle* end,
                                 size_t offset, bool prox = true, std::set<BoundaryLocation> locations = {})
            : cur(cur + offset),
              end(end),
              radius(radius),
              center_or_domain(center_or_domain),
              prox(prox),
              locations(std::move(locations)) {  // NOLINT
            satisfy();
        }

        reference operator*() const { return *cur; }
        pointer operator->() const { return cur; }

        const_proximity_iterator& operator++() {
            ++cur;  // NOLINT
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
     * @brief Iterator that iterates over all particles that apply a force on a given particle.
     *
     * Satisfies the forward iterator concept.
     * Only iterates over particles within the given radius of the center, or all if the radius is infinite.
     *
     * Enables cutoff radius, BUT still iterates over all particles when checking proximity.
     * To optimize that behavior and enable optimizations, use the LinkedCellContainer.
     *
     * This class is mainly used to fulfill the interface with the calculateF() function in Simulation.
     */
    class proximity_iterator {
        Particle* cur;
        Particle* end;
        double radius;
        R3 center_or_domain;
        bool prox;
        std::set<BoundaryLocation> locations;

        bool fitBoundary() {  // NOLINT
            if (cur->getX()[0] > center_or_domain[0] || cur->getX()[1] > center_or_domain[1] ||
                cur->getX()[2] > center_or_domain[2] || cur->getX()[0] < 0.0 || cur->getX()[1] < 0.0 ||
                cur->getX()[2] < 0.0) {
                return false;
            }
            for (auto location : locations) {
                switch (location) {
                    case BoundaryLocation::UPPER: {
                        if (cur->getX()[2] >= center_or_domain[2] + radius && cur->getX()[2] <= center_or_domain[2]) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::LOWER: {
                        if (cur->getX()[2] >= 0.0 && cur->getX()[2] <= -radius) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::FRONT: {
                        if (cur->getX()[1] >= 0.0 && cur->getX()[1] <= -radius) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::BACK: {
                        if (cur->getX()[1] >= center_or_domain[1] + radius && cur->getX()[1] <= center_or_domain[1]) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::LEFT: {
                        if (cur->getX()[0] >= 0.0 && cur->getX()[0] <= -radius) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::RIGHT: {
                        if (cur->getX()[0] >= center_or_domain[0] + radius && cur->getX()[0] <= center_or_domain[0]) {
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
                        if (cur->getX()[2] >= center_or_domain[2] && cur->getX()[2] <= center_or_domain[2] + radius) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::LOWER: {
                        if (cur->getX()[2] >= -radius && cur->getX()[2] <= 0.0) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::FRONT: {
                        if (cur->getX()[1] >= -radius && cur->getX()[1] <= 0.0) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::BACK: {
                        if (cur->getX()[1] >= center_or_domain[1] && cur->getX()[1] <= center_or_domain[1] + radius) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::LEFT: {
                        if (cur->getX()[0] >= -radius && cur->getX()[0] <= 0.0) {
                            return true;
                        }
                        break;
                    }
                    case BoundaryLocation::RIGHT: {
                        if (cur->getX()[0] >= center_or_domain[0] && cur->getX()[0] <= center_or_domain[0] + radius) {
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
        using value_type = Particle;
        using difference_type = std::ptrdiff_t;
        using pointer = Particle*;
        using reference = Particle&;

        proximity_iterator() noexcept : cur(nullptr), end(nullptr), radius(0.0), prox(true) {}
        proximity_iterator(R3 center_or_domain, double radius, Particle* cur, Particle* end, size_t offset,
                           bool prox = true, std::set<BoundaryLocation> locations = {})
            : cur(cur + offset),
              end(end),
              radius(radius),
              center_or_domain(center_or_domain),
              prox(prox),
              locations(std::move(locations)) {
            satisfy();
        }

        reference operator*() const { return *cur; }
        pointer operator->() const { return cur; }

        proximity_iterator& operator++() {
            ++cur;  // NOLINT
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
     *
     * @return Const iterator after the last particle within the given radius of the center.
     */
    [[nodiscard]] const_proximity_iterator proximityEnd(R3 center, double radius) const;

    // boundary and halo iterators

    /**
     * @brief Iterator over particles in halo cells.
     *
     * @param locations Boundary types to specify which halo cells to iterate over. Defaults to all sides.
     *
     * @return Iterator to the first particle within the given halo cells.
     */
    [[nodiscard]] proximity_iterator haloBegin(const std::set<BoundaryLocation>& locations = {
                                                   BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                   BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                   BoundaryLocation::LEFT, BoundaryLocation::RIGHT});

    /**
     * @brief Const Iterator over particles in halo cells.
     *
     * @param locations Boundary types to specify which halo cells to iterate over. Defaults to all sides.
     *
     * @return Const iterator to the first particle within the given halo cells.
     */
    [[nodiscard]] const_proximity_iterator haloBegin(const std::set<BoundaryLocation>& locations = {
                                                         BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                         BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                         BoundaryLocation::LEFT, BoundaryLocation::RIGHT}) const;

    /**
     * @brief Iterator over particles in halo cells.
     *
     * @param locations Boundary types to specify which halo cells to iterate over. Defaults to all sides.
     *
     * @return Iterator after the last particle within the given halo cells.
     */
    [[nodiscard]] proximity_iterator haloEnd(const std::set<BoundaryLocation>& locations = {
                                                 BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                 BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                 BoundaryLocation::LEFT, BoundaryLocation::RIGHT});

    /**
     * @brief Const Iterator over particles in halo cells.
     *
     * @param locations Boundary types to specify which halo cells to iterate over. Defaults to all sides.
     *
     * @return Const iterator after the last particle within the given halo cells.
     */
    [[nodiscard]] const_proximity_iterator haloEnd(const std::set<BoundaryLocation>& locations = {
                                                       BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                       BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                       BoundaryLocation::LEFT, BoundaryLocation::RIGHT}) const;

    /**
     * @brief Iterator over particles in boundary cells.
     *
     * @param locations Boundary types to specify which boundary cells to iterate over. Defaults to all sides.
     *
     * @return Iterator to the first particle within the given boundary cells.
     */
    [[nodiscard]] proximity_iterator boundaryBegin(const std::set<BoundaryLocation>& locations = {
                                                       BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                       BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                       BoundaryLocation::LEFT, BoundaryLocation::RIGHT});

    /**
     * @brief Const Iterator over particles in boundary cells.
     *
     * @param locations Boundary types to specify which boundary cells to iterate over. Defaults to all sides.
     *
     * @return Const iterator to the first particle within the given boundary cells.
     */
    [[nodiscard]] const_proximity_iterator boundaryBegin(const std::set<BoundaryLocation>& locations = {
                                                             BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                             BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                             BoundaryLocation::LEFT, BoundaryLocation::RIGHT}) const;

    /**
     * @brief Iterator over particles in boundary cells.
     *
     * @param locations Boundary types to specify which boundary cells to iterate over. Defaults to all sides.
     *
     * @return Iterator after the last particle within the given boundary cells.
     */
    [[nodiscard]] proximity_iterator boundaryEnd(const std::set<BoundaryLocation>& locations = {
                                                     BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                     BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                     BoundaryLocation::LEFT, BoundaryLocation::RIGHT});

    /**
     * @brief Const Iterator over particles in boundary cells.
     *
     * @param locations Boundary types to specify which boundary cells to iterate over. Defaults to all sides.
     *
     * @return Const iterator after the last particle within the given boundary cells.
     */
    [[nodiscard]] const_proximity_iterator boundaryEnd(const std::set<BoundaryLocation>& locations = {
                                                           BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                           BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                           BoundaryLocation::LEFT, BoundaryLocation::RIGHT}) const;
};

}  // namespace mol_sim

#endif