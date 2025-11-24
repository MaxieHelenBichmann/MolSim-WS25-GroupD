#ifndef SIMPLE_CONTAINER_H
#define SIMPLE_CONTAINER_H

#include <vector>

#include "particles/ParticleContainer.h"

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
   public:
    // constructors
    using std::vector<Particle>::vector;

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
        R3 center;

        void satisfy() {
            if (std::isinf(radius)) {
                return;
            }
            while (cur != end && !((center - cur->getX()).euclidNorm() <= radius)) {
                ++cur;  // NOLINT
            }
        }

       public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = Particle;
        using difference_type = std::ptrdiff_t;
        using pointer = const Particle*;
        using reference = const Particle&;

        const_proximity_iterator() noexcept : cur(nullptr), end(nullptr), radius(0.0) {}
        const_proximity_iterator(R3 center, double radius, const Particle* cur, const Particle* end, size_t offset)
            : cur(cur + offset), end(end), radius(radius), center(center) {  // NOLINT
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
        R3 center;

        void satisfy() {
            if (std::isinf(radius)) {
                return;
            }
            while (cur != end && !((center - cur->getX()).euclidNorm() <= radius)) {
                ++cur;  // NOLINT
            }
        }

       public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = Particle;
        using difference_type = std::ptrdiff_t;
        using pointer = Particle*;
        using reference = Particle&;

        proximity_iterator() noexcept : cur(nullptr), end(nullptr), radius(0.0) {}
        proximity_iterator(R3 center, double radius, Particle* cur, Particle* end, size_t offset)
            : cur(cur + offset), end(end), radius(radius), center(center) {  // NOLINT
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
};

}  // namespace mol_sim

#endif
