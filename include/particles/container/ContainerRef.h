#ifndef CONTAINER_REF_H
#define CONTAINER_REF_H

#include <variant>

#include "particles/container/SimpleContainer.h"

namespace mol_sim {

/**
 * @brief Container Reference for Particles
 *
 * This container implements the concept ParticleContainer.
 * It essentially stores a pointer (not null) to an arbitrary ParticleContainer - currently only SimpleContainer
 * possible - but does not gain ownership over the elements in the Container. It calls the methods of the actual
 * SimpleContainer from which it was constructed. It is virtually a Generic Type for a ParticleContainer, which can be
 * instantiated.
 * Not being based on templates, it does not increase compile time as much, but checks for which function to call are
 * made at run time.
 *
 * Only use for UNCRITICAL functions/operations (e.g., IO). DO NOT USE IN PERFORMANCE-RELEVANT FUNCTIONS!
 *
 * Do not take References of a ContainerRef object, as it is a reference itself, and trivially copyable.
 * Alway pass-by-value (copy) - similar to std::span.
 *
 * @note maybe will get erased later, if IO functions will be templated as well and compile time and code size
 * acceptable
 *
 */
class ContainerRef {
    /**
     * Pointer to the current instance of the ContainerRef.
     * More containers can be added in the std::variant later.
     */
    std::variant<SimpleContainer*> instance;  // NOLINT

   public:
    // constructors

    /**
     * @brief Constructor, initializing ContainerRef with a reference to a SimpleContainer.
     */
    ContainerRef(SimpleContainer& c);

    // retrieve data
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
     * @brief Directly constructing a Particle with its required parameters in-place.
     *
     * @param x_arg Initial coordinates of the Particle.
     * @param v_arg Initial velocities of the Particle.
     * @param m_arg Mass of the Particle.
     */
    void addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg);

    /**
     * @brief Directly constructing a Particle with its required parameters, including type, in-place.
     *
     * @param x_arg Initial coordinates of the Particle.
     * @param v_arg Initial velocities of the Particle.
     * @param m_arg Mass of the Particle.
     * @param type Type of the Particle.
     */
    void addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg, int type);
    /**
     * @brief Directly constructing a Particle with its required parameters, including type, in-place.
     *
     * @param x_arg Initial coordinates of the Particle.
     * @param v_arg Initial velocities of the Particle.
     * @param m_arg Mass of the Particle.
     * @param type Type of the Particle.
     */
    void addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg, double epsilon_arg,
                     double sigma_arg);
    void addParticle(Vector<double, 3> x_arg, Vector<double, 3> v_arg, double m_arg, double epsilon_arg,
                     double sigma_arg, int type);

    // iterators

    /**
     * @brief Mutable iterator.
     *
     * @return Mutable iterator to the first element of the container.
     */
    std::vector<Particle>::iterator begin();

    /**
     * @brief Const iterator.
     *
     * @return Const iterator to the first element of the container.
     */
    [[nodiscard]] std::vector<Particle>::const_iterator begin() const;

    /**
     * @brief Const iterator.
     *
     * @return Const iterator to the first element of the container.
     */
    [[nodiscard]] std::vector<Particle>::const_iterator cbegin() const;

    /**
     * @brief Mutable iterator.
     *
     * @return Mutable past-the-end iterator of the container.
     */
    std::vector<Particle>::iterator end();

    /**
     * @brief Const iterator.
     *
     * @return Const past-the-end iterator of the container.
     */
    [[nodiscard]] std::vector<Particle>::const_iterator end() const;

    /**
     * @brief Const iterator.
     *
     * @return Const past-the-end iterator of the container
     */
    [[nodiscard]] std::vector<Particle>::const_iterator cend() const;
};

}  // namespace mol_sim

#endif
