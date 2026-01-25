#ifndef CONTAINER_REF_H
#define CONTAINER_REF_H

#include <set>
#include <type_traits>
#include <variant>
#include <vector>

#include "particles/Particle.h"
#include "particles/container/LinkedCellContainer.h"
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
using CONTAINER_REF = std::variant<SimpleContainer*, LinkedCellContainer*>;
class ContainerRef {
    /**
     * Pointer to the current instance of the ContainerRef.
     * More containers can be added in the std::variant later.
     */
    CONTAINER_REF instance;  // NOLINT

   public:
    // constructors
    ContainerRef();
    /**
     * @brief Constructor, initializing ContainerRef with a reference to a SimpleContainer.
     */
    ContainerRef(SimpleContainer& c);

    /**
     * @brief Constructor, initializing ContainerRef with a reference to a LinkedCellContainer.
     */
    ContainerRef(LinkedCellContainer& c);

    /**
     * @brief Constructor, initializing ContainerRef with another CONTAINER_REF.
     */
    ContainerRef(CONTAINER_REF& c);

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
     * @brief Directly constructing a Particle with its required parameters, including type, in-place
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
     * @brief Removes a given particle from the container.
     *
     * @param p Iterator to the Particle to be removed.
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

    /**
     * @brief Get the Instance object.
     *
     */
    CONTAINER_REF getInstance();

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

    /**
     * @brief Iterator that iterates over all particles that apply a force on a given particle.
     *
     * Satisfies the forward iterator concept.
     * Can be instantiated as mutable and const iterator with P = (const) Particle and C = (const) Cell.
     * @tparam P Either Particle or const Particle
     * @tparam C Either Cell or const Cell

     * Only iterates over particles within the given radius of the center, or all if the radius is infinite.
     * Or is used for the iteration over halo and boundary particles.
     */
    template <typename P, typename C>
        requires((std::is_same_v<P, Particle> && std::is_same_v<C, Cell>) ||
                 (std::is_same_v<P, const Particle> && std::is_same_v<C, const Cell>))
    class proximity_iterator {
        std::variant<SimpleContainer::proximity_iterator<P>, LinkedCellContainer::proximity_iterator<P, C>> cur;

       public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = P;
        using difference_type = std::ptrdiff_t;
        using pointer = P*;
        using reference = P&;

        proximity_iterator() noexcept = default;
        proximity_iterator(SimpleContainer::proximity_iterator<P> it) : cur(it) {}         // NOLINT
        proximity_iterator(LinkedCellContainer::proximity_iterator<P, C> it) : cur(it) {}  // NOLINT

        reference operator*() const {
            return std::visit([](auto& c) -> reference { return c.operator*(); }, cur);
        }
        pointer operator->() const {
            return std::visit([](auto& c) { return c.operator->(); }, cur);
        }
        proximity_iterator& operator++() {
            std::visit([](auto& c) { ++c; }, cur);
            return *this;
        }

        proximity_iterator<P, C> operator++(int) {
            proximity_iterator<P, C> tmp = *this;
            ++(*this);
            return tmp;
        }

        friend bool operator==(const proximity_iterator<P, C>& a, const proximity_iterator<P, C>& b) {
            return std::visit(
                [](const auto& lhs, const auto& rhs) {
                    using L = std::decay_t<decltype(lhs)>;
                    using R = std::decay_t<decltype(rhs)>;
                    if constexpr (!std::is_same_v<L, R>) {
                        return false;
                    } else {
                        return lhs == rhs;
                    }
                },
                a.cur, b.cur);
        }
        friend bool operator!=(const proximity_iterator<P, C>& a, const proximity_iterator<P, C>& b) {
            return !(a == b);
        }

        operator LinkedCellContainer::proximity_iterator<P, C>() {
            return std::get<LinkedCellContainer::proximity_iterator<P, C>>(cur);
        }

        operator SimpleContainer::proximity_iterator<P>() {
            return std::get<SimpleContainer::proximity_iterator<P>>(cur);
        }
    };
    static_assert(std::forward_iterator<proximity_iterator<Particle, Cell>>);
    static_assert(std::forward_iterator<proximity_iterator<const Particle, const Cell>>);

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

    /**
     * @brief Mutable Iterator over particles in proximity, but does not use the N3L optimization.
     *
     * @param center Center point to check proximity from (position of the particle).
     *
     * @return Mutable iterator to the first particle within the given radius of the center.
     */
    [[nodiscard]] proximity_iterator<Particle, Cell> proximityBegin_no_N3L(R3 center);  // NOLINT

    /**
     * @brief Mutable Iterator over particles in proximity, but does not use the N3L optimization.
     *
     * @param center Center point to check proximity from (position of the particle).
     *
     * @return Mutable iterator after the last particle within the given radius of the center.
     */
    [[nodiscard]] proximity_iterator<Particle, Cell> proximityEnd_no_N3L(R3 center);  // NOLINT

    /**
     * @brief Const Iterator over particles in proximity, but does not use the N3L optimization.
     *
     * @param center Center point to check proximity from (position of the particle).
     *
     * @return Const iterator to the first particle within the given radius of the center.
     */
    [[nodiscard]] proximity_iterator<const Particle, const Cell> proximityBegin_no_N3L(R3 center) const;  // NOLINT

    /**
     * @brief Const Iterator over particles in proximity, but does not use the N3L optimization.
     *
     * @param center Center point to check proximity from (position of the particle).
     *
     * @return Const iterator after the last particle within the given radius of the center.
     */
    [[nodiscard]] proximity_iterator<const Particle, const Cell> proximityEnd_no_N3L(R3 center) const;  // NOLINT

    // boundary and halo iterators

    /**
     * @brief Iterator over particles in halo.
     *
     * @param locations Boundary types to specify which halo boundary to iterate over. Defaults to all sides.
     *
     * @return Iterator to the first particle within the given halo.
     */
    [[nodiscard]] proximity_iterator<Particle, Cell> haloBegin(const std::set<BoundaryLocation>& locations = {
                                                                   BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                                   BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                                   BoundaryLocation::LEFT, BoundaryLocation::RIGHT});

    /**
     * @brief Const Iterator over particles in halo.
     *
     * @param locations Boundary types to specify which halo boundary to iterate over. Defaults to all sides.
     *
     * @return Const iterator to the first particle within the given halo.
     */
    [[nodiscard]] proximity_iterator<const Particle, const Cell> haloBegin(
        const std::set<BoundaryLocation>& locations = {BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                       BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                       BoundaryLocation::LEFT, BoundaryLocation::RIGHT}) const;

    /**
     * @brief Iterator over particles in halo.
     *
     * @param locations Boundary types to specify which halo boundary to iterate over. Defaults to all sides.
     *
     * @return Iterator after the last particle within the given halo.
     */
    [[nodiscard]] proximity_iterator<Particle, Cell> haloEnd(const std::set<BoundaryLocation>& locations = {
                                                                 BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                                 BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                                 BoundaryLocation::LEFT, BoundaryLocation::RIGHT});

    /**
     * @brief Const Iterator over particles in halo.
     *
     * @param locations Boundary types to specify which halo boundary to iterate over. Defaults to all sides.
     *
     * @return Const iterator after the last particle within the given halo.
     */
    [[nodiscard]] proximity_iterator<const Particle, const Cell> haloEnd(
        const std::set<BoundaryLocation>& locations = {BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                       BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                       BoundaryLocation::LEFT, BoundaryLocation::RIGHT}) const;

    /**
     * @brief Iterator over particles in boundary.
     *
     * @param locations Boundary types to specify which boundary to iterate over. Defaults to all sides.
     *
     * @return Iterator to the first particle within the given boundary.
     */
    [[nodiscard]] proximity_iterator<Particle, Cell> boundaryBegin(const std::set<BoundaryLocation>& locations = {
                                                                       BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                                       BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                                       BoundaryLocation::LEFT,
                                                                       BoundaryLocation::RIGHT});

    /**
     * @brief Const Iterator over particles in boundary.
     *
     * @param locations Boundary types to specify which boundary to iterate over. Defaults to all sides.
     *
     * @return Const iterator to the first particle within the given boundary.
     */
    [[nodiscard]] proximity_iterator<const Particle, const Cell> boundaryBegin(
        const std::set<BoundaryLocation>& locations = {BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                       BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                       BoundaryLocation::LEFT, BoundaryLocation::RIGHT}) const;

    /**
     * @brief Iterator over particles in boundary.
     *
     * @param locations Boundary types to specify which boundary to iterate over. Defaults to all sides.
     *
     * @return Iterator after the last particle within the given boundary.
     */
    [[nodiscard]] proximity_iterator<Particle, Cell> boundaryEnd(const std::set<BoundaryLocation>& locations = {
                                                                     BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                                     BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                                     BoundaryLocation::LEFT, BoundaryLocation::RIGHT});

    /**
     * @brief Const Iterator over particles in boundary.
     *
     * @param locations Boundary types to specify which boundary to iterate over. Defaults to all sides.
     *
     * @return Const iterator after the last particle within the given boundary.
     */
    [[nodiscard]] proximity_iterator<const Particle, const Cell> boundaryEnd(
        const std::set<BoundaryLocation>& locations = {BoundaryLocation::UPPER, BoundaryLocation::LOWER,
                                                       BoundaryLocation::FRONT, BoundaryLocation::BACK,
                                                       BoundaryLocation::LEFT, BoundaryLocation::RIGHT}) const;
};
static_assert(ParticleContainer<ContainerRef>);

}  // namespace mol_sim

#endif
