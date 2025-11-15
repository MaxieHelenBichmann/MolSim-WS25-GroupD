#ifndef LINKEDCELL_CONTAINER_H
#define LINKEDCELL_CONTAINER_H

#include <particles/ParticleContainer.h>

#include <vector>

namespace mol_sim {

/**
 * @brief Linked-Cell Container for Particles
 *
 * This container implements the concept ParticleContainer.
 * It stored the Particles in linked cells to optimize proximity queries. However,
 * it is also possible to simply iterate over all Particles, because they are stored in a linearized vector.
 *
 */
class LinkedCellContainer {
    /**
     * @brief Cell with specific boundaries and type in the Linked-Cell container.
     *
     * Unites methods for access and modification of Particles in the cell.
     */
    struct Cell {
        enum class CellType : std::uint8_t { INNER, BOUNDARY, HALO };

        Cell(CellType cell_type, std::array<double, 6> bounds);

        void addParticle(Particle* value);
        void removeParticle(size_t idx);
        void removeParticle(Particle* p);
        void clear();

        Particle* operator[](size_t idx);
        bool fits(Particle* p) const;

       private:
        std::vector<Particle*> data;
        [[maybe_unused]] CellType type = CellType::INNER;
        std::array<double, 6> bounds;  // xmin, xmax, ymin, ymax, zmin, zmax
    };

    void switchCell(Particle& p, size_t old_cell_idx, size_t new_cell_idx);

    void switchCell(Particle* p, size_t old_cell_idx, size_t new_cell_idx);

    void switchCell(size_t p, size_t old_cell_idx, size_t new_cell_idx);

    size_t findCellIndex(Particle* p);

    /**
     * Vector storing all Particles in the container.
     */
    std::vector<Particle> data;

    /**
     * Vector storing all cells of the container.
     */
    std::vector<Cell> cells;

    R3 domain_size;
    std::array<size_t, 3> num_cells;
    double cutoff_radius;

   public:
    // constructors

    /**
     * @brief Constructor, initializing a LinkedCellContainer.
     */
    LinkedCellContainer(R3 domain_size, double cutoff_radius);

    // retrieve data
    [[nodiscard]] bool fits(R3 v) const;

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
     * @param epsilon_arg Epsilopn of the Particle.
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

    class proximity_iterator {
       public:
        proximity_iterator() = default;
    };
    // static_assert(std::forward_iterator<proximity_iterator>);

    class const_proximity_iterator {
       public:
        const_proximity_iterator() = default;
    };
    // static_assert(std::forward_iterator<const_proximity_iterator>);

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
     * @param offset Offset from the beginning of the container to start the iteration (used for N3L optimization).
     *
     * @return Const iterator to the first particle within the given radius of the center.
     */
    [[nodiscard]] const_proximity_iterator proximityBegin(R3 center, double radius, size_t offset = 0) const;

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
