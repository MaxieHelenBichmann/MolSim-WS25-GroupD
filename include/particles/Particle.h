/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#ifndef PARTICLE_H
#define PARTICLE_H

#include <cstdint>
#include <string>
#include <vector>

#include "utils/Vector.h"

namespace mol_sim {

/**
 * @brief Stores all important attributes of a Particle.
 *
 */
class Particle {
   private:
    // ========== HOT DATA ==========

    /**
     * @brief Position of the particle
     */
    R3 x;

    /**
     * @brief Force effective on this particle
     */
    R3 f;

    /**
     * @brief Epsilon of this particle.
     *
     */
    double epsilon;
    /**
     * @brief sigma of this particle.
     *
     */
    double sigma;

    // ========== WARM DATA ==========

    /**
     * @brief Velocity of the particle
     */
    R3 v;

    /**
     * @brief Mass of this particle
     */
    double m;

    // ========== COLD DATA ==========

    /**
     * @brief New position of the particle
     */
    R3 old_x;

    /**
     * @brief Force which was effective on this particle
     */
    R3 old_f;

    /**
     * @brief Type of the particle.
     * A particle with a negative type should be created and destroyed within the same
     * simulation iteration. Positive particle types may also exist beyond one iteration.
     * -1 if particle is GHOST particle
     *  0 default
     *  1 if particle is a mirrored particle (in periodic boundaries)
     */
    int type;

    // NOLINTBEGIN
    /**
     * @brief A bitmap indicating the locations the particle has been mirrored to.
     * This is relevant for periodic boundaries. See also Assignment4/Periodic slides
     * for a visualization.
     *
     * n-th bit = n-th mirror location. The coordinates are *relative* to particle.getX()
     * d0 = domain_size[0], d1 = domain_size[1], d2 = domain_size[2]
     * --------------------------------------------------------------------------
     * n-th bit         Location
     * --------------------------------------------------------------------------
     * 0                (-d0, -d1, -d2)
     * 1                (-d0, -d1,   0)
     * 2                (-d0, -d1,  d2)
     * 3                (  0, -d0, -d2)
     * 4                (  0, -d1,   0)
     * 5                (  0, -d1,  d2)
     * 6                ( d0, -d1, -d2)
     * 7                ( d0, -d1,   0)
     * 8                ( d0, -d1,  d2)
     * 9                (-d0,   0, -d2)
     * 10               (-d0,   0,   0)
     * 11               (-d0,   0,  d2)
     * 12               (  0,   0, -d2)
     * 13               (  0,   0,   0)         this is where the particle is
     * 14               (  0,   0,  d2)
     * 15               ( d0,   0, -d2)
     * 16               ( d0,   0,   0)
     * 17               ( d0,   0,  d2)
     * 18               (-d0,  d1, -d2)
     * 19               (-d0,  d1,   0)
     * 20               (-d0,  d1,  d2)
     * 21               (  0,  d1, -d2)
     * 22               (  0,  d1,   0)
     * 23               (  0,  d1,  d2)
     * 24               ( d0,  d1, -d2)
     * 25               ( d0,  d1,   0)
     * 26               ( d0,  d1,  d2)
     */
    uint32_t mirror_locations = 0;

    // Positions Particle is mirrored to
    std::vector<R3> mirror_positions;
    // NOLINTEND

   public:
    explicit Particle(int type = 0);

    Particle(const Particle& other);

    Particle& operator=(const Particle& other);

    Particle(
        // for visualization, we need always 3 coordinates
        // -> in case of 2d, we use only the first and the second
        R3 x_arg, R3 v_arg, double m_arg, double epsilon_arg, double sigma_arg, int type_arg = 0);
    Particle(
        // for visualization, we need always 3 coordinates
        // -> in case of 2d, we use only the first and the second
        R3 x_arg, R3 v_arg, R3 f_arg, double m_arg, double epsilon_arg, double sigma_arg, int type_arg = 0);

    Particle(
        // for visualization, we need always 3 coordinates
        // -> in case of 2d, we use only the first and the second
        R3 x_arg, R3 old_x_arg, R3 v_arg, R3 f_arg, R3 old_f_arg, double m_arg, double epsilon_arg, double sigma_arg,
        int type_arg);

    ~Particle();

    /**
     * @brief Access the current coordinates of a Particle.
     *
     * @return Const reference to the coordinates of the Particle.
     */
    [[nodiscard]] const R3& getX() const noexcept { return x; };
    /**
     * @brief Access the current coordinates of a Particle.
     *
     * @return Reference to the coordinates of the Particle.
     */
    R3& getX() noexcept { return x; };
    /**
     * @brief Access the previous coordinates of a Particle.
     *
     * @return Const reference to the previous coordinates of the Particle.
     */
    [[nodiscard]] const R3& getOldX() const noexcept { return old_x; };
    /**
     * @brief Access the previous coordinates of a Particle.
     *
     * @return Reference to the previous coordinates of the Particle.
     */
    R3& getOldX() noexcept { return old_x; };

    /**
     * @brief Access the current velocity Vector of a Particle.
     *
     * @return Const reference to the velocity of the Particle.
     */
    [[nodiscard]] const R3& getV() const noexcept { return v; };
    /**
     * @brief Access the current velocity Vector of a Particle.
     *
     * @return Reference to the velocity of the Particle.
     */
    R3& getV() noexcept { return v; };

    /**
     * @brief Access the current force acting on a Particle.
     *
     * @return Const reference to the force on the Particle.
     */
    [[nodiscard]] const R3& getF() const noexcept { return f; };
    /**
     * @brief Access the current force acting on a Particle.
     *
     * @return Reference to the force on the Particle.
     */
    R3& getF() noexcept { return f; };

    /**
     * @brief Access the previous force acting on a Particle.
     *
     * @return Const reference to the previous force on the Particle.
     */
    [[nodiscard]] const R3& getOldF() const noexcept { return old_f; };
    /**
     * @brief Access the previous force acting on a Particle.
     *
     * @return Reference to the previous force on the Particle.
     */
    R3& getOldF() noexcept { return old_f; };

    /**
     * @brief Access the mass of a Particle.
     *
     * @return Double of the mass of the Particle.
     */
    [[nodiscard]] double getM() const noexcept { return m; };

    /**
     * @brief Access the type of a Particle.
     *
     * @return Int of the type of the Particle.
     */
    [[nodiscard]] int& getType() noexcept { return type; };
    /**
     * @brief Access the type of a Particle.
     *
     * @return Int of the type of the Particle.
     */
    [[nodiscard]] const int& getType() const noexcept { return type; };

    /**
     * @brief Access to the sigma of this Particle
     *
     * @return double of the sigma of the Particle
     */
    [[nodiscard]] double getSigma() const noexcept { return sigma; };

    /**
     * @brief Access to the sigma of this Particle
     *
     * @return Reference to the double of the sigma of the Particle
     */
    double& getSigma() noexcept { return sigma; };

    /**
     * @brief Access to the epsilon of this Particle
     *
     * @return double of the epsilon of the Particle
     */
    [[nodiscard]] double getEpsilon() const noexcept { return epsilon; };

    /**
     * @brief Access to the epsilon of this Particle
     *
     * @return Reference to the double of the epsilon of the Particle
     */
    double& getEpsilon() noexcept { return epsilon; };

    /**
     * @brief Get the Mirror Locations bitmap of the particle
     *
     * @return uint32_t& A reference to the bitmap indicating where the particle has already been mirrored
     */
    uint32_t& getMirrorLocations() noexcept { return mirror_locations; }

    /**
     * @brief Get the Mirror Locations bitmap of the particle
     *
     * @return uint32_t& A const reference to the bitmap indicating where the particle has already been mirrored
     */
    [[nodiscard]] const uint32_t& getMirrorLocations() const noexcept { return mirror_locations; }

    bool operator==(const Particle& other) const noexcept;

    /**
     * @brief Conversion to a std::string.
     *
     * @return String representation of the Particle.
     */
    [[nodiscard]] std::string toString() const;

    [[nodiscard]] std::vector<R3>& getMirrorPositions() noexcept { return mirror_positions; }
    [[nodiscard]] const std::vector<R3>& getMirrorPositions() const noexcept { return mirror_positions; }
};

std::ostream& operator<<(std::ostream& stream, const Particle& p);

}  // namespace mol_sim

#endif
