/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#ifndef PARTICLE_H
#define PARTICLE_H

#include <string>
#include <cstdint>

#include "utils/Vector.h"

namespace mol_sim {

/**
 * @brief Stores all important attributes of a Particle.
 *
 */
class Particle {
   private:
    /**
     * @brief Position of the particle
     */
    R3 x;
    /**
     * @brief New position of the particle
     */
    R3 old_x;

    /**
     * @brief Velocity of the particle
     */
    R3 v;

    /**
     * @brief Force effective on this particle
     */
    R3 f;

    /**
     * @brief Force which was effective on this particle
     */
    R3 old_f;

    /**
     * @brief Mass of this particle
     */
    double m;
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
    /**
     * @brief Type of the particle. 
     * A particle with a negative type should be created and destroyed within the same 
     * simulation iteration. Positive particle types may also exist beyond one iteration.
     * -1 if particle is GHOST particle
     *  0 default
     *  1 if particle is a mirrored particle (in periodic boundaries)
     */
    int type;
    /**
     * @brief A bitmap indicating the locations the particle has been mirrored to.
     * This is relevant for Periodic boundaries.
     * 
     * nth bit = 1 means the particle has been mirrored to the nth mirror domain.
     * -------------------------------------------------------------------------- 
     * nth bit          Location
     * -------------------------------------------------------------------------- 
     * 0                (0,0,0)
     * ... TODO: finish writing this.
     * 26               (max,max,max)
     */
    uint32_t mirror_locations = 0;

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
    [[nodiscard]] const R3& getX() const noexcept;
    /**
     * @brief Access the current coordinates of a Particle.
     *
     * @return Reference to the coordinates of the Particle.
     */
    R3& getX() noexcept;
    /**
     * @brief Access the previous coordinates of a Particle.
     *
     * @return Const reference to the previous coordinates of the Particle.
     */
    [[nodiscard]] const R3& getOldX() const noexcept;
    /**
     * @brief Access the previous coordinates of a Particle.
     *
     * @return Reference to the previous coordinates of the Particle.
     */
    R3& getOldX() noexcept;

    /**
     * @brief Access the current velocity Vector of a Particle.
     *
     * @return Const reference to the velocity of the Particle.
     */
    [[nodiscard]] const R3& getV() const noexcept;
    /**
     * @brief Access the current velocity Vector of a Particle.
     *
     * @return Reference to the velocity of the Particle.
     */
    R3& getV() noexcept;

    /**
     * @brief Access the current force acting on a Particle.
     *
     * @return Const reference to the force on the Particle.
     */
    [[nodiscard]] const R3& getF() const noexcept;
    /**
     * @brief Access the current force acting on a Particle.
     *
     * @return Reference to the force on the Particle.
     */
    R3& getF() noexcept;

    /**
     * @brief Access the previous force acting on a Particle.
     *
     * @return Const reference to the previous force on the Particle.
     */
    [[nodiscard]] const R3& getOldF() const noexcept;
    /**
     * @brief Access the previous force acting on a Particle.
     *
     * @return Reference to the previous force on the Particle.
     */
    R3& getOldF() noexcept;

    /**
     * @brief Access the mass of a Particle.
     *
     * @return Double of the mass of the Particle.
     */
    [[nodiscard]] double getM() const noexcept;

    /**
     * @brief Access the type of a Particle.
     *
     * @return Int of the type of the Particle.
     */
    [[nodiscard]] int& getType() noexcept;
    /**
     * @brief Access the type of a Particle.
     *
     * @return Int of the type of the Particle.
     */
    [[nodiscard]] const int& getType() const noexcept;

    /**
     * @brief Access to the sigma of this Particle
     *
     * @return double of the sigma of the Particle
     */
    [[nodiscard]] double getSigma() const noexcept;

    /**
     * @brief Access to the sigma of this Particle
     *
     * @return Reference to the double of the sigma of the Particle
     */
    double& getSigma() noexcept;

    /**
     * @brief Access to the epsilon of this Particle
     *
     * @return double of the epsilon of the Particle
     */
    [[nodiscard]] double getEpsilon() const noexcept;

    /**
     * @brief Access to the epsilon of this Particle
     *
     * @return Reference to the double of the epsilon of the Particle
     */
    double& getEpsilon() noexcept;

    /**
     * @brief Get the Mirror Locations bitmap of the particle
     * 
     * @return uint32_t& A reference to the bitmap indicating where the particle has already been mirrored
     */
    uint32_t& getMirrorLocations() noexcept;
    
    /**
     * @brief Get the Mirror Locations bitmap of the particle
     * 
     * @return uint32_t& A const reference to the bitmap indicating where the particle has already been mirrored
     */
    [[nodiscard]] const uint32_t& getMirrorLocations() const noexcept;

    bool operator==(const Particle& other) const noexcept;

    /**
     * @brief Conversion to a std::string.
     *
     * @return String representation of the Particle.
     */
    [[nodiscard]] std::string toString() const;
};

std::ostream& operator<<(std::ostream& stream, const Particle& p);

}  // namespace mol_sim

#endif
