/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#ifndef PARTICLE_H
#define PARTICLE_H

#include <string>

#include "utils/Vector.h"

namespace mol_sim {

/**
 * @brief Stores all important attributes of a Particle.
 *
 */
class Particle {
   private:
    /**
     * Position of the particle
     */
    R3 x;

    /**
     * Velocity of the particle
     */
    R3 v;

    /**
     * Force effective on this particle
     */
    R3 f;

    /**
     * Force which was effective on this particle
     */
    R3 old_f;

    /**
     * Mass of this particle
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
     * Type of the particle. Use it for whatever you want (e.g. to separate
     * molecules belonging to different bodies, matters, and so on)
     */
    int type;

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

    virtual ~Particle();

    /**
     * @brief Access the current coordinates of a Particle.
     *
     * @return Const reference to the coordinates of the Particle.
     */
    [[nodiscard]] const R3& getX() const;
    /**
     * @brief Access the current coordinates of a Particle.
     *
     * @return Reference to the coordinates of the Particle.
     */
    R3& getX();

    /**
     * @brief Access the current velocity Vector of a Particle.
     *
     * @return Const reference to the velocity of the Particle.
     */
    [[nodiscard]] const R3& getV() const;
    /**
     * @brief Access the current velocity Vector of a Particle.
     *
     * @return Reference to the velocity of the Particle.
     */
    R3& getV();

    /**
     * @brief Access the current force acting on a Particle.
     *
     * @return Const reference to the force on the Particle.
     */
    [[nodiscard]] const R3& getF() const;
    /**
     * @brief Access the current force acting on a Particle.
     *
     * @return Reference to the force on the Particle.
     */
    R3& getF();

    /**
     * @brief Access the previous force acting on a Particle.
     *
     * @return Const reference to the previous force on the Particle.
     */
    [[nodiscard]] const R3& getOldF() const;
    /**
     * @brief Access the previous force acting on a Particle.
     *
     * @return Reference to the previous force on the Particle.
     */
    R3& getOldF();

    /**
     * @brief Access the mass of a Particle.
     *
     * @return Double of the mass of the Particle.
     */
    [[nodiscard]] double getM() const;

    /**
     * @brief Access the type of a Particle.
     *
     * @return Int of the type of the Particle.
     */
    [[nodiscard]] int getType() const;

    /**
     * @brief Access to the sigma of this Particle
     *
     * @return double of the sigma of the Particle
     */
    [[nodiscard]] double getSigma() const;

    /**
     * @brief Access to the epsilon of this Particle
     *
     * @return double of the epsilon of the Particle
     */
    [[nodiscard]] double getEpsilon() const;

    /**
     * @brief Set the Sigma of this Particle.
     *
     * @param sigma_arg new sigma value.
     */
    void setSigma(double sigma_arg);

    /**
     * @brief Set the Epsilon of this Particle.
     *
     * @param epsilon_arg new epsilon value.
     */
    void setEpsilon(double epsilon_arg);

    bool operator==(const Particle& other) const;

    /**
     * @brief Conversion to a std::string.
     *
     * @return String representation of the Particle.
     */
    [[nodiscard]] std::string toString() const;
};

std::ostream& operator<<(std::ostream& stream, Particle& p);

}  // namespace mol_sim

#endif
