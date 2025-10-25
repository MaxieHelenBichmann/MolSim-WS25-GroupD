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

using R3 = Vector<double, 3>;

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
   * Type of the particle. Use it for whatever you want (e.g. to separate
   * molecules belonging to different bodies, matters, and so on)
   */
  int type;

 public:
  explicit Particle(int type = 0);

  Particle(const Particle& other);

  Particle(
      // for visualization, we need always 3 coordinates
      // -> in case of 2d, we use only the first and the second
      R3 x_arg, R3 v_arg, double m_arg, int type = 0);

  virtual ~Particle();

  [[nodiscard]] const R3& getX() const;
  R3& getX();

  [[nodiscard]] const R3& getV() const;
  R3& getV();

  [[nodiscard]] const R3& getF() const;
  R3& getF();

  [[nodiscard]] const R3& getOldF() const;
  R3& getOldF();

  [[nodiscard]] double getM() const;

  [[nodiscard]] int getType() const;

  bool operator==(const Particle& other) const;

  [[nodiscard]] std::string toString() const;
};

std::ostream& operator<<(std::ostream& stream, Particle& p);

}  // namespace mol_sim

#endif