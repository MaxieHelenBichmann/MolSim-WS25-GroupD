/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include <array>
#include <string>

#include "../utils/Vector.h"

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

  const R3& getX() const;
  R3& getX();

  const R3& getV() const;
  R3& getV();

  const R3& getF() const;
  R3& getF();

  const R3& getOldF() const;
  R3& getOldF();

  double getM() const;

  int getType() const;

  bool operator==(const Particle& other) const;

  std::string toString() const;
};

std::ostream& operator<<(std::ostream& stream, Particle& p);
