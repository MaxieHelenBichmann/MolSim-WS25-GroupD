#pragma once

#include "io/FileReader.h"
#include "particles/ParticleContainer.h"

/**
 * @brief Reads the particles from a file that contains the initial position, force, and mass of particles.
 *
 * Reads the particles from a file that is formatted like the following regex:
 *
 * (#comment\n || empty line\n)*
 * (number of particles\n)
 * ((x y z)(whitespace)(x y z)(whitespace)(m)\n)*
 *
 * where the first (x y z) tuple represents the initial position of the particle
 * the second (x y z) tuple represents the initial velocity of the particle
 * and m represents the mass of the particle.
 * All these values are floats.
 */
class XVMReader : public FileReader {
 public:
  XVMReader();
  ~XVMReader() override;
  void readFile(particle_containers::ParticleContainer& particles, char* filename) override;
};