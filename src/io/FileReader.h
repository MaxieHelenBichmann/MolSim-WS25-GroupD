/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include "particles/Particle.h"
#include "particles/ParticleContainer.h"

/**
 * @brief Interface for file readers.
 *
 * Interface for file readers.
 * These classes all implement the readFile method.
 */
class FileReader {
 public:
  virtual ~FileReader() = default;
  /**
   * @brief Reads the particles encoded in an appropriately formatted file into a ParticleContainer.
   *
   * @param particles The ParticleContainer the particles contained in the file will be stored in.
   * @param filename The path to the file containing the particles to be stored in 'particles'.
   */
  virtual void readFile(particle_containers::ParticleContainer& particles, char* filename) = 0;
};
