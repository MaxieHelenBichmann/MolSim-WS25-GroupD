/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#ifndef FILE_READER_H
#define FILE_READER_H

#include "particles/container/ContainerRef.h"

namespace mol_sim {

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
  virtual void readFile(ContainerRef particles, const std::string& filename) = 0;
};

}  // namespace mol_sim

#endif
