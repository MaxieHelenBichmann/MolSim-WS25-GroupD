/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include <list>

#include "particles/Particle.h"
#include "particles/ParticleContainer.h"

class FileReader {
 public:
  FileReader();
  virtual ~FileReader();

  void readFile(particle_containers::ParticleContainer& particles, char* filename);
};
