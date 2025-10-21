/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include <list>

#include "particles/Particle.h"

class FileReader {
 public:
  FileReader();
  virtual ~FileReader();

  void readFile(std::list<Particle> &particles, char *filename);
};
