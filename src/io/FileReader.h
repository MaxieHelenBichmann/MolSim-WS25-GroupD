/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include "particles/container/ContainerRef.h"

class FileReader {
 public:
  FileReader();
  virtual ~FileReader();

  static void readFile(particle_containers::ContainerRef particles, char* filename);
};
