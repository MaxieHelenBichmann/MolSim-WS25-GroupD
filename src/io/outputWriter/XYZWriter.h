/*
 * XYZWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#pragma once

#include <fstream>
#include <list>

#include "particles/Particle.h"
#include "particles/container/ContainerRef.h"

namespace outputWriter {

class XYZWriter {
 public:
  XYZWriter();

  virtual ~XYZWriter();

  static void plotParticles(particle_containers::ContainerRef particles, const std::string& filename, int iteration);
};

}  // namespace outputWriter
