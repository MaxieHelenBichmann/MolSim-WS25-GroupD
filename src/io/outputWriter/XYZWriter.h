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
#include "particles/container/SimpleContainer.h"
#include "io/OutputWriter.h"

namespace outputWriter {

/**
 * @brief Writes the (x,y,z) position of the given set of particles into a .xyz file
 *
 * Writes the (x,y,z) position of the given set of particles into a .xyz file. The filename
 * will consist of the provided filename + the current iteration.
 */
class XYZWriter : public OutputWriter {
 public:
  XYZWriter();
  ~XYZWriter() override;
  void plotParticles(particle_containers::ParticleContainer& particles, const std::string& filename, int iteration) override;
};

}  // namespace outputWriter
