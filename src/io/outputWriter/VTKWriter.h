/*
 * VTKWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#pragma once
/**
 * @brief Writes the (x,y,z) position of the given set of particles into a .vtu file
 *
 * Writes the (x,y,z) position of the given set of particles into a .vtu file.The filename
 * will consist of the provided filename + the current iteration.
 */
#ifdef ENABLE_VTK_OUTPUT

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <list>
#include <string>

#include "particles/Particle.h"
#include "particles/container/SimpleContainer.h"
#include "io/OutputWriter.h"

namespace outputWriter {

/**
 * This class implements the functionality to generate vtk output from
 * particles using the official VTK library.
 */
class VTKWriter : public OutputWriter {
 public:
  VTKWriter() = default;
  ~VTKWriter() override = default;

  // Delete copy constructor and assignment operator
  VTKWriter(const VTKWriter&) = delete;
  VTKWriter& operator=(const VTKWriter&) = delete;

  /**
   * Write VTK output of particles.
   * @param particles Particles to add to the output
   * @param filename Output filename
   * @param iteration Current iteration number
   */
  void plotParticles(particle_containers::ParticleContainer& particles, const std::string& filename, int iteration);
};

}  // namespace outputWriter
#endif
