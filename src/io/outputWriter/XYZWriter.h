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
/**
 * @brief Writes the (x,y,z) position of the given set of particles into a .xyz file
 *
 * Writes the (x,y,z) position of the given set of particles into a .xyz file. The filename
 * will consist of the provided filename + the current iteration.
 */
class XYZWriter {
public:
 XYZWriter();

 virtual ~XYZWriter();

 static void plotParticles(particle_containers::ContainerRef particles, const std::string& filename, int iteration);
};

}  // namespace outputWriter
