/*
 * XYZWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#ifndef XYZ_WRITER_H
#define XYZ_WRITER_H

#include "io/OutputWriter.h"
#include "particles/container/ContainerRef.h"

namespace mol_sim {

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

    void plotParticles(ContainerRef particles, const std::string& filename, int iteration) const override;
};

}  // namespace mol_sim

#endif
