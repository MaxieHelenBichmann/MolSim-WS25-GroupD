/*
 * VTKWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#ifndef VTK_WRITER_H
#define VTK_WRITER_H
/**
 * @brief Writes the (x,y,z) position of the given set of particles into a .vtu file
 *
 * Writes the (x,y,z) position of the given set of particles into a .vtu file.The filename
 * will consist of the provided filename + the current iteration.
 */
#include "io/OutputWriter.h"
#ifdef ENABLE_VTK_OUTPUT

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <string>

#include "particles/container/ContainerRef.h"

namespace mol_sim {

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
    void plotParticles(ContainerRef particles, const std::string& filename, int iteration) override;
};

}  // namespace mol_sim

#endif
#endif
