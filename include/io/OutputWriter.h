#ifndef OUTPUT_WRITER_H
#define OUTPUT_WRITER_H

#include "particles/container/ContainerRef.h"

namespace mol_sim {

/**
 * @brief Interface for output writers.
 *
 * Interface for output writers.
 * These classes all implement the plotParticles method.
 */
class OutputWriter {
   public:
    virtual ~OutputWriter() = default;

    /**
     * @brief Plots the particles to a file to be able to visualize the simulation state.
     *
     * @param particles The particles to plot.
     * @param filename The base filename to use.
     * @param iteration The current iteration of the simulation.
     */
    virtual void plotParticles(ContainerRef particles, const std::string& filename, int iteration) const = 0;
};

}  // namespace mol_sim

#endif
