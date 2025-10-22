#pragma once

#include "particles/container/SimpleContainer.h"

/**
 * @brief Interface for output writers.
 *
 * Interface for output writers.
 * These classes all implement the plotParticles method.
 */
class OutputWriter {
 public:
  virtual ~OutputWriter() = default;
  virtual void plotParticles(particle_containers::ParticleContainer& particles, const std::string& filename, int iteration) = 0;
};