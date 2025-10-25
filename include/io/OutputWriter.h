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
  virtual void plotParticles(ContainerRef particles, const std::string& filename, int iteration) = 0;
};
}  // namespace mol_sim
#endif
