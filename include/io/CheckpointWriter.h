#ifndef CHECKPOINT_WRITER_H
#define CHECKPOINT_WRITER_H

#include "particles/container/ContainerRef.h"
#include "particles/container/domain/Domain.h"
#include "physics/ForceSource.h"

namespace mol_sim {

/**
 * @brief Interface for checkpoint writers.
 *
 * Creates a checkpoint of the current simulation state.
 * This checkpoint will be named cp_<iteration>.{<file_extension>}
 */
class CheckpointWriter {
   public:
    virtual ~CheckpointWriter() = default;
    virtual void createCheckpoint(const Domain& domain, ContainerRef particles, int iteration, Force force,
                                  double delta_t, double current_time, double end_time, size_t frequency_output,
                                  size_t frequency_checkpoint, const std::string& base_name, double cutoff_radius,
                                  double target_temp, double delta_temp, size_t thermostat_freq, size_t N) const = 0;
};

}  // namespace mol_sim

#endif