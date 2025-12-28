#ifndef CHECKPOINT_WRITER_H
#define CHECKPOINT_WRITER_H

#include "particles/container/ContainerRef.h"
#include "particles/container/domain/Domain.h"
#include "utils/Settings.h"

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
    virtual void createCheckpoint(SettingsParam& settings, const Domain& domain, ContainerRef particles, int iteration,
                                  size_t N) const = 0;
};

}  // namespace mol_sim

#endif