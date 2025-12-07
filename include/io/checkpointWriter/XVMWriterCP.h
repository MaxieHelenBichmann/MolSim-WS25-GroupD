#ifndef XVM_WRITER_CP_H
#define XVM_WRITER_CP_H

#include "io/CheckpointWriter.h"
#include "particles/container/ContainerRef.h"

namespace mol_sim {

/**
 * @brief Writes all particle data into a .txt file that can be read by the XVMReader, but misses other state of a
simulation (as the XVM format does not support this). Also misses certain components of the particle data, like
Forces, Sigma, Epsilon or type.
 * Basically very useless and only included for backward compatibility with the XVMReader.
 *
 * The filename will consist of cp_<iteration>.txt
 */
class XVMWriterCP : public CheckpointWriter {
   public:
    XVMWriterCP();

    ~XVMWriterCP() override;
    void createCheckpoint(const Domain& domain, ContainerRef particles, int iteration, Force force, double delta_t,
                          double start_time, double end_time, size_t frequency_output, size_t frequency_checkpoint,
                          const std::string& base_name, double cutoff_radius, size_t N) const override;
};

}  // namespace mol_sim

#endif
