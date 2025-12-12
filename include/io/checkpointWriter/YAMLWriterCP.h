#ifndef YAML_WRITER_CP_H
#define YAML_WRITER_CP_H

#include "io/CheckpointWriter.h"
#include "particles/container/ContainerRef.h"

namespace mol_sim {

/**
 * @brief Writes the whole state of a simulation with all particle data into a .txt file
 *
 * The filename will consist of cp_<iteration>.yaml
 */
class YAMLWriterCP : public CheckpointWriter {
   public:
    YAMLWriterCP();

    ~YAMLWriterCP() override;
    void createCheckpoint(const Domain& domain, ContainerRef particles, int iteration, Force force, double delta_t,
                          double current_time, double end_time, size_t frequency_output, size_t frequency_checkpoint,
                          const std::string& base_name, double cutoff_radius, double target_temp, double delta_temp,
                          size_t thermostat_freq, size_t N) const override;
};

}  // namespace mol_sim

#endif