#include "io/checkpointWriter/YAMLWriterCP.h"

using namespace mol_sim;

YAMLWriterCP::YAMLWriterCP() = default;
YAMLWriterCP::~YAMLWriterCP() = default;

void YAMLWriterCP::createCheckpoint(const Domain& domain, ContainerRef particles, int iteration, Force force,
                                    double delta_t, double start_time, double end_time, size_t frequency_output,
                                    size_t frequency_checkpoint, const std::string& base_name, double cutoff_radius,
                                    size_t N) const {
    (void)domain;
    (void)particles;
    (void)iteration;
    (void)force;
    (void)delta_t;
    (void)start_time;
    (void)end_time;
    (void)frequency_output;
    (void)frequency_checkpoint;
    (void)base_name;
    (void)cutoff_radius;
    (void)N;
}