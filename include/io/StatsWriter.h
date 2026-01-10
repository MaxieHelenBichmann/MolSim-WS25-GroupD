#ifndef STATS_WRITER_H
#define STATS_WRITER_H

#include "particles/container/ContainerRef.h"

namespace mol_sim {

/**
 * @brief Writes the diffusion and RDF of the given set of particles into a .csv file
 *
 * Writes the diffusion of the particles into the file diffusion.csv in the format: <iteration>,<diffusion>
 * Writes the Radial Distribution Function (RDF) into the file rdf.csv in the format: <distance>,<rdf_value>
 */
class StatsWriter {
    bool compute_rdf = false;
    bool compute_diffusion = false;
    double sample_radius = 1.0;

    [[nodiscard]] double computeDiffusion(ContainerRef particles) const;
    [[nodiscard]] double computeRDF(ContainerRef particles, std::vector<double>& rdf_bins, double bin_width,
                                    size_t n_bins) const;

   public:
    StatsWriter() = default;
    ~StatsWriter();

    void initStats(bool compute_rdf, bool compute_diffusion, double sample_radius) {
        this->compute_rdf = compute_rdf;
        this->compute_diffusion = compute_diffusion;
        this->sample_radius = sample_radius;
    };
    void plotStatistics(ContainerRef particles, int iteration) const;
};

}  // namespace mol_sim

#endif