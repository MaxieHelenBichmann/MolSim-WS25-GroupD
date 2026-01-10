#ifndef STATS_WRITER_H
#define STATS_WRITER_H

#include <string>
#include <vector>

#include "particles/container/ContainerRef.h"

namespace mol_sim {

/**
 * @brief Writes the diffusion and RDF of the given set of particles into a .csv file
 *
 * Writes the diffusion of the particles into the file diffusion.csv in the format: <iteration>,<diffusion>
 * Writes the Radial Distribution Function (RDF) into the file rdf.csv in the format: <iteration>,<distance>,<amount>
 */
class StatsWriter {
    bool compute_rdf = false;
    bool compute_diffusion = false;
    double sample_radius = 1.0;
    double window_size = 10.0;

    std::string filename_diffusion = "diffusion.csv";
    std::string filename_rdf = "rdf.csv";

    [[nodiscard]] double computeDiffusion(ContainerRef particles);
    void computeRDF(ContainerRef particles, std::vector<size_t>& results) const;

   public:
    StatsWriter() = default;

    StatsWriter(bool compute_rdf, bool compute_diffusion, double sample_radius, double window_size)
        : compute_rdf(compute_rdf),
          compute_diffusion(compute_diffusion),
          sample_radius(sample_radius),
          window_size(window_size) {};

    ~StatsWriter();

    void plotDiffusion(ContainerRef particles, int iteration);
    void plotRDF(ContainerRef particles, int iteration) const;
};

}  // namespace mol_sim

#endif