#ifndef STATS_WRITER_H
#define STATS_WRITER_H

#include <string>
#include <vector>

#include "particles/container/ContainerRef.h"

namespace mol_sim {

/**
 * @brief Writes the diffusion and RDF of the given set of particles into a corresponding .csv file
 */
class StatsWriter {
    bool compute_rdf = false;
    bool compute_diffusion = false;
    double sample_radius = 1.0;
    double window_size = 10.0;

    std::string filename_diffusion = "diffusion.csv";
    std::string filename_rdf = "rdf.csv";

   public:
    StatsWriter() = default;

    StatsWriter(bool compute_rdf, bool compute_diffusion, double sample_radius, double window_size);

    ~StatsWriter() = default;

    /**
     * @brief Computes the variance of the movement of particles since the last call, using the Particle member
     * reference_position, which takes care of Periodic Boundary Conditions as well.
     * Public method so it can be tested.
     *
     * @param particles ContainerRef to the particles.
     * @return double Computed diffusion value.
     */
    [[nodiscard]] double computeDiffusion(ContainerRef particles) const;
    /**
     * @brief Computes the local densities of the Radial Distribution Function (RDF) of the given particles and stores
     * the results in the provided vector.
     * Public method so it can be tested.
     *
     * @param particles ContainerRef to the particles.
     * @param results Reference to a vector where the results will be stored. The index corresponds to distance /
     * sample_radius.
     * @param mirrored Vector of mirrored particles to consider for RDF calculation.
     */
    void computeRDF(ContainerRef particles, std::vector<double>& results, const std::vector<Particle>& mirrored) const;

    /**
     * @brief Writes the diffusion data into the diffusion.csv file. Each row has the format: <iteration>,<diffusion>
     *
     * @param particles ContainerRef to the particles.
     * @param iteration Current iteration of the simulation.
     */
    void plotDiffusion(ContainerRef particles, int iteration) const;
    /**
     * @brief Writes the RDF data into the rdf.csv file. Each row has the format: <iteration>,<distance>,<density>
     *
     * @param particles ContainerRef to the particles.
     * @param iteration Current iteration of the simulation.
     * @param mirrored Vector of mirrored particles to consider for RDF calculation.
     */
    void plotRDF(ContainerRef particles, int iteration, const std::vector<Particle>& mirrored) const;
};

}  // namespace mol_sim

#endif