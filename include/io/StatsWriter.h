#ifndef STATS_WRITER_H
#define STATS_WRITER_H

#include <cstddef>
#include <string>
#include <vector>

#include "particles/container/ContainerRef.h"

namespace mol_sim {

/**
 * @brief Writes the diffusion and RDF of the given set of particles into a corresponding .csv file
 */
class StatsWriter {
    /**
     * @brief Flag whether to compute the RDF or not.
     */
    bool compute_rdf = false;

    /**
     * @brief Flag whether to compute the Diffusion and Temperature or not.
     */
    bool compute_diff_temp = false;

    /**
     * @brief Sample radius for the RDF computation.
     */
    double sample_radius = 1.0;

    /**
     * @brief Maximal radius that is considered for the RDF computation.
     */
    double window_size = 10.0;

    /**
     * @brief Filename for diffusion output.
     */
    std::string filename_diffusion = "diffusion.csv";
    /**
     * @brief Filename for RDF output.
     */
    std::string filename_rdf = "rdf.csv";
    /**
     * @brief Filename for temperature output.
     */
    std::string filename_temp = "temp.csv";

   public:
    StatsWriter() = default;

    StatsWriter(bool compute_rdf, bool compute_diff_temp, double sample_radius, double window_size);

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
     */
    void computeRDF(ContainerRef particles, std::vector<double>& results) const;

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
     */
    void plotRDF(ContainerRef particles, int iteration) const;
    /**
     * @brief Writes the temperature data into the temp.csv file. Each row has the format: <iteration>,<temperature>
     * Mainly used for debugging purposes.
     *
     * @param total_energy Total energy of the system.
     * @param dimension Dimension of the system.
     * @param N Number of particles.
     * @param iteration Current iteration of the simulation.
     */
    void plotTemp(double total_energy, size_t dimension, size_t N, int iteration) const;
};

}  // namespace mol_sim

#endif