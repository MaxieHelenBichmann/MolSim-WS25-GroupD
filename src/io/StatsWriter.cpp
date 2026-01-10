#include "io/StatsWriter.h"

#include <fstream>

using namespace mol_sim;

double StatsWriter::computeDiffusion([[maybe_unused]] ContainerRef particles) const { return 0.0; }

double StatsWriter::computeRDF([[maybe_unused]] ContainerRef particles) const { return 0.0; }

void StatsWriter::plotDiffusion(ContainerRef particles, int iteration) const {
    if (compute_diffusion) {
        double diffusion = computeDiffusion(particles);

        std::ofstream diffusion_file(filename_diffusion, std::ios::app);
        if (diffusion_file.is_open()) {
            diffusion_file << iteration << "," << diffusion << "\n";
            diffusion_file.close();
        } else {
            SPDLOG_ERROR("Failed to open diffusion.csv for writing.");
        }
    }
}

void StatsWriter::plotRDF([[maybe_unused]] ContainerRef particles, [[maybe_unused]] int iteration) const {
    if (compute_rdf) {
        double rdf = computeRDF(particles);

        std::ofstream rdf_file(filename_rdf, std::ios::app);
        if (rdf_file.is_open()) {
            rdf_file << iteration << "," << rdf << "\n";
            rdf_file.close();
        } else {
            SPDLOG_ERROR("Failed to open rdf.csv for writing.");
        }
    }
}