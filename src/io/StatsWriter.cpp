#include "io/StatsWriter.h"

#include <cstddef>
#include <fstream>

using namespace mol_sim;

double StatsWriter::computeDiffusion([[maybe_unused]] ContainerRef particles) {
    double result = 0.0;
    for (auto& p : particles) {
        result += (p.getX() - p.getRefX()).sqrEuclidNorm();
        p.getRefX() = p.getX();
    }

    return result / static_cast<double>(particles.size());
}

void StatsWriter::computeRDF([[maybe_unused]] ContainerRef particles, std::vector<size_t>& results) const {
    for (auto& p1 : particles) {
        for (auto& p2 : particles) {
            if (&p1 != &p2) {
                auto index = static_cast<size_t>((p1.getX() - p2.getX()).euclidNorm() / sample_radius);
                if (index < results.size()) {
                    results[index]++;
                }
            }
        }
    }
}

void StatsWriter::plotDiffusion(ContainerRef particles, int iteration) {
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
        std::vector<size_t> results;
        results.reserve(static_cast<size_t>(std::ceil(window_size / sample_radius)));

        computeRDF(particles, results);

        std::ofstream rdf_file(filename_rdf, std::ios::app);
        if (rdf_file.is_open()) {
            for (size_t i = 0; i < results.size(); ++i) {
                rdf_file << iteration << "," << (static_cast<double>(i) * sample_radius) << "," << results[i] << "\n";
            }
            rdf_file.close();
        } else {
            SPDLOG_ERROR("Failed to open rdf.csv for writing.");
        }
    }
}