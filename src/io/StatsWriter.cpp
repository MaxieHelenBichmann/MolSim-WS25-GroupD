#include "io/StatsWriter.h"

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <numbers>

using namespace mol_sim;

StatsWriter::StatsWriter(bool compute_rdf, bool compute_diff_temp, double sample_radius, double window_size)
    : compute_rdf(compute_rdf),
      compute_diff_temp(compute_diff_temp),
      sample_radius(sample_radius),
      window_size(window_size) {
    auto remove_if_exists = [](const std::string& filename) {
        std::error_code ec;
        if (std::filesystem::exists(filename, ec) && !ec) {
            std::filesystem::remove(filename, ec);
            if (ec) {
                SPDLOG_ERROR("Failed to remove {}: {}", filename, ec.message());
            }
        } else if (ec) {
            SPDLOG_ERROR("Failed to check existence of {}: {}", filename, ec.message());
        }
    };

    if (this->compute_diff_temp) {
        remove_if_exists(filename_diffusion);
        remove_if_exists(filename_temp);
    }
    if (this->compute_rdf) {
        remove_if_exists(filename_rdf);
    }
}

double StatsWriter::computeDiffusion(ContainerRef particles) const {
    double result = 0.0;
    for (auto& p : particles) {
        result += (p.getX() - p.getRefX()).sqrEuclidNorm();
        p.getRefX() = p.getX();
    }

    return !particles.empty() ? result / static_cast<double>(particles.size()) : 0.0;
}

void StatsWriter::computeRDF(ContainerRef particles, std::vector<double>& results) const {  // NOLINT
    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            auto index = static_cast<size_t>((particles[i].getX() - particles[j].getX()).euclidNorm() / sample_radius);
            for (const auto& mirrored_pos : particles[j].getMirrorPositions()) {
                size_t idx = static_cast<size_t>((particles[i].getX() - mirrored_pos).euclidNorm() / sample_radius);
                index = idx < index ? idx : index;
            }
            if (index < results.size()) {
                results[index]++;
            }
        }
    }
    for (size_t i = 0; i < results.size(); ++i) {
        const double end_interval = (static_cast<double>(i + 1)) * sample_radius;
        const double start_interval = static_cast<double>(i) * sample_radius;
        double vol = (end_interval * end_interval * end_interval) - (start_interval * start_interval * start_interval);
        results[i] = 0.75 * results[i] / (vol * std::numbers::pi);
    }
}

void StatsWriter::plotDiffusion(ContainerRef particles, int iteration) const {
    if (compute_diff_temp) {
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

void StatsWriter::plotRDF(ContainerRef particles, int iteration) const {
    if (compute_rdf) {
        const auto bin_count = static_cast<size_t>(std::ceil(window_size / sample_radius));
        std::vector<double> results(bin_count, 0.0);

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

void StatsWriter::plotTemp(double total_energy, size_t dimension, size_t N, int iteration) const {
    if (compute_diff_temp) {
        double temperature = 0.0;
        if (N != 0 && dimension != 0) {
            temperature = (2.0 * total_energy) / (static_cast<double>(dimension) * static_cast<double>(N));
        }
        std::ofstream temp_file(filename_temp, std::ios::app);
        if (temp_file.is_open()) {
            temp_file << iteration << "," << temperature << "\n";
            temp_file.close();
        } else {
            SPDLOG_ERROR("Failed to open temp.csv for writing.");
        }
    }
}