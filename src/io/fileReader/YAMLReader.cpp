#include "io/fileReader/YAMLReader.h"

#include <spdlog/spdlog.h>
#include <yaml-cpp/exceptions.h>
#include <yaml-cpp/node/node.h>

#include <array>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <string>
#include <utility>

#include "io/fileReader/YAMLReaderException.h"
#include "particles/boundaries/Boundary.h"
#include "particles/boundaries/Outflow.h"
#include "particles/boundaries/Reflecting.h"
#include "particles/container/domain/Domain.h"
#include "particles/generators/CuboidGenerator.h"
#include "particles/generators/DiscGenerator.h"
#include "physics/ForceSource.h"
#include "utils/Settings.h"

namespace mol_sim {
YAMLReader::YAMLReader() = default;

YAMLReader::~YAMLReader() = default;

void YAMLReader::readSettings(SettingsParam& settings, const std::string& filename) {
    try {
        YAML::Node root = YAML::LoadFile(filename);

        if (root.size() == 0) {
            throw YAMLReaderException("Empty YAML file");
        }

        // First block must be settings
        auto first_it = root.begin();
        YAML::Node node = first_it->second;

        if (!node["format"] || node["format"].as<std::string>() != "Settings") {
            throw YAMLReaderException("First block must be 'Settings' format");
        }

        if (node["delta_t"] && !settings.delta_t.has_value()) {
            settings.delta_t = node["delta_t"].as<double>();
        }
        if (node["end_time"] && !settings.end_time.has_value()) {
            settings.end_time = node["end_time"].as<double>();
        }
        if (node["start_time"] && !settings.start_time.has_value()) {
            settings.start_time = node["start_time"].as<double>();
        }
        if (node["base_name"] && !settings.base_name.has_value()) {
            settings.base_name = node["base_name"].as<std::string>();
        }
        if (node["force"] && !settings.force.has_value()) {
            auto force_str = node["force"].as<std::string>();
            if (force_str == "Lennard Jones") {
                settings.force = LENNARDJONES;
            } else if (force_str == "Gravitational") {
                settings.force = GRAVITATIONAL;
            } else {
                SPDLOG_WARN("Unknown Force Type, defaulting to Lennard Jones!");
            }
        }
        if (node["container"] && !settings.container_type.has_value()) {
            settings.container_type = node["container"].as<std::string>();
        }
        if (node["frequency"] && !settings.frequency.has_value()) {
            settings.frequency = node["frequency"].as<size_t>();
        }
        if (node["cutoff"] && !settings.cutoff.has_value()) {
            settings.cutoff = node["cutoff"].as<double>();
        }
        if (node["domain"]) {
            parseDomain(settings, node["domain"]);
        }
    } catch (const YAML::Exception& e) {
        SPDLOG_ERROR("Error parsing YAML settings: {}", e.what());
        throw YAMLReaderException(e.what());
    }
}

void YAMLReader::readParticles(ContainerRef particles, const std::string& filename) {
    try {
        YAML::Node root = YAML::LoadFile(filename);

        if (root.size() == 0) {
            throw YAMLReaderException("Empty YAML file");
        }

        // Process all blocks looking for particle definitions
        bool has_particle_definition = false;
        for (auto it = root.begin(); it != root.end(); ++it) {
            const auto& name = it->first.as<std::string>();
            YAML::Node node = it->second;

            if (!node["format"]) {
                SPDLOG_WARN("Block '{}' has no format field, skipping", name);
                continue;
            }

            auto format = node["format"].as<std::string>();
            if (format == "XVM") {
                has_particle_definition = true;
                readXVM(particles, node);
            } else if (format == "Cuboid") {
                has_particle_definition = true;
                readCube(particles, node);
            } else if (format == "Disc") {
                has_particle_definition = true;
                readDisc(particles, node);
            }
            // Skip Settings and unknown formats silently in phase 2
        }

        if (!has_particle_definition) {
            throw YAMLReaderException("No particle definitions (XVM, Cuboid, or Disc) found in the input file");
        }
    } catch (const YAML::Exception& e) {
        SPDLOG_ERROR("Error parsing YAML particles: {}", e.what());
        throw YAMLReaderException(e.what());
    }
}

void YAMLReader::readXVM(ContainerRef particles, const YAML::Node& node) {
    try {
        if (node["particles"] && node["particles"].IsSequence()) {
            auto num_part = node["num_particles"].as<size_t>();
            particles.reserve(num_part);
            for (const auto& curr : node["particles"]) {
                R3 position;
                const YAML::Node& coordinates = curr["coordinates"];
                position[0] = coordinates["x"].as<double>();
                position[1] = coordinates["y"].as<double>();
                position[2] = coordinates["z"].as<double>();

                R3 velocity;
                const YAML::Node& velocity_node = curr["velocity"];
                velocity[0] = velocity_node["vx"].as<double>();
                velocity[1] = velocity_node["vy"].as<double>();
                velocity[2] = velocity_node["vz"].as<double>();

                auto mass = curr["mass"].as<double>();

                double epsilon = SettingsParam::EPSILON_DEFAULT;
                double sigma = SettingsParam::SIGMA_DEFAULT;
                if (curr["epsilon"]) {
                    epsilon = curr["epsilon"].as<double>();
                }
                if (curr["sigma"]) {
                    sigma = curr["sigma"].as<double>();
                }
                particles.addParticle(position, velocity, mass, epsilon, sigma);
            }
        }
    } catch (const YAML::Exception& e) {
        SPDLOG_ERROR("Error parsing YAML: {}", e.what());
        throw YAMLReaderException(e.what());
    }
}

std::vector<YAMLReader::CuboidData> YAMLReader::parseCuboids(const YAML::Node& node) {
    std::vector<CuboidData> cuboids;
    try {
        if (node["cuboids"] && node["cuboids"].IsSequence()) {
            for (const auto& curr : node["cuboids"]) {
                CuboidData data;

                const YAML::Node& coordinates = curr["coordinates"];
                data.position[0] = coordinates["x"].as<double>();
                data.position[1] = coordinates["y"].as<double>();
                data.position[2] = coordinates["z"].as<double>();

                const YAML::Node& velocity_node = curr["velocity"];
                data.velocity[0] = velocity_node["vx"].as<double>();
                data.velocity[1] = velocity_node["vy"].as<double>();
                data.velocity[2] = velocity_node["vz"].as<double>();

                const YAML::Node& count_node = curr["particleNum"];
                data.num_particles[0] = count_node["nx"].as<size_t>();
                data.num_particles[1] = count_node["ny"].as<size_t>();
                data.num_particles[2] = count_node["nz"].as<size_t>();

                data.mass = curr["mass"].as<double>();
                data.distance = curr["distance"].as<double>();
                data.avg_velo = curr["mean_velo"].as<double>();
                data.epsilon = curr["epsilon"].as<double>();
                data.sigma = curr["sigma"].as<double>();

                cuboids.push_back(data);
            }
        }
    } catch (const YAML::Exception& e) {
        SPDLOG_ERROR("Error parsing cuboids: {}", e.what());
        throw YAMLReaderException(e.what());
    }
    return cuboids;
}

std::vector<YAMLReader::DiscData> YAMLReader::parseDiscs(const YAML::Node& node) {
    std::vector<DiscData> discs;
    try {
        if (node["discs"] && node["discs"].IsSequence()) {
            for (const auto& curr : node["discs"]) {
                DiscData data;

                const YAML::Node& coordinates = curr["coordinates"];
                data.position[0] = coordinates["x"].as<double>();
                data.position[1] = coordinates["y"].as<double>();
                data.position[2] = coordinates["z"].as<double>();

                const YAML::Node& velocity_node = curr["velocity"];
                data.velocity[0] = velocity_node["vx"].as<double>();
                data.velocity[1] = velocity_node["vy"].as<double>();
                data.velocity[2] = velocity_node["vz"].as<double>();

                data.radius = curr["radius"].as<size_t>();
                data.mass = curr["mass"].as<double>();
                data.distance = curr["distance"].as<double>();
                data.avg_velo = curr["mean_velo"].as<double>();
                data.epsilon = curr["epsilon"].as<double>();
                data.sigma = curr["sigma"].as<double>();

                discs.push_back(data);
            }
        }
    } catch (const YAML::Exception& e) {
        SPDLOG_ERROR("Error parsing discs: {}", e.what());
        throw YAMLReaderException(e.what());
    }
    return discs;
}

void YAMLReader::readCube(ContainerRef particles, const YAML::Node& node) {
    auto cuboids = parseCuboids(node);
    for (const auto& data : cuboids) {
        CuboidGenerator generator(data.position, data.velocity, data.num_particles, data.mass, data.distance,
                                  data.avg_velo, data.epsilon, data.sigma);
        generator.generateParticles(particles);
    }
}

void YAMLReader::readDisc(ContainerRef particles, const YAML::Node& node) {
    auto discs = parseDiscs(node);
    for (const auto& data : discs) {
        DiscGenerator generator(data.position, data.velocity, data.radius, data.mass, data.distance, data.avg_velo,
                                data.epsilon, data.sigma);
        generator.generateParticles(particles);
    }
}

void YAMLReader::parseDomain(SettingsParam& settings, const YAML::Node& node) {
    try {
        R3 dimension = {node["x"].as<double>(), node["y"].as<double>(), node["z"].as<double>()};

        // Define boundary locations and their YAML keys
        static const std::array<std::pair<BoundaryLocation, std::string>, 6> boundary_mappings = {{
            {BoundaryLocation::LEFT, "left"},
            {BoundaryLocation::RIGHT, "right"},
            {BoundaryLocation::FRONT, "front"},
            {BoundaryLocation::BACK, "back"},
            {BoundaryLocation::UPPER, "upper"},
            {BoundaryLocation::LOWER, "lower"},
        }};

        std::array<std::unique_ptr<Boundary>, 6> boundaries;
        const YAML::Node& bounds_node = node["boundaries"];

        for (size_t i = 0; i < boundary_mappings.size(); ++i) {
            const auto& [location, key] = boundary_mappings[i];
            if (bounds_node && bounds_node[key]) {
                std::unique_ptr<Boundary> boundary;
                const YAML::Node& curr_node = bounds_node[key];
                auto type_str = curr_node["type"].as<std::string>("OUTFLOW");
                BoundaryType boundary_type = mol_sim::parseBoundaryType(type_str);

                switch (boundary_type) {
                    case BoundaryType::REFLECTING: {
                        std::optional<double> sigma =
                            curr_node["sigma"] ? std::optional<double>(curr_node["sigma"].as<double>()) : std::nullopt;
                        std::optional<double> epsilon = curr_node["epsilon"]
                                                            ? std::optional<double>(curr_node["epsilon"].as<double>())
                                                            : std::nullopt;

                        boundary = std::make_unique<Reflecting>(location, dimension, sigma, epsilon);
                    }
                    case BoundaryType::OUTFLOW:
                    default:
                        boundary = std::make_unique<Outflow>(location);
                }
                boundaries[i] = std::move(boundary);
            } else {
                boundaries[i] = std::make_unique<Outflow>(location);
            }
        }
        settings.domain.emplace(dimension, std::move(boundaries));
    } catch (YAML::Exception& e) {
        SPDLOG_ERROR("Error parsing domain: {}", e.what());
        throw YAMLReaderException(e.what());
    }
}

}  // namespace mol_sim
