#include "io/fileReader/YAMLReader.h"

#include <spdlog/spdlog.h>
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

void YAMLReader::readFile(ContainerRef particles, SettingsParam& settings, const std::string& filename) {
    try {
        YAML::Node root = YAML::LoadFile(filename);

        if (root.size() == 0) {
            throw YAMLReaderException("Empty YAML file");
        }

        // First block must be settings
        auto first_it = root.begin();
        YAML::Node first_node = first_it->second;

        if (!first_node["format"] || first_node["format"].as<std::string>() != "Settings") {
            throw YAMLReaderException("First block must be 'Settings' format");
        }

        readSettings(settings, first_node);
        settings.setDefaults();

        // Process remaining blocks (particle definitions)
        bool has_particle_definition = false;
        for (auto it = ++root.begin(); it != root.end(); ++it) {
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
            } else if (format == "Settings") {
                SPDLOG_WARN("Additional Settings block '{}' ignored (only first is used)", name);
            } else {
                SPDLOG_WARN("Unknown format '{}' for entry '{}'", format, name);
            }
        }

        if (!has_particle_definition) {
            throw YAMLReaderException("No particle definitions (XVM, Cuboid, or Disc) found in the input file");
        }
    } catch (const YAML::Exception& e) {
        SPDLOG_ERROR("Error parsing YAML: {}", e.what());
        throw YAMLReaderException(e.what());
    }
}

void YAMLReader::readSettings(SettingsParam& settings, const YAML::Node& node) {
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
    if (node["frequency"] && !settings.frequency.has_value()) {
        settings.frequency = node["frequency"].as<size_t>();
    }
    if (node["cutoff"] && !settings.cutoff.has_value()) {
        settings.cutoff = node["cutoff"].as<double>();
    }
    if (node["domain"]) {
        parseDomain(settings, node["domain"]);
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

// Add these new methods that parse without generating
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

// Now refactor readCube to use parseCuboids
void YAMLReader::readCube(ContainerRef particles, const YAML::Node& node) {
    auto cuboids = parseCuboids(node);
    for (const auto& data : cuboids) {
        CuboidGenerator generator(data.position, data.velocity, data.num_particles, data.mass, data.distance,
                                  data.avg_velo, data.epsilon, data.sigma);
        generator.generateParticles(particles);
    }
}

// Refactor readDisc similarly
void YAMLReader::readDisc(ContainerRef particles, const YAML::Node& node) {
    auto discs = parseDiscs(node);
    for (const auto& data : discs) {
        DiscGenerator generator(data.position, data.velocity, data.radius, data.mass, data.distance, data.avg_velo,
                                data.epsilon, data.sigma);
        generator.generateParticles(particles);
    }
}

double YAMLReader::getBoundaryPosition(BoundaryLocation location, const R3& dimension) {
    switch (location) {
        case BoundaryLocation::LEFT:
            return 0.0;
        case BoundaryLocation::RIGHT:
            return dimension[0];
        case BoundaryLocation::LOWER:
            return 0.0;
        case BoundaryLocation::UPPER:
            return dimension[1];
        case BoundaryLocation::FRONT:
            return 0.0;
        case BoundaryLocation::BACK:
            return dimension[2];
        default:
            return 0.0;
    }
}

std::unique_ptr<Boundary> YAMLReader::parseBoundary(BoundaryLocation location, const R3& dimension,
                                                    const YAML::Node& node) {
    auto type_str = node["type"].as<std::string>("OUTFLOW");
    BoundaryType boundary_type = mol_sim::parseBoundaryType(type_str);

    switch (boundary_type) {
        case BoundaryType::REFLECTING: {
            std::optional<double> sigma =
                node["sigma"] ? std::optional<double>(node["sigma"].as<double>()) : std::nullopt;
            std::optional<double> epsilon =
                node["epsilon"] ? std::optional<double>(node["epsilon"].as<double>()) : std::nullopt;
            double position = getBoundaryPosition(location, dimension);
            return std::make_unique<Reflecting>(location, position, sigma, epsilon);
        }
        case BoundaryType::OUTFLOW:
        default:
            return std::make_unique<Outflow>(location);
    }
}

void YAMLReader::parseDomain(SettingsParam& settings, const YAML::Node& node) {
    R3 dimension = {node["x"].as<double>(1.0), node["y"].as<double>(1.0), node["z"].as<double>(1.0)};

    auto domain_type = node["domain_type"].as<std::string>(LINKED);
    settings.domain_type = domain_type;

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
            boundaries[i] = parseBoundary(location, dimension, bounds_node[key]);
        } else {
            boundaries[i] = std::make_unique<Outflow>(location);
        }
    }

    settings.domain.emplace(dimension, std::move(boundaries));
}

}  // namespace mol_sim
