#include "io/fileReader/YAMLReader.h"

#include <spdlog/spdlog.h>
#include <yaml-cpp/node/node.h>

#include <cstddef>
#include <cstdlib>
#include <string>

#include "exceptions/YAMLReaderException.h"
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
        bool has_particle_definition = false;
        for (YAML::const_iterator it = root.begin(); it != root.end(); ++it) {
            const auto& name = it->first.as<std::string>();
            YAML::Node node = it->second;
            if (node["format"]) {
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
                    readSettings(settings, node);
                } else {
                    SPDLOG_WARN("Unknown Format '{}' for entry '{}'", format, name);
                }
            }
        }

        if (!has_particle_definition) {
            throw YAMLReaderException("No particle definitions (XVM or Cuboid) found in the input file.");
        }
    } catch (const YAML::Exception& e) {
        SPDLOG_ERROR("Error parsing YAML: {}", e.what());
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
    if (node["domain"] && !settings.domain.has_value()) {
        const YAML::Node& domain_node = node["domain"];
        auto x = domain_node["x"].as<double>();
        auto y = domain_node["y"].as<double>();
        auto z = domain_node["z"].as<double>();
        Domain domain;
        domain.dimension = {x, y, z};
        settings.domain = domain;
    }
}
// Add these new methods that parse without generating
std::vector<YAMLReader::CuboidData> YAMLReader::parseCuboids(const YAML::Node& node) {
    std::vector<CuboidData> cuboids;
    try {
        CuboidData data;

        const YAML::Node& coordinates = node["coordinates"];
        data.position[0] = coordinates["x"].as<double>();
        data.position[1] = coordinates["y"].as<double>();
        data.position[2] = coordinates["z"].as<double>();

        const YAML::Node& velocity_node = node["velocity"];
        data.velocity[0] = velocity_node["vx"].as<double>();
        data.velocity[1] = velocity_node["vy"].as<double>();
        data.velocity[2] = velocity_node["vz"].as<double>();

        const YAML::Node& count_node = node["particleNum"];
        data.num_particles[0] = count_node["nx"].as<size_t>();
        data.num_particles[1] = count_node["ny"].as<size_t>();
        data.num_particles[2] = count_node["nz"].as<size_t>();

        data.mass = node["mass"].as<double>();
        data.distance = node["distance"].as<double>();
        data.avg_velo = node["mean_velo"].as<double>();
        data.epsilon = node["epsilon"].as<double>();
        data.sigma = node["sigma"].as<double>();

        cuboids.push_back(data);
    }

    catch (const YAML::Exception& e) {
        SPDLOG_ERROR("Error parsing cuboids: {}", e.what());
        throw YAMLReaderException(e.what());
    }
    return cuboids;
}

std::vector<YAMLReader::DiscData> YAMLReader::parseDiscs(const YAML::Node& node) {
    std::vector<DiscData> discs;
    try {
        DiscData data;

        const YAML::Node& coordinates = node["coordinates"];
        data.position[0] = coordinates["x"].as<double>();
        data.position[1] = coordinates["y"].as<double>();
        data.position[2] = coordinates["z"].as<double>();

        const YAML::Node& velocity_node = node["velocity"];
        data.velocity[0] = velocity_node["vx"].as<double>();
        data.velocity[1] = velocity_node["vy"].as<double>();
        data.velocity[2] = velocity_node["vz"].as<double>();

        data.radius = node["radius"].as<size_t>();
        data.mass = node["mass"].as<double>();
        data.distance = node["distance"].as<double>();
        data.avg_velo = node["mean_velo"].as<double>();
        data.epsilon = node["epsilon"].as<double>();
        data.sigma = node["sigma"].as<double>();

        discs.push_back(data);

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

}  // namespace mol_sim
