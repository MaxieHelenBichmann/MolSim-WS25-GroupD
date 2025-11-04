#include "io/fileReader/YAMLReader.h"

#include <spdlog/spdlog.h>

#include <cstddef>
#include <cstdlib>
#include <string>

#include "io/fileReader/YAMLReaderException.h"
#include "particles/generators/CuboidGenerator.h"
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

                particles.addParticle(position, velocity, mass);
            }
        }
    } catch (const YAML::Exception& e) {
        SPDLOG_ERROR("Error parsing YAML: {}", e.what());
        throw YAMLReaderException(e.what());
    }
}
void YAMLReader::readCube(ContainerRef particles, const YAML::Node& node) {
    try {
        if (node["cuboids"] && node["cuboids"].IsSequence()) {
            auto num_cube = node["num_cuboids"].as<size_t>();
            particles.reserve(num_cube);
            for (const auto& curr : node["cuboids"]) {
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

                N3 num_particles;
                const YAML::Node& count_node = curr["particleNum"];
                num_particles[0] = count_node["nx"].as<size_t>();
                num_particles[1] = count_node["ny"].as<size_t>();
                num_particles[2] = count_node["nz"].as<size_t>();

                auto mass = curr["mass"].as<double>();
                auto distance = curr["distance"].as<double>();
                auto avg_velo = curr["mean_velo"].as<double>();

                CuboidGenerator generator(position, velocity, num_particles, mass, distance, avg_velo);
                generator.generateParticles(particles);
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
}
}  // namespace mol_sim
