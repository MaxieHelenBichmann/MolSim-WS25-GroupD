#include "io/fileReader/YAMLReader.h"

#include <cstddef>
#include <cstdlib>
#include <string>

#include "io/fileReader/YAMLReaderException.h"
#include "particles/generators/CuboidGenerator.h"
#include "utils/Logging.h"

namespace mol_sim {
YAMLReader::YAMLReader() = default;

YAMLReader::~YAMLReader() = default;
void YAMLReader::readFile(ContainerRef particles, const std::string& filename) {
    try {
        YAML::Node root = YAML::LoadFile(filename);
        auto format = root["format"].as<std::string>();
        if (format == "XVM") {
            readXVM(particles, root);
        } else if (format == "Cuboid") {
            readCube(particles, root);
        } else {
            SPDLOG_ERROR("Unknown YAML Format");
            throw YAMLReaderException("Unknown YAML Format");
        }
    } catch (const YAML::Exception& e) {
        SPDLOG_ERROR("Error parsing YAML: {}", e.what());
        throw YAMLReaderException(e.what());
    }
}
void YAMLReader::readXVM(ContainerRef particles, YAML::Node& node) {
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
void YAMLReader::readCube(ContainerRef particles, YAML::Node& node) {
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
}  // namespace mol_sim
