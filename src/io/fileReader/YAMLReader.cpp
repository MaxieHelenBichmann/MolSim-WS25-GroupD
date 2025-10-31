#include "io/fileReader/YAMLReader.h"

#include <cstdlib>
#include <iostream>
#include <string>

namespace mol_sim {
void YAMLReader::readFile(ContainerRef particles, const std::string& filename) {
    try {
        YAML::Node root = YAML::LoadFile(filename);
        auto format = root["format"].as<std::string>();
        if (format == "XVM") {
            readXVM(particles, root);
        } else if (format == "Cuboid") {
            readCube(particles, root);
        } else {
            std::cerr << "Unknown YAML Format" << '\n';
        }
    } catch (const YAML::Exception& e) {
        std::cerr << "Error parsing YAML: " << e.what() << '\n';
        exit(EXIT_FAILURE);
    }
}
void YAMLReader::readXVM(ContainerRef particles, YAML::Node& node) {
    try {
        if (node["particles"] && node["particles"].IsSequence()) {
            for (const auto& curr : node["particles"]) {
                R3 position;
                const YAML::Node& coordinates = curr["coordinates"];
                position[0] = coordinates[0].as<double>();
                position[1] = coordinates[1].as<double>();
                position[2] = coordinates[2].as<double>();

                R3 velocity;
                const YAML::Node& velocity_node = curr["velocity"];
                velocity[0] = velocity_node[0].as<double>();
                velocity[1] = velocity_node[1].as<double>();
                velocity[2] = velocity_node[2].as<double>();

                auto mass = curr["mass"].as<double>();

                particles.addParticle(position, velocity, mass);
            };
        }
    } catch (const YAML::Exception& e) {
        std::cerr << "Error parsing YAML: " << e.what() << '\n';
        exit(EXIT_FAILURE);
    }
}
void YAMLReader::readCube(ContainerRef particles, YAML::Node& node) {}
}  // namespace mol_sim
