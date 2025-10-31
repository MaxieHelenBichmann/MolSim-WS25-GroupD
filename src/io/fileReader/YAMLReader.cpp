#include "io/fileReader/YAMLReader.h"

#include <string>

namespace mol_sim {
void YAMLReader::readFile(ContainerRef particles, const std::string& filename) {
    YAML::Node root = YAML::LoadFile(filename);
    if (root["format"].as<std::string>() == "XVM") {
    }
}
}  // namespace mol_sim
