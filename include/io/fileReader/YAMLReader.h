#ifndef YAML_READER_H
#define YAML_READER_H

#include <yaml-cpp/node/node.h>

#include <cstddef>

#include "io/FileReader.h"
#include "particles/container/ContainerRef.h"
#include "utils/Vector.h"
#include "yaml-cpp/yaml.h"

namespace mol_sim {
using R3 = Vector<double, 3>;
using N3 = Vector<size_t, 3>;
class YAMLReader : public FileReader {
   public:
    YAMLReader();
    ~YAMLReader() override;
    void readFile(ContainerRef particles, const std::string& filename) override;

   private:
    void readXVM(ContainerRef particles, YAML::Node& node);
    void readCube(ContainerRef particles, YAML::Node& node);
};

}  // namespace mol_sim

#endif
