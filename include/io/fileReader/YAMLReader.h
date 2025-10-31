#ifndef YAML_READER_H
#define YAML_READER_H

#include "io/FileReader.h"
#include "particles/container/ContainerRef.h"
#include "yaml-cpp/yaml.h"

namespace mol_sim {

class YAMLReader : public FileReader {
   public:
    YAMLReader();
    ~YAMLReader() override;
    void readFile(ContainerRef particles, const std::string& filename) override;
};

}  // namespace mol_sim

#endif
