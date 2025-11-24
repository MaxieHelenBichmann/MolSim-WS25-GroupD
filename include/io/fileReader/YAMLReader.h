#ifndef YAML_READER_H
#define YAML_READER_H

#include <yaml-cpp/node/node.h>

#include <cstddef>

#include "io/FileReader.h"
#include "particles/container/ContainerRef.h"
#include "utils/Settings.h"
#include "utils/Vector.h"
#include "yaml-cpp/yaml.h"

namespace mol_sim {
using R3 = Vector<double, 3>;
using N3 = Vector<size_t, 3>;
/**
 * @brief Reads in .yaml files that follow supported Structure
 * Supported formats: XVM, Cuboid
 */
class YAMLReader : public FileReader {
   public:
    YAMLReader();
    ~YAMLReader() override;
    /**
     * @brief function that reads in parameters from .yaml files
     *
     * @param particles Container to place particles in.
     * @param settings SettingsParam where Simulation parameters are stored.
     * @param filename Path to input file.
     */
    void readFile(ContainerRef particles, SettingsParam& settings, const std::string& filename) override;

   private:
    static void readXVM(ContainerRef particles, const YAML::Node& node);
    static void readCube(ContainerRef particles, const YAML::Node& node);
    static void readSettings(SettingsParam& Settings, const YAML::Node& node);
};

}  // namespace mol_sim

#endif
