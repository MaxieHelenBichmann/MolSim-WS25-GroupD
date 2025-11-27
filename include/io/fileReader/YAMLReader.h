#ifndef YAML_READER_H
#define YAML_READER_H

#include <sys/stat.h>
#include <yaml-cpp/node/node.h>
#include <cstddef>
#include <variant>

#include "io/FileReader.h"
#include "particles/container/ContainerRef.h"
#include "utils/Settings.h"
#include "utils/Vector.h"
#include "yaml-cpp/yaml.h"

namespace mol_sim {
using N3 = Vector<size_t, 3>;

/**
  * @brief Names of the supported ParticleContainers.
  * @note May be deprecated soon-ish
  * 
  * SIMPLE = SimpleContainer 
  * LINKED = LinkedCellContainer
  */
const std::string SIMPLE = "SIMPLE";
const std::string LINKED = "LINKED";

/**
 * @brief Reads in .yaml files that follow supported Structure
 * Supported formats: XVM, Cuboid, Disc
 */
class YAMLReader : public FileReader {
    public:
    /**
     * @brief Struct to hold Cuboid Generation Data
    *
    */
    struct CuboidData {
        R3 position;
        R3 velocity;
        N3 num_particles;
        double mass;
        double distance;
        double avg_velo;
        double epsilon;
        double sigma;
    };
    /**
     * @brief Struct to hold Disc Generation Data
     *
     */
    struct DiscData {
        R3 position;
        R3 velocity;
        size_t radius;
        double mass;
        double distance;
        double avg_velo;
        double epsilon;
        double sigma;
    };
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
    std::vector<YAMLReader::CuboidData> parseCuboids(const YAML::Node& node);
    std::vector<YAMLReader::DiscData> parseDiscs(const YAML::Node& node);

   private:
    static void readSettings(SettingsParam& settings, const YAML::Node& node);
    static void parseDomain(SettingsParam& settings, const YAML::Node& node);
    static void readBoundaryCondition(std::optional<BoundaryCondition&>& boundary, BoundaryLocation location, R3 dimension, const YAML::Node& node);
    void readXVM(ContainerRef particles, const YAML::Node& node);
    void readCube(ContainerRef particles, const YAML::Node& node);
    void readDisc(ContainerRef particles, const YAML::Node& node);
};

}  // namespace mol_sim

#endif
