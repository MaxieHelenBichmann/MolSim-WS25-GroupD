#ifndef YAML_READER_H
#define YAML_READER_H

#include <sys/stat.h>
#include <yaml-cpp/node/node.h>

#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "io/FileReader.h"
#include "particles/boundaries/Boundary.h"
#include "particles/container/ContainerRef.h"
#include "utils/Settings.h"
#include "utils/Vector.h"
#include "yaml-cpp/yaml.h"

namespace mol_sim {
using N3 = Vector<size_t, 3>;

/**
 * @brief Names of the supported ParticleContainers.
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
     * @brief Phase 1: Parses settings from the YAML file.
     *
     * @param settings SettingsParam where Simulation parameters are stored.
     * @param filename Path to input file.
     */
    void readSettings(SettingsParam& settings, const std::string& filename) override;

    /**
     * @brief Phase 2: Reads particles from the YAML file into a container.
     *
     * @param particles Container to place particles in.
     * @param filename Path to input file.
     */
    void readParticles(ContainerRef particles, const std::string& filename) override;

    std::vector<CuboidData> parseCuboids(const YAML::Node& node);
    std::vector<DiscData> parseDiscs(const YAML::Node& node);

   private:
    void parseDomain(SettingsParam& settings, const YAML::Node& node);
    void readXVM(ContainerRef particles, const YAML::Node& node);
    void readCube(ContainerRef particles, const YAML::Node& node);
    void readDisc(ContainerRef particles, const YAML::Node& node);
};

}  // namespace mol_sim

#endif
