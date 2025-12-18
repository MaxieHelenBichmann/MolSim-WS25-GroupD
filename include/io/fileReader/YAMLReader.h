#ifndef YAML_READER_H
#define YAML_READER_H

#include <sys/stat.h>
#include <yaml-cpp/node/node.h>

#include <cstddef>
#include <string>
#include <vector>

#include "io/FileReader.h"
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
     * @throws YAMLReaderException if the file cannot be parsed or is malformed.
     * @throws ValidationException if settings values are invalid.
     */
    void readSettings(SettingsParam& settings, const std::string& filename) override;

    /**
     * @brief Phase 2: Reads particles from the YAML file into a container.
     *
     * @param particles Container to place particles in.
     * @param filename Path to input file.
     * @throws YAMLReaderException if the file cannot be parsed or has no particle definitions.
     * @throws ValidationException if particle parameters are invalid.
     */
    void readParticles(ContainerRef particles, const std::string& filename) override;
    /**
     * @brief      Helper function to parse Cuboid Format
     *
     * @param[in]  node YAML::Node of the start of the cuboid block
     *
     * @return     Vector of read in Cuboids
     * @throws YAMLReaderException if parsing fails.
     * @throws ValidationException if cuboid parameters are invalid.
     */
    std::vector<CuboidData> parseCuboids(const YAML::Node& node);
    /**
     * @brief      Helper function to parse Disc Format
     *
     * @param[in]  node YAML::Node of the start of the disc block
     *
     * @return     Vector of read in Discs
     * @throws YAMLReaderException if parsing fails.
     * @throws ValidationException if disc parameters are invalid.
     */
    std::vector<DiscData> parseDiscs(const YAML::Node& node);

   private:
    /**
     * @brief      Helper function to parse Domains
     *
     * @param[in]  node YAML::Node of the start of the domain block
     */
    void parseDomain(SettingsParam& settings, const YAML::Node& node);
    /**
     * @brief      Helper function to create XVM particles
     *
     * @param[in]  node YAML::Node of the start of the particles block
     */
    void readXVM(ContainerRef particles, const YAML::Node& node);
    /**
     * @brief      Helper function to create cuboid particles
     *
     * @param[in]  node YAML::Node of the start of the cuboid block
     */
    void readCube(ContainerRef particles, const YAML::Node& node);
    /**
     * @brief      Helper function to create disc particles
     *
     * @param[in]  node YAML::Node of the start of the disc block
     */
    void readDisc(ContainerRef particles, const YAML::Node& node);
};

}  // namespace mol_sim

#endif
