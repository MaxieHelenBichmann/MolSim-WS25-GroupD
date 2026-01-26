#include "io/fileReader/YAMLReader.h"

#include <spdlog/spdlog.h>
#include <yaml-cpp/exceptions.h>
#include <yaml-cpp/node/node.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <string>
#include <utility>

#include "exceptions/ValidationException.h"
#include "exceptions/YAMLReaderException.h"
#include "particles/boundaries/Boundary.h"
#include "particles/boundaries/Outflow.h"
#include "particles/boundaries/Periodic.h"
#include "particles/boundaries/Reflecting.h"
#include "particles/boundaries/VelocityReflect.h"
#include "particles/container/domain/Domain.h"
#include "particles/generators/CuboidGenerator.h"
#include "particles/generators/DiscGenerator.h"
#include "particles/generators/MembraneGenerator.h"
#include "physics/pairwiseforces/PairwiseForceSource.h"
#include "utils/Settings.h"
#include "utils/Vector.h"

namespace mol_sim {

namespace {
/**
 * @brief      Used for validating general Simulation Settings.
 *
 * @param[in]  settings
 * @throws ValidationException when bad values are read.
 */
void validateSettings(SettingsParam& settings) {
    if (settings.delta_t <= 0) {
        SPDLOG_ERROR("delta_t must be positive, got: " + std::to_string(settings.delta_t));
        throw ValidationException("delta_t must be positive, got: " + std::to_string(settings.delta_t));
    }
    if (settings.end_time < settings.start_time) {
        SPDLOG_ERROR("end_time must be greater than start_time");
        throw ValidationException("end_time must be greater than start_time");
    }
    if (settings.cutoff < 0) {
        SPDLOG_ERROR("cutoff must be non-negative, got: " + std::to_string(settings.cutoff));
        throw ValidationException("cutoff must be non-negative, got: " + std::to_string(settings.cutoff));
    }
    if (settings.smoothing < 0) {
        SPDLOG_ERROR("smoothing must be non-negative, got: " + std::to_string(settings.smoothing));
        throw ValidationException("smoothing must be non-negative, got: " + std::to_string(settings.smoothing));
    }
    // Check if smooth LJ is used, then smoothing must be <= cutoff
    if (std::ranges::find(settings.pairwise_forces, S_LENNARDJONES) != settings.pairwise_forces.end()) {
        if (settings.smoothing > settings.cutoff) {
            SPDLOG_ERROR(
                "For SMOOTHLENNARDJONES: smoothing radius must be <= cutoff radius, got smoothing={}, cutoff={}",
                settings.smoothing, settings.cutoff);
            throw ValidationException("For SMOOTHLENNARDJONES: smoothing radius must be <= cutoff radius");
        }
    }
    if (settings.delta_temp <= 0) {
        SPDLOG_ERROR("delta_temp must be greater then 0, got: " + std::to_string(settings.delta_temp));
        throw ValidationException("delta_temp must be greater then 0, got: " + std::to_string(settings.delta_temp));
    }
    if (settings.frequency_output <= 0) {
        SPDLOG_ERROR("frequency must be non-negative, got: " + std::to_string(settings.frequency_output));
        throw ValidationException("frequency must be non-negative, got: " + std::to_string(settings.frequency_output));
    }
    if (settings.frequency_checkpoint <= 0) {
        SPDLOG_ERROR("checkpoint must be non-negative, got: " + std::to_string(settings.frequency_checkpoint));
        throw ValidationException("checkpoint must be non-negative, got: " +
                                  std::to_string(settings.frequency_checkpoint));
    }
    if (settings.thermostat_freq <= 0) {
        SPDLOG_ERROR("thermostat_freq must be non-negative, got: " + std::to_string(settings.thermostat_freq));
        throw ValidationException("thermostat_freq must be non-negative, got: " +
                                  std::to_string(settings.thermostat_freq));
    }
    if (settings.pairwise_forces.empty()) {
        SPDLOG_ERROR("At least one pairwise force has to be specified");
        throw ValidationException("At least one pairwise force has to be specified");
    }
    if (settings.strategy == ParallelizationStrategy::COLORING) {
#ifndef _OPENMP
        SPDLOG_WARN("COLORING strategy requires OpenMP, falling back to NAIVE");
        settings.strategy = ParallelizationStrategy::NAIVE;
#else
        if (settings.container_type != "LINKED") {
            SPDLOG_WARN("COLORING strategy only works with LINKED container, falling back to NAIVE");
            settings.strategy = ParallelizationStrategy::NAIVE;
        }
#endif
    }
}
/**
 * @brief      Used for validating particle related Settings.
 *
 * @param[in]  settings
 * @throws ValidationException when bad values are read.
 */
void validateParticleParams(double mass, double epsilon, double sigma, const std::string& context) {
    if (mass <= 0) {
        SPDLOG_ERROR(context + ": mass must be positive, got: " + std::to_string(mass));
        throw ValidationException(context + ": mass must be positive, got: " + std::to_string(mass));
    }
    if (epsilon <= 0) {
        SPDLOG_ERROR(context + ": epsilon must be positive, got: " + std::to_string(epsilon));
        throw ValidationException(context + ": epsilon must be positive, got: " + std::to_string(epsilon));
    }
    if (sigma <= 0) {
        SPDLOG_ERROR(context + ": sigma must be positive, got: " + std::to_string(sigma));
        throw ValidationException(context + ": sigma must be positive, got: " + std::to_string(sigma));
    }
}
}  // namespace

YAMLReader::YAMLReader() = default;

YAMLReader::~YAMLReader() = default;

void YAMLReader::readSettings(SettingsParam& settings, const std::string& filename) {
    try {
        YAML::Node root = YAML::LoadFile(filename);

        if (root.size() == 0) {
            SPDLOG_ERROR("Empty YAML file");
            throw YAMLReaderException("Empty YAML file");
        }

        // First block must be settings
        auto first_it = root.begin();
        YAML::Node node = first_it->second;

        if (!node["format"] || node["format"].as<std::string>() != "Settings") {
            SPDLOG_ERROR("First block must be 'Settings' format");
            throw YAMLReaderException("First block must be 'Settings' format");
        }

        if (node["delta_t"]) {
            settings.delta_t = node["delta_t"].as<double>();
        }
        if (node["end_time"]) {
            settings.end_time = node["end_time"].as<double>();
        }
        if (node["start_time"]) {
            settings.start_time = node["start_time"].as<double>();
        }
        if (node["base_name"]) {
            settings.base_name = node["base_name"].as<std::string>();
        }
        if (node["container"]) {
            auto container_str = node["container"].as<std::string>();
            if (container_str != "SIMPLE" && container_str != "LINKED") {
                throw ValidationException("Unknown container type: " + container_str + " (expected SIMPLE or LINKED)");
            }
            settings.container_type = container_str;
        }
        if (node["strategy"]) {
            auto strategy_str = node["strategy"].as<std::string>();
            if (strategy_str == "NAIVE") {
                settings.strategy = ParallelizationStrategy::NAIVE;
            } else if (strategy_str == "COLORING") {
                settings.strategy = ParallelizationStrategy::COLORING;
            } else {
                throw ValidationException("Unknown parallelization strategy: " + strategy_str +
                                          " (expected NAIVE or COLORING)");
            }
        }
        if (node["frequency"]) {
            settings.frequency_output = node["frequency"].as<size_t>();
        }
        if (node["checkpoint"]) {
            settings.frequency_checkpoint = node["checkpoint"].as<size_t>();
        }
        if (node["cutoff"]) {
            settings.cutoff = node["cutoff"].as<double>();
        }
        if (node["smooth"]) {
            settings.smoothing = node["smooth"].as<double>();
        }
        if (node["enable_brownian"]) {
            settings.brownian = node["enable_brownian"].as<bool>();
        }
        if (node["pairwise_forces"]) {
            settings.pairwise_forces.clear();
            YAML::Node pairwise_node = node["pairwise_forces"];
            if (pairwise_node.IsSequence()) {
                for (const auto& force_item : pairwise_node) {
                    auto force_str = force_item.as<std::string>();
                    if (force_str == "GRAVITATIONAL") {
                        settings.pairwise_forces.push_back(PairwiseForce::GRAVITATIONAL);
                    } else if (force_str == "LENNARDJONES") {
                        settings.pairwise_forces.push_back(PairwiseForce::LENNARDJONES);
                    } else if (force_str == "TRUNCLENNARDJONES") {
                        settings.pairwise_forces.push_back(PairwiseForce::TRUNCLENNARDJONES);
                    } else if (force_str == "SMOOTHLENNARDJONES") {
                        settings.pairwise_forces.push_back(PairwiseForce::S_LENNARDJONES);
                    } else {
                        SPDLOG_ERROR("Unknown pairwise force type: " + force_str);
                        throw ValidationException("Unknown pairwise force type: " + force_str);
                    }
                }
            } else {
                SPDLOG_ERROR("pairwise_forces must be a sequence");
                throw ValidationException("pairwise_forces must be a sequence");
            }
        }
        if (node["single_forces"]) {
            settings.single_forces.clear();
            YAML::Node single_node = node["single_forces"];
            if (single_node.IsSequence()) {
                for (const auto& force_item : single_node) {
                    if (force_item.IsMap() && force_item["type"]) {
                        // New object-based format
                        auto force_type = force_item["type"].as<std::string>();
                        if (force_type == "GRAV") {
                            settings.single_forces.push_back(SingleForce::GRAV);
                            if (force_item["g_grav"]) {
                                const YAML::Node& g_grav_node = force_item["g_grav"];
                                if (g_grav_node.IsSequence() && g_grav_node.size() == 3) {
                                    settings.g_grav_vec[0] = g_grav_node[0].as<double>();
                                    settings.g_grav_vec[1] = g_grav_node[1].as<double>();
                                    settings.g_grav_vec[2] = g_grav_node[2].as<double>();
                                } else {
                                    SPDLOG_ERROR("GRAV force: g_grav must be [gx, gy, gz]");
                                    throw ValidationException("GRAV force: g_grav must be [gx, gy, gz]");
                                }
                            }
                        } else if (force_type == "HARMONIC") {
                            settings.single_forces.push_back(SingleForce::HARMONIC);
                            if (force_item["k"]) {
                                settings.k = force_item["k"].as<double>();
                            }
                            if (force_item["r_0"]) {
                                settings.r_0 = force_item["r_0"].as<double>();
                            }
                        } else {
                            SPDLOG_ERROR("Unknown single force type: " + force_type);
                            throw ValidationException("Unknown single force type: " + force_type);
                        }
                    } else {
                        SPDLOG_ERROR("single_forces items must be objects with 'type' field or strings");
                        throw ValidationException("single_forces items must be objects with 'type' field or strings");
                    }
                }
            } else {
                SPDLOG_ERROR("single_forces must be a sequence");
                throw ValidationException("single_forces must be a sequence");
            }
        }

        if (node["target_force"]) {
            settings.target_force_enabled = true;
            YAML::Node target_node = node["target_force"];
            if (target_node["direction"]) {
                const YAML::Node& dir_node = target_node["direction"];
                if (dir_node.IsSequence() && dir_node.size() == 3) {
                    settings.target_force_direction[0] = dir_node[0].as<double>();
                    settings.target_force_direction[1] = dir_node[1].as<double>();
                    settings.target_force_direction[2] = dir_node[2].as<double>();
                } else {
                    SPDLOG_ERROR("target_force: direction must be [dx, dy, dz]");
                    throw ValidationException("target_force: direction must be [dx, dy, dz]");
                }
            }
            if (target_node["magnitude"]) {
                settings.target_force_magnitude = target_node["magnitude"].as<double>();
            }
            if (target_node["max_iterations"]) {
                settings.target_force_max_iterations = target_node["max_iterations"].as<size_t>();
            }
        } else {
            settings.target_force_enabled = false;
        }

        if (node["thermostat"]) {
            settings.thermo = true;
            YAML::Node t_node = node["thermostat"];
            if (t_node["initial_temp"]) {
                settings.init_temp = t_node["initial_temp"].as<double>();
            }
            if (t_node["target_temp"]) {
                settings.target_temp = t_node["target_temp"].as<double>();
            } else {
                settings.target_temp = settings.init_temp;
            }
            if (t_node["n_thermostat"]) {
                settings.thermostat_freq = t_node["n_thermostat"].as<size_t>();
            }
            if (t_node["delta_temp"]) {
                settings.delta_temp = t_node["delta_temp"].as<double>();
            }
        } else {
            settings.thermo = false;
        }
        if (node["statistics"]) {
            YAML::Node t_node = node["statistics"];
            if (t_node["diffusion"]) {
                settings.stats_freq_diffusion = t_node["diffusion"].as<size_t>();
                settings.diff = true;
            } else {
                settings.diff = false;
            }
            if (t_node["rdf"]) {
                settings.stats_freq_rdf = t_node["rdf"].as<size_t>();
                settings.rdf = true;
            } else {
                settings.rdf = false;
            }
            if (t_node["sample_r"]) {
                if (!settings.rdf) {
                    throw ValidationException(
                        "Sample radius can only be set if data collection for the RDF is enabled");
                }
                settings.sample_radius = t_node["sample_r"].as<double>();
            }
            if (t_node["window_size"]) {
                if (!settings.rdf) {
                    throw ValidationException("Window size can only be set if data collection for the RDF is enabled");
                }
                settings.window_size = t_node["window_size"].as<double>();
            }
        } else {
            settings.rdf = false;
            settings.diff = false;
        }

        if (node["domain"]) {
            parseDomain(settings, node["domain"]);
        }

        validateSettings(settings);
        SPDLOG_DEBUG("Settings validated: delta_t={}, t=[{}, {}], cutoff={}", settings.delta_t, settings.start_time,
                     settings.end_time, settings.cutoff);
    } catch (const YAML::Exception& e) {
        SPDLOG_ERROR("Error parsing YAML settings: {}", e.what());
        throw YAMLReaderException(e.what());
    }
}

void YAMLReader::readParticles(ContainerRef particles, const SettingsParam& settings, const std::string& filename) {
    try {
        YAML::Node root = YAML::LoadFile(filename);

        if (root.size() == 0) {
            SPDLOG_ERROR("Empty YAML file");
            throw YAMLReaderException("Empty YAML file");
        }

        // Process all blocks looking for particle definitions
        bool has_particle_definition = false;
        for (auto it = root.begin(); it != root.end(); ++it) {
            const auto& name = it->first.as<std::string>();
            YAML::Node node = it->second;

            if (!node["format"]) {
                SPDLOG_WARN("Block '{}' has no format field, skipping", name);
                continue;
            }

            auto format = node["format"].as<std::string>();
            if (format == "XVM") {
                has_particle_definition = true;
                readXVM(particles, node);
            } else if (format == "Cuboid") {
                has_particle_definition = true;
                readCube(particles, settings, node);
            } else if (format == "Disc") {
                has_particle_definition = true;
                readDisc(particles, settings, node);
            } else if (format == "Membrane") {
                has_particle_definition = true;
                readMembrane(particles, settings, node);
            }
            // Skip Settings and unknown formats silently in phase 2
        }

        //this assumes that all files that contain particles (namely in the case of checkpoint files) have already been
        //parsed before this one.
        if (!has_particle_definition && particles.empty()) {
            throw YAMLReaderException("No particle definitions (XVM, Cuboid, or Disc) found in the input file");
        }
    } catch (const YAML::Exception& e) {
        SPDLOG_ERROR("Error parsing YAML particles: {}", e.what());
        throw YAMLReaderException(e.what());
    }
}

void YAMLReader::readXVM(ContainerRef particles, const YAML::Node& node) {
    try {
        if (node["particles"] && node["particles"].IsSequence()) {
            auto num_part = node["num_particles"].as<size_t>();
            particles.reserve(num_part);
            size_t particle_idx = 0;
            for (const auto& curr : node["particles"]) {
                R3 position;
                const YAML::Node& coordinates = curr["coordinates"];
                if (coordinates.IsSequence() && coordinates.size() == 3) {
                    position[0] = coordinates[0].as<double>();
                    position[1] = coordinates[1].as<double>();
                    position[2] = coordinates[2].as<double>();
                } else {
                    SPDLOG_ERROR("XVM: coordinates must be [x, y, z]");
                    throw ValidationException("XVM: coordinates must be [x, y, z]");
                }

                R3 old_position = R3{0., 0., 0.};
                if (curr["old_coordinates"]) {
                    const YAML::Node& old_coordinates_node = curr["old_coordinates"];
                    if (old_coordinates_node.IsSequence() && old_coordinates_node.size() == 3) {
                        old_position[0] = old_coordinates_node[0].as<double>();
                        old_position[1] = old_coordinates_node[1].as<double>();
                        old_position[2] = old_coordinates_node[2].as<double>();
                    }
                }

                R3 velocity;
                const YAML::Node& velocity_node = curr["velocity"];
                if (velocity_node.IsSequence() && velocity_node.size() == 3) {
                    velocity[0] = velocity_node[0].as<double>();
                    velocity[1] = velocity_node[1].as<double>();
                    velocity[2] = velocity_node[2].as<double>();
                } else {
                    SPDLOG_ERROR("XVM: velocity must be [vx, vy, vz]");
                    throw ValidationException("XVM: velocity must be [vx, vy, vz]");
                }

                R3 force = R3{0., 0., 0.};
                if (curr["force"]) {
                    const YAML::Node& force_node = curr["force"];
                    if (force_node.IsSequence() && force_node.size() == 3) {
                        force[0] = force_node[0].as<double>();
                        force[1] = force_node[1].as<double>();
                        force[2] = force_node[2].as<double>();
                    }
                }

                R3 old_force = R3{0., 0., 0.};
                if (curr["old_force"]) {
                    const YAML::Node& old_force_node = curr["old_force"];
                    if (old_force_node.IsSequence() && old_force_node.size() == 3) {
                        old_force[0] = old_force_node[0].as<double>();
                        old_force[1] = old_force_node[1].as<double>();
                        old_force[2] = old_force_node[2].as<double>();
                    }
                }

                auto mass = curr["mass"].as<double>();

                double epsilon = SettingsParam::EPSILON_DEFAULT;
                double sigma = SettingsParam::SIGMA_DEFAULT;
                if (curr["epsilon"]) {
                    epsilon = curr["epsilon"].as<double>();
                }
                if (curr["sigma"]) {
                    sigma = curr["sigma"].as<double>();
                }
                int type = 0;
                if (curr["type"]) {
                    type = curr["type"].as<int>();
                }
                validateParticleParams(mass, epsilon, sigma, "XVM particle " + std::to_string(particle_idx));
                particles.addParticle(position, old_position, velocity, force, old_force, mass, epsilon, sigma, type);
                particle_idx++;
            }
            SPDLOG_DEBUG("Parsed {} XVM particles", particle_idx);
        }
    } catch (const YAML::Exception& e) {
        SPDLOG_ERROR("Error parsing XVM particles: {}", e.what());
        throw YAMLReaderException(e.what());
    }
}

std::vector<YAMLReader::CuboidData> YAMLReader::parseCuboids(const YAML::Node& node) {
    std::vector<CuboidData> cuboids;
    try {
        CuboidData data;

        const YAML::Node& coordinates = node["coordinates"];
        if (coordinates.IsSequence() && coordinates.size() == 3) {
            data.position[0] = coordinates[0].as<double>();
            data.position[1] = coordinates[1].as<double>();
            data.position[2] = coordinates[2].as<double>();
        } else {
            SPDLOG_ERROR("Cuboid: coordinates must be [x, y, z]");
            throw ValidationException("Cuboid: coordinates must be [x, y, z]");
        }

        const YAML::Node& velocity_node = node["velocity"];
        if (velocity_node.IsSequence() && velocity_node.size() == 3) {
            data.velocity[0] = velocity_node[0].as<double>();
            data.velocity[1] = velocity_node[1].as<double>();
            data.velocity[2] = velocity_node[2].as<double>();
        } else {
            SPDLOG_ERROR("Cuboid: velocity must be [vx, vy, vz]");
            throw ValidationException("Cuboid: velocity must be [vx, vy, vz]");
        }

        const YAML::Node& count_node = node["particleNum"];
        if (count_node.IsSequence() && count_node.size() == 3) {
            data.num_particles[0] = count_node[0].as<size_t>();
            data.num_particles[1] = count_node[1].as<size_t>();
            data.num_particles[2] = count_node[2].as<size_t>();
        } else {
            SPDLOG_ERROR("Cuboid: particleNum must be [nx, ny, nz]");
            throw ValidationException("Cuboid: particleNum must be [nx, ny, nz]");
        }

        data.mass = node["mass"].as<double>();
        data.distance = node["distance"].as<double>();
        if (node["mean_velo"]) {
            data.avg_velo = node["mean_velo"].as<double>();
        }
        data.epsilon = node["epsilon"].as<double>();
        data.sigma = node["sigma"].as<double>();

        if (node["targets"] && node["targets"].IsSequence()) {
            for (const auto& target : node["targets"]) {
                if (target.IsSequence() && target.size() == 3) {
                    N3 target_pos;
                    target_pos[0] = target[0].as<size_t>();
                    target_pos[1] = target[1].as<size_t>();
                    target_pos[2] = target[2].as<size_t>();
                    data.targets.push_back(target_pos);
                } else {
                    SPDLOG_ERROR("Cuboid: target must be a sequence of 3 integers [x, y, z]");
                    throw ValidationException("Cuboid: target must be a sequence of 3 integers [x, y, z]");
                }
            }
            SPDLOG_DEBUG("Parsed {} target particles for cuboid", data.targets.size());
        }

        validateParticleParams(data.mass, data.epsilon, data.sigma, "Cuboid");
        if (data.distance <= 0) {
            SPDLOG_ERROR("Cuboid: distance must be positive, got: " + std::to_string(data.distance));
            throw ValidationException("Cuboid: distance must be positive, got: " + std::to_string(data.distance));
        }

        SPDLOG_DEBUG("Parsed cuboid: {}x{}x{} particles at ({}, {}, {})", data.num_particles[0], data.num_particles[1],
                     data.num_particles[2], data.position[0], data.position[1], data.position[2]);

        cuboids.push_back(data);
    }

    catch (const YAML::Exception& e) {
        SPDLOG_ERROR("Error parsing cuboids: {}", e.what());
        throw YAMLReaderException(e.what());
    }
    return cuboids;
}

std::vector<YAMLReader::DiscData> YAMLReader::parseDiscs(const YAML::Node& node) {
    std::vector<DiscData> discs;
    try {
        DiscData data;

        const YAML::Node& coordinates = node["coordinates"];
        if (coordinates.IsSequence() && coordinates.size() == 3) {
            data.position[0] = coordinates[0].as<double>();
            data.position[1] = coordinates[1].as<double>();
            data.position[2] = coordinates[2].as<double>();
        } else {
            SPDLOG_ERROR("Disc: coordinates must be [x, y, z]");
            throw ValidationException("Disc: coordinates must be [x, y, z]");
        }

        const YAML::Node& velocity_node = node["velocity"];
        if (velocity_node.IsSequence() && velocity_node.size() == 3) {
            data.velocity[0] = velocity_node[0].as<double>();
            data.velocity[1] = velocity_node[1].as<double>();
            data.velocity[2] = velocity_node[2].as<double>();
        } else {
            SPDLOG_ERROR("Disc: velocity must be [vx, vy, vz]");
            throw ValidationException("Disc: velocity must be [vx, vy, vz]");
        }

        data.radius = node["radius"].as<size_t>();
        data.mass = node["mass"].as<double>();
        data.distance = node["distance"].as<double>();
        if (node["mean_velo"]) {
            data.avg_velo = node["mean_velo"].as<double>();
        }
        data.epsilon = node["epsilon"].as<double>();
        data.sigma = node["sigma"].as<double>();

        if (node["targets"] && node["targets"].IsSequence()) {
            for (const auto& target : node["targets"]) {
                if (target.IsSequence() && target.size() == 3) {
                    N3 target_pos;
                    target_pos[0] = target[0].as<size_t>();
                    target_pos[1] = target[1].as<size_t>();
                    target_pos[2] = target[2].as<size_t>();
                    data.targets.push_back(target_pos);
                } else {
                    SPDLOG_ERROR("Disc: target must be a sequence of 3 integers [x, y, z]");
                    throw ValidationException("Disc: target must be a sequence of 3 integers [x, y, z]");
                }
            }
            SPDLOG_DEBUG("Parsed {} target particles for disc", data.targets.size());
        }

        validateParticleParams(data.mass, data.epsilon, data.sigma, "Disc");
        if (data.distance <= 0) {
            SPDLOG_ERROR("Disc: distance must be positive, got: " + std::to_string(data.distance));
            throw ValidationException("Disc: distance must be positive, got: " + std::to_string(data.distance));
        }

        SPDLOG_DEBUG("Parsed disc: radius={} at ({}, {}, {})", data.radius, data.position[0], data.position[1],
                     data.position[2]);

        discs.push_back(data);

    } catch (const YAML::Exception& e) {
        SPDLOG_ERROR("Error parsing discs: {}", e.what());
        throw YAMLReaderException(e.what());
    }
    return discs;
}

void YAMLReader::readCube(ContainerRef particles, const SettingsParam& settings, const YAML::Node& node) {
    auto cuboids = parseCuboids(node);
    for (const auto& data : cuboids) {
        CuboidGenerator generator(data.position, data.velocity, data.num_particles, data.targets, data.mass,
                                  data.distance, data.avg_velo, data.epsilon, data.sigma, settings.init_temp);
        generator.generateParticles(particles, settings.brownian, settings.thermo);
    }
}

void YAMLReader::readDisc(ContainerRef particles, const SettingsParam& settings, const YAML::Node& node) {
    auto discs = parseDiscs(node);
    for (const auto& data : discs) {
        DiscGenerator generator(data.position, data.velocity, data.radius, data.mass, data.distance, data.avg_velo,
                                data.epsilon, data.sigma, settings.init_temp);
        generator.generateParticles(particles, settings.brownian, settings.thermo);
    }
}

void YAMLReader::parseDomain(SettingsParam& settings, const YAML::Node& node) {
    try {
        R3 dimension;
        auto coordinates = node["coordinates"];
        if (coordinates && coordinates.IsSequence() && coordinates.size() == 3) {
            dimension = {coordinates[0].as<double>(), coordinates[1].as<double>(), coordinates[2].as<double>()};
        } else {
            SPDLOG_ERROR("Domain: coordinates must be [x, y, z]");
            throw ValidationException("Domain: coordinates must be [x, y, z]");
        }

        const YAML::Node& dimensions_node = node["dimensions"];  // NOLINT
        if (dimensions_node) {
            settings.dimensions = dimensions_node.as<size_t>();  // NOLINT
        }

        // Define boundary locations and their YAML keys
        static const std::array<std::pair<BoundaryLocation, std::string>, 6> boundary_mappings = {{
            {BoundaryLocation::LEFT, "left"},
            {BoundaryLocation::RIGHT, "right"},
            {BoundaryLocation::FRONT, "front"},
            {BoundaryLocation::BACK, "back"},
            {BoundaryLocation::UPPER, "upper"},
            {BoundaryLocation::LOWER, "lower"},
        }};

        std::array<std::unique_ptr<Boundary>, 6> boundaries;
        const YAML::Node& bounds_node = node["boundaries"];

        for (size_t i = 0; i < boundary_mappings.size(); ++i) {
            const auto& [location, key] = boundary_mappings[i];
            if (bounds_node && bounds_node[key]) {
                std::unique_ptr<Boundary> boundary;
                const YAML::Node& curr_node = bounds_node[key];
                auto type_str = curr_node["type"].as<std::string>("OUTFLOW");
                BoundaryType boundary_type = mol_sim::parseBoundaryType(type_str);

                switch (boundary_type) {
                    case BoundaryType::REFLECTING: {
                        std::optional<double> sigma =
                            curr_node["sigma"] ? std::optional<double>(curr_node["sigma"].as<double>()) : std::nullopt;
                        std::optional<double> epsilon = curr_node["epsilon"]
                                                            ? std::optional<double>(curr_node["epsilon"].as<double>())
                                                            : std::nullopt;
                        bool ghost_on_boundary =
                            curr_node["ghost_on_boundary"] ? curr_node["ghost_on_boundary"].as<bool>() : false;
                        boundary = std::make_unique<Reflecting>(location, dimension, ghost_on_boundary, sigma, epsilon);
                        break;
                    }
                    case BoundaryType::VELOCITYREFLECT:
                        boundary = std::make_unique<VelocityReflect>(location, dimension);
                        break;
                    case BoundaryType::PERIODIC: {
                        boundary =
                            std::make_unique<Periodic>(location, dimension, settings.cutoff, settings.dimensions);
                        break;
                    }
                    case BoundaryType::OUTFLOW:
                    default:
                        boundary = std::make_unique<Outflow>(location, dimension);
                }
                boundaries[i] = std::move(boundary);
            } else {
                boundaries[i] = std::make_unique<Outflow>(location, dimension);
            }
        }
        settings.domain = Domain(dimension, std::move(boundaries));
    } catch (YAML::Exception& e) {
        SPDLOG_ERROR("Error parsing domain: {}", e.what());
        throw YAMLReaderException(e.what());
    }
}

std::vector<YAMLReader::MembraneData> YAMLReader::parseMembranes(const YAML::Node& node) {
    std::vector<MembraneData> membranes;
    try {
        MembraneData data;

        const YAML::Node& coordinates = node["coordinates"];
        if (!coordinates.IsSequence() || coordinates.size() != 3) {
            SPDLOG_ERROR("Membrane: coordinates must be [x, y, z]");
            throw ValidationException("Membrane: coordinates must be [x, y, z]");
        }
        data.position[0] = coordinates[0].as<double>();
        data.position[1] = coordinates[1].as<double>();
        data.position[2] = coordinates[2].as<double>();

        const YAML::Node& velocity_node = node["velocity"];
        if (!velocity_node.IsSequence() || velocity_node.size() != 3) {
            SPDLOG_ERROR("Membrane: velocity must be [vx, vy, vz]");
            throw ValidationException("Membrane: velocity must be [vx, vy, vz]");
        }
        data.velocity[0] = velocity_node[0].as<double>();
        data.velocity[1] = velocity_node[1].as<double>();
        data.velocity[2] = velocity_node[2].as<double>();

        const YAML::Node& count_node = node["particleNum"];
        if (!count_node.IsSequence() || count_node.size() != 2) {
            SPDLOG_ERROR("Membrane: particleNum must be [nx, ny]");
            throw ValidationException("Membrane: particleNum must be [nx, ny]");
        }
        data.num_particles[0] = count_node[0].as<size_t>();
        data.num_particles[1] = count_node[1].as<size_t>();

        data.mass = node["mass"].as<double>();
        data.distance = node["distance"].as<double>();
        if (node["mean_velo"]) {
            data.avg_velo = node["mean_velo"].as<double>();
        }
        data.epsilon = node["epsilon"].as<double>();
        data.sigma = node["sigma"].as<double>();

        if (node["targets"] && node["targets"].IsSequence()) {
            for (const auto& target : node["targets"]) {
                if (target.IsSequence() && target.size() == 2) {
                    Vector<size_t, 2> target_pos;
                    target_pos[0] = target[0].as<size_t>();
                    target_pos[1] = target[1].as<size_t>();
                    data.targets.push_back(target_pos);
                } else {
                    SPDLOG_ERROR("Membrane: target must be a sequence of 2 integers [x, y]");
                    throw ValidationException("Membrane: target must be a sequence of 2 integers [x, y]");
                }
            }
            SPDLOG_DEBUG("Parsed {} target particles for membrane", data.targets.size());
        }

        validateParticleParams(data.mass, data.epsilon, data.sigma, "Membrane");
        if (data.distance <= 0) {
            SPDLOG_ERROR("Membrane: distance must be positive, got: " + std::to_string(data.distance));
            throw ValidationException("Membrane: distance must be positive, got: " + std::to_string(data.distance));
        }

        SPDLOG_DEBUG("Parsed membrane: {}x{} particles at ({}, {}, {})", data.num_particles[0], data.num_particles[1],
                     data.position[0], data.position[1], data.position[2]);

        membranes.push_back(data);

    } catch (const YAML::Exception& e) {
        SPDLOG_ERROR("Error parsing membranes: {}", e.what());
        throw YAMLReaderException(e.what());
    }
    return membranes;
}

void YAMLReader::readMembrane(ContainerRef particles, const SettingsParam& settings, const YAML::Node& node) {
    auto membranes = parseMembranes(node);
    for (const auto& data : membranes) {
        MembraneGenerator generator(data.position, data.velocity, data.num_particles, data.targets, data.mass,
                                    data.distance, data.avg_velo, data.epsilon, data.sigma, settings.init_temp);
        generator.generateParticles(particles, settings.brownian, settings.thermo);
    }
}

}  // namespace mol_sim
