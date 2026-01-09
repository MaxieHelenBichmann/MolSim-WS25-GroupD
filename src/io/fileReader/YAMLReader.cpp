#include "io/fileReader/YAMLReader.h"

#include <spdlog/spdlog.h>
#include <yaml-cpp/exceptions.h>
#include <yaml-cpp/node/node.h>

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
#include "physics/ForceSource.h"
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
void validateSettings(const SettingsParam& settings) {
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
// NOLINTNEXTLINE
void YAMLReader::readSettings(SettingsParam& settings, const std::string& filename) {
    try {
        YAML::Node root = YAML::LoadFile(filename);

        if (root.size() == 0) {
            throw YAMLReaderException("Empty YAML file");
        }

        // First block must be settings
        auto first_it = root.begin();
        YAML::Node node = first_it->second;

        if (!node["format"] || node["format"].as<std::string>() != "Settings") {
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
        if (node["force"]) {
            auto force_str = node["force"].as<std::string>();
            if (force_str == "Lennard Jones") {
                settings.force = LENNARDJONES;
            } else if (force_str == "Smooth Lennard Jones") {
                settings.force = S_LENNARDJONES;
            } else if (force_str == "Gravitational") {
                settings.force = GRAVITATIONAL;
            } else {
                throw ValidationException("Unknown force type: " + force_str);
            }
        }
        if (node["container"]) {
            auto container_str = node["container"].as<std::string>();
            if (container_str != "SIMPLE" && container_str != "LINKED") {
                throw ValidationException("Unknown container type: " + container_str + " (expected SIMPLE or LINKED)");
            }
            settings.container_type = container_str;
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
            if (settings.force != S_LENNARDJONES) {
                throw ValidationException("Smoothing radius can only be set for Smooth Lennard-Jones force");
            }
            settings.smoothing = node["smooth"].as<double>();
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
            }
            // Skip Settings and unknown formats silently in phase 2
        }

        if (!has_particle_definition) {
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
                position[0] = coordinates["x"].as<double>();
                position[1] = coordinates["y"].as<double>();
                position[2] = coordinates["z"].as<double>();

                R3 old_position = R3{0., 0., 0.};
                if (curr["old_coordinates"]) {
                    const YAML::Node& old_coordinates_node = curr["old_coordinates"];
                    old_position[0] = old_coordinates_node["ox"].as<double>();
                    old_position[1] = old_coordinates_node["oy"].as<double>();
                    old_position[2] = old_coordinates_node["oz"].as<double>();
                }

                R3 velocity;
                const YAML::Node& velocity_node = curr["velocity"];
                velocity[0] = velocity_node["vx"].as<double>();
                velocity[1] = velocity_node["vy"].as<double>();
                velocity[2] = velocity_node["vz"].as<double>();

                R3 force = R3{0., 0., 0.};
                if (curr["force"]) {
                    const YAML::Node& force_node = curr["force"];
                    force[0] = force_node["fx"].as<double>();
                    force[1] = force_node["fy"].as<double>();
                    force[2] = force_node["fz"].as<double>();
                }

                R3 old_force = R3{0., 0., 0.};
                if (curr["old_force"]) {
                    const YAML::Node& old_force_node = curr["old_force"];
                    old_force[0] = old_force_node["ofx"].as<double>();
                    old_force[1] = old_force_node["ofy"].as<double>();
                    old_force[2] = old_force_node["ofz"].as<double>();
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
        data.position[0] = coordinates["x"].as<double>();
        data.position[1] = coordinates["y"].as<double>();
        data.position[2] = coordinates["z"].as<double>();

        const YAML::Node& velocity_node = node["velocity"];
        data.velocity[0] = velocity_node["vx"].as<double>();
        data.velocity[1] = velocity_node["vy"].as<double>();
        data.velocity[2] = velocity_node["vz"].as<double>();

        const YAML::Node& count_node = node["particleNum"];
        data.num_particles[0] = count_node["nx"].as<size_t>();
        data.num_particles[1] = count_node["ny"].as<size_t>();
        data.num_particles[2] = count_node["nz"].as<size_t>();

        data.mass = node["mass"].as<double>();
        data.distance = node["distance"].as<double>();
        if (node["mean_velo"]) {
            data.avg_velo = node["mean_velo"].as<double>();
        }
        data.epsilon = node["epsilon"].as<double>();
        data.sigma = node["sigma"].as<double>();

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
        data.position[0] = coordinates["x"].as<double>();
        data.position[1] = coordinates["y"].as<double>();
        data.position[2] = coordinates["z"].as<double>();

        const YAML::Node& velocity_node = node["velocity"];
        data.velocity[0] = velocity_node["vx"].as<double>();
        data.velocity[1] = velocity_node["vy"].as<double>();
        data.velocity[2] = velocity_node["vz"].as<double>();

        data.radius = node["radius"].as<size_t>();
        data.mass = node["mass"].as<double>();
        data.distance = node["distance"].as<double>();
        if (node["mean_velo"]) {
            data.avg_velo = node["mean_velo"].as<double>();
        }
        data.epsilon = node["epsilon"].as<double>();
        data.sigma = node["sigma"].as<double>();

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
        CuboidGenerator generator(data.position, data.velocity, data.num_particles, data.mass, data.distance,
                                  data.avg_velo, data.epsilon, data.sigma, settings.init_temp);
        generator.generateParticles(particles);
    }
}

void YAMLReader::readDisc(ContainerRef particles, const SettingsParam& settings, const YAML::Node& node) {
    auto discs = parseDiscs(node);
    for (const auto& data : discs) {
        DiscGenerator generator(data.position, data.velocity, data.radius, data.mass, data.distance, data.avg_velo,
                                data.epsilon, data.sigma, settings.init_temp);
        generator.generateParticles(particles);
    }
}

void YAMLReader::parseDomain(SettingsParam& settings, const YAML::Node& node) {
    try {
        R3 dimension = {node["x"].as<double>(), node["y"].as<double>(), node["z"].as<double>()};
        const YAML::Node& g_grav_node = node["g_grav"];
        if (g_grav_node) {
            settings.g_grav = g_grav_node.as<double>();
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

}  // namespace mol_sim
