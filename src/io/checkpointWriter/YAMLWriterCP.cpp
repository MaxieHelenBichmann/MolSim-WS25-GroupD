#include "io/checkpointWriter/YAMLWriterCP.h"

#include <spdlog/spdlog.h>
#include <yaml-cpp/emitter.h>
#include <yaml-cpp/exceptions.h>
#include <yaml-cpp/node/node.h>
#include <yaml-cpp/yaml.h>

#include <cstddef>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "particles/boundaries/Boundary.h"
#include "particles/boundaries/Reflecting.h"
#include "physics/ForceSource.h"

using namespace mol_sim;

YAMLWriterCP::YAMLWriterCP() = default;
YAMLWriterCP::~YAMLWriterCP() = default;

void YAMLWriterCP::createCheckpoint(SettingsParam& settings, const Domain& domain, ContainerRef particles,
                                    int iteration, size_t N) const {
    int decimal_places = 0;
    while (N >= 10) {
        N /= 10;
        decimal_places++;
    }
    std::ofstream file;
    std::stringstream strstr;
    strstr << "cp_" << std::setfill('0') << std::setw(decimal_places) << iteration << ".yaml";

    file.open(strstr.str().c_str());

    YAML::Emitter out;
    out << YAML::BeginMap;  // open root

    // settings block
    out << YAML::Key << "settings" << YAML::Value << YAML::BeginMap;

    out << YAML::Key << "format" << YAML::Value << "Settings";
    out << YAML::Key << "delta_t" << YAML::Value << settings.delta_t;
    out << YAML::Key << "end_time" << YAML::Value << settings.end_time;
    out << YAML::Key << "start_time" << YAML::Value
        << (settings.start_time > settings.end_time ? settings.end_time : settings.start_time);
    ;
    out << YAML::Key << "base_name" << YAML::Value << settings.base_name;
    switch (settings.force) {
        case LENNARDJONES: {
            out << YAML::Key << "force" << YAML::Value << "Lennard Jones";
            break;
        }
        case S_LENNARDJONES: {
            out << YAML::Key << "force" << YAML::Value << "Smooth Lennard Jones";
            break;
        }
        case GRAVITATIONAL: {
            out << YAML::Key << "force" << YAML::Value << "Gravitational";
            break;
        }
    };
    out << YAML::Key << "container" << YAML::Value << settings.container_type;
    out << YAML::Key << "frequency" << YAML::Value << settings.frequency_output;
    out << YAML::Key << "checkpoint" << YAML::Value << settings.frequency_checkpoint;
    out << YAML::Key << "cutoff" << YAML::Value << settings.cutoff;
    if (settings.force == S_LENNARDJONES) {
        out << YAML::Key << "smooth" << YAML::Value << settings.smoothing;
    }

    // Write thermostat settings in nested format
    if (settings.thermo) {
        out << YAML::Key << "thermostat" << YAML::Value << YAML::BeginMap;
        out << YAML::Key << "target_temp" << YAML::Value << settings.target_temp;
        out << YAML::Key << "n_thermostat" << YAML::Value << settings.thermostat_freq;
        out << YAML::Key << "delta_temp" << YAML::Value << settings.delta_temp;
        out << YAML::EndMap;  // close thermostat
    }

    out << YAML::Key << "domain" << YAML::Value << YAML::BeginMap;  // open domain
    out << YAML::Key << "x" << YAML::Value << domain.getDimension()[0];
    out << YAML::Key << "y" << YAML::Value << domain.getDimension()[1];
    out << YAML::Key << "z" << YAML::Value << domain.getDimension()[2];
    out << YAML::Key << "g_grav" << YAML::Value << settings.g_grav;
    out << YAML::Key << "dimensions" << YAML::Value << settings.dimensions;

    const auto boundary_type_string = [](BoundaryType type) -> std::string {
        switch (type) {
            case BoundaryType::OUTFLOW:
                return "OUTFLOW";
            case BoundaryType::REFLECTING:
                return "REFLECTING";
            case BoundaryType::VELOCITYREFLECT:
                return "VELOCITYREFLECT";
            case BoundaryType::PERIODIC:
                return "PERIODIC";
        }
        return "OUTFLOW";
    };

    const auto boundary_extra = [](const Boundary& b, YAML::Emitter& out) {
        if (b.getType() == BoundaryType::REFLECTING) {
            const auto& refl_boundary = dynamic_cast<const Reflecting&>(b);
            if (refl_boundary.getBoundarySigma().has_value()) {
                out << YAML::Key << "sigma" << YAML::Value << refl_boundary.getBoundarySigma().value();
            }
            if (refl_boundary.getBoundaryEpsilon().has_value()) {
                out << YAML::Key << "epsilon" << YAML::Value << refl_boundary.getBoundaryEpsilon().value();
            }
            out << YAML::Key << "ghost_on_boundary" << YAML::Value << refl_boundary.isGhostOnBoundary();
        }
    };

    out << YAML::Key << "boundaries" << YAML::Value << YAML::BeginMap;  // open boundaries

    out << YAML::Key << "upper" << YAML::Value << YAML::BeginMap << YAML::Key << "type" << YAML::Value
        << boundary_type_string(domain.getBoundary(BoundaryLocation::UPPER).getType());
    boundary_extra(domain.getBoundary(BoundaryLocation::UPPER), out);
    out << YAML::EndMap;

    out << YAML::Key << "lower" << YAML::Value << YAML::BeginMap << YAML::Key << "type" << YAML::Value
        << boundary_type_string(domain.getBoundary(BoundaryLocation::LOWER).getType());
    boundary_extra(domain.getBoundary(BoundaryLocation::LOWER), out);
    out << YAML::EndMap;

    out << YAML::Key << "left" << YAML::Value << YAML::BeginMap << YAML::Key << "type" << YAML::Value
        << boundary_type_string(domain.getBoundary(BoundaryLocation::LEFT).getType());
    boundary_extra(domain.getBoundary(BoundaryLocation::LEFT), out);
    out << YAML::EndMap;

    out << YAML::Key << "right" << YAML::Value << YAML::BeginMap << YAML::Key << "type" << YAML::Value
        << boundary_type_string(domain.getBoundary(BoundaryLocation::RIGHT).getType());
    boundary_extra(domain.getBoundary(BoundaryLocation::RIGHT), out);
    out << YAML::EndMap;

    out << YAML::Key << "front" << YAML::Value << YAML::BeginMap << YAML::Key << "type" << YAML::Value
        << boundary_type_string(domain.getBoundary(BoundaryLocation::FRONT).getType());
    boundary_extra(domain.getBoundary(BoundaryLocation::FRONT), out);
    out << YAML::EndMap;

    out << YAML::Key << "back" << YAML::Value << YAML::BeginMap << YAML::Key << "type" << YAML::Value
        << boundary_type_string(domain.getBoundary(BoundaryLocation::BACK).getType());
    boundary_extra(domain.getBoundary(BoundaryLocation::BACK), out);
    out << YAML::EndMap;

    out << YAML::EndMap;  // close boundaries
    out << YAML::EndMap;  // close domain
    out << YAML::EndMap;  // close settings

    // xvm section
    out << YAML::Key << "xvm_data" << YAML::Value << YAML::BeginMap;
    out << YAML::Key << "format" << YAML::Value << "XVM";
    out << YAML::Key << "num_particles" << YAML::Value << particles.size();
    out << YAML::Key << "particles" << YAML::Value << YAML::BeginSeq;  // open particles sequence
    for (const auto& p : particles) {
        out << YAML::BeginMap;

        out << YAML::Key << "coordinates" << YAML::Value << YAML::BeginMap;
        out << YAML::Key << "x" << YAML::Value << p.getX()[0];
        out << YAML::Key << "y" << YAML::Value << p.getX()[1];
        out << YAML::Key << "z" << YAML::Value << p.getX()[2];
        out << YAML::EndMap;

        out << YAML::Key << "old_coordinates" << YAML::Value << YAML::BeginMap;
        out << YAML::Key << "ox" << YAML::Value << p.getOldX()[0];
        out << YAML::Key << "oy" << YAML::Value << p.getOldX()[1];
        out << YAML::Key << "oz" << YAML::Value << p.getOldX()[2];
        out << YAML::EndMap;

        out << YAML::Key << "velocity" << YAML::Value << YAML::BeginMap;
        out << YAML::Key << "vx" << YAML::Value << p.getV()[0];
        out << YAML::Key << "vy" << YAML::Value << p.getV()[1];
        out << YAML::Key << "vz" << YAML::Value << p.getV()[2];
        out << YAML::EndMap;

        out << YAML::Key << "force" << YAML::Value << YAML::BeginMap;
        out << YAML::Key << "fx" << YAML::Value << p.getF()[0];
        out << YAML::Key << "fy" << YAML::Value << p.getF()[1];
        out << YAML::Key << "fz" << YAML::Value << p.getF()[2];
        out << YAML::EndMap;

        out << YAML::Key << "old_force" << YAML::Value << YAML::BeginMap;
        out << YAML::Key << "ofx" << YAML::Value << p.getOldF()[0];
        out << YAML::Key << "ofy" << YAML::Value << p.getOldF()[1];
        out << YAML::Key << "ofz" << YAML::Value << p.getOldF()[2];
        out << YAML::EndMap;

        out << YAML::Key << "mass" << YAML::Value << p.getM();
        out << YAML::Key << "epsilon" << YAML::Value << p.getEpsilon();
        out << YAML::Key << "sigma" << YAML::Value << p.getSigma();
        out << YAML::Key << "type" << YAML::Value << p.getType();

        out << YAML::EndMap;  // close particle
    }
    out << YAML::EndSeq;  // close particles sequence
    out << YAML::EndMap;  // close xvm_data

    out << YAML::EndMap;  // close root

    file << out.c_str();
    file.close();
}
