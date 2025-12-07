#include "io/checkpointWriter/YAMLWriterCP.h"

#include <spdlog/spdlog.h>
#include <yaml-cpp/exceptions.h>
#include <yaml-cpp/node/node.h>
#include <yaml-cpp/yaml.h>

#include <cstddef>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "particles/boundaries/Boundary.h"

using namespace mol_sim;

YAMLWriterCP::YAMLWriterCP() = default;
YAMLWriterCP::~YAMLWriterCP() = default;

void YAMLWriterCP::createCheckpoint(const Domain& domain, ContainerRef particles, int iteration, Force force,
                                    double delta_t, double start_time, double end_time, size_t frequency_output,
                                    size_t frequency_checkpoint, const std::string& base_name, double cutoff_radius,
                                    size_t N) const {
    int decimal_places = 0;
    while (N >= 10) {
        N /= 10;
        decimal_places++;
    }
    std::ofstream file;
    std::stringstream strstr;
    strstr << "cp_" << std::setfill('0') << std::setw(decimal_places) << iteration << ".txt";

    file.open(strstr.str().c_str());

    YAML::Emitter out;
    out << YAML::BeginMap;  // open root

    // settings block
    out << YAML::Key << "settings" << YAML::Value << YAML::BeginMap;

    out << YAML::Key << "format" << YAML::Value << "Settings";
    out << YAML::Key << "delta_t" << YAML::Value << delta_t;
    out << YAML::Key << "end_time" << YAML::Value << end_time;
    out << YAML::Key << "start_time" << YAML::Value << start_time;
    out << YAML::Key << "base_name" << YAML::Value << base_name;
    out << YAML::Key << "force" << YAML::Value << (force == LENNARDJONES ? "Lennard Jones" : "Gravitational");
    out << YAML::Key << "frequency" << YAML::Value << frequency_output;
    out << YAML::Key << "checkpoint" << YAML::Value << frequency_checkpoint;
    out << YAML::Key << "cutoff" << YAML::Value << cutoff_radius;

    out << YAML::Key << "domain" << YAML::Value << YAML::BeginMap;  // open domain
    out << YAML::Key << "x" << YAML::Value << domain.getDimension()[0];
    out << YAML::Key << "y" << YAML::Value << domain.getDimension()[1];
    out << YAML::Key << "z" << YAML::Value << domain.getDimension()[2];

    const auto boundary_type_string = [](BoundaryType type) -> std::string {
        switch (type) {
            case BoundaryType::OUTFLOW:
                return "OUTFLOW";
            case BoundaryType::REFLECTING:
                return "REFLECTING";
            case BoundaryType::VELOCITYREFLECT:
                return "VELOCITYREFLECT";
        }
        return "OUTFLOW";
    };

    out << YAML::Key << "boundaries" << YAML::Value << YAML::BeginMap;  // open boundaries

    out << YAML::Key << "UPPER" << YAML::Value << YAML::BeginMap << YAML::Key << "type" << YAML::Value
        << boundary_type_string(domain.getBoundary(BoundaryLocation::UPPER).getType()) << YAML::EndMap;
    out << YAML::Key << "LOWER" << YAML::Value << YAML::BeginMap << YAML::Key << "type" << YAML::Value
        << boundary_type_string(domain.getBoundary(BoundaryLocation::LOWER).getType()) << YAML::EndMap;
    out << YAML::Key << "LEFT" << YAML::Value << YAML::BeginMap << YAML::Key << "type" << YAML::Value
        << boundary_type_string(domain.getBoundary(BoundaryLocation::LEFT).getType()) << YAML::EndMap;
    out << YAML::Key << "RIGHT" << YAML::Value << YAML::BeginMap << YAML::Key << "type" << YAML::Value
        << boundary_type_string(domain.getBoundary(BoundaryLocation::RIGHT).getType()) << YAML::EndMap;
    out << YAML::Key << "FRONT" << YAML::Value << YAML::BeginMap << YAML::Key << "type" << YAML::Value
        << boundary_type_string(domain.getBoundary(BoundaryLocation::FRONT).getType()) << YAML::EndMap;
    out << YAML::Key << "BACK" << YAML::Value << YAML::BeginMap << YAML::Key << "type" << YAML::Value
        << boundary_type_string(domain.getBoundary(BoundaryLocation::BACK).getType()) << YAML::EndMap;

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

        out << YAML::Key << "oldCoordinates" << YAML::Value << YAML::BeginMap;
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

        out << YAML::Key << "oldForce" << YAML::Value << YAML::BeginMap;
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