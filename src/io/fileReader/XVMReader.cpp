/*
 * FileReader.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "io/fileReader/XVMReader.h"

#include <spdlog/spdlog.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#include "exceptions/ValidationException.h"
#include "exceptions/XVMReaderException.h"
#include "utils/Settings.h"

using namespace mol_sim;

namespace {

int readParticleCount(std::ifstream& input_file, std::string& line, const std::string& filename) {
    if (!std::getline(input_file, line)) {
        SPDLOG_ERROR("Error reading file {}: missing particle count", filename);
        throw XVMReaderException("Missing particle count in file: " + filename);
    }
    SPDLOG_TRACE("Read line: {}", line);

    while (line.empty() || line[0] == '#') {
        if (!std::getline(input_file, line)) {
            SPDLOG_ERROR("Error reading file {}: missing particle count after header", filename);
            throw XVMReaderException("Missing particle count in file: " + filename);
        }
        SPDLOG_TRACE("Read line: {}", line);
    }

    std::istringstream numstream(line);
    int num_particles = 0;
    numstream >> num_particles;

    if (!numstream || num_particles <= 0) {
        SPDLOG_ERROR("Invalid particle count: {}", line);
        throw XVMReaderException("Invalid particle count: " + line);
    }

    return num_particles;
}

bool readNextDataLine(std::ifstream& input_file, std::string& line) {
    while (std::getline(input_file, line)) {
        SPDLOG_TRACE("Read line: {}", line);
        if (!line.empty() && line[0] != '#') {
            return true;
        }
    }
    return false;
}

void parseParticleLine(int index, const std::string& line, std::array<double, 3>& x, std::array<double, 3>& v,
                       double& m) {
    std::istringstream datastream(line);

    for (auto& xj : x) {
        if (!(datastream >> xj)) {
            SPDLOG_ERROR("Malformed position data for particle {}", index);
            throw XVMReaderException("Malformed position data for particle " + std::to_string(index));
        }
    }
    for (auto& vj : v) {
        if (!(datastream >> vj)) {
            SPDLOG_ERROR("Malformed velocity data for particle {}", index);
            throw XVMReaderException("Malformed velocity data for particle " + std::to_string(index));
        }
    }
    if (!(datastream >> m)) {
        SPDLOG_ERROR("Malformed mass entry for particle {}", index);
        throw XVMReaderException("Malformed mass entry for particle " + std::to_string(index));
    }
}

}  // namespace

XVMReader::XVMReader() = default;

XVMReader::~XVMReader() = default;

void XVMReader::readSettings([[maybe_unused]] SettingsParam& settings, [[maybe_unused]] const std::string& filename) {
    SPDLOG_DEBUG("XVM files do not contain settings, using defaults");
}

void XVMReader::readParticles(ContainerRef particles, const SettingsParam& settings, const std::string& filename) {
    (void)settings;
    std::array<double, 3> x;
    std::array<double, 3> v;
    double m;
    int num_particles = 0;

    std::ifstream input_file(filename);
    std::string tmp_string;

    if (!input_file.is_open()) {
        SPDLOG_ERROR("Error: could not open file {}", filename);
        throw XVMReaderException("Could not open file: " + filename);
    }

    const int num_particles = readParticleCount(input_file, tmp_string, filename);

    SPDLOG_DEBUG("Reading {} particles from XVM file", num_particles);
    particles.reserve(num_particles);

    if (!readNextDataLine(input_file, tmp_string)) {
        SPDLOG_ERROR("Error reading file: no particle data after header in {}", filename);
        throw XVMReaderException("Error reading file: no particle data after header in " + filename);
    }

    for (int i = 0; i < num_particles; ++i) {
        parseParticleLine(i, tmp_string, x, v, m);

        if (m <= 0) {
            SPDLOG_ERROR("XVM particle {}: mass must be positive, got: {}", i, m);
            throw ValidationException("XVM particle " + std::to_string(i) +
                                      ": mass must be positive, got: " + std::to_string(m));
        }

        particles.addParticle(x, v, m, SettingsParam::EPSILON_DEFAULT, SettingsParam::SIGMA_DEFAULT);

        if (i < num_particles - 1 && !readNextDataLine(input_file, tmp_string)) {
            SPDLOG_ERROR("Error reading file: eof reached unexpectedly at particle {}", i + 1);
            throw XVMReaderException("Error reading file: eof reached unexpectedly at particle " +
                                     std::to_string(i + 1));
        }
    }

    if (particles.size() != static_cast<size_t>(num_particles)) {
        SPDLOG_WARN("Expected {} particles, read {}", num_particles, particles.size());
    }
    SPDLOG_DEBUG("Successfully parsed {} XVM particles", particles.size());
}
