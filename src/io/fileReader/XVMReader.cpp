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

XVMReader::XVMReader() = default;

XVMReader::~XVMReader() = default;

void XVMReader::readSettings([[maybe_unused]] SettingsParam& settings, [[maybe_unused]] const std::string& filename) {
    SPDLOG_DEBUG("XVM files do not contain settings, using defaults");
}

void XVMReader::readParticles(ContainerRef particles, const std::string& filename) {
    std::array<double, 3> x;
    std::array<double, 3> v;
    double m;
    int num_particles = 0;

    std::ifstream input_file(filename);
    std::string tmp_string;

    if (input_file.is_open()) {
        getline(input_file, tmp_string);

        SPDLOG_TRACE("Read line: {}", tmp_string);

        while (tmp_string.empty() or tmp_string[0] == '#') {
            getline(input_file, tmp_string);
            SPDLOG_TRACE("Read line: {}", tmp_string);
        }

        std::istringstream numstream(tmp_string);
        numstream >> num_particles;

        if (num_particles <= 0) {
            SPDLOG_ERROR("Invalid particle count: " + std::to_string(num_particles));
            throw XVMReaderException("Invalid particle count: " + std::to_string(num_particles));
        }

        SPDLOG_DEBUG("Reading {} particles from XVM file", num_particles);
        getline(input_file, tmp_string);
        SPDLOG_TRACE("Read line: {}", tmp_string);

        particles.reserve(num_particles);

        for (int i = 0; i < num_particles; i++) {
            std::istringstream datastream(tmp_string);

            for (auto& xj : x) {
                datastream >> xj;
            }
            for (auto& vj : v) {
                datastream >> vj;
            }
            if (datastream.eof()) {
                SPDLOG_ERROR("Error reading file: eof reached unexpectedly reading from line {}", i);
                throw XVMReaderException("Error reading file: eof reached unexpectedly at line " + std::to_string(i));
            }
            datastream >> m;

            if (m <= 0) {
                SPDLOG_ERROR("XVM particle " + std::to_string(i) +
                             ": mass must be positive, got: " + std::to_string(m));
                throw ValidationException("XVM particle " + std::to_string(i) +
                                          ": mass must be positive, got: " + std::to_string(m));
            }

            particles.addParticle(x, v, m, SettingsParam::EPSILON_DEFAULT, SettingsParam::SIGMA_DEFAULT);

            getline(input_file, tmp_string);
            SPDLOG_TRACE("Read line: {}", tmp_string);
        }

        if (particles.size() != static_cast<size_t>(num_particles)) {
            SPDLOG_WARN("Expected {} particles, read {}", num_particles, particles.size());
        }
        SPDLOG_DEBUG("Successfully parsed {} XVM particles", particles.size());
    } else {
        SPDLOG_ERROR("Error: could not open file {}", filename);
        throw XVMReaderException("Could not open file: " + filename);
    }
}
