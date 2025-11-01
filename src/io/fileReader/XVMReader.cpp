/*
 * FileReader.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "io/fileReader/XVMReader.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#include "utils/Logging.h"

using namespace mol_sim;

XVMReader::XVMReader() = default;

XVMReader::~XVMReader() = default;

void XVMReader::readFile(ContainerRef particles, const std::string& filename) {
    std::array<double, 3> x;
    std::array<double, 3> v;
    double m;
    int num_particles = 0;

    std::ifstream input_file(filename);
    std::string tmp_string;

    if (input_file.is_open()) {
        getline(input_file, tmp_string);

        LOG_TRACE("Read line: {}", tmp_string);

        while (tmp_string.empty() or tmp_string[0] == '#') {
            getline(input_file, tmp_string);
            LOG_TRACE("Read line: {}", tmp_string);
        }

        std::istringstream numstream(tmp_string);
        numstream >> num_particles;
        LOG_TRACE("Reading {} particles", num_particles);
        getline(input_file, tmp_string);
        LOG_TRACE("Read line: {}", tmp_string);

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
                LOG_ERROR("Error reading file: eof reached unexpectedly reading from line {}", i);
                exit(-1);
            }
            datastream >> m;
            particles.addParticle(x, v, m);

            getline(input_file, tmp_string);
            LOG_TRACE("Read line: {}", tmp_string);
        }
    } else {
        LOG_ERROR("Error: could not open file {}", filename);
        exit(-1);
    }
}
