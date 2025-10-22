/*
 * FileReader.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#include "XVMReader.h"

XVMReader::XVMReader() = default;

XVMReader::~XVMReader() = default;

void XVMReader::readFile(particle_containers::ParticleContainer& particles, char* filename) {
  std::array<double, 3> x = {.0,.0,.0};
  std::array<double, 3> v = {.0,.0,.0};
  double m;
  int num_particles = 0;

  std::ifstream input_file(filename);
  std::string tmp_string;

  if (input_file.is_open()) {
    getline(input_file, tmp_string);
    std::cout << "Read line: " << tmp_string << '\n';

    while (tmp_string.empty() or tmp_string[0] == '#') {
      getline(input_file, tmp_string);
      std::cout << "Read line: " << tmp_string << '\n';
    }

    std::istringstream numstream(tmp_string);
    numstream >> num_particles;
    std::cout << "Reading " << num_particles << "." << '\n';
    getline(input_file, tmp_string);
    std::cout << "Read line: " << tmp_string << '\n';

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
        std::cout << "Error reading file: eof reached unexpectedly reading from line " << i << '\n';
        exit(-1);
      }
      datastream >> m;
      particles.addParticle(x, v, m);

      getline(input_file, tmp_string);
      std::cout << "Read line: " << tmp_string << '\n';
    }
  } else {
    std::cout << "Error: could not open file " << filename << '\n';
    exit(-1);
  }
}
