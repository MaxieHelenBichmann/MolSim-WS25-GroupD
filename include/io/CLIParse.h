#ifndef CLI_PARSE_H
#define CLI_PARSE_H

#include <iostream>
#include <string>

#include "io/FileReader.h"

namespace mol_sim {

/**
 * @brief Parses command line input.
 *
 * Parses command line input.
 * @param argc Number of arguments including the program name
 * @param argsv Array of arguments
 * @param fileReader The FileReader that will read the input file
 * @param delta_t The variable that will store the provided timestep interval
 * @param end_time The variable that will store the provided end time of the simulation
 * @param particles A ParticleContainer wherein the particles read from the input file will be stored
 */
void cliParse(int argc, char* argsv[], FileReader& fileReader, double& delta_t, double& end_time,
              SimpleContainer& particles) {
  std::cout << "Hello from MolSim for PSE!" << '\n';
  if (argc != 4) {
    std::cout << "Erroneous programme call! " << '\n' << "./MolSim filename delta_t end_time" << '\n';
    exit(-1);
  }
  fileReader.readFile(particles, argsv[1]);
  try {
    delta_t = std::stod(argsv[2]);
    end_time = std::stod(argsv[3]);
  } catch (std::invalid_argument& e) {
    std::cout << "Erroneous programme call! " << '\n'
              << "delta_t and end_time must be valid floating point numbers" << "\n";
    exit(-1);
  } catch (std::out_of_range& e) {
    std::cout << "Erroneous programme call! " << '\n' << "delta_t or end_time out of floating point range" << "\n";
    exit(-1);
  }
}

}  // namespace mol_sim

#endif