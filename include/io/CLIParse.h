#ifndef CLI_PARSE_H
#define CLI_PARSE_H

#include <iostream>
#include <string>
#include <algorithm>

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

const std::string HELP_MSG =
    "Usage: ./MolSim path/to/input/file <ARGS>\n"
    "| -d <DOUBLE>              : sets DELTA_T (default = 0.014)\n"
    "| -t <DOUBLE>              : sets END_TIME (default = 1000)\n"
    "| -h                       : displays this message";
const std::string& h = "-h";
const std::string& d = "-d";
const std::string& t = "-t";

void cliParse(int argc, char** argsv, FileReader& fileReader, double& delta_t, double& end_time,
              SimpleContainer& particles) {
  std::cout << "Hello from MolSim for PSE!" << '\n';
  char** help = std::find(argsv, argsv+argc, h);
  if (help != &argsv[argc]) {
    std::cout << HELP_MSG << '\n';
    exit(0);
  }
  if (argc == 1) {
    std::cout << "Erroneous programme call! " << '\n'
              << HELP_MSG << '\n';
    exit(-1);
  }
  fileReader.readFile(particles, argsv[1]);
  char** delta_t_opt = std::find(argsv, &argsv[argc], d);
  char** end_time_opt = std::find(argsv, &argsv[argc], t);

  try {
    if (delta_t_opt != &argsv[argc]) {
      delta_t = std::stod(*(++delta_t_opt));
    }
    if (end_time_opt != &argsv[argc]) {
      end_time = std::stod(*(++end_time_opt));
    }
  } catch (std::invalid_argument& e) {
    std::cout << "Erroneous programme call! " << '\n'
              << "delta_t and end_time must be valid floating point numbers" << "\n"
              << HELP_MSG << '\n';
    exit(-1);
  } catch (std::out_of_range& e) {
    std::cout << "Erroneous programme call! " << '\n'
              << "delta_t or end_time out of floating point range" << "\n"
              << HELP_MSG << '\n';
    exit(-1);
  }
}

}  // namespace mol_sim

#endif