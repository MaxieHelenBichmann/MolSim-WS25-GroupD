#ifndef CLI_PARSE_H
#define CLI_PARSE_H

#include <algorithm>
#include <iostream>
#include <string>
#include <cstring>

#include "io/FileReader.h"
#include "utils/Logging.h"
#include "physics/ForceSource.h"

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
    "| -f {GRAV, LJ}            : sets the force-type (default = GRAV)\n"
    "| -e <DOUBLE>              : sets epsilon for Lennard-Jones force (default = 5)\n"
    "| -o <DOUBLE>              : sets sigma for Lennard-Jones force (default = 1)\n"
#ifdef LOG_LEVEL_DEBUG
    "| -l {Info, Debug, Trace}  : sets Log Level (default = Info)\n"
#endif
    "| -h                       : displays this message";
const std::string& d = "-d";
const std::string& t = "-t";
const std::string& f = "-f";
const std::string& e = "-e";
const std::string& o = "-o";
#ifdef LOG_LEVEL_DEBUG
const std::string& l = "-l";
#endif
const std::string& h = "-h";

void cliParse(int argc, char** argsv, FileReader& fileReader, double& delta_t, double& end_time, 
              Force& force, std::vector<double>& params, SimpleContainer& particles) {
    std::cout << "Hello from MolSim for PSE!" << '\n';
    char** help = std::find(argsv, argsv + argc, h);
    if (help != &argsv[argc]) {
        std::cout << HELP_MSG << '\n';
        exit(0);
    }
    if (argc == 1) {
        std::cout << "Erroneous programme call! " << '\n' << HELP_MSG << '\n';
        exit(-1);
    }

    int parsed_args = 2;  // program name + assume file name is OK (bad files handled in fileReader.readfile)
    char** delta_t_opt = std::find(argsv, &argsv[argc], d);
    char** end_time_opt = std::find(argsv, &argsv[argc], t);
    char** force_opt = std::find(argsv, &argsv[argc], f);
    char** epsilon_opt = std::find(argsv, &argsv[argc], e);
    char** sigma_opt = std::find(argsv, &argsv[argc], o);
#ifdef LOG_LEVEL_DEBUG
    char** log_opt = std::find(argsv, &argsv[argc], l);
#endif
    try {
        if (delta_t_opt != &argsv[argc]) {
            delta_t = std::stod(*(++delta_t_opt));
            parsed_args += 2;
        }
        if (end_time_opt != &argsv[argc]) {
            end_time = std::stod(*(++end_time_opt));
            parsed_args += 2;
        }
        if (force_opt != &argsv[argc]) {
            char* force_string = *(++force_opt);
            bool grav = !static_cast<bool>(std::strcmp(force_string, "GRAV"));
            bool lj   = !static_cast<bool>(std::strcmp(force_string, "LJ"));
            
            if (grav) { 
                force = GRAVITATIONAL; 
                parsed_args += 2;
            }
            else if (lj) { 
                force = LENNARDJONES;
                parsed_args += 2;
                double epsilon = 5;
                double sigma = 1;
                if (epsilon_opt != &argsv[argc]) {
                    epsilon = std::stod(*(++epsilon_opt));
                    parsed_args += 2;
                } else if (sigma_opt == &argsv[argc]) { 
                    sigma = std::stod(*(++sigma_opt));
                    parsed_args += 2;
                }
                params.push_back(epsilon);
                params.push_back(sigma);
            }
        }
#ifdef LOG_LEVEL_DEBUG
        if (log_opt != &argsv[argc]) {
            logInit(*(++log_opt));
            parsed_args += 2;
        } else {
            logInit("Default");
        }
#endif
    } catch (std::invalid_argument& e) {
        std::cout << "Erroneous programme call! " << '\n'
                  << "Given options must be valid floating point numbers or strings!" << "\n"
                  << HELP_MSG << '\n';
        exit(-1);
    } catch (std::out_of_range& e) {
        std::cout << "Erroneous programme call! " << '\n'
                  << "At least one floating point out of double precision floating point range!" << "\n"
                  << HELP_MSG << '\n';
        exit(-1);
    } catch (std::logic_error& e) {
        std::cout << "Erroneous programme call! " << '\n'
                  << "Flag set but no value provided!" << "\n"
                  << HELP_MSG << '\n';
        exit(-1);
    }

    if (parsed_args - argc != 0) {
        std::cout << "Erroneous programme call! " << '\n'
                  << "Unrecognized arguments in programme call!" << "\n"
                  << HELP_MSG << '\n';
        exit(-1);
    }

    fileReader.readFile(particles, argsv[1]);
}

}  // namespace mol_sim

#endif
