#ifndef CLI_PARSE_H
#define CLI_PARSE_H

#include <spdlog/spdlog.h>

#include <algorithm>
#include <iostream>
#include <string>

#include "io/FileReader.h"
#include "utils/Logging.h"
#include "utils/Settings.h"


namespace mol_sim {

/**
 * @brief Parses command line input.
 *
 * Parses command line input.
 * @param argc Number of arguments including the program name
 * @param argsv Array of arguments
 * @param fileReader The FileReader that will read the input file
 * @param particles Container to place the generated particles in.
 * @param settings SettingsParam where options for the simulation are stored.
 */

const std::string HELP_MSG =
    "Usage: ./MolSim path/to/input/file <ARGS>\n"
    "| -d <DOUBLE>              : sets DELTA_T (default = 0.014)\n"
    "| -t <DOUBLE>              : sets END_TIME (default = 1000)\n"
#if SPDLOG_ACTIVE_LEVEL == SPDLOG_LEVEL_TRACE
    "| -l {Info, Debug, Trace}  : sets Log Level (default = Info)\n"
#endif
    "| -h                       : displays this message";
const std::string& h = "-h";
const std::string& d = "-d";
const std::string& t = "-t";
#if SPDLOG_ACTIVE_LEVEL == SPDLOG_LEVEL_TRACE
const std::string& l = "-l";
#endif
void cliParse(int argc, char** argsv, FileReader& fileReader, SimpleContainer& particles, SettingsParam& settings) {
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
#if SPDLOG_ACTIVE_LEVEL == SPDLOG_LEVEL_TRACE
    char** log_opt = std::find(argsv, &argsv[argc], l);
#endif
    try {
        if (delta_t_opt != &argsv[argc]) {
            settings.delta_t = std::stod(*(++delta_t_opt));
            parsed_args += 2;
        }
        if (end_time_opt != &argsv[argc]) {
            settings.end_time = std::stod(*(++end_time_opt));
            parsed_args += 2;
        }
#if SPDLOG_ACTIVE_LEVEL == SPDLOG_LEVEL_TRACE
        if (log_opt != &argsv[argc]) {
            logInit(*(++log_opt));
            parsed_args += 2;
        } else {
            logInit("Default");
        }
#endif
    } catch (std::invalid_argument& e) {
        std::cout << "Erroneous programme call! " << '\n'
                  << "given options must be valid floating point numbers or strings!" << "\n"
                  << HELP_MSG << '\n';
        exit(-1);
    } catch (std::out_of_range& e) {
        std::cout << "Erroneous programme call! " << '\n'
                  << "delta_t or end_time out of floating point range!" << "\n"
                  << HELP_MSG << '\n';
        exit(-1);
    } catch (std::logic_error& e) {
        std::cout << "Erroneous programme call! " << '\n'
                  << "flag -d or -t set but no value provided!" << "\n"
                  << HELP_MSG << '\n';
        exit(-1);
    }

    if (parsed_args - argc != 0) {
        std::cout << "Erroneous programme call! " << '\n'
                  << "Unrecognized arguments in programme call!" << "\n"
                  << HELP_MSG << '\n';
        exit(-1);
    }

    fileReader.readFile(particles, settings, argsv[1]);
}

}  // namespace mol_sim

#endif
