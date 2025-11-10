#ifndef CLI_PARSE_H
#define CLI_PARSE_H

#include <spdlog/spdlog.h>

#include <algorithm>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>

#include "io/FileReader.h"
#include "io/fileReader/XVMReader.h"
#include "io/fileReader/YAMLReader.h"
#include "utils/Logging.h"
#include "utils/Settings.h"

namespace mol_sim {

/**
 * @brief Parses command line input.
 *
 * Parses command line input.
 * @param argc Number of arguments including the program name
 * @param argsv Array of arguments
 * @param particles Container to place the generated particles in.
 * @param settings SettingsParam where options for the simulation are stored.
 */

const std::string HELP_MSG =
    "Usage: ./MolSim path/to/input/file <ARGS>\n"
    "| -d <DOUBLE>              : sets DELTA_T (default = 0.014)\n"
    "| -t <DOUBLE>              : sets END_TIME (default = 1000)\n"
    "| -f {GRAV, LJ}            : sets the force-type (default = LJ)\n"
    "| -e <DOUBLE>              : sets epsilon for Lennard-Jones force (default = 5)\n"
    "| -o <DOUBLE>              : sets sigma for Lennard-Jones force (default = 1)\n"
#if SPDLOG_ACTIVE_LEVEL == SPDLOG_LEVEL_TRACE
    "| -l {Info, Debug, Trace}  : sets Log Level (default = Info)\n"
#endif
    "| -h                       : displays this message";
const std::string& d = "-d";
const std::string& t = "-t";
const std::string& f = "-f";
const std::string& e = "-e";
const std::string& o = "-o";
#if SPDLOG_ACTIVE_LEVEL == SPDLOG_LEVEL_TRACE
const std::string& l = "-l";
#endif
const std::string& h = "-h";
const char* GRAV = "GRAV";
const char* LJ = "LJ";

void cliParse(int argc, char** argsv, SimpleContainer& particles, SettingsParam& settings) {
    std::cout << "Hello from MolSim for PSE!" << '\n';
    char** help = std::find(argsv, argsv + argc, h);
    if (help != &argsv[argc]) {
        std::cout << HELP_MSG << '\n';
        exit(-1);
    }
    if (argc == 0) {
        std::cout << "Erroneous programme call! " << '\n' << HELP_MSG << '\n';
        exit(-2);
    }

    int parsed_args = 2;  // program name + assume file name is OK (bad files handled in fileReader.readfile)
    char** delta_t_opt = std::find(argsv, &argsv[argc], d);
    char** end_time_opt = std::find(argsv, &argsv[argc], t);
    char** force_opt = std::find(argsv, &argsv[argc], f);
    char** epsilon_opt = std::find(argsv, &argsv[argc], e);
    char** sigma_opt = std::find(argsv, &argsv[argc], o);
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
        if (force_opt != &argsv[argc]) {
            char* force_string = *(++force_opt);
            bool grav = !static_cast<bool>(std::strcmp(force_string, GRAV));
            bool lj = !static_cast<bool>(std::strcmp(force_string, LJ));

            if (grav) {
                parsed_args += 2;
            } else if (lj) {
                parsed_args += 2;
                if (epsilon_opt != &argsv[argc]) {
                    settings.epsilon = std::stod(*(++epsilon_opt));
                    parsed_args += 2;
                }
                if (sigma_opt != &argsv[argc]) {
                    settings.sigma = std::stod(*(++sigma_opt));
                    parsed_args += 2;
                }
            }
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
        SPDLOG_ERROR("Erroneous programme call!: Given options must be valid floating point numbers or strings!");
        std::cout << HELP_MSG << '\n';
        exit(-1);
    } catch (std::out_of_range& e) {
        SPDLOG_ERROR(
            "Erroneous programme call!: At least one floating point out of double precision floating point range!");
        std::cout << HELP_MSG << '\n';
        exit(-1);
    } catch (std::logic_error& e) {
        SPDLOG_ERROR("Erroneous programme call!: Flag set but no value provided!");
        std::cout << HELP_MSG << '\n';
        exit(-1);
    }

    if (parsed_args - argc != 0) {
        SPDLOG_WARN("Erroneous programme call!: Unrecognized arguments in programme call!");
    }
    std::unique_ptr<FileReader> file_reader;

    std::filesystem::path path = argsv[1];
    if (path.extension() == ".txt") {
        file_reader = std::make_unique<XVMReader>();
    } else if (path.extension() == ".yaml") {
        file_reader = std::make_unique<YAMLReader>();
    } else {
        SPDLOG_ERROR("Unsupported File Extension!");
        exit(-1);
    }

    file_reader->readFile(particles, settings, argsv[1]);
}

}  // namespace mol_sim

#endif
