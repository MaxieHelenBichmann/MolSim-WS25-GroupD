#ifndef CLI_PARSE_H
#define CLI_PARSE_H

#include <spdlog/spdlog.h>

#include <CLI/CLI.hpp>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "exceptions/CLIException.h"
#include "io/FileReader.h"
#include "utils/Logging.h"

namespace mol_sim {

/**
 * @brief Parses command line input.
 *
 * Parses command line input and validates parameters.
 * @param argc Number of arguments including the program name
 * @param argsv Array of arguments
 * @throws CLIException if CLI parsing or file reading fails
 */
std::vector<std::string> cliParse(int argc, char** argv) {
    SPDLOG_INFO("Hello from MolSim for PSE!");
    CLI::App app{"MolSim - Molecular Dynamics Simulator"};
    argv = app.ensure_utf8(argv);
    std::vector<std::string> files;

    // Customize help formatting
    app.get_formatter()->column_width(40);
    app.get_formatter()->label("REQUIRED", "");
    app.get_formatter()->label("TEXT", "");
    std::unique_ptr<FileReader> file_reader;
    std::filesystem::path filepath;
    std::string log_level = "Default";  // NOLINT

    // Define custom validator for file extensions
    auto file_ext_validator = [](const std::string& filename) -> std::string {
        std::filesystem::path p(filename);
        auto ext = p.extension();
        if (ext != ".txt" && ext != ".yaml") {
            return "File must have .txt or .yaml extension";
        }
        return "";
    };

    // Add files with least settings priority

    app.add_option("filepath,-f,--file", filepath, "Input file path (.txt or .yaml)")
        ->required()
        ->check(CLI::ExistingFile.description(""))
        ->check(CLI::Validator(file_ext_validator, ""));
    files.insert(files.begin(), filepath);

    app.add_option("checkpoint,-c,--checkpoint", filepath, "Input file path for checkpoint file (.txt or .yaml)")
        ->check(CLI::ExistingFile.description(""))
        ->check(CLI::Validator(file_ext_validator, ""));
    files.insert(files.begin(), filepath);

#if SPDLOG_ACTIVE_LEVEL == SPDLOG_LEVEL_TRACE
    app.add_option("-l,--log_level", log_level, "Logging verbosity level")
        ->check(CLI::IsMember({"Info", "Trace", "Debug"}).description("{Info,Debug,Trace}"));
#endif

    // Parse with error handling
    try {
        app.parse(argc, argv);
    } catch (const CLI::Success& e) {
        // User requested help or version - exit cleanly
        exit(0);
    } catch (const CLI::RequiredError& e) {
        SPDLOG_ERROR("Missing required argument: {}", e.what());
        std::cout << app.help() << '\n';
        throw CLIException("Required argument missing: " + std::string(e.what()));
    } catch (const CLI::ValidationError& e) {
        SPDLOG_ERROR("Validation failed: {}", e.what());
        std::cout << app.help() << '\n';
        throw CLIException("Validation error: " + std::string(e.what()));
    } catch (const CLI::FileError& e) {
        SPDLOG_ERROR("File error: {}", e.what());
        throw CLIException("File error: " + std::string(e.what()));
    } catch (const CLI::ConversionError& e) {
        SPDLOG_ERROR("Invalid value provided: {}", e.what());
        std::cout << app.help() << '\n';
        throw CLIException("Invalid value: " + std::string(e.what()));
    } catch (const CLI::ParseError& e) {
        SPDLOG_ERROR("Command line parsing failed: {}", e.what());
        std::cout << app.help() << '\n';
        throw CLIException("CLI parsing error: " + std::string(e.what()));
    }

#if SPDLOG_ACTIVE_LEVEL == SPDLOG_LEVEL_TRACE
    logInit(log_level);
#endif

    return files;
}
}  // namespace mol_sim

#endif
