#ifndef LOGGING_H
#define LOGGING_H

#ifndef LOG_LEVEL_DEBUG
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_ERROR  // NOLINT
#else
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE  // NOLINT
#endif
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

namespace mol_sim {

/**
 * @brief Initialize the global logger.
 *
 * Ensures that a colored stdout logger exists, installs it as the default
 * logger if necessary, and applies the runtime-time logging level configured via CLI
 * when not building in Release or MinSizeRel.
 */
inline void logInit(const std::string& option) {
    if (option == "Info" || option == "Default") {
        spdlog::set_level(spdlog::level::info);
    } else if (option == "Debug") {
        spdlog::set_level(spdlog::level::debug);
    } else if (option == "Trace") {
        spdlog::set_level(spdlog::level::trace);
    } else {
        throw std::invalid_argument{"Invalid Log Level"};
    }
};

}  // namespace mol_sim

#endif  // LOGGING_H
