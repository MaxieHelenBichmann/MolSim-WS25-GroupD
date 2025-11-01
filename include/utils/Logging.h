#ifndef LOGGING_H
#define LOGGING_H

#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

namespace mol_sim {

inline void logInit() {
    auto logger = spdlog::default_logger();
    if (!logger) {
        logger = spdlog::stdout_color_mt("default");
        spdlog::set_default_logger(logger);
    }

    spdlog::level::level_enum level =
#ifdef LOG_LEVEL_TRACE
        spdlog::level::trace;
#elif defined(LOG_LEVEL_DEBUG)
        spdlog::level::debug;
#elif defined(LOG_LEVEL_INFO)
        spdlog::level::info;
#elif defined(LOG_LEVEL_ERROR)
        spdlog::level::err;
#else
        spdlog::level::info;
#endif

    logger->set_level(level);
}

}  // namespace mol_sim

#define LOG_TRACE(...) spdlog::trace(__VA_ARGS__)  // NOLINT(cppcoreguidelines-macro-usage)
#define LOG_DEBUG(...) spdlog::debug(__VA_ARGS__)  // NOLINT(cppcoreguidelines-macro-usage)
#define LOG_INFO(...) spdlog::info(__VA_ARGS__)    // NOLINT(cppcoreguidelines-macro-usage)
#define LOG_ERROR(...) spdlog::error(__VA_ARGS__)  // NOLINT(cppcoreguidelines-macro-usage)

#endif  // LOGGING_H
