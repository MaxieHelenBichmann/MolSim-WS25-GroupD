#ifndef LOGGING_H
#define LOGGING_H

#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

/**
 * @brief Log a trace-level message.
 *
 * Accepts a format string followed by optional arguments.
 */
#define LOG_TRACE(...) spdlog::trace(__VA_ARGS__)  // NOLINT(cppcoreguidelines-macro-usage)

/**
 * @brief Log a debug-level message.
 *
 * Accepts a format string followed by optional arguments.
 */
#define LOG_DEBUG(...) spdlog::debug(__VA_ARGS__)  // NOLINT(cppcoreguidelines-macro-usage)

/**
 * @brief Log a info-level message.
 *
 * Accepts a format string followed by optional arguments.
 */
#define LOG_INFO(...) spdlog::info(__VA_ARGS__)  // NOLINT(cppcoreguidelines-macro-usage)

/**
 * @brief Log a error-level message.
 *
 * Accepts a format string followed by optional arguments.
 */
#define LOG_ERROR(...) spdlog::error(__VA_ARGS__)  // NOLINT(cppcoreguidelines-macro-usage)

namespace mol_sim {

/**
 * @brief Initialize the global logger.
 *
 * Ensures that a colored stdout logger exists, installs it as the default
 * logger if necessary, and applies the compile-time logging level configured via
 * `LOG_LEVEL_DEBUG` (Debug), `LOG_LEVEL_TRACE` (RelWithDebInfo), `LOG_LEVEL_INFO` (Default), or `LOG_LEVEL_ERROR`
 * (Release, MinSizeRel).
 */
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
};

}  // namespace mol_sim

#endif  // LOGGING_H
