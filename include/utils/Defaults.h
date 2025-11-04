#ifndef DEFAULTS_H
#define DEFAULTS_H

#include "physics/ForceSource.h"
#include "utils/Settings.h"

namespace mol_sim {
/**
 * @brief Default value for the delta_t parameter of the simulation.
 *
 */
constexpr double DELTA_T_DEFAULT = 0.014;
/**
 * @brief Default value for the end_time parameter of the simulation.
 *
 */
constexpr double END_TIME_DEFAULT = 1000;
/**
 * @brief Default value for the start_time parameter of the simulation.
 *
 */
constexpr double START_TIME_DEFAULT = 0;
/**
 * @brief Default value for the start_time parameter of the simulation.
 *
 */
constexpr Force FORCE_TYPE_DEFAULT = GRAVITATIONAL;
/**
 * @brief Default value for the epsilon parameter of the simulation if 
 * the force_type is LJ (Lennard-Jones).
 *
 */
constexpr double EPSILON_DEFAULT = 5;
/**
 * @brief Default value for the sigma parameter of the simulation if 
 * the force_type is LJ (Lennard-Jones).
 *
 */
constexpr double SIGMA_DEFAULT = 1;


/**
 * @brief Provide default values for settings that have not been set.
 * @param settings A SettingsParam of the settings that were made.
 *
 * */
void setDefaults(SettingsParam& settings) {
    if (!settings.delta_t.has_value()) {
        settings.delta_t = DELTA_T_DEFAULT;
    }
    if (!settings.end_time.has_value()) {
        settings.end_time = END_TIME_DEFAULT;
    }
    if (!settings.start_time.has_value()) {
        settings.start_time = START_TIME_DEFAULT;
    }
    if (!settings.force_type.has_value()) {
        settings.force_type = FORCE_TYPE_DEFAULT;
    }
    if (!settings.start_time.has_value()) {
        settings.epsilon = EPSILON_DEFAULT;
    }
    if (!settings.start_time.has_value()) {
        settings.sigma = SIGMA_DEFAULT;
    }
}

} // namespace mol_sim

#endif
