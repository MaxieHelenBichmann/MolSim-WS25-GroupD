#ifndef SETTINGS_H
#define SETTINGS_H

#include <optional>

#include "physics/ForceSource.h"

namespace mol_sim {
/**
 * @brief Provides a wrapper for settings of the simulation set during config
 *
 */
class SettingsParam {
    /**
     * @brief Default value for the delta_t parameter of the simulation.
     *
     */
    constexpr static double DELTA_T_DEFAULT = 0.014;
    /**
     * @brief Default value for the end_time parameter of the simulation.
     *
     */
    constexpr static double END_TIME_DEFAULT = 1000;
    /**
     * @brief Default value for the start_time parameter of the simulation.
     *
     */
    constexpr static double START_TIME_DEFAULT = 0;
    /**
     * @brief Default value for the start_time parameter of the simulation.
     *
     */
    constexpr static Force FORCE_TYPE_DEFAULT = LENNARDJONES;
    /**
     * @brief Default value for the epsilon parameter of the simulation if
     * the force_type is LJ (Lennard-Jones).
     *
     */
    constexpr static double EPSILON_DEFAULT = 5;
    /**
     * @brief Default value for the sigma parameter of the simulation if
     * the force_type is LJ (Lennard-Jones).
     *
     */
    constexpr static double SIGMA_DEFAULT = 1;

   public:
    /**
     * @brief delta_t of the simulation.
     *
     */
    std::optional<double> delta_t;
    /**
     * @brief Start time of the simulation.
     *
     */
    std::optional<double> start_time;
    /**
     * @brief End time of the simulation.
     *
     */
    std::optional<double> end_time;
    /**
     * @brief Force type used in the simulation.
     *
     */
    std::optional<Force> force_type;
    /**
     * @brief Force type used in the simulation.
     *
     */
    std::optional<double> epsilon;
    /**
     * @brief Force type used in the simulation.
     *
     */
    std::optional<double> sigma;
    /**
     * @brief Construct new SettingsParam.
     * All values will be set to null_opt if not specified otherwise.
     * Default Values will be set in MolSim.cpp
     * @param delta_t
     * @param start_time
     * @param end_time
     */
    SettingsParam(std::optional<double> delta_t = std::nullopt, std::optional<double> start_time = std::nullopt,
                  std::optional<double> end_time = std::nullopt, std::optional<Force> force_type = std::nullopt,
                  std::optional<double> epsilon = std::nullopt, std::optional<double> sigma = std::nullopt)
        : delta_t(delta_t),
          start_time(start_time),
          end_time(end_time),
          force_type(force_type),
          epsilon(epsilon),
          sigma(sigma) {}

    /**
     * @brief Provide default values for settings that have not been set.
     * @param settings A SettingsParam of the settings that were made.
     *
     * */

    void setDefaults() {
        if (!delta_t.has_value()) {
            delta_t = DELTA_T_DEFAULT;
        }
        if (!end_time.has_value()) {
            end_time = END_TIME_DEFAULT;
        }
        if (!start_time.has_value()) {
            start_time = START_TIME_DEFAULT;
        }
        if (!force_type.has_value()) {
            force_type = FORCE_TYPE_DEFAULT;
        }
        if (!start_time.has_value()) {
            epsilon = EPSILON_DEFAULT;
        }
        if (!start_time.has_value()) {
            sigma = SIGMA_DEFAULT;
        }
    }
};
}  // namespace mol_sim
#endif
