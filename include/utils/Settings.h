#ifndef SETTINGS_H
#define SETTINGS_H

#include <optional>

namespace mol_sim {
/**
 * @brief Provides a wrapper for settings of the simulation set during config
 *
 */
class SettingsParam {
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
     * @brief Construct new SettingsParam.
     * All values will be set to null_opt if not specified otherwise.
     * Default Values will be set in MolSim.cpp
     * @param delta_t
     * @param start_time
     * @param end_time
     */
    SettingsParam(std::optional<double> delta_t = std::nullopt, std::optional<double> start_time = std::nullopt,
                  std::optional<double> end_time = std::nullopt)
        : delta_t(delta_t), start_time(start_time), end_time(end_time) {}
};
}  // namespace mol_sim
#endif
