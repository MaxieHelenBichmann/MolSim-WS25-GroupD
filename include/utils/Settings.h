#ifndef SETTINGS_H
#define SETTINGS_H

#include <cstddef>
#include <limits>
#include <string>

#include "particles/container/domain/Domain.h"
#include "physics/ForceSource.h"

namespace mol_sim {
/**
 * @brief Provides a wrapper for settings of the simulation set during config
 *
 */
class SettingsParam {
   public:
    /**
     * @brief Default value for the delta_t parameter of the simulation.
     */
    static constexpr double DELTA_T_DEFAULT = 0.014;
    /**
     * @brief Default value for the end_time parameter of the simulation.
     */
    static constexpr double END_TIME_DEFAULT = 1000;
    /**
     * @brief Default value for the start_time parameter of the simulation.
     */
    static constexpr double START_TIME_DEFAULT = 0;
    /**
     * @brief Default value for the epsilon parameter of the simulation if
     * the force_type is LJ (Lennard-Jones).
     */
    static constexpr double EPSILON_DEFAULT = 5;
    /**
     * @brief Default value for the sigma parameter of the simulation if
     * the force_type is LJ (Lennard-Jones).
     */
    static constexpr double SIGMA_DEFAULT = 1;
    /**
     * @brief Default force type of the simulation.
     */
    static constexpr Force FORCE_DEFAULT = LENNARDJONES;
    /**
     * @brief Default write frequency of the simulation.
     */
    static constexpr size_t FREQUENCY_OUTPUT_DEFAULT = 10;
    /**
     * @brief Default checkpoint frequency of the simulation.
     */
    static constexpr size_t FREQUENCY_CHECKPOINT_DEFAULT = 5000;
    /**
     * @brief Default cutoff radius for the linked cells.
     */
    static constexpr double CUTOFF_DEFAULT = std::numeric_limits<double>::infinity();
    /**
     * @brief Default initial temperature for the system.
     */
    static constexpr double INIT_TEMP_DEFAULT = 40.;
    /**
     * @brief Default delta temperature of the simulation.
     */
    static constexpr double DELTA_TEMP_DEFAULT = std::numeric_limits<double>::infinity();
    /**
     * @brief Default thermostat frequency of the simulation.
     */
    static constexpr size_t THERMOSTAT_FREQ_DEFAULT = 1000;
    /**
     * @brief Default dimensions of the simulation..
     */
    static constexpr size_t DIMENSIONS_DEFAULT = 3;
    /**
     * @brief Default value for g_grav used to simulate gravitational pull in simulation.
     */
    static constexpr double G_GRAV_DEFAULT = 0.0;
    /**
     * @brief delta_t of the simulation.
     */
    double delta_t = DELTA_T_DEFAULT;
    /**
     * @brief Start time of the simulation.
     */
    double start_time = START_TIME_DEFAULT;
    /**
     * @brief End time of the simulation.
     */
    double end_time = END_TIME_DEFAULT;
    /**
     * @brief Epsilon parameter for Lennard-Jones force.
     */
    double epsilon = EPSILON_DEFAULT;
    /**
     * @brief Sigma parameter for Lennard-Jones force.
     */
    double sigma = SIGMA_DEFAULT;
    /**
     * @brief Base name of the output files.
     */
    std::string base_name = "MD";
    /**
     * @brief Type of force used in the simulation.
     */
    Force force = FORCE_DEFAULT;
    /**
     * @brief Frequency of output files being written.
     * Output is written every *frequency* iterations.
     */
    size_t frequency_output = FREQUENCY_OUTPUT_DEFAULT;
    /**
     * @brief Frequency of checkpoint files being written.
     * Checkpoints are written every *frequency* iterations.
     */
    size_t frequency_checkpoint = FREQUENCY_CHECKPOINT_DEFAULT;
    /**
     * @brief Dimensions of the simulation.
     * TODO: use this for 2D optimizations.
     */
    size_t dimensions = DIMENSIONS_DEFAULT;
    /**
     * @brief Cutoff radius for the linked cells algorithm.
     */
    double cutoff = CUTOFF_DEFAULT;
    /**
     * @brief The type in string format of the particle container.
     *
     * SIMPLE = SimpleContainer
     * LINKED = LinkedCellContainer
     */
    std::string container_type = "SIMPLE";
    /**
     * @brief The domain of the simulation.
     */
    Domain domain;

    /**
     * @brief Toggles if thermostat is enabled for this simulation.
     *
     */
    bool thermo = false;

    /**
     * @brief Target initial temperature of the system.
     */
    double init_temp = INIT_TEMP_DEFAULT;
    /**
     * @brief General target temperature of the system.
     */
    double target_temp = init_temp;
    /**
     * @brief Frequency with which the thermostat gets applied.
     */
    size_t thermostat_freq = THERMOSTAT_FREQ_DEFAULT;
    /**
     * @brief Maximum allowed temperature change of the system with one thermostat application.
     */
    double delta_temp = DELTA_TEMP_DEFAULT;
    /**
     * @brief Used to simulate of gravitational pull in simulation.
     */
    double g_grav = G_GRAV_DEFAULT;

    /**
     * @brief Construct new SettingsParam with default values.
     */
    SettingsParam() = default;
};
}  // namespace mol_sim
#endif
