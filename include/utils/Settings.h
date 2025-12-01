#ifndef SETTINGS_H
#define SETTINGS_H

#include <cmath>
#include <cstddef>
#include <limits>
#include <string>

#include "particles/boundaries/Boundary.h"
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
    static constexpr size_t FREQUENCY_DEFAULT = 10;
    /**
     * @brief Default cutoff radius for the linked cells.
     */
    static constexpr double CUTOFF_DEFAULT = std::numeric_limits<double>::infinity();

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
    size_t frequency = FREQUENCY_DEFAULT;
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
     * @brief Construct new SettingsParam with default values.
     */
    SettingsParam() = default;
};
}  // namespace mol_sim
#endif
