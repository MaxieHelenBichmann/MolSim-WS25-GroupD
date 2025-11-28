#ifndef SETTINGS_H
#define SETTINGS_H

#include <optional>
#include <string>
#include <utility>
#include <variant>

#include "particles/Particle.h"
#include "particles/boundaries/Boundary.h"
#include "particles/container/LinkedCellContainer.h"
#include "particles/container/SimpleContainer.h"
#include "particles/container/domain/Domain.h"
#include "physics/ForceSource.h"
#include "physics/LennardJonesForce.h"

namespace mol_sim {
/**
 * @brief Provides a wrapper for settings of the simulation set during config
 *
 */
class SettingsParam {
   public:
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
    /**
     * @brief Default base string of the output files of the simulation.
     *
     */
    constexpr static std::string NAME_DEFAULT = "MD";
    /**
     * @brief Default force type of the simulation.
     *
     */
    constexpr static Force FORCE_DEFAULT = LENNARDJONES;
    /**
     * @brief Default write frequency of the simulation.
     *
     */
    constexpr static size_t FREQUENCY_DEFAULT = 10;
    /**
     * @brief Default cutoff radius for the linked cells.
     *
     */
    constexpr static double CUTOFF_DEFAULT = 0.5;
    /**
     * @brief Default container type. SIMPLE = SimpleContainer
     *
     */
    constexpr static std::string CONTAINER_TYPE_DEFAULT = "SIMPLE";
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
    std::optional<double> epsilon;
    /**
     * @brief Force type used in the simulation.
     *
     */
    std::optional<double> sigma;
    /**
     * @brief Base name of the output files
     *
     */
    std::optional<std::string> base_name;
    /**
     * @brief Type of force used in the simulation
     *
     */
    std::optional<Force> force;
    /**
     * @brief Frequency of output files being written
     * All *frequency* iterations file is written
     */
    std::optional<size_t> frequency;
    /**
     * @brief Cutoff of the linked cells algorithm
     *
     */
    std::optional<double> cutoff;
    /**
     * @brief The type in string format (as it would be expected in the .yaml files) of the particle container.
     *
     * SIMPLE = SimpleContainer
     * LINKED = LinkedCellContainer
     *
     * @note May be deprecated depending on future changes with typing.
     */
    std::optional<std::string> container_type;
    /**
     * @brief The domain of the simulation
     */
    std::optional<Domain> domain;
    /**
     * @brief Construct new SettingsParam.
     * All values will be set to nullopt if not specified otherwise.
     * Default values will be set in YAMLReader.cpp::readFile
     */
    SettingsParam() = default;
    /**
     * @brief Provide default values for settings that have not been set.
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
        if (!epsilon.has_value()) {
            epsilon = EPSILON_DEFAULT;
        }
        if (!sigma.has_value()) {
            sigma = SIGMA_DEFAULT;
        }
        if (!base_name.has_value()) {
            base_name = NAME_DEFAULT;
        }
        if (!force.has_value()) {
            force = FORCE_DEFAULT;
        }
        if (!frequency.has_value()) {
            frequency = FREQUENCY_DEFAULT;
        }
        if (!cutoff.has_value()) {
            cutoff = CUTOFF_DEFAULT;
        }
        if (!container_type.has_value()) {
            container_type = CONTAINER_TYPE_DEFAULT;
        }
        if (!domain.has_value()) {
            domain.emplace();
        }
    }
};
}  // namespace mol_sim
#endif
