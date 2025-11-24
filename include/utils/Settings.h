#ifndef SETTINGS_H
#define SETTINGS_H

#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <variant>

#include "particles/Particle.h"
#include "particles/boundaries/BoundaryCondition.h"
#include "particles/container/LinkedCellContainer.h"
#include "particles/container/domain/Domain.h"
#include "physics/ForceSource.h"
#include "physics/LennardJonesForce.h"

namespace mol_sim {
/**
 * @brief Provides a wrapper for settings of the simulation set during config
 *
 */
using DomainVariant =
    std::variant<Domain<SimpleContainer>, Domain<LinkedCellContainer>>;  // if more container types are needed, just add
                                                                         // another Domain<containerType> here.

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

    DomainVariant domain;
    /**
     * @brief Construct new SettingsParam.
     * All values will be set to nullopt if not specified otherwise.
     * Default Values will be set in MolSim.cpp
     * @param delta_t
     * @param start_time
     * @param end_time
     * @param epsilon
     * @param sigma
     * @param base_name
     * @param force
     * @param frequency
     * @param cutoff
     */
    SettingsParam(std::optional<double> delta_t = std::nullopt, std::optional<double> start_time = std::nullopt,
                  std::optional<double> end_time = std::nullopt, std::optional<double> epsilon = std::nullopt,
                  std::optional<double> sigma = std::nullopt, std::optional<std::string> base_name = std::nullopt,
                  std::optional<Force> force = std::nullopt, std::optional<size_t> frequency = std::nullopt,
                  std::optional<double> cutoff = std::nullopt, DomainVariant domain = DomainVariant())
        : delta_t(delta_t),
          start_time(start_time),
          end_time(end_time),
          epsilon(epsilon),
          sigma(sigma),
          base_name(std::move(base_name)),
          force(force),
          frequency(frequency),
          cutoff(cutoff),
          domain(std::move(domain)) {}
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
    }
};
}  // namespace mol_sim
#endif
