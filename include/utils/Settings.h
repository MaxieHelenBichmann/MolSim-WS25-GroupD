#ifndef SETTINGS_H
#define SETTINGS_H

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "particles/container/domain/Domain.h"
#include "physics/pairwiseforces/PairwiseForceSource.h"
#include "physics/singleforces/SingleForceSource.h"

namespace mol_sim {

/**
 * @brief Parallelization strategy for force calculation.
 */
enum class ParallelizationStrategy : uint8_t { NAIVE, COLORING };

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
    inline static const std::vector<PairwiseForce> PAIRWISE_FORCE_DEFAULT = {LENNARDJONES};
    inline static const std::vector<SingleForce> SINGLE_FORCE_DEFAULT = {};
    /**
     * @brief Default write frequency of the simulation.
     */
    static constexpr size_t FREQUENCY_OUTPUT_DEFAULT = 10;
    /**
     * @brief Default checkpoint frequency of the simulation.
     */
    static constexpr size_t FREQUENCY_CHECKPOINT_DEFAULT = 5000;
    /**
     * @brief Default frequency for thermodynamical statistics (Diffusion and Temperature).
     */
    static constexpr size_t FREQUENCY_STATS_DEFAULT_DIFF = 1000;
    /**
     * @brief Default frequency for thermodynamical statistics (RDF).
     */
    static constexpr size_t FREQUENCY_STATS_DEFAULT_RDF = 7000;
    /**
     * @brief Default sample width for the radial distribution function.
     */
    static constexpr double RDF_DEFAULT = 1.0;
    /**
     * @brief Default window size for the radial distribution function.
     */
    static constexpr double RDF_SIZE_DEFAULT = 10.0;
    /**
     * @brief Default cutoff radius for the linked cells.
     */
    static constexpr double CUTOFF_DEFAULT = std::numeric_limits<double>::infinity();
    /**
     * @brief Default smoothing radius for the linked cells.
     */
    static constexpr double SMOOTHING_DEFAULT = std::numeric_limits<double>::infinity();
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
     * @brief Default stiffness constant of the simulation..
     */
    static constexpr double K_DEFAULT = 300.0;

    /**
     * @brief Default average bond length of the simulation..
     */
    static constexpr double R_0_DEFAULT = 2.2;

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
     * @note: could use this for 2D optimizations.
     */
    size_t dimensions = DIMENSIONS_DEFAULT;

    /**
     * @brief Cutoff radius for the linked cells algorithm.
     */
    double cutoff = CUTOFF_DEFAULT;

    /**
     * @brief Smoothing radius for the linked cells algorithm.
     */
    double smoothing = SMOOTHING_DEFAULT;

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
     * @brief Frequency for thermodynamical statistics (Diffusion).
     */
    size_t stats_freq_diffusion = FREQUENCY_STATS_DEFAULT_DIFF;

    /**
     * @brief Frequency for thermodynamical statistics (RDF).
     */
    size_t stats_freq_rdf = FREQUENCY_STATS_DEFAULT_RDF;

    /**
     * @brief Sample width for the radial distribution function.
     */
    double sample_radius = RDF_DEFAULT;

    /**
     * @brief Size of the window for the radial distribution function.
     */
    double window_size = RDF_SIZE_DEFAULT;

    /**
     * @brief Maximum allowed temperature change of the system with one thermostat application.
     */
    double delta_temp = DELTA_TEMP_DEFAULT;

    /**
     * @brief Stiffness constant K for Harmonic forces in particle membranes.
     *
     */
    double k = K_DEFAULT;

    /**
     * @brief Average bond length of a molecule pair in particle membranes.
     *
     */
    double r_0 = R_0_DEFAULT;

    /**
     * @brief Target force magnitude.
     */
    double target_force_magnitude = 0.0;

    /**
     * @brief Maximum iterations for target force.
     */
    size_t target_force_max_iterations = 0;

    /**
     * @brief Types of pairwise forces used this simulation
     */
    std::vector<PairwiseForce> pairwise_forces;

    /**
     * @brief Types of single forces used this simulation
     */
    std::vector<SingleForce> single_forces = SINGLE_FORCE_DEFAULT;

    /**
     * @brief g_grav vector for GRAV single force (gravitational acceleration vector).
     */
    R3 g_grav_vec = {0.0, -9.81, 0.0};

    /**
     * @brief Target force direction vector.
     */
    R3 target_force_direction = {0.0, 1.0, 0.0};

    /**
     * @brief The type in string format of the particle container.
     *
     * SIMPLE = SimpleContainer
     * LINKED = LinkedCellContainer
     */
    std::string container_type = "LINKED";

    /**
     * @brief Base name of the output files of the Simulation.
     *
     */
    std::string base_name = "MD";

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
     * @brief Toggles if data collection for the radial distribution function is enabled for this simulation.
     *
     */
    bool rdf = false;

    /**
     * @brief Toggles if data collection for the diffusion is enabled for this simulation.
     *
     */
    bool diff = false;

    /**
     * @brief Toggles particle generation with brownian motion
     *
     */
    bool brownian = true;

    /**
     * @brief Enable target force.
     */
    bool target_force_enabled = false;

    /**
     * @brief Parallelization strategy for force calculation.
     */
    ParallelizationStrategy strategy = ParallelizationStrategy::NAIVE;

    /**
     * @brief Construct new SettingsParam with default values.
     */
    SettingsParam() = default;
};
}  // namespace mol_sim
#endif
