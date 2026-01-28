#ifndef SIMULATION_EXCEPTION_H
#define SIMULATION_EXCEPTION_H

#include <string>

#include "exceptions/MolSimException.h"

namespace mol_sim {
/**
 * @brief Exception thrown by the Program.
 *
 */
class SimulationException : public MolSimException {
   public:
    explicit SimulationException(const std::string& message) : MolSimException(message) {}
};
}  // namespace mol_sim
#endif
