#ifndef BOUNDARY_EXCEPTION_H
#define BOUNDARY_EXCEPTION_H

#include <string>

#include "exceptions/SimulationException.h"

namespace mol_sim {

/**
 * @brief Exception thrown for boundary-related errors during simulation.
 *
 */
class BoundaryException : public SimulationException {
   public:
    explicit BoundaryException(const std::string& message) : SimulationException(message) {}
};
}  // namespace mol_sim

#endif
