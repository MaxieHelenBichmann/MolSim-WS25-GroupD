#ifndef CONTAINER_EXCEPTION_H
#define CONTAINER_EXCEPTION_H

#include <string>

#include "exceptions/MolSimException.h"

namespace mol_sim {

/**
 * @brief Exception thrown by particle containers for container-related errors.
 *
 */
class ContainerException : public MolSimException {
   public:
    explicit ContainerException(const std::string& message) : MolSimException(message) {}
};
}  // namespace mol_sim

#endif
