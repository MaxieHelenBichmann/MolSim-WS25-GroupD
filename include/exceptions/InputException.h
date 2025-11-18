
#ifndef INPUT_EXCEPTION_H
#define INPUT_EXCEPTION_H

#include <string>

#include "exceptions/MolSimException.h"
namespace mol_sim {

/**
 * @brief Exception thrown by Input handling;
 *
 */
class InputException : public MolSimException {
   public:
    InputException(const std::string& message) : MolSimException(message) {}
};
}  // namespace mol_sim
#endif
