#ifndef CLI_EXCEPTION_H
#define CLI_EXCEPTION_H

#include <string>

#include "exceptions/InputException.h"

namespace mol_sim {

/**
 * @brief Exception thrown by the CLI;
 *
 */
class CLIException : public InputException {
   public:
    explicit CLIException(const std::string& message) : InputException(message) {}
};
}  // namespace mol_sim
#endif
