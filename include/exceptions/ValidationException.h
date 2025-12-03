#ifndef VALIDATION_EXCEPTION_H
#define VALIDATION_EXCEPTION_H

#include <string>

#include "exceptions/InputException.h"

namespace mol_sim {

/**
 * @brief Exception thrown for input validation errors.
 *
 */
class ValidationException : public InputException {
   public:
    explicit ValidationException(const std::string& message) : InputException(message) {}
};
}  // namespace mol_sim

#endif
