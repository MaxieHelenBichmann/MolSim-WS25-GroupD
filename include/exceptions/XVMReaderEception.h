#ifndef XVM_READER_EXCEPTION_H
#define XVM_READER_EXCEPTION_H

#include <string>

#include "exceptions/InputException.h"

namespace mol_sim {

/**
 * @brief Exception thrown by the YAMLReader;
 *
 */
class XVMReaderException : public InputException {
   public:
    explicit XVMReaderException(const std::string& message) : InputException(message) {};
};
}  // namespace mol_sim

#endif
