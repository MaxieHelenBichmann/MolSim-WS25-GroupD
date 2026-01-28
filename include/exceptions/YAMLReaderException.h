#ifndef YAML_READER_EXCEPTION_H
#define YAML_READER_EXCEPTION_H

#include <string>

#include "exceptions/InputException.h"

namespace mol_sim {

/**
 * @brief Exception thrown by the YAMLReader.
 *
 */
class YAMLReaderException : public InputException {
   public:
    explicit YAMLReaderException(const std::string& message) : InputException(message) {};
};
}  // namespace mol_sim

#endif
