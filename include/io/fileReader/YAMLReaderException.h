#ifndef YAMLREADEREXCEPTION_H
#define YAMLREADEREXCEPTION_H

#include <stdexcept>
#include <string>

namespace mol_sim {

class YAMLReaderException : public std::runtime_error {
public:
    explicit YAMLReaderException(const std::string& message)
        : std::runtime_error(message) {}
};

} // namespace mol_sim

#endif // YAMLREADEREXCEPTION_H
