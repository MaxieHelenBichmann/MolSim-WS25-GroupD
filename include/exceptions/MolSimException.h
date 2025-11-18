#ifndef MOLSIM_EXCEPTION_H
#define MOLSIM_EXCEPTION_H

#include <stdexcept>
#include <string>
namespace mol_sim {
/**
 * @brief Exception thrown by the Program;
 *
 */
class MolSimException : public std::runtime_error {
   public:
    MolSimException(const std::string& message) : std::runtime_error(message) {}
};
}  // namespace mol_sim
#endif
