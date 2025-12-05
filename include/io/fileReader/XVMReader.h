#ifndef XVM_READER_H
#define XVM_READER_H

#include "io/FileReader.h"
#include "particles/container/ContainerRef.h"
#include "utils/Settings.h"

namespace mol_sim {

/**
 * @brief Reads the particles from a file that contains the initial position, force, and mass of particles.
 *
 * Reads the particles from a file that is formatted like the following regex:
 *
 * (#comment\n || empty line\n)*
 * (number of particles\n)
 * ((x y z)(whitespace)(x y z)(whitespace)(m)\n)*
 *
 * where the first (x y z) tuple represents the initial position of the particle
 * the second (x y z) tuple represents the initial velocity of the particle
 * and m represents the mass of the particle.
 * All these values are floats.
 */
class XVMReader : public FileReader {
   public:
    XVMReader();
    ~XVMReader() override;

    /**
     * @brief Phase 1: XVM files have no settings, just sets defaults.
     * @param settings SettingsParam to populate (unchanged for XVM files).
     * @param filename Path to input file (unused).
     */
    void readSettings(SettingsParam& settings, const std::string& filename) override;

    /**
     * @brief Phase 2: Reads particles from the XVM file.
     * @param particles Container to place particles in.
     * @param filename Path to input file.
     * @throws XVMReaderException if file cannot be opened or is malformed.
     * @throws ValidationException if particle mass is invalid.
     */
    void readParticles(ContainerRef particles, const SettingsParam& settings, const std::string& filename) override;
};

}  // namespace mol_sim

#endif
