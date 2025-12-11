/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#ifndef FILE_READER_H
#define FILE_READER_H

#include <string>

#include "particles/container/ContainerRef.h"
#include "utils/Settings.h"

namespace mol_sim {

/**
 * @brief Interface for file readers.
 *
 * Interface for file readers.
 * These classes implement two-phase reading:
 * 1. parseSettings - reads configuration/settings only
 * 2. readParticles - reads particles into a container
 */
class FileReader {
   public:
    virtual ~FileReader() = default;

    /**
     * @brief Phase 1: Reads settings/configuration from the file.
     *
     * @param settings SettingsParam where read settings are stored.
     * @param filename The path to the input file.
     */
    virtual void readSettings(SettingsParam& settings, const std::string& filename) = 0;

    /**
     * @brief Phase 2: Reads particles from the file into a container.
     *
     * @param particles The ParticleContainer to store particles in.
     * @param filename The path to the input file.
     */
    virtual void readParticles(ContainerRef particles, const SettingsParam& settings, const std::string& filename) = 0;
};

}  // namespace mol_sim

#endif
