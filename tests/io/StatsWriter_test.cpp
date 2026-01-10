#include "io/StatsWriter.h"

#include <gtest/gtest.h>
#include <particles/Particle.h>

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

#include "io/fileReader/XVMReader.h"
#include "particles/container/ContainerRef.h"
#include "particles/container/SimpleContainer.h"
#include "utils/Settings.h"

namespace mol_sim {

/**
 * @brief Tests collecting and writing thermostatistical data with the StatsWriter.
 */
class StatsWriterTest : public testing::Test {
   protected:
    StatsWriter writer;
    Domain domain;
    SimpleContainer container;
    std::vector<std::filesystem::path> created_files;
    Particle p1;
    Particle p2;
    Particle p3;
    Particle p4;
    Particle p5;
    size_t n = 0;

    StatsWriterTest() {
        container.reserve(5);
        container.addParticle(R3{1.0, 2.0, 3.0}, R3{0.1, 0.2, 0.3}, 1.5, SettingsParam::EPSILON_DEFAULT,
                              SettingsParam::SIGMA_DEFAULT);
        container.addParticle(R3{-4.0, 5.5, 0.0}, R3{1.0, -0.4, 2.5}, 2.0, SettingsParam::EPSILON_DEFAULT,
                              SettingsParam::SIGMA_DEFAULT);
    }

    void TearDown() override {
        for (const auto& path : created_files) {
            std::error_code ec;
            std::filesystem::remove(path, ec);
        }
        created_files.clear();
    }
};

// ----------------------------------------------------------------------------------------------------

/**
 * @brief Tests collecting and writing RDF data.
 */
class DiffusionWriterTtest : public StatsWriterTest {
   protected:
    void SetUp() override { StatsWriterTest::SetUp(); }
    void TearDown() override { StatsWriterTest::TearDown(); }
    /**
     * @brief Writes the data for diffusion file using the StatsWriter and returns the path to the created file, as well
     * as tracks the created files for cleanup.
     */
    std::filesystem::path writeDiff(SimpleContainer& container, int iteration) {
        ContainerRef particles(container);
        writer.plotDiffusion(particles, iteration);

        std::filesystem::path path = std::filesystem::current_path() / "diff.csv";
        created_files.push_back(path);

        return path;
    }
};

TEST_F(DiffusionWriterTtest, testComputeNormal) {}
TEST_F(DiffusionWriterTtest, testComputePeriodic) {}
TEST_F(DiffusionWriterTtest, testOutput) {}

// ----------------------------------------------------------------------------------------------------

/**
 * @brief Tests collecting and writing RDF data.
 */
class RDFWriterTtest : public StatsWriterTest {
   protected:
    void SetUp() override { StatsWriterTest::SetUp(); }
    void TearDown() override { StatsWriterTest::TearDown(); }
    /**
     * @brief Writes a the data for RDF file using the StatsWriter and returns the path to the created file, as well as
     * tracks the created files for cleanup.
     */
    std::filesystem::path writeRDF(SimpleContainer& container, int iteration) {
        ContainerRef particles(container);
        writer.plotRDF(particles, iteration);

        std::filesystem::path path = std::filesystem::current_path() / "rdf.csv";
        created_files.push_back(path);

        return path;
    }
};

TEST_F(RDFWriterTtest, testComputeNormal) {}
TEST_F(RDFWriterTtest, testComputePeriodic) {}
TEST_F(RDFWriterTtest, testOutput) {}

}  // namespace mol_sim