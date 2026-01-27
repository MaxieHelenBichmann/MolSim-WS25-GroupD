#include "io/StatsWriter.h"

#include <gtest/gtest.h>

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

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
    SimpleContainer container;
    std::vector<std::filesystem::path> created_files;
    Particle p1;
    Particle p2;
    Particle p3;
    Particle p4;
    Particle p5;

    StatsWriterTest() {
        container.reserve(5);
        p1 = Particle(R3{1.0, 2.0, 3.0}, R3{0.1, 0.2, 0.3}, 1.5, SettingsParam::EPSILON_DEFAULT,
                      SettingsParam::SIGMA_DEFAULT);
        p2 = Particle(R3{-4.0, 5.5, 0.0}, R3{1.0, -0.4, 2.5}, 2.0, SettingsParam::EPSILON_DEFAULT,
                      SettingsParam::SIGMA_DEFAULT);
        p3 = Particle(R3{1.0, 0.0, 0.0}, R3{10.0, 0., 0.}, 1.5, SettingsParam::EPSILON_DEFAULT,
                      SettingsParam::SIGMA_DEFAULT);
        p4 = Particle(R3{0.0, 1.0, 0.0}, R3{0.0, -4., 0.0}, 2.0, SettingsParam::EPSILON_DEFAULT,
                      SettingsParam::SIGMA_DEFAULT);
        p5 = Particle(R3{0.0, 0.0, 1.0}, R3{0.0, 0.0, 6.0}, 1.0, SettingsParam::EPSILON_DEFAULT,
                      SettingsParam::SIGMA_DEFAULT);
        container.push_back(p1);
        container.push_back(p2);
        container.push_back(p3);
        container.push_back(p4);
        container.push_back(p5);
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
class DiffusionWriterTest : public StatsWriterTest {
   protected:
    void SetUp() override { StatsWriterTest::SetUp(); }
    void TearDown() override { StatsWriterTest::TearDown(); }
    /**
     * @brief Writes the data for diffusion file using the StatsWriter and returns the path to the created file, as well
     * as tracks the created files for cleanup.
     */
    std::filesystem::path writeDiff(const StatsWriter& w, SimpleContainer& container, int iteration) {
        std::filesystem::path path = std::filesystem::current_path() / "diffusion.csv";
        {
            std::error_code ec;
            std::filesystem::remove(path, ec);
        }

        ContainerRef particles(container);
        w.plotDiffusion(particles, iteration);
        created_files.push_back(path);

        return path;
    }
};

namespace {

/**
 * @brief Helper function to read all lines from a file.
 */
std::vector<std::string> readAllLines(const std::filesystem::path& path) {
    std::ifstream in(path);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (!line.empty()) {
            lines.push_back(line);
        }
    }
    return lines;
}

/**
 * @brief Helper function to split a CSV line into its components.
 */
std::vector<std::string> splitCsvLine(const std::string& line) {
    std::vector<std::string> parts;
    std::stringstream ss(line);
    std::string item;
    while (std::getline(ss, item, ',')) {
        parts.push_back(item);
    }
    return parts;
}

}  // namespace

/**
 * @brief Tests computing and writing diffusion data.
 */
TEST_F(DiffusionWriterTest, testComputeDisplacement) {
    ContainerRef particles(container);

    container[0].getX() = container[0].getX() + R3{1.0, 0.0, 0.0};  // sq = 1
    container[1].getX() = container[1].getX() + R3{3.0, 0.0, 0.0};  // sq = 9

    const double diffusion = writer.computeDiffusion(particles);
    EXPECT_DOUBLE_EQ(diffusion, 2.0);

    const double diffusion2 = writer.computeDiffusion(particles);
    EXPECT_DOUBLE_EQ(diffusion2, 0.0);
}

/**
 * @brief Tests computing and writing diffusion data with explicit reference positions.
 */
TEST_F(DiffusionWriterTest, testComputeExplicit) {
    ContainerRef particles(container);

    container[0].getRefX() = R3{10.0, 0.0, 0.0};
    container[0].getX() = R3{-10.0, 0.0, 0.0};
    container[1].getRefX() = R3{0.0, 0.0, 0.0};
    container[1].getX() = R3{0.0, 0.0, 0.0};

    const double diffusion = writer.computeDiffusion(particles);
    EXPECT_DOUBLE_EQ(diffusion, 80.0);
}

/**
 * @brief Tests computing and writing diffusion data with a bigger container.
 */
TEST_F(DiffusionWriterTest, testComputeBigger) {
    ContainerRef particles(container);

    for (Particle& p : container) {
        p.getX() += p.getV();
    }

    const double diffusion = writer.computeDiffusion(particles);

    EXPECT_DOUBLE_EQ(diffusion, 31.91);
}

/**
 * @brief Tests writing diffusion data to file.
 */
TEST_F(DiffusionWriterTest, testOutput) {
    writer = StatsWriter(false, true, 1.0, 10.0);

    container[0].getX() = container[0].getX() + R3{1.0, 0.0, 0.0};
    container[1].getX() = container[1].getX() + R3{3.0, 0.0, 0.0};

    const std::filesystem::path path = writeDiff(writer, container, 7);

    const auto lines = readAllLines(path);
    ASSERT_EQ(lines.size(), 1U);
    const auto parts = splitCsvLine(lines[0]);
    ASSERT_EQ(parts.size(), 2U);
    EXPECT_EQ(std::stoi(parts[0]), 7);
    EXPECT_NEAR(std::stod(parts[1]), 2.0, 1e-12);
}

// ----------------------------------------------------------------------------------------------------

/**
 * @brief Tests collecting and writing Temperature data.
 */
class TemperatureWriterTest : public StatsWriterTest {
   protected:
    void SetUp() override { StatsWriterTest::SetUp(); }
    void TearDown() override { StatsWriterTest::TearDown(); }
    /**
     * @brief Writes the data for temperature file using the StatsWriter and returns the path to the created file, as
     * well as tracks the created files for cleanup.
     */
    std::filesystem::path writeTemp(const StatsWriter& w, double total_energy, size_t dimension, size_t N,
                                    int iteration) {
        std::filesystem::path path = std::filesystem::current_path() / "temp.csv";
        {
            std::error_code ec;
            std::filesystem::remove(path, ec);
        }

        ContainerRef particles(container);
        w.plotTemp(total_energy, dimension, N, iteration);
        created_files.push_back(path);

        return path;
    }
};

/**
 * @brief Tests writing temperature data to file.
 */
TEST_F(TemperatureWriterTest, testOutput) {
    writer = StatsWriter(false, true, 1.0, 10.0);

    const std::filesystem::path path = writeTemp(writer, 100.0, 3, 2, 7);

    const auto lines = readAllLines(path);
    ASSERT_EQ(lines.size(), 1U);
    const auto parts = splitCsvLine(lines[0]);
    ASSERT_EQ(parts.size(), 2U);
    EXPECT_EQ(std::stoi(parts[0]), 7);
    EXPECT_NEAR(std::stod(parts[1]), (2 * 100.0) / (3.0 * 2.0), 1e-4);
}

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
    std::filesystem::path writeRDF(const StatsWriter& w, SimpleContainer& container, int iteration) {
        std::filesystem::path path = std::filesystem::current_path() / "rdf.csv";
        {
            std::error_code ec;
            std::filesystem::remove(path, ec);
        }

        ContainerRef particles(container);
        w.plotRDF(particles, iteration);
        created_files.push_back(path);

        return path;
    }
};

/**
 * @brief Tests computing RDF data.
 */
TEST_F(RDFWriterTtest, testComputeNormal) {
    // Enable RDF computation with a known bin width.
    writer = StatsWriter(true, false, 1.0, 10.0);
    SimpleContainer container;
    container.reserve(2);
    container.push_back(Particle(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, SettingsParam::EPSILON_DEFAULT,
                                 SettingsParam::SIGMA_DEFAULT));
    container.push_back(Particle(R3{5.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, SettingsParam::EPSILON_DEFAULT,
                                 SettingsParam::SIGMA_DEFAULT));

    std::vector<double> results(10, 0.0);

    writer.computeRDF(container, results);

    const double expected = 0.75 / (91.0 * std::numbers::pi);
    EXPECT_NEAR(results[5], expected, 1e-12);
}

/**
 * @brief Tests computing and writing RDF data with simulated periodic boundary conditions, meaning that mirrored
 * particles contribute to the RDF bins.
 */
TEST_F(RDFWriterTtest, testComputePeriodic) {
    writer = StatsWriter(true, false, 1.0, 10.0);

    SimpleContainer container;
    container.reserve(2);
    container.push_back(Particle(R3{0.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, SettingsParam::EPSILON_DEFAULT,
                                 SettingsParam::SIGMA_DEFAULT));
    container.push_back(Particle(R3{99.5, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, SettingsParam::EPSILON_DEFAULT,
                                 SettingsParam::SIGMA_DEFAULT));

    ContainerRef particles(container);

    std::vector<double> results(10, 0.0);

    container[0].getMirrorPositions().push_back(R3{100, 0.0, 0.0});   // NOLINT
    container[1].getMirrorPositions().push_back(R3{-0.5, 0.0, 0.0});  // NOLINT

    writer.computeRDF(container, results);

    const double expected = 1.0 / ((4.0 * std::numbers::pi) / 3.0);
    EXPECT_NEAR(results[0], expected, 1e-12);
}

/**
 * @brief Tests writing RDF data to file.
 */
TEST_F(RDFWriterTtest, testOutput) {
    writer = StatsWriter(true, false, 1.0, 3.0);

    const std::vector<Particle> mirrored;
    const std::filesystem::path path = writeRDF(writer, container, 42);

    const auto lines = readAllLines(path);

    ASSERT_EQ(lines.size(), 3U);

    for (size_t i = 0; i < lines.size(); ++i) {
        const auto parts = splitCsvLine(lines[i]);
        ASSERT_EQ(parts.size(), 3U);
        EXPECT_EQ(std::stoi(parts[0]), 42);
        EXPECT_NEAR(std::stod(parts[1]), static_cast<double>(i) * 1.0, 1e-12);
    }
}

}  // namespace mol_sim
