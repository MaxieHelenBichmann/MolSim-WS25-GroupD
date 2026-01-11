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

namespace {

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

TEST_F(DiffusionWriterTtest, testComputeNormal) {
    ContainerRef particles(container);

    // Move particles by known displacements.
    container[0].getX() = container[0].getX() + R3{1.0, 0.0, 0.0};  // sq = 1
    container[1].getX() = container[1].getX() + R3{3.0, 0.0, 0.0};  // sq = 9

    const double diffusion = writer.computeDiffusion(particles);
    EXPECT_DOUBLE_EQ(diffusion, 5.0);

    // Calling again without movement should yield 0 (reference positions updated).
    const double diffusion2 = writer.computeDiffusion(particles);
    EXPECT_DOUBLE_EQ(diffusion2, 0.0);
}

TEST_F(DiffusionWriterTtest, testComputePeriodic) {
    // This mainly checks that reference positions are honored even if coordinates are negative
    // (the PBC handling is encapsulated in how ref positions are managed elsewhere).
    ContainerRef particles(container);

    container[0].getRefX() = R3{10.0, 0.0, 0.0};
    container[0].getX() = R3{-10.0, 0.0, 0.0};
    container[1].getRefX() = R3{0.0, 0.0, 0.0};
    container[1].getX() = R3{0.0, 0.0, 0.0};

    const double diffusion = writer.computeDiffusion(particles);
    // ((-20)^2 + 0) / 2 = 200
    EXPECT_DOUBLE_EQ(diffusion, 200.0);
}

TEST_F(DiffusionWriterTtest, testOutput) {
    // Enable diffusion output.
    writer = StatsWriter(false, true, 1.0, 10.0);
    ContainerRef particles(container);

    const std::filesystem::path path = std::filesystem::current_path() / "diffusion.csv";
    {
        std::error_code ec;
        std::filesystem::remove(path, ec);
    }
    created_files.push_back(path);

    container[0].getX() = container[0].getX() + R3{1.0, 0.0, 0.0};
    container[1].getX() = container[1].getX() + R3{3.0, 0.0, 0.0};

    writer.plotDiffusion(particles, 7);

    const auto lines = readAllLines(path);
    ASSERT_EQ(lines.size(), 1U);
    const auto parts = splitCsvLine(lines[0]);
    ASSERT_EQ(parts.size(), 2U);
    EXPECT_EQ(std::stoi(parts[0]), 7);
    EXPECT_NEAR(std::stod(parts[1]), 5.0, 1e-12);
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
    std::filesystem::path writeRDF(SimpleContainer& container, int iteration, const std::vector<Particle>& mirrored) {
        ContainerRef particles(container);
        writer.plotRDF(particles, iteration, mirrored);

        std::filesystem::path path = std::filesystem::current_path() / "rdf.csv";
        created_files.push_back(path);

        return path;
    }
};

TEST_F(RDFWriterTtest, testComputeNormal) {
    // Enable RDF computation with a known bin width.
    writer = StatsWriter(true, false, 1.0, 10.0);
    ContainerRef particles(container);

    container[0].getX() = R3{0.0, 0.0, 0.0};
    container[1].getX() = R3{5.0, 0.0, 0.0};

    std::vector<double> results(10, 0.0);
    const std::vector<Particle> mirrored;

    writer.computeRDF(particles, results, mirrored);

    // Two particles -> pair counted twice (p0->p1, p1->p0) into bin index 5.
    // density = 0.375 * 2 / ( (6^3-5^3) * pi )
    const double expected = 0.75 / (91.0 * std::numbers::pi);
    EXPECT_NEAR(results[5], expected, 1e-12);
}

TEST_F(RDFWriterTtest, testComputePeriodic) {
    // This checks that mirrored particles contribute to the RDF bins.
    writer = StatsWriter(true, false, 1.0, 10.0);
    ContainerRef particles(container);

    container[0].getX() = R3{0.0, 0.0, 0.0};
    container[1].getX() = R3{100.0, 0.0, 0.0};  // out of RDF window for bin_count=10

    std::vector<double> results(10, 0.0);
    std::vector<Particle> mirrored;
    mirrored.emplace_back(R3{1.0, 0.0, 0.0}, R3{0.0, 0.0, 0.0}, 1.0, SettingsParam::EPSILON_DEFAULT,
                          SettingsParam::SIGMA_DEFAULT);

    writer.computeRDF(particles, results, mirrored);

    // Only p0 sees the mirrored particle at distance 1 -> count 1 into bin index 1.
    const double expected = 0.375 / (7.0 * std::numbers::pi);  // vol = 2^3 - 1^3 = 7
    EXPECT_NEAR(results[1], expected, 1e-12);
}

TEST_F(RDFWriterTtest, testOutput) {
    writer = StatsWriter(true, false, 1.0, 3.0);
    ContainerRef particles(container);

    const std::filesystem::path path = std::filesystem::current_path() / "rdf.csv";
    {
        std::error_code ec;
        std::filesystem::remove(path, ec);
    }
    created_files.push_back(path);

    const std::vector<Particle> mirrored;
    writer.plotRDF(particles, 42, mirrored);

    const auto lines = readAllLines(path);
    // bin_count = ceil(3/1) = 3
    ASSERT_EQ(lines.size(), 3U);

    for (size_t i = 0; i < lines.size(); ++i) {
        const auto parts = splitCsvLine(lines[i]);
        ASSERT_EQ(parts.size(), 3U);
        EXPECT_EQ(std::stoi(parts[0]), 42);
        EXPECT_NEAR(std::stod(parts[1]), static_cast<double>(i) * 1.0, 1e-12);
    }
}

}  // namespace mol_sim