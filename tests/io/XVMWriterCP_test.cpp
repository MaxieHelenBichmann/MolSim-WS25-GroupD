#include "io/checkpointWriter/XVMWriterCP.h"

#include <gtest/gtest.h>

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

class XVMWriterCPTest : public testing::Test {
   protected:
    XVMWriterCP writer;
    Domain domain;
    SimpleContainer container;
    std::vector<std::filesystem::path> created_files;

    XVMWriterCPTest() : domain(R3{10.0, 10.0, 10.0}) {
        container.reserve(2);
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

    /**
     * @brief Writes a checkpoint file using the XVMWriterCP and returns the path to the created file, as well as tracks
     * the created files for cleanup.
     */
    std::filesystem::path writeCheckpoint(SimpleContainer& container, int iteration, size_t N) {
        ContainerRef particles(container);
        writer.createCheckpoint(domain, particles, iteration, Force::GRAVITATIONAL, 0.01, 0.0, 1.0, 1, 1, "test_xvm",
                                1.0, N);

        int decimal_places = 0;
        while (N >= 10) {
            N /= 10;
            decimal_places++;
        }

        std::ostringstream oss;
        oss << "cp_" << std::setfill('0') << std::setw(decimal_places) << iteration << ".txt";
        std::filesystem::path path = std::filesystem::current_path() / oss.str();
        created_files.push_back(path);

        return path;
    }
};

TEST_F(XVMWriterCPTest, WritesParticleSectionWithCorrectValues) {  // NOLINT
    constexpr int iteration = 7;
    constexpr size_t iteration_cap = 100;
    auto file_path = writeCheckpoint(container, iteration, iteration_cap);

    ASSERT_TRUE(std::filesystem::exists(file_path));
    std::ifstream file(file_path);
    ASSERT_TRUE(file.is_open());

    std::string line;
    ASSERT_TRUE(std::getline(file, line));
    EXPECT_EQ("# Number Particles", line);

    ASSERT_TRUE(std::getline(file, line));
    EXPECT_EQ(std::to_string(container.size()), line);

    ASSERT_TRUE(std::getline(file, line));
    EXPECT_EQ("# Position | Velocity | Mass", line);

    for (const auto& expected : container) {
        ASSERT_TRUE(std::getline(file, line));
        std::istringstream iss(line);
        double x0 = 0.0;
        double x1 = 0.0;
        double x2 = 0.0;
        double v0 = 0.0;
        double v1 = 0.0;
        double v2 = 0.0;
        double mass = 0.0;
        iss >> x0 >> x1 >> x2 >> v0 >> v1 >> v2 >> mass;
        EXPECT_DOUBLE_EQ(expected.getX()[0], x0);
        EXPECT_DOUBLE_EQ(expected.getX()[1], x1);
        EXPECT_DOUBLE_EQ(expected.getX()[2], x2);
        EXPECT_DOUBLE_EQ(expected.getV()[0], v0);
        EXPECT_DOUBLE_EQ(expected.getV()[1], v1);
        EXPECT_DOUBLE_EQ(expected.getV()[2], v2);
        EXPECT_DOUBLE_EQ(expected.getM(), mass);
    }
}

TEST_F(XVMWriterCPTest, ReadBackCheckpointWithXVMReader) {
    constexpr int iteration = 42;
    constexpr size_t iteration_cap = 1000;
    auto file_path = writeCheckpoint(container, iteration, iteration_cap);

    SimpleContainer restored;
    ContainerRef restored_ref(restored);
    XVMReader reader;
    SettingsParam dummy_settings;
    reader.readParticles(restored_ref, dummy_settings, file_path.string());

    const auto contains = [](const Particle& particle, SimpleContainer& cont) -> bool {
        for (auto p : cont) {
            if (p.getX() == particle.getX() && p.getV() == particle.getV() && p.getM() == particle.getM()) {
                return true;
            }
        }
        return false;
    };

    ASSERT_EQ(restored.size(), container.size());
    for (auto& p : container) {
        EXPECT_TRUE(contains(p, restored));
    }
}

}  // namespace mol_sim