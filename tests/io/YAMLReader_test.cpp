#include "io/fileReader/YAMLReader.h"

#include <gtest/gtest.h>
#include <physics/ForceSource.h>
#include <spdlog/spdlog.h>

#include <mutex>
#include <string>
#include <variant>

#include "exceptions/YAMLReaderException.h"
#include "particles/Particle.h"
#include "particles/boundaries/Boundary.h"
#include "particles/boundaries/Reflecting.h"
#include "particles/container/ContainerRef.h"
#include "particles/container/SimpleContainer.h"
#include "spdlog/sinks/ostream_sink.h"
#include "testingUtils.h"
#include "utils/Settings.h"

namespace mol_sim {
// The TEST_DATA_DIR macro is defined in CMakeLists.txt
const std::string test_data_dir = TEST_DATA_DIR;  // NOLINT
/**
 * @brief Test Fixture for YAMLReader tests.
 * This is used to be able to easily reuse the data throughout the tests.
 */
class YAMLReaderTest : public testing::Test {
   protected:
    SimpleContainer part_container;
    ContainerRef particles;
    SettingsParam settings;
    YAMLReaderTest() : part_container(), particles(part_container) {}
    // This string stream will capture all log output during a test.
    std::shared_ptr<std::ostringstream> log_stream;

    // Store the original logger to restore it after the test is done.
    std::shared_ptr<spdlog::logger> original_logger;

    static std::mutex mtx;

    void SetUp() override {
        std::lock_guard<std::mutex> lock(mtx);  // NOLINT
        // 1. Save the existing default logger so we can restore it later.
        original_logger = spdlog::default_logger();

        // 2. Create an ostringstream to capture the logs.
        log_stream = std::make_shared<std::ostringstream>();

        // 3. Create a sink that writes to our stringstream.
        auto ostream_sink = std::make_shared<spdlog::sinks::ostream_sink_st>(*log_stream, true);

        // 4. Create a new logger with our sink.
        auto test_logger = std::make_shared<spdlog::logger>("test_logger", ostream_sink);

        // 5. Set the new logger as the default for the duration of this test.
        spdlog::set_default_logger(test_logger);

        // 6. Set the log level to the most permissive level to ensure we capture everything.
        spdlog::set_level(spdlog::level::trace);
    }

    void TearDown() override {
        std::lock_guard<std::mutex> lock(mtx);  // NOLINT
        // Restore the original logger to avoid side-effects between tests.
        spdlog::set_default_logger(original_logger);
    }
};

std::mutex YAMLReaderTest::mtx;

/**
 * @brief Tests the YAML Readers XVM Format support
 * Tests this by reading in simple_XVM.yaml
 */
TEST_F(YAMLReaderTest, ReadSimpleXVM) {
    YAMLReader reader;
    reader.readParticles(particles, test_data_dir + "/simple_XVM.yaml");
    std::string output = log_stream->str();
    R3 expected_pos = {0., 0., 0.};
    R3 expected_velo = {0., 0., 0.};
    ASSERT_EQ(particles.size(), 1);
    EXPECT_R3_EQ(particles[0].getX(), expected_pos);
    EXPECT_R3_EQ(particles[0].getV(), expected_velo);
    EXPECT_EQ(particles[0].getM(), 1.0);
    EXPECT_EQ(output.find("Error"), std::string::npos);
}
/**
 * @brief Tests the YAML Readers Cuboid Format support.
 * Tests this by reading in simple_cuboid.yaml
 */
TEST_F(YAMLReaderTest, ReadSimpleCuboid) {
    YAMLReader reader;
    YAML::Node root = YAML::LoadFile(test_data_dir + "/simple_cuboid.yaml");
    auto cuboids = reader.parseCuboids(root["simple_cuboid"]);
    ASSERT_EQ(cuboids.size(), 1);
    R3 expected_pos = {0., 0., 0.};
    R3 expected_velo = {0., 0., 0.};
    auto cube = cuboids[0];
    EXPECT_EQ(cube.position, expected_pos);
    EXPECT_EQ(cube.velocity, expected_pos);
    EXPECT_EQ(cube.num_particles[0], 2);
    EXPECT_EQ(cube.num_particles[1], 2);
    EXPECT_EQ(cube.num_particles[2], 2);
    EXPECT_DOUBLE_EQ(cube.avg_velo, 0.);
    EXPECT_DOUBLE_EQ(cube.mass, 1.);
    EXPECT_DOUBLE_EQ(cube.distance, 1.);
    EXPECT_DOUBLE_EQ(cube.epsilon, 5.);
    EXPECT_DOUBLE_EQ(cube.sigma, 1.);
}
/**
 * @brief Tests the YAML Readers Disc Format support.
 * Tests this by reading in simple_disc.yaml
 */
TEST_F(YAMLReaderTest, ReadSimpleDisc) {
    YAMLReader reader;
    YAML::Node root = YAML::LoadFile(test_data_dir + "/simple_disc.yaml");
    auto discs = reader.parseDiscs(root["simple_disc"]);
    ASSERT_EQ(discs.size(), 1);
    R3 expected_pos = {0., 0., 0.};
    R3 expected_velo = {0., 0., 0.};
    auto disc = discs[0];
    EXPECT_EQ(disc.position, expected_pos);
    EXPECT_EQ(disc.velocity, expected_pos);
    EXPECT_EQ(disc.radius, 2);

    EXPECT_DOUBLE_EQ(disc.avg_velo, 0.);
    EXPECT_DOUBLE_EQ(disc.mass, 1.);
    EXPECT_DOUBLE_EQ(disc.distance, 1.);
    EXPECT_DOUBLE_EQ(disc.epsilon, 5.);
    EXPECT_DOUBLE_EQ(disc.sigma, 1.);
}
/**
 * @brief Tests the behaviour of YAMLReader when an invalid filepath is given
 *
 */
TEST_F(YAMLReaderTest, ReadSettingsNonExistentFile) {
    YAMLReader reader;
    EXPECT_THROW(reader.readSettings(settings, test_data_dir + "/bogus_file.yaml"), YAMLReaderException);
}
/**
 * @brief Tests the behaviour of YAMLReader when an unknown format option is parsed
 *
 */
TEST_F(YAMLReaderTest, ReadWrongFileFormat) {
    YAMLReader reader;
    // File with unknown format has no valid particle definitions
    EXPECT_THROW(reader.readParticles(particles, test_data_dir + "/unknown_format.yaml"), YAMLReaderException);
}

/**
 * @brief Tests the behaviour of YAMLReader when a malformed file is given
 *
 */
TEST_F(YAMLReaderTest, ReadSettingsMalformedFile) {
    YAMLReader reader;
    EXPECT_THROW(reader.readSettings(settings, test_data_dir + "/malformed.yaml"), YAMLReaderException);
}

/**
 * @brief Tests the behaviour of YAMLReader when a file with missing fields is given
 *
 */
TEST_F(YAMLReaderTest, ReadParticleMissingFields) {
    YAMLReader reader;
    EXPECT_THROW(reader.readParticles(particles, test_data_dir + "/missing_fields.yaml"), YAMLReaderException);
}

/**
 * @brief Tests the behaviour of YAMLReader when an empty file is given
 *
 */
TEST_F(YAMLReaderTest, ReadSettingsEmptyFile) {
    YAMLReader reader;
    EXPECT_THROW(reader.readSettings(settings, test_data_dir + "/empty.yaml"), YAMLReaderException);
}

/**
 * @brief Tests the YAML Readers ability to read multiple objects from a single file
 *
 */
TEST_F(YAMLReaderTest, ReadMultipleObjects) {
    YAMLReader reader;
    reader.readSettings(settings, test_data_dir + "/multiple_objects.yaml");
    reader.readParticles(particles, test_data_dir + "/multiple_objects.yaml");
    std::string output = log_stream->str();
    EXPECT_EQ(particles.size(), 5);
    EXPECT_EQ(output.find("Error"), std::string::npos);
}

/**
 * @brief Tests the behaviour of YAMLReader when a file with only settings is given
 *
 */
TEST_F(YAMLReaderTest, ReadFileWithOnlySettings) {
    YAMLReader reader;
    EXPECT_THROW(reader.readParticles(particles, test_data_dir + "/only_settings.yaml"), YAMLReaderException);
}

/**
 * @brief  Tests the behaviour of YAMLReader when a file with no settings is given
 */
TEST_F(YAMLReaderTest, ReadFileWithNoSettings) {
    YAMLReader reader;
    EXPECT_THROW(reader.readSettings(settings, test_data_dir + "/no_settings.yaml"), YAMLReaderException);
}

/**
 * @brief Tests the YAML Readers ability to read a full config file with settings (except domain) and multiple objects
 *
 */
TEST_F(YAMLReaderTest, ReadFullConfigFile) {
    YAMLReader reader;
    reader.readSettings(settings, test_data_dir + "/full_config.yaml");
    reader.readParticles(particles, test_data_dir + "/full_config.yaml");
    std::string output = log_stream->str();
    EXPECT_EQ(particles.size(), 1);

    ASSERT_TRUE(settings.delta_t.has_value());
    EXPECT_EQ(settings.delta_t.value(), 0.005);

    ASSERT_TRUE(settings.end_time.has_value());
    EXPECT_DOUBLE_EQ(settings.end_time.value(), 500.0);

    ASSERT_TRUE(settings.start_time.has_value());
    EXPECT_DOUBLE_EQ(settings.start_time.value(), 0.0);

    ASSERT_TRUE(settings.base_name.has_value());
    EXPECT_EQ(settings.base_name.value(), "MD");

    ASSERT_TRUE(settings.force.has_value());
    EXPECT_EQ(settings.force.value(), LENNARDJONES);

    ASSERT_TRUE(settings.frequency.has_value());
    EXPECT_EQ(settings.frequency.value(), 10);

    ASSERT_TRUE(settings.cutoff.has_value());
    EXPECT_DOUBLE_EQ(settings.cutoff.value(), 1.);

    EXPECT_EQ(output.find("Error"), std::string::npos);
}

/**
 * @brief Tests that readSettings only reads settings and doesn't affect particle container
 */
TEST_F(YAMLReaderTest, ReadSettingsOnly) {
    YAMLReader reader;
    reader.readSettings(settings, test_data_dir + "/full_config.yaml");

    // Settings should be populated
    ASSERT_TRUE(settings.delta_t.has_value());
    EXPECT_DOUBLE_EQ(settings.delta_t.value(), 0.005);
    ASSERT_TRUE(settings.end_time.has_value());
    EXPECT_DOUBLE_EQ(settings.end_time.value(), 500.0);

    // Particles should remain empty (we only called readSettings)
    EXPECT_EQ(particles.size(), 0);
}

/**
 * @brief Tests that readParticles only reads particles and doesn't affect settings
 */
TEST_F(YAMLReaderTest, ReadParticlesOnly) {
    YAMLReader reader;
    reader.readParticles(particles, test_data_dir + "/full_config.yaml");

    // Particles should be populated
    EXPECT_EQ(particles.size(), 1);

    // Settings should remain unset (we only called readParticles)
    EXPECT_FALSE(settings.delta_t.has_value());
    EXPECT_FALSE(settings.end_time.has_value());
}

/**
 * @brief Tests the two-step workflow: settings first, then particles
 */
TEST_F(YAMLReaderTest, TwoStepReading) {
    YAMLReader reader;

    // Step 1: Read settings
    reader.readSettings(settings, test_data_dir + "/full_config.yaml");
    ASSERT_TRUE(settings.delta_t.has_value());
    EXPECT_EQ(particles.size(), 0);

    // Step 2: Read particles
    reader.readParticles(particles, test_data_dir + "/full_config.yaml");
    EXPECT_EQ(particles.size(), 1);

    // Both should now be populated
    EXPECT_DOUBLE_EQ(settings.delta_t.value(), 0.005);
    EXPECT_EQ(particles.size(), 1);
}

/**
 * @brief Tests reading settings from only_settings file succeeds
 */
TEST_F(YAMLReaderTest, ReadSettingsFromOnlySettingsFile) {
    YAMLReader reader;
    reader.readSettings(settings, test_data_dir + "/only_settings.yaml");

    ASSERT_TRUE(settings.delta_t.has_value());
    EXPECT_DOUBLE_EQ(settings.delta_t.value(), 0.005);
    ASSERT_TRUE(settings.end_time.has_value());
    EXPECT_DOUBLE_EQ(settings.end_time.value(), 500.0);
}

/**
 * @brief Tests reading particles from no_settings file succeeds
 */
TEST_F(YAMLReaderTest, ReadParticlesFromNoSettingsFile) {
    YAMLReader reader;
    reader.readParticles(particles, test_data_dir + "/no_settings.yaml");

    ASSERT_EQ(particles.size(), 1);
    R3 expected_pos = {0., 0., 0.};
    R3 expected_velo = {0., 0., 0.};
    EXPECT_R3_EQ(particles[0].getX(), expected_pos);
    EXPECT_R3_EQ(particles[0].getV(), expected_velo);
    EXPECT_EQ(particles[0].getM(), 1.0);
}

/**
 * @brief Tests that readParticles on empty file throws
 */
TEST_F(YAMLReaderTest, ReadParticlesEmptyFile) {
    YAMLReader reader;
    EXPECT_THROW(reader.readParticles(particles, test_data_dir + "/empty.yaml"), YAMLReaderException);
}

/**
 * @brief Tests that readParticles on non-existent file throws
 */
TEST_F(YAMLReaderTest, ReadParticlesNonExistentFile) {
    YAMLReader reader;
    EXPECT_THROW(reader.readParticles(particles, test_data_dir + "/bogus_file.yaml"), YAMLReaderException);
}

/**
 * @brief Tests the YAML Readers ability to read domain and boundary configurations
 */
TEST_F(YAMLReaderTest, ReadDomainAndBoundaries) {
    YAMLReader reader;
    reader.readSettings(settings, test_data_dir + "/domain_boundaries.yaml");
    std::string output = log_stream->str();

    // Verify domain was read
    ASSERT_TRUE(settings.domain.has_value());
    R3 expected_domain = {180.0, 90.0, 50.0};
    EXPECT_R3_EQ(settings.domain->getDimension(), expected_domain);

    // Verify LEFT boundary (REFLECTING with sigma=1.2, epsilon=5.0)
    const Boundary* left = settings.domain->getBoundary(BoundaryLocation::LEFT);
    ASSERT_NE(left, nullptr);
    EXPECT_EQ(left->getType(), BoundaryType::REFLECTING);
    const auto* left_reflecting = dynamic_cast<const Reflecting*>(left);
    ASSERT_NE(left_reflecting, nullptr);
    ASSERT_TRUE(left_reflecting->getBoundarySigma().has_value());
    EXPECT_DOUBLE_EQ(left_reflecting->getBoundarySigma().value(), 1.2);
    ASSERT_TRUE(left_reflecting->getBoundaryEpsilon().has_value());
    EXPECT_DOUBLE_EQ(left_reflecting->getBoundaryEpsilon().value(), 5.0);

    // Verify RIGHT boundary (REFLECTING without custom sigma/epsilon)
    const Boundary* right = settings.domain->getBoundary(BoundaryLocation::RIGHT);
    ASSERT_NE(right, nullptr);
    EXPECT_EQ(right->getType(), BoundaryType::REFLECTING);

    // Verify FRONT boundary (OUTFLOW)
    const Boundary* front = settings.domain->getBoundary(BoundaryLocation::FRONT);
    ASSERT_NE(front, nullptr);
    EXPECT_EQ(front->getType(), BoundaryType::OUTFLOW);

    // Verify BACK boundary (OUTFLOW)
    const Boundary* back = settings.domain->getBoundary(BoundaryLocation::BACK);
    ASSERT_NE(back, nullptr);
    EXPECT_EQ(back->getType(), BoundaryType::OUTFLOW);

    // Verify UPPER boundary (REFLECTING with sigma=2.0, epsilon=10.0)
    const Boundary* upper = settings.domain->getBoundary(BoundaryLocation::UPPER);
    ASSERT_NE(upper, nullptr);
    EXPECT_EQ(upper->getType(), BoundaryType::REFLECTING);
    const auto* upper_reflecting = dynamic_cast<const Reflecting*>(upper);
    ASSERT_NE(upper_reflecting, nullptr);
    ASSERT_TRUE(upper_reflecting->getBoundarySigma().has_value());
    EXPECT_DOUBLE_EQ(upper_reflecting->getBoundarySigma().value(), 2.0);
    ASSERT_TRUE(upper_reflecting->getBoundaryEpsilon().has_value());
    EXPECT_DOUBLE_EQ(upper_reflecting->getBoundaryEpsilon().value(), 10.0);

    // Verify LOWER boundary (OUTFLOW)
    const Boundary* lower = settings.domain->getBoundary(BoundaryLocation::LOWER);
    ASSERT_NE(lower, nullptr);
    EXPECT_EQ(lower->getType(), BoundaryType::OUTFLOW);

    EXPECT_EQ(output.find("Error"), std::string::npos);
}

}  // namespace mol_sim
