#include "io/fileReader/YAMLReader.h"

#include <gtest/gtest.h>
#include <spdlog/spdlog.h>

#include <mutex>
#include <string>

#include "io/fileReader/YAMLReaderException.h"
#include "particles/Particle.h"
#include "particles/container/ContainerRef.h"
#include "particles/container/SimpleContainer.h"
#include "spdlog/sinks/ostream_sink.h"
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
        std::lock_guard<std::mutex> lock(mtx);
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
        std::lock_guard<std::mutex> lock(mtx);
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
    reader.readFile(particles, settings, test_data_dir + "/simple_XVM.yaml");
    std::string output = log_stream->str();
    EXPECT_EQ(particles.size(), 1);
    EXPECT_EQ(particles[0].getM(), 1.0);
    EXPECT_EQ(output.find("Error"), std::string::npos);
}
/**
 * @brief Tests the YAML Readers Cuboid Format support.
 * Tests this by reading in simple_cuboid.yaml
 */
TEST_F(YAMLReaderTest, ReadSimpleCuboid) {
    YAMLReader reader;
    reader.readFile(particles, settings, test_data_dir + "/simple_cuboid.yaml");
    std::string output = log_stream->str();
    EXPECT_EQ(particles.size(), 8);
    for (auto& p : particles) {
        EXPECT_EQ(p.getM(), 1.0);
    }
    EXPECT_EQ(output.find("Error"), std::string::npos);
}
/**
 * @brief Tests the behaviour of YAMLReader when an invalid filepath is given
 *
 */
TEST_F(YAMLReaderTest, ReadNonExistentFile) {
    YAMLReader reader;
    EXPECT_THROW(reader.readFile(particles, settings, test_data_dir + "/bogus_file.yaml"), YAMLReaderException);
}
/**
 * @brief Tests the behaviour of YAMLReader when an unknown format option is parsed
 *
 */
TEST_F(YAMLReaderTest, ReadWrongFileFormat) {
    YAMLReader reader;
    EXPECT_THROW(reader.readFile(particles, settings, test_data_dir + "/unknown_format.yaml"), YAMLReaderException);
    std::string output = log_stream->str();
    EXPECT_NE(output.find("Unknown Format"), std::string::npos);
}

/**
 * @brief Tests the behaviour of YAMLReader when a malformed file is given
 *
 */
TEST_F(YAMLReaderTest, ReadMalformedFile) {
    YAMLReader reader;
    EXPECT_THROW(reader.readFile(particles, settings, test_data_dir + "/malformed.yaml"), YAMLReaderException);
}

/**
 * @brief Tests the behaviour of YAMLReader when a file with missing fields is given
 *
 */
TEST_F(YAMLReaderTest, ReadMissingFields) {
    YAMLReader reader;
    EXPECT_THROW(reader.readFile(particles, settings, test_data_dir + "/missing_fields.yaml"), YAMLReaderException);
}

/**
 * @brief Tests the behaviour of YAMLReader when an empty file is given
 *
 */
TEST_F(YAMLReaderTest, ReadEmptyFile) {
    YAMLReader reader;
    EXPECT_THROW(reader.readFile(particles, settings, test_data_dir + "/empty.yaml"), YAMLReaderException);
}

/**
 * @brief Tests the YAML Readers ability to read multiple objects from a single file
 *
 */
TEST_F(YAMLReaderTest, ReadMultipleObjects) {
    YAMLReader reader;
    reader.readFile(particles, settings, test_data_dir + "/multiple_objects.yaml");
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
    EXPECT_THROW(reader.readFile(particles, settings, test_data_dir + "/only_settings.yaml"), YAMLReaderException);
}

/**
 * @brief Tests the YAML Readers ability to read a full config file with settings and multiple objects
 *
 */
TEST_F(YAMLReaderTest, ReadFullConfigFile) {
    YAMLReader reader;
    reader.readFile(particles, settings, test_data_dir + "/full_config.yaml");
    std::string output = log_stream->str();
    EXPECT_EQ(particles.size(), 5);

    ASSERT_TRUE(settings.delta_t.has_value());
    EXPECT_EQ(settings.delta_t.value(), 0.005);

    ASSERT_TRUE(settings.end_time.has_value());
    EXPECT_EQ(settings.end_time.value(), 500.0);
    EXPECT_EQ(output.find("Error"), std::string::npos);
}

}  // namespace mol_sim
