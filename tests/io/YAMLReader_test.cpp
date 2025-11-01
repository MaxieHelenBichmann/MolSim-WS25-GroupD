#include "io/fileReader/YAMLReader.h"

#include <gtest/gtest.h>

#include <string>

#include "io/fileReader/YAMLReaderException.h"
#include "particles/Particle.h"
#include "particles/container/ContainerRef.h"
#include "particles/container/SimpleContainer.h"
#include "spdlog/sinks/ostream_sink.h"
#include "utils/Logging.h"

namespace mol_sim {
// The TEST_DATA_DIR macro is defined in CMakeLists.txt
const std::string test_data_dir = TEST_DATA_DIR;

class YAMLReaderTest : public testing::Test {
   protected:
    SimpleContainer part_container;
    ContainerRef particles;
    YAMLReaderTest() : part_container(), particles(part_container) {}
    // This string stream will capture all log output during a test.
    std::shared_ptr<std::ostringstream> log_stream;

    // Store the original logger to restore it after the test is done.
    std::shared_ptr<spdlog::logger> original_logger;

    void SetUp() override {
        // 1. Save the existing default logger so we can restore it later.
        original_logger = spdlog::default_logger();

        // 2. Create an ostringstream to capture the logs.
        log_stream = std::make_shared<std::ostringstream>();

        // 3. Create a sink that writes to our stringstream.
        auto ostream_sink = std::make_shared<spdlog::sinks::ostream_sink_st>(*log_stream);

        // 4. Create a new logger with our sink.
        auto test_logger = std::make_shared<spdlog::logger>("test_logger", ostream_sink);

        // 5. Set the new logger as the default for the duration of this test.
        spdlog::set_default_logger(test_logger);

        // 6. Set the log level to the most permissive level to ensure we capture everything.
        spdlog::set_level(spdlog::level::trace);
    }

    void TearDown() override {
        // Restore the original logger to avoid side-effects between tests.
        spdlog::set_default_logger(original_logger);
    }
};
TEST_F(YAMLReaderTest, ReadSimpleXVM) {
    YAMLReader reader;
    reader.readFile(particles, test_data_dir + "/simple_XVM.yaml");
    std::string output = log_stream->str();
    EXPECT_EQ(particles.size(), 1);
    EXPECT_EQ(particles[0].getM(), 1.0);
    EXPECT_EQ(output.find("Error"), std::string::npos);
}
TEST_F(YAMLReaderTest, ReadSimpleCuboid) {
    YAMLReader reader;
    reader.readFile(particles, test_data_dir + "/simple_cuboid.yaml");
    std::string output = log_stream->str();
    EXPECT_EQ(particles.size(), 8);
    for (auto& p : particles) {
        EXPECT_EQ(p.getM(), 1.0);
    }
    EXPECT_EQ(output.find("Error"), std::string::npos);
}

TEST_F(YAMLReaderTest, ReadNonExistentFile) {
    YAMLReader reader;
    EXPECT_THROW(reader.readFile(particles, test_data_dir + "/bogus_file.yaml"), YAMLReaderException);
}

TEST_F(YAMLReaderTest, ReadWrongFileFormat) {
    YAMLReader reader;
    EXPECT_THROW(reader.readFile(particles, test_data_dir + "/unknown_format.yaml"), YAMLReaderException);
    std::string output = log_stream->str();
    EXPECT_NE(output.find("Unknown YAML Format"), std::string::npos);
}

}  // namespace mol_sim
