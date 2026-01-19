#include <gtest/gtest.h>
#include <particles/Particle.h>

#include <memory>
#include <string>
#include <vector>

#include "../testingUtils.h"
#include "io/fileReader/YAMLReader.h"
#include "particles/container/SimpleContainer.h"
#include "testingUtils.h"
#include "utils/Settings.h"

namespace mol_sim {

const std::string test_data_dir = TEST_DATA_DIR;  // NOLINT

/**
 * @brief Test fixture for multi-file reading scenarios
 */
class MultiFileReaderTest : public testing::Test {
   protected:
    SimpleContainer part_container;
    SettingsParam settings;
    std::unique_ptr<YAMLReader> reader;

    void SetUp() override {
        reader = std::make_unique<YAMLReader>();
        part_container.clear();
        settings = SettingsParam();
    }
};

/**
 * @brief Test that settings from the last file override earlier files
 */
TEST_F(MultiFileReaderTest, SettingsOverrideFromMultipleFiles) {
    std::vector<std::string> files = {test_data_dir + "/checkpoint_input.yaml", test_data_dir + "/override_input.yaml"};

    for (const auto& file : files) {
        reader->readSettings(settings, file);
    }

    EXPECT_DOUBLE_EQ(settings.delta_t, 0.005);
    EXPECT_DOUBLE_EQ(settings.end_time, 200.0);
    EXPECT_EQ(settings.base_name, "OVERRIDE");
    EXPECT_EQ(settings.pairwise_forces.size(), 1);
    EXPECT_EQ(settings.pairwise_forces[0], LENNARDJONES);

    EXPECT_DOUBLE_EQ(settings.start_time, 0.0);
    EXPECT_EQ(settings.frequency_output, 5);
    EXPECT_EQ(settings.frequency_checkpoint, 50);
}

/**
 * @brief Test that particles accumulate from multiple files
 */
TEST_F(MultiFileReaderTest, ParticlesAccumulateFromMultipleFiles) {
    std::vector<std::string> files = {test_data_dir + "/checkpoint_input.yaml", test_data_dir + "/override_input.yaml"};

    for (const auto& file : files) {
        reader->readParticles(part_container, settings, file);
    }

    ASSERT_EQ(part_container.size(), 3);
    Particle expected_1 = {
        {1.0, 2.0, 3.0}, {0.1, 0.2, 0.3}, 1, SettingsParam::EPSILON_DEFAULT, SettingsParam::SIGMA_DEFAULT};

    EXPECT_PARTICLE_EQ(part_container[0], expected_1);

    Particle expected_2 = {
        {4.0, 5.0, 6.0}, {0.4, 0.5, 0.6}, 2., SettingsParam::EPSILON_DEFAULT, SettingsParam::SIGMA_DEFAULT};

    EXPECT_PARTICLE_EQ(part_container[1], expected_2);

    Particle expected_3 = {
        {10.0, 20.0, 30.0}, {1.0, 2.0, 3.0}, 5., SettingsParam::EPSILON_DEFAULT, SettingsParam::SIGMA_DEFAULT};

    EXPECT_PARTICLE_EQ(part_container[2], expected_3);
}  // namespace mol_sim

/**
 * @brief Test complete checkpoint + input workflow
 */
TEST_F(MultiFileReaderTest, CheckpointAndInputFileWorkflow) {
    std::vector<std::string> files = {test_data_dir + "/checkpoint_input.yaml", test_data_dir + "/override_input.yaml"};

    for (const auto& file : files) {
        reader->readSettings(settings, file);
    }

    for (const auto& file : files) {
        reader->readParticles(part_container, settings, file);
    }

    EXPECT_DOUBLE_EQ(settings.delta_t, 0.005);
    EXPECT_DOUBLE_EQ(settings.end_time, 200.0);
    EXPECT_EQ(settings.base_name, "OVERRIDE");
    EXPECT_EQ(settings.pairwise_forces.size(), 1);
    EXPECT_EQ(settings.pairwise_forces[0], LENNARDJONES);
    EXPECT_EQ(part_container.size(), 3);
}

/**
 * @brief Test that only non-conflicting settings from checkpoint persist
 */
TEST_F(MultiFileReaderTest, PartialOverrideKeepsNonConflictingSettings) {
    std::vector<std::string> files = {test_data_dir + "/checkpoint_input.yaml", test_data_dir + "/override_input.yaml"};

    for (const auto& file : files) {
        reader->readSettings(settings, file);
    }

    EXPECT_DOUBLE_EQ(settings.delta_t, 0.005);
    EXPECT_DOUBLE_EQ(settings.end_time, 200.0);
    EXPECT_EQ(settings.base_name, "OVERRIDE");

    EXPECT_DOUBLE_EQ(settings.start_time, 0.0);
    EXPECT_EQ(settings.frequency_output, 5);
    EXPECT_EQ(settings.frequency_checkpoint, 50);
}

/**
 * @brief Test reading from single file still works
 */
TEST_F(MultiFileReaderTest, SingleFileReading) {
    std::vector<std::string> files = {test_data_dir + "/full_config.yaml"};

    for (const auto& file : files) {
        reader->readSettings(settings, file);
        reader->readParticles(part_container, settings, file);
    }

    EXPECT_DOUBLE_EQ(settings.delta_t, 0.005);
    EXPECT_DOUBLE_EQ(settings.end_time, 500.0);
    EXPECT_EQ(part_container.size(), 1);
}

}  // namespace mol_sim
