#include "io/checkpointWriter/YAMLWriterCP.h"

#include <gtest/gtest.h>
#include <spdlog/sinks/ostream_sink.h>
#include <yaml-cpp/yaml.h>

#include <array>
#include <filesystem>
#include <iomanip>
#include <memory>
#include <sstream>
#include <vector>

#include "io/fileReader/YAMLReader.h"
#include "particles/boundaries/Outflow.h"
#include "particles/boundaries/Reflecting.h"
#include "particles/container/ContainerRef.h"
#include "particles/container/SimpleContainer.h"
#include "utils/Settings.h"

namespace mol_sim {
class YAMLWriterCPTest : public testing::Test {
   protected:
    YAMLWriterCP writer;
    SimpleContainer container;
    std::vector<std::filesystem::path> created_files;

    YAMLWriterCPTest() {
        container.reserve(2);
        container.addParticle(R3{1.0, 2.0, 3.0}, R3{0.1, 0.2, 0.3}, 1.5, 2.5, 1.7);
        container.addParticle(R3{-4.0, 5.5, 0.0}, R3{1.0, -0.4, 2.5}, 2.0, 5.0, 0.3);
    }

    void TearDown() override {
        for (const auto& path : created_files) {
            std::error_code ec;
            std::filesystem::remove(path, ec);
        }
        created_files.clear();
    }

    /**
     * @brief Writes a checkpoint file using the YAMLWriterCP and returns the path to the created file, as well as
     * tracks the created files for cleanup.
     */
    std::filesystem::path writeCheckpoint(const Domain& domain, SimpleContainer& container, int iteration, PairwiseForce force,
                                          double delta_t, double current_time, double end_time, size_t frequency_output,
                                          size_t frequency_checkpoint, const std::string& base_name,
                                          double cutoff_radius, double target_temp, double delta_temp,
                                          size_t thermostat_freq, size_t dimensions, size_t N) {
        ContainerRef particles(container);
        SettingsParam settings;
        settings.delta_t = delta_t;
        settings.end_time = end_time;
        settings.start_time = current_time;
        settings.base_name = base_name;
        settings.pairwise_forces = {force};
        settings.container_type = "SIMPLE";
        settings.frequency_output = frequency_output;
        settings.frequency_checkpoint = frequency_checkpoint;
        settings.cutoff = cutoff_radius;
        settings.target_temp = target_temp;
        settings.delta_temp = delta_temp;
        settings.thermostat_freq = thermostat_freq;
        settings.dimensions = dimensions;
        settings.thermo = true;
        writer.createCheckpoint(settings, domain, particles, iteration, N);

        int decimal_places = 0;
        while (N >= 10) {
            N /= 10;
            decimal_places++;
        }

        std::ostringstream oss;
        oss << "cp_" << std::setfill('0') << std::setw(decimal_places) << iteration << ".yaml";
        std::filesystem::path path = std::filesystem::current_path() / oss.str();
        created_files.push_back(path);

        return path;
    }
};

TEST_F(YAMLWriterCPTest, testWritesSimpleSettingsAndParticleData) {  // NOLINT
    constexpr int iteration = 7;
    constexpr size_t iteration_cap = 100;
    constexpr double delta_t = 0.02;
    constexpr double start_time = 1.0;
    constexpr double end_time = 5.0;
    constexpr size_t frequency_output = 5;
    constexpr size_t frequency_checkpoint = 25;
    const std::string base_name = "test_simple";
    constexpr double cutoff_radius = 1.25;
    constexpr double target_temp = 300.0;
    constexpr double delta_temp = 0.5;
    constexpr size_t thermostat_freq = 10;
    constexpr size_t dimensions = 3;
    const Domain domain{R3{10.0, 10.0, 10.0}};

    auto file_path = writeCheckpoint(domain, container, iteration, PairwiseForce::LENNARDJONES, delta_t, start_time, end_time,
                                     frequency_output, frequency_checkpoint, base_name, cutoff_radius, target_temp,
                                     delta_temp, thermostat_freq, dimensions, iteration_cap);

    ASSERT_TRUE(std::filesystem::exists(file_path));

    YAML::Node root = YAML::LoadFile(file_path.string());
    ASSERT_TRUE(root["settings"]);
    const YAML::Node settings_node = root["settings"];
    ASSERT_TRUE(settings_node);

    EXPECT_EQ("Settings", settings_node["format"].as<std::string>());
    EXPECT_DOUBLE_EQ(delta_t, settings_node["delta_t"].as<double>());
    EXPECT_DOUBLE_EQ(end_time, settings_node["end_time"].as<double>());
    EXPECT_DOUBLE_EQ(start_time, settings_node["start_time"].as<double>());
    EXPECT_EQ(base_name, settings_node["base_name"].as<std::string>());
    
    // Check pairwise_forces instead of deprecated "force"
    ASSERT_TRUE(settings_node["pairwise_forces"]);
    const YAML::Node pairwise_forces = settings_node["pairwise_forces"];
    ASSERT_TRUE(pairwise_forces.IsSequence());
    ASSERT_EQ(1, pairwise_forces.size());
    EXPECT_EQ("LENNARDJONES", pairwise_forces[0].as<std::string>());
    
    EXPECT_EQ("SIMPLE", settings_node["container"].as<std::string>());
    EXPECT_EQ("NAIVE", settings_node["strategy"].as<std::string>());
    EXPECT_EQ(frequency_output, settings_node["frequency"].as<size_t>());
    EXPECT_EQ(frequency_checkpoint, settings_node["checkpoint"].as<size_t>());
    EXPECT_DOUBLE_EQ(cutoff_radius, settings_node["cutoff"].as<double>());

    // Check thermostat settings (now nested)
    ASSERT_TRUE(settings_node["thermostat"]);
    const YAML::Node thermostat_node = settings_node["thermostat"];
    EXPECT_DOUBLE_EQ(target_temp, thermostat_node["target_temp"].as<double>());
    EXPECT_DOUBLE_EQ(delta_temp, thermostat_node["delta_temp"].as<double>());
    EXPECT_EQ(thermostat_freq, thermostat_node["n_thermostat"].as<size_t>());

    const YAML::Node domain_node = settings_node["domain"];
    ASSERT_TRUE(domain_node);
    const YAML::Node coords_node = domain_node["coordinates"];
    ASSERT_TRUE(coords_node);
    ASSERT_TRUE(coords_node.IsSequence());
    EXPECT_EQ(3, coords_node.size());
    EXPECT_DOUBLE_EQ(domain.getDimension()[0], coords_node[0].as<double>());
    EXPECT_DOUBLE_EQ(domain.getDimension()[1], coords_node[1].as<double>());
    EXPECT_DOUBLE_EQ(domain.getDimension()[2], coords_node[2].as<double>());
    EXPECT_EQ(dimensions, domain_node["dimensions"].as<size_t>());

    const YAML::Node boundaries = domain_node["boundaries"];
    ASSERT_TRUE(boundaries);
    EXPECT_EQ("OUTFLOW", boundaries["left"]["type"].as<std::string>());
    EXPECT_EQ("OUTFLOW", boundaries["right"]["type"].as<std::string>());
    EXPECT_EQ("OUTFLOW", boundaries["upper"]["type"].as<std::string>());
    EXPECT_EQ("OUTFLOW", boundaries["lower"]["type"].as<std::string>());
    EXPECT_EQ("OUTFLOW", boundaries["front"]["type"].as<std::string>());
    EXPECT_EQ("OUTFLOW", boundaries["back"]["type"].as<std::string>());

    ASSERT_TRUE(root["xvm_data"]);
    const YAML::Node xvm_data = root["xvm_data"];
    EXPECT_EQ("XVM", xvm_data["format"].as<std::string>());
    EXPECT_EQ(container.size(), xvm_data["num_particles"].as<size_t>());

    const YAML::Node particles_node = xvm_data["particles"];
    ASSERT_TRUE(particles_node.IsSequence());
    ASSERT_EQ(container.size(), particles_node.size());
    size_t idx = 0;
    for (const auto& expected : container) {
        const YAML::Node particle_node = particles_node[idx++];
        ASSERT_TRUE(particle_node);
        
        const YAML::Node coords = particle_node["coordinates"];
        ASSERT_TRUE(coords.IsSequence());
        EXPECT_EQ(3, coords.size());
        EXPECT_DOUBLE_EQ(expected.getX()[0], coords[0].as<double>());
        EXPECT_DOUBLE_EQ(expected.getX()[1], coords[1].as<double>());
        EXPECT_DOUBLE_EQ(expected.getX()[2], coords[2].as<double>());
        
        const YAML::Node old_coords = particle_node["old_coordinates"];
        ASSERT_TRUE(old_coords.IsSequence());
        EXPECT_EQ(3, old_coords.size());
        EXPECT_DOUBLE_EQ(expected.getOldX()[0], old_coords[0].as<double>());
        EXPECT_DOUBLE_EQ(expected.getOldX()[1], old_coords[1].as<double>());
        EXPECT_DOUBLE_EQ(expected.getOldX()[2], old_coords[2].as<double>());
        
        const YAML::Node velocity = particle_node["velocity"];
        ASSERT_TRUE(velocity.IsSequence());
        EXPECT_EQ(3, velocity.size());
        EXPECT_DOUBLE_EQ(expected.getV()[0], velocity[0].as<double>());
        EXPECT_DOUBLE_EQ(expected.getV()[1], velocity[1].as<double>());
        EXPECT_DOUBLE_EQ(expected.getV()[2], velocity[2].as<double>());
        
        const YAML::Node force = particle_node["force"];
        ASSERT_TRUE(force.IsSequence());
        EXPECT_EQ(3, force.size());
        EXPECT_DOUBLE_EQ(expected.getF()[0], force[0].as<double>());
        EXPECT_DOUBLE_EQ(expected.getF()[1], force[1].as<double>());
        EXPECT_DOUBLE_EQ(expected.getF()[2], force[2].as<double>());
        
        const YAML::Node old_force = particle_node["old_force"];
        ASSERT_TRUE(old_force.IsSequence());
        EXPECT_EQ(3, old_force.size());
        EXPECT_DOUBLE_EQ(expected.getOldF()[0], old_force[0].as<double>());
        EXPECT_DOUBLE_EQ(expected.getOldF()[1], old_force[1].as<double>());
        EXPECT_DOUBLE_EQ(expected.getOldF()[2], old_force[2].as<double>());
        EXPECT_DOUBLE_EQ(expected.getM(), particle_node["mass"].as<double>());
        EXPECT_DOUBLE_EQ(expected.getEpsilon(), particle_node["epsilon"].as<double>());
        EXPECT_DOUBLE_EQ(expected.getSigma(), particle_node["sigma"].as<double>());
        EXPECT_EQ(expected.getType(), particle_node["type"].as<int>());
    }
}

TEST_F(YAMLWriterCPTest, testEmptyContainer) {  // NOLINT
    SimpleContainer empty;
    constexpr int iteration = 0;
    constexpr size_t iteration_cap = 1;
    const Domain domain{R3{10.0, 10.0, 10.0}};
    auto file_path = writeCheckpoint(domain, empty, iteration, PairwiseForce::GRAVITATIONAL, 0.01, 0.0, 1.0, 1, 1, "test_base",
                                     1.0, 1.0, 1.0, 10, 3, iteration_cap);

    ASSERT_TRUE(std::filesystem::exists(file_path));
    YAML::Node root = YAML::LoadFile(file_path.string());
    const YAML::Node xvm_data = root["xvm_data"];
    ASSERT_TRUE(xvm_data);
    EXPECT_EQ(0U, xvm_data["num_particles"].as<size_t>());
    const YAML::Node particles_node = xvm_data["particles"];
    ASSERT_TRUE(particles_node.IsSequence());
    EXPECT_EQ(0U, particles_node.size());
}

TEST_F(YAMLWriterCPTest, testWritesReflectingBoundaryMetadata) {  // NOLINT
    const R3 dimensions{15.0, 12.0, 6.0};
    std::array<std::unique_ptr<Boundary>, 6> bds{};
    bds[0] = std::make_unique<Reflecting>(BoundaryLocation::LEFT, dimensions, true, 1.1, 5.5);
    bds[1] = std::make_unique<Reflecting>(BoundaryLocation::RIGHT, dimensions, false);
    bds[2] = std::make_unique<Outflow>(BoundaryLocation::FRONT, dimensions);
    bds[3] = std::make_unique<Outflow>(BoundaryLocation::BACK, dimensions);
    bds[4] = std::make_unique<Outflow>(BoundaryLocation::UPPER, dimensions);
    bds[5] = std::make_unique<Outflow>(BoundaryLocation::LOWER, dimensions);
    Domain reflecting_domain = {dimensions, std::move(bds)};

    constexpr int iteration = 5;
    constexpr size_t iteration_cap = 10;
    auto file_path = writeCheckpoint(reflecting_domain, container, iteration, PairwiseForce::GRAVITATIONAL, 0.02, 0.0, 2.0, 2,
                                     4, "test_reflect", 2.0, 1.0, 1.0, 10, 3, iteration_cap);

    YAML::Node root = YAML::LoadFile(file_path.string());
    const YAML::Node boundaries = root["settings"]["domain"]["boundaries"];
    ASSERT_TRUE(boundaries);

    const YAML::Node left = boundaries["left"];
    ASSERT_TRUE(left);
    EXPECT_EQ("REFLECTING", left["type"].as<std::string>());
    EXPECT_TRUE(left["ghost_on_boundary"].as<bool>());
    EXPECT_DOUBLE_EQ(1.1, left["sigma"].as<double>());
    EXPECT_DOUBLE_EQ(5.5, left["epsilon"].as<double>());

    const YAML::Node right = boundaries["right"];
    ASSERT_TRUE(right);
    EXPECT_EQ("REFLECTING", right["type"].as<std::string>());
    EXPECT_FALSE(right["ghost_on_boundary"].as<bool>());
    EXPECT_FALSE(right["sigma"]);
    EXPECT_FALSE(right["epsilon"]);

    const YAML::Node front = boundaries["front"];
    ASSERT_TRUE(front);
    EXPECT_EQ("OUTFLOW", front["type"].as<std::string>());
}

TEST_F(YAMLWriterCPTest, testReadBackCheckpointWithYAMLReader) {  // NOLINT
    constexpr int iteration = 12;
    constexpr size_t iteration_cap = 1000;
    constexpr double delta_t = 0.005;
    constexpr double start_time = 2.0;
    constexpr double end_time = 20.0;
    constexpr size_t frequency_output = 4;
    constexpr size_t frequency_checkpoint = 12;
    const std::string base_name = "test_restart";
    constexpr double cutoff_radius = 0.8;
    constexpr double target_temp = 300.0;
    constexpr double delta_temp = 0.5;
    constexpr size_t thermostat_freq = 10;
    constexpr size_t dimensions = 3;
    const Domain domain{R3{10.0, 10.0, 10.0}};

    auto file_path = writeCheckpoint(domain, container, iteration, PairwiseForce::LENNARDJONES, delta_t, start_time, end_time,
                                     frequency_output, frequency_checkpoint, base_name, cutoff_radius, target_temp,
                                     delta_temp, thermostat_freq, dimensions, iteration_cap);

    // Settings
    SettingsParam settings;
    YAMLReader reader;
    reader.readSettings(settings, file_path.string());
    EXPECT_DOUBLE_EQ(delta_t, settings.delta_t);
    EXPECT_DOUBLE_EQ(start_time, settings.start_time);
    EXPECT_DOUBLE_EQ(end_time, settings.end_time);
    EXPECT_EQ(base_name, settings.base_name);
    EXPECT_EQ(1, settings.pairwise_forces.size());
    EXPECT_EQ(PairwiseForce::LENNARDJONES, settings.pairwise_forces[0]);
    EXPECT_EQ("SIMPLE", settings.container_type);
    EXPECT_EQ(frequency_output, settings.frequency_output);
    EXPECT_EQ(frequency_checkpoint, settings.frequency_checkpoint);
    EXPECT_DOUBLE_EQ(cutoff_radius, settings.cutoff);
    EXPECT_DOUBLE_EQ(target_temp, settings.target_temp);
    EXPECT_DOUBLE_EQ(delta_temp, settings.delta_temp);
    EXPECT_EQ(thermostat_freq, settings.thermostat_freq);
    EXPECT_DOUBLE_EQ(domain.getDimension()[0], settings.domain.getDimension()[0]);
    EXPECT_DOUBLE_EQ(domain.getDimension()[1], settings.domain.getDimension()[1]);
    EXPECT_DOUBLE_EQ(domain.getDimension()[2], settings.domain.getDimension()[2]);

    // Particle data
    SimpleContainer restored;
    ContainerRef restored_ref(restored);
    SettingsParam dummy_settings;
    reader.readParticles(restored_ref, dummy_settings, file_path.string());

    ASSERT_EQ(restored.size(), container.size());

    for (auto& particle : container) {
        EXPECT_TRUE(std::ranges::any_of(restored, [&particle](const auto& p) { return p == particle; }));
    }
}

TEST_F(YAMLWriterCPTest, testWritesSingleForces) {  // NOLINT
    constexpr int iteration = 3;
    constexpr size_t iteration_cap = 10;
    const Domain domain{R3{10.0, 10.0, 10.0}};
    
    ContainerRef particles(container);
    SettingsParam settings;
    settings.delta_t = 0.01;
    settings.end_time = 1.0;
    settings.start_time = 0.0;
    settings.base_name = "test_single";
    settings.pairwise_forces = {PairwiseForce::LENNARDJONES};
    settings.single_forces = {SingleForce::GRAV, SingleForce::HARMONIC};
    settings.g_grav_vec = {0.0, -9.81, 0.0};
    settings.k = 300.0;
    settings.r_0 = 2.2;
    settings.container_type = "SIMPLE";
    settings.frequency_output = 1;
    settings.frequency_checkpoint = 1;
    settings.cutoff = 1.0;
    settings.dimensions = 3;
    settings.thermo = false;
    
    writer.createCheckpoint(settings, domain, particles, iteration, iteration_cap);
    
    std::ostringstream oss;
    oss << "cp_" << std::setfill('0') << std::setw(1) << iteration << ".yaml";
    std::filesystem::path path = std::filesystem::current_path() / oss.str();
    created_files.push_back(path);
    
    ASSERT_TRUE(std::filesystem::exists(path));
    YAML::Node root = YAML::LoadFile(path.string());
    
    const YAML::Node settings_node = root["settings"];
    ASSERT_TRUE(settings_node["single_forces"]);
    const YAML::Node single_forces = settings_node["single_forces"];
    ASSERT_TRUE(single_forces.IsSequence());
    ASSERT_EQ(2, single_forces.size());
    
    const YAML::Node grav_force = single_forces[0];
    EXPECT_EQ("GRAV", grav_force["type"].as<std::string>());
    ASSERT_TRUE(grav_force["g_grav"]);
    const YAML::Node g_grav = grav_force["g_grav"];
    ASSERT_TRUE(g_grav.IsSequence());
    EXPECT_DOUBLE_EQ(0.0, g_grav[0].as<double>());
    EXPECT_DOUBLE_EQ(-9.81, g_grav[1].as<double>());
    EXPECT_DOUBLE_EQ(0.0, g_grav[2].as<double>());
    
    const YAML::Node harmonic_force = single_forces[1];
    EXPECT_EQ("HARMONIC", harmonic_force["type"].as<std::string>());
    EXPECT_DOUBLE_EQ(300.0, harmonic_force["k"].as<double>());
    EXPECT_DOUBLE_EQ(2.2, harmonic_force["r_0"].as<double>());
}

TEST_F(YAMLWriterCPTest, testWritesTargetForce) {  // NOLINT
    constexpr int iteration = 2;
    constexpr size_t iteration_cap = 5;
    const Domain domain{R3{10.0, 10.0, 10.0}};
    
    ContainerRef particles(container);
    SettingsParam settings;
    settings.delta_t = 0.01;
    settings.end_time = 1.0;
    settings.start_time = 0.0;
    settings.base_name = "test_target";
    settings.pairwise_forces = {PairwiseForce::LENNARDJONES};
    settings.target_force_enabled = true;
    settings.target_force_direction = {1.0, 0.0, 0.0};
    settings.target_force_magnitude = 50.0;
    settings.target_force_max_iterations = 1000;
    settings.container_type = "SIMPLE";
    settings.frequency_output = 1;
    settings.frequency_checkpoint = 1;
    settings.cutoff = 1.0;
    settings.dimensions = 3;
    settings.thermo = false;
    
    writer.createCheckpoint(settings, domain, particles, iteration, iteration_cap);
    
    std::ostringstream oss;
    oss << "cp_" << std::setfill('0') << std::setw(1) << iteration << ".yaml";
    std::filesystem::path path = std::filesystem::current_path() / oss.str();
    created_files.push_back(path);
    
    ASSERT_TRUE(std::filesystem::exists(path));
    YAML::Node root = YAML::LoadFile(path.string());
    
    const YAML::Node settings_node = root["settings"];
    ASSERT_TRUE(settings_node["target_force"]);
    const YAML::Node target_force = settings_node["target_force"];
    
    ASSERT_TRUE(target_force["direction"]);
    const YAML::Node direction = target_force["direction"];
    ASSERT_TRUE(direction.IsSequence());
    EXPECT_DOUBLE_EQ(1.0, direction[0].as<double>());
    EXPECT_DOUBLE_EQ(0.0, direction[1].as<double>());
    EXPECT_DOUBLE_EQ(0.0, direction[2].as<double>());
    
    EXPECT_DOUBLE_EQ(50.0, target_force["magnitude"].as<double>());
    EXPECT_EQ(1000, target_force["max_iterations"].as<size_t>());
}

TEST_F(YAMLWriterCPTest, testWritesInitialTemp) {  // NOLINT
    constexpr int iteration = 1;
    constexpr size_t iteration_cap = 5;
    const Domain domain{R3{10.0, 10.0, 10.0}};
    
    ContainerRef particles(container);
    SettingsParam settings;
    settings.delta_t = 0.01;
    settings.end_time = 1.0;
    settings.start_time = 0.0;
    settings.base_name = "test_init_temp";
    settings.pairwise_forces = {PairwiseForce::LENNARDJONES};
    settings.container_type = "SIMPLE";
    settings.frequency_output = 1;
    settings.frequency_checkpoint = 1;
    settings.cutoff = 1.0;
    settings.dimensions = 3;
    settings.thermo = true;
    settings.init_temp = 200.0;
    settings.target_temp = 300.0;
    settings.delta_temp = 0.5;
    settings.thermostat_freq = 10;
    
    writer.createCheckpoint(settings, domain, particles, iteration, iteration_cap);
    
    std::ostringstream oss;
    oss << "cp_" << std::setfill('0') << std::setw(1) << iteration << ".yaml";
    std::filesystem::path path = std::filesystem::current_path() / oss.str();
    created_files.push_back(path);
    
    ASSERT_TRUE(std::filesystem::exists(path));
    YAML::Node root = YAML::LoadFile(path.string());
    
    const YAML::Node settings_node = root["settings"];
    ASSERT_TRUE(settings_node["thermostat"]);
    const YAML::Node thermostat = settings_node["thermostat"];
    
    EXPECT_DOUBLE_EQ(200.0, thermostat["initial_temp"].as<double>());
    EXPECT_DOUBLE_EQ(300.0, thermostat["target_temp"].as<double>());
}

TEST_F(YAMLWriterCPTest, testReadBackCheckpointWithSingleForces) {  // NOLINT
    constexpr int iteration = 5;
    constexpr size_t iteration_cap = 10;
    const Domain domain{R3{10.0, 10.0, 10.0}};
    
    ContainerRef particles(container);
    SettingsParam settings;
    settings.delta_t = 0.01;
    settings.end_time = 1.0;
    settings.start_time = 0.0;
    settings.base_name = "test_readback";
    settings.pairwise_forces = {PairwiseForce::LENNARDJONES};
    settings.single_forces = {SingleForce::GRAV, SingleForce::HARMONIC};
    settings.g_grav_vec = {0.0, -9.81, 0.0};
    settings.k = 300.0;
    settings.r_0 = 2.2;
    settings.target_force_enabled = true;
    settings.target_force_direction = {1.0, 0.0, 0.0};
    settings.target_force_magnitude = 50.0;
    settings.target_force_max_iterations = 1000;
    settings.container_type = "SIMPLE";
    settings.frequency_output = 1;
    settings.frequency_checkpoint = 1;
    settings.cutoff = 1.0;
    settings.dimensions = 3;
    settings.thermo = false;
    
    writer.createCheckpoint(settings, domain, particles, iteration, iteration_cap);
    
    std::ostringstream oss;
    oss << "cp_" << std::setfill('0') << std::setw(1) << iteration << ".yaml";
    std::filesystem::path path = std::filesystem::current_path() / oss.str();
    created_files.push_back(path);
    
    // Read back
    SettingsParam read_settings;
    YAMLReader reader;
    reader.readSettings(read_settings, path.string());
    
    ASSERT_EQ(2, read_settings.single_forces.size());
    EXPECT_EQ(SingleForce::GRAV, read_settings.single_forces[0]);
    EXPECT_EQ(SingleForce::HARMONIC, read_settings.single_forces[1]);
    EXPECT_DOUBLE_EQ(0.0, read_settings.g_grav_vec[0]);
    EXPECT_DOUBLE_EQ(-9.81, read_settings.g_grav_vec[1]);
    EXPECT_DOUBLE_EQ(0.0, read_settings.g_grav_vec[2]);
    EXPECT_DOUBLE_EQ(300.0, read_settings.k);
    EXPECT_DOUBLE_EQ(2.2, read_settings.r_0);
    
    EXPECT_TRUE(read_settings.target_force_enabled);
    EXPECT_DOUBLE_EQ(1.0, read_settings.target_force_direction[0]);
    EXPECT_DOUBLE_EQ(0.0, read_settings.target_force_direction[1]);
    EXPECT_DOUBLE_EQ(0.0, read_settings.target_force_direction[2]);
    EXPECT_DOUBLE_EQ(50.0, read_settings.target_force_magnitude);
    EXPECT_EQ(1000, read_settings.target_force_max_iterations);
}

}  // namespace mol_sim