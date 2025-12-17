#include "particles/boundaries/Periodic.h"

#include <gtest/gtest.h>

#include <limits>
#include <memory>

#include "io/checkpointWriter/XVMWriterCP.h"
#include "io/outputWriter/XYZWriter.h"
#include "particles/Particle.h"
#include "particles/boundaries/Boundary.h"
#include "physics/LennardJonesForce.h"
#include "utils/Vector.h"
#include "utils/Settings.h"
#include "utils/Simulation.h"
#include "particles/boundaries/Reflecting.h"

namespace mol_sim {
/**
 * @brief Test Fixture for testing the Outflow boundary condition
 * Base data/setup:
 * settings: default values below, other values set in tests
 * p:
 *  position: (5,5,5)
 *  velocity: (0,0,0)
 *  (mass: 1)
 * dimension: (10,10,10)
 * delta_t: 0.5
 * end_time = 1
 */
class PeriodicTest : public testing::Test {
   protected:
    R3 dimension;
    Periodic boundary;

    PeriodicTest() 
    : dimension{10.0, 10.0, 10.0}, boundary{BoundaryLocation::LEFT, dimension, 1.0, false} {};
};

TEST_F(PeriodicTest, ParticleExactlyOnBoundaryCopyTestLimitedCutoff) {
    boundary = Periodic(BoundaryLocation::LEFT, dimension, 1.0, false);
    R3 x = {.0, 5.0, 5.0};
    R3 v = {.0, .0, .0};
    Particle p = Particle(x, v, 1.0, 1.0, 1.0, 0);
    const ForceSource& force_source = LennardJonesForce();
    auto new_particles = boundary.applyBoundary(p, force_source);
   
    Particle p1(p);
    p1.getX() = {10.0, 5.0, 5.0};
    p1.getType() = 1;
    EXPECT_EQ(new_particles.value().size(), 1);
    EXPECT_EQ(p1, new_particles.value()[0]);
}

/* TEST_F(PeriodicTest, SimulationParticleExactlyOnBoundaryCopyTestLimitedCutoff) {
    double cutoff = 3.0;
    SettingsParam settings;
    std::array<std::unique_ptr<Boundary>, 6> boundaries;
        boundaries[0] = std::make_unique<Periodic>(BoundaryLocation::LEFT, dimension, cutoff);
        boundaries[1] = std::make_unique<Periodic>(BoundaryLocation::RIGHT, dimension, cutoff);
        boundaries[2] = std::make_unique<Reflecting>(BoundaryLocation::FRONT, dimension, false);
        boundaries[3] = std::make_unique<Reflecting>(BoundaryLocation::BACK, dimension, false);
        settings.domain = Domain(dimension, std::move(boundaries));
        settings.delta_t = 0.0005;
        settings.end_time = 0.0005;


  domain:
    x: 63
    y: 36
    z: 1
    g_grav: -12.44
    boundaries:
      left:
        type: Periodic
      right:
        type: Periodic
      front:
        type: Reflecting
      back:
        type: Reflectingd_time = 1;
    settings.container_type = "LINKED";
    LinkedCellContainer particles(dimension, cutoff);
    R3 x = {.5, 5.0, 5.0};
    R3 v = {.0, .0, .0};
    Particle p = Particle(x, v, 1.0, 1.0, 1.0, 0);
    particles.addParticle(p);

    auto force_source = std::make_unique<LennardJonesForce>();
    auto writer = std::make_unique<XYZWriter>();
    auto cp_writer = std::make_unique<XVMWriterCP>();
    Simulation<LinkedCellContainer> simulation(particles, *force_source, settings, *writer, *cp_writer);
    simulation.run();
   
    Particle p1(p);
    p1.getX() = {10.0, 5.0, 5.0};
    p1.getType() = 1;
}
 */
TEST_F(PeriodicTest, ParticleExactlyOnBorderInnerBoundaryCopyTestLimitedCutoff) {
    boundary = Periodic(BoundaryLocation::LEFT, dimension, 1.0, false);
    R3 x = {1.0, 5.0, 5.0};
    R3 v = {.0, .0, .0};
    Particle p = Particle(x, v, 1.0, 1.0, 1.0, 0);
    const ForceSource& force_source = LennardJonesForce();
    auto new_particles = boundary.applyBoundary(p, force_source);
   
    Particle p1(p);
    p1.getX() = {11.0, 5.0, 5.0};
    p1.getType() = 1;
    EXPECT_EQ(new_particles.value().size(), 1);
    EXPECT_EQ(p1, new_particles.value()[0]);
}

TEST_F(PeriodicTest, ParticleInInnerCellNoCopyTestLimitedCutoff) {
    boundary = Periodic(BoundaryLocation::LEFT, dimension, 1.0, false);
    R3 x = {5.0, 5.0, 5.0};
    R3 v = {.0, .0, .0};
    Particle p = Particle(x, v, 1.0, 1.0, 1.0, 0);
    const ForceSource& force_source = LennardJonesForce();
    auto new_particles = boundary.applyBoundary(p, force_source);
   
    EXPECT_TRUE(new_particles == std::nullopt);
}

TEST_F(PeriodicTest, ParticleOnLeftBoundaryCopyTestLimitedCutoff) {
    boundary = Periodic(BoundaryLocation::LEFT, dimension, 1.0, false);
    R3 x = {.5, 5.0, 5.0};
    R3 v = {.0, .0, .0};
    Particle p = Particle(x, v, 1.0, 1.0, 1.0, 0);
    const ForceSource& force_source = LennardJonesForce();
    auto new_particles = boundary.applyBoundary(p, force_source);
   
    Particle p1(p);
    p1.getX() = {10.5, 5.0, 5.0};
    p1.getType() = 1;
    EXPECT_EQ(new_particles.value().size(), 1);
    EXPECT_EQ(p1, new_particles.value()[0]);
}

TEST_F(PeriodicTest, ParticleOnEdgeCopyTestLimitedCutoff) {
    boundary = Periodic(BoundaryLocation::LEFT, dimension, 1.0, false);
    R3 x = {.0, .0, 5.0};
    R3 v = {.0, .0, .0};
    Particle p = Particle(x, v, 1.0, 1.0, 1.0, 0);
    const ForceSource& force_source = LennardJonesForce();
    auto new_particles = boundary.applyBoundary(p, force_source);
   
    Particle p1(p);
    Particle p2(p);
    Particle p3(p);
    p1.getX() = {10.0, .0, 5.0};
    p2.getX() = {10.0, 10.0, 5.0};
    p3.getX() = {.0, 10.0, 5.0};
    p1.getType() = 1;
    p2.getType() = 1;
    p3.getType() = 1;
    EXPECT_TRUE(new_particles.value().size() == 3);
    EXPECT_TRUE(std::find(new_particles.value().begin(), new_particles.value().end(), p1) != new_particles.value().end());
    EXPECT_TRUE(std::find(new_particles.value().begin(), new_particles.value().end(), p2) != new_particles.value().end());
    EXPECT_TRUE(std::find(new_particles.value().begin(), new_particles.value().end(), p3) != new_particles.value().end());
}

TEST_F(PeriodicTest, ParticleOnCornerCopyTestLimitedCutoff) {
    boundary = Periodic(BoundaryLocation::LEFT, dimension, 1.0, false);
    R3 x = {.0, .0, .0};
    R3 v = {.0, .0, .0};
    Particle p = Particle(x, v, 1.0, 1.0, 1.0, 0);
    const ForceSource& force_source = LennardJonesForce();
    auto new_particles = boundary.applyBoundary(p, force_source);
   
    Particle p1(p);
    Particle p2(p);
    Particle p3(p);
    Particle p4(p);
    Particle p5(p);
    Particle p6(p);
    Particle p7(p);
    p1.getX() = {10.0, .0, .0};
    p2.getX() = {10.0, .0, 10.0};
    p3.getX() = {10.0, 10.0, .0};
    p4.getX() = {10.0, 10.0, 10.0};
    p5.getX() = {.0, 10.0, .0};
    p6.getX() = {.0, 10.0, 10.0};
    p7.getX() = {.0, .0, 10.0};
    p1.getType() = 1;
    p2.getType() = 1;
    p3.getType() = 1;
    p4.getType() = 1;
    p5.getType() = 1;
    p6.getType() = 1;
    p7.getType() = 1;
    EXPECT_TRUE(new_particles.value().size() == 7);
    EXPECT_TRUE(std::find(new_particles.value().begin(), new_particles.value().end(), p1) != new_particles.value().end());
    EXPECT_TRUE(std::find(new_particles.value().begin(), new_particles.value().end(), p2) != new_particles.value().end());
    EXPECT_TRUE(std::find(new_particles.value().begin(), new_particles.value().end(), p3) != new_particles.value().end());
    EXPECT_TRUE(std::find(new_particles.value().begin(), new_particles.value().end(), p4) != new_particles.value().end());
    EXPECT_TRUE(std::find(new_particles.value().begin(), new_particles.value().end(), p5) != new_particles.value().end());
    EXPECT_TRUE(std::find(new_particles.value().begin(), new_particles.value().end(), p6) != new_particles.value().end());
    EXPECT_TRUE(std::find(new_particles.value().begin(), new_particles.value().end(), p7) != new_particles.value().end());
}

TEST_F(PeriodicTest, ParticleMoveTestLimitedCutoff) {
    boundary = Periodic(BoundaryLocation::LEFT, dimension, 1.0, false);
    R3 x = {-.5, 5.0, 5.0};
    R3 v = {.0, .0, .0};
    Particle p = Particle(x, v, 1.0, 1.0, 1.0, 0);
    const ForceSource& force_source = LennardJonesForce();
    auto new_particles = boundary.applyBoundary(p, force_source);

    Particle p_exp(p);
    p_exp.getX() = {9.5, 5.0, 5.0};
    EXPECT_EQ(p, p_exp);
}

TEST_F(PeriodicTest, _2DtestTeleport) {
    boundary = Periodic(BoundaryLocation::LEFT, dimension, 1.0, true);
    R3 x = {-.5, 5.0, 0.0};
    R3 v = {.0, .0, .0};
    Particle p = Particle(x, v, 1.0, 1.0, 1.0, 0);
    const ForceSource& force_source = LennardJonesForce();
    auto new_particles = boundary.applyBoundary(p, force_source);

    Particle p_exp(p);
    p_exp.getX() = {9.5, 5.0, 0.0};
    EXPECT_EQ(p, p_exp);
}

TEST_F(PeriodicTest, _2DtestMirrorSimple) {
    dimension = {10.0, 20.0, 1.0};
    boundary = Periodic(BoundaryLocation::LEFT, dimension, 1.0, true);
    R3 x = {.5, 5.0, 0.0};
    R3 v = {.0, .0, .0};
    Particle p = Particle(x, v, 1.0, 1.0, 1.0, 0);
    const ForceSource& force_source = LennardJonesForce();
    auto new_particles = boundary.applyBoundary(p, force_source);

    Particle p_exp(p);
    p_exp.getX() = {10.5, 5.0, 0.0};
    p_exp.getType() = 1;
    EXPECT_TRUE(new_particles.value().size() == 1);
    EXPECT_EQ(new_particles.value()[0], p_exp);
}

TEST_F(PeriodicTest, _2DtestMirrorOnBoundary) {
    dimension = {10.0, 20.0, 1.0};
    boundary = Periodic(BoundaryLocation::LEFT, dimension, 1.0, true);
    R3 x = {.0, 5.0, 0.0};
    R3 v = {.0, .0, .0};
    Particle p = Particle(x, v, 1.0, 1.0, 1.0, 0);
    const ForceSource& force_source = LennardJonesForce();
    auto new_particles = boundary.applyBoundary(p, force_source);

    Particle p_exp(p);
    p_exp.getX() = {10.0, 5.0, 0.0};
    p_exp.getType() = 1;
    EXPECT_TRUE(new_particles.value().size() == 1);
    EXPECT_EQ(new_particles.value()[0], p_exp);
}

}  // namespace mol_sim
