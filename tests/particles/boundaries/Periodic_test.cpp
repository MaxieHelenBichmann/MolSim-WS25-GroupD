#include "particles/boundaries/Periodic.h"

#include <gtest/gtest.h>

#include <algorithm>

#include "particles/Particle.h"
#include "particles/boundaries/Boundary.h"
#include "physics/LennardJonesForce.h"
#include "utils/Vector.h"

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

    PeriodicTest() : dimension{10.0, 10.0, 10.0}, boundary{BoundaryLocation::LEFT, dimension, 1.0, false} {};
};

TEST_F(PeriodicTest, ParticleExactlyOnBoundaryCopyTestLimitedCutoff) {
    boundary = Periodic(BoundaryLocation::LEFT, dimension, 1.0, 3);
    R3 x = {.0, 5.0, 5.0};
    R3 v = {.0, .0, .0};
    Particle p = Particle(x, v, 1.0, 1.0, 1.0, 0);
    const ForceSource& force_source = LennardJonesForce();
    auto new_particles = boundary.applyBoundary(p, force_source);

    Particle p1(p);
    p1.getX() = {10.0, 5.0, 5.0};
    p1.getType() = 1;
    auto np = new_particles.value_or(std::vector<Particle>());
    EXPECT_EQ(np.size(), 1);
    EXPECT_EQ(p1, np[0]);
}

TEST_F(PeriodicTest, ParticleExactlyOnBorderInnerBoundaryCopyTestLimitedCutoff) {
    boundary = Periodic(BoundaryLocation::LEFT, dimension, 1.0, 3);
    R3 x = {1.0, 5.0, 5.0};
    R3 v = {.0, .0, .0};
    Particle p = Particle(x, v, 1.0, 1.0, 1.0, 0);
    const ForceSource& force_source = LennardJonesForce();
    auto new_particles = boundary.applyBoundary(p, force_source);

    Particle p1(p);
    p1.getX() = {11.0, 5.0, 5.0};
    p1.getType() = 1;
    auto np = new_particles.value_or(std::vector<Particle>());
    EXPECT_EQ(np.size(), 1);
    EXPECT_EQ(p1, np[0]);
}

TEST_F(PeriodicTest, ParticleInInnerCellNoCopyTestLimitedCutoff) {
    boundary = Periodic(BoundaryLocation::LEFT, dimension, 1.0, 3);
    R3 x = {5.0, 5.0, 5.0};
    R3 v = {.0, .0, .0};
    Particle p = Particle(x, v, 1.0, 1.0, 1.0, 0);
    const ForceSource& force_source = LennardJonesForce();
    auto new_particles = boundary.applyBoundary(p, force_source);

    EXPECT_TRUE(new_particles == std::nullopt);
}

TEST_F(PeriodicTest, ParticleOnLeftBoundaryCopyTestLimitedCutoff) {
    boundary = Periodic(BoundaryLocation::LEFT, dimension, 1.0, 3);
    R3 x = {.5, 5.0, 5.0};
    R3 v = {.0, .0, .0};
    Particle p = Particle(x, v, 1.0, 1.0, 1.0, 0);
    const ForceSource& force_source = LennardJonesForce();
    auto new_particles = boundary.applyBoundary(p, force_source);

    Particle p1(p);
    p1.getX() = {10.5, 5.0, 5.0};
    p1.getType() = 1;
    auto np = new_particles.value_or(std::vector<Particle>());
    EXPECT_EQ(np.size(), 1);
    EXPECT_EQ(p1, np[0]);
}

TEST_F(PeriodicTest, ParticleOnEdgeCopyTestLimitedCutoff) {
    boundary = Periodic(BoundaryLocation::LEFT, dimension, 1.0, 3);
    R3 x = {.0, .0, 5.0};
    R3 v = {.0, .0, .0};
    Particle p = Particle(x, v, 1.0, 1.0, 1.0, 0);
    const ForceSource& force_source = LennardJonesForce();
    auto new_particles = boundary.applyBoundary(p, force_source);

    Particle p1(p);
    Particle p2(p);
    p1.getX() = {10.0, .0, 5.0};
    p2.getX() = {10.0, 10.0, 5.0};
    p1.getType() = 1;
    p2.getType() = 1;
    auto np = new_particles.value_or(std::vector<Particle>());
    EXPECT_TRUE(np.size() == 2);
    EXPECT_TRUE(std::ranges::find(np.begin(), np.end(), p1) != np.end());
    EXPECT_TRUE(std::ranges::find(np.begin(), np.end(), p2) != np.end());
}

TEST_F(PeriodicTest, ParticleOnCornerCopyTestLimitedCutoff) {
    boundary = Periodic(BoundaryLocation::LEFT, dimension, 1.0, 3);
    R3 x = {.0, .0, .0};
    R3 v = {.0, .0, .0};
    Particle p = Particle(x, v, 1.0, 1.0, 1.0, 0);
    const ForceSource& force_source = LennardJonesForce();
    auto new_particles = boundary.applyBoundary(p, force_source);

    Particle p1(p);
    Particle p2(p);
    Particle p3(p);
    Particle p4(p);
    p1.getX() = {10.0, .0, .0};
    p2.getX() = {10.0, .0, 10.0};
    p3.getX() = {10.0, 10.0, .0};
    p4.getX() = {10.0, 10.0, 10.0};
    p1.getType() = 1;
    p2.getType() = 1;
    p3.getType() = 1;
    p4.getType() = 1;
    auto np = new_particles.value_or(std::vector<Particle>());
    EXPECT_TRUE(np.size() == 4);
    EXPECT_TRUE(std::ranges::find(np.begin(), np.end(), p1) != np.end());
    EXPECT_TRUE(std::ranges::find(np.begin(), np.end(), p2) != np.end());
    EXPECT_TRUE(std::ranges::find(np.begin(), np.end(), p3) != np.end());
    EXPECT_TRUE(std::ranges::find(np.begin(), np.end(), p4) != np.end());
}

TEST_F(PeriodicTest, ParticleMoveTestLimitedCutoff) {
    boundary = Periodic(BoundaryLocation::LEFT, dimension, 1.0, 3);
    R3 x = {-.5, 5.0, 5.0};
    R3 v = {.0, .0, .0};
    Particle p = Particle(x, v, 1.0, 1.0, 1.0, 0);
    const ForceSource& force_source = LennardJonesForce();
    auto new_particles = boundary.applyBoundary(p, force_source);

    Particle p_exp(p);
    p_exp.getX() = {9.5, 5.0, 5.0};
    EXPECT_EQ(p, p_exp);
}

//--------------------------------------------2D--------------------------------------------------
TEST_F(PeriodicTest, _2DtestTeleport) {
    boundary = Periodic(BoundaryLocation::LEFT, dimension, 1.0, 2);
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
    boundary = Periodic(BoundaryLocation::LEFT, dimension, 1.0, 2);
    R3 x = {.5, 5.0, 0.0};
    R3 v = {.0, .0, .0};
    Particle p = Particle(x, v, 1.0, 1.0, 1.0, 0);
    const ForceSource& force_source = LennardJonesForce();
    auto new_particles = boundary.applyBoundary(p, force_source);

    Particle p_exp(p);
    p_exp.getX() = {10.5, 5.0, 0.0};
    p_exp.getType() = 1;
    auto np = new_particles.value_or(std::vector<Particle>());
    EXPECT_TRUE(np.size() == 1);
    EXPECT_EQ(np[0], p_exp);
}

TEST_F(PeriodicTest, _2DtestMirrorOnBoundary) {
    dimension = {10.0, 20.0, 1.0};
    boundary = Periodic(BoundaryLocation::LEFT, dimension, 1.0, 2);
    R3 x = {.0, 5.0, 0.0};
    R3 v = {.0, .0, .0};
    Particle p = Particle(x, v, 1.0, 1.0, 1.0, 0);
    const ForceSource& force_source = LennardJonesForce();
    auto new_particles = boundary.applyBoundary(p, force_source);

    Particle p_exp(p);
    p_exp.getX() = {10.0, 5.0, 0.0};
    p_exp.getType() = 1;
    auto np = new_particles.value_or(std::vector<Particle>());
    EXPECT_TRUE(np.size() == 1);
    EXPECT_EQ(np[0], p_exp);
}

TEST_F(PeriodicTest, _2DParticleOnCorner0CopyTestLimitedCutoff) {
    boundary = Periodic(BoundaryLocation::LEFT, dimension, 1.0, 2);
    R3 x = {.0, .0, .0};
    R3 v = {.0, .0, .0};
    Particle p = Particle(x, v, 1.0, 1.0, 1.0, 0);
    const ForceSource& force_source = LennardJonesForce();
    auto new_particles = boundary.applyBoundary(p, force_source);
   
    Particle p1(p);
    Particle p2(p);
    p1.getX() = {10.0, .0, .0};
    p2.getX() = {10.0, 10.0, .0};
    p1.getType() = 1;
    p2.getType() = 1;
    auto np = new_particles.value_or(std::vector<Particle>());
    EXPECT_TRUE(np.size() == 2);
    EXPECT_TRUE(std::find(np.begin(), np.end(), p1) != np.end());
    EXPECT_TRUE(std::find(np.begin(), np.end(), p2) != np.end());
}

TEST_F(PeriodicTest, _2DParticleOnCorner1CopyTestLimitedCutoff) {
    boundary = Periodic(BoundaryLocation::RIGHT, dimension, 1.0, 2);
    R3 x = {9.5, .5, .0};
    R3 v = {.0, .0, .0};
    Particle p = Particle(x, v, 1.0, 1.0, 1.0, 0);
    const ForceSource& force_source = LennardJonesForce();
    auto new_particles = boundary.applyBoundary(p, force_source);
   
    Particle p1(p);
    Particle p2(p);
    p1.getX() = {-.5, .5, .0};
    p2.getX() = {-.5, 10.5, .0};
    p1.getType() = 1;
    p2.getType() = 1;
    auto np = new_particles.value_or(std::vector<Particle>());
    EXPECT_TRUE(np.size() == 2);
    EXPECT_TRUE(std::find(np.begin(), np.end(), p1) != np.end());
    EXPECT_TRUE(std::find(np.begin(), np.end(), p2) != np.end());
}

TEST_F(PeriodicTest, _2DParticleOnCorner2CopyTestLimitedCutoff) {
    boundary = Periodic(BoundaryLocation::RIGHT, dimension, 1.0, 2);
    R3 x = {9.0, 9.0, .0};
    R3 v = {.0, .0, .0};
    Particle p = Particle(x, v, 1.0, 1.0, 1.0, 0);
    const ForceSource& force_source = LennardJonesForce();
    auto new_particles = boundary.applyBoundary(p, force_source);
   
    Particle p1(p);
    Particle p2(p);
    p1.getX() = {-1.0, 9.0, .0};
    p2.getX() = {-1.0, -1.0, .0};
    p1.getType() = 1;
    p2.getType() = 1;
    auto np = new_particles.value_or(std::vector<Particle>());
    EXPECT_TRUE(np.size() == 2);
    EXPECT_TRUE(std::find(np.begin(), np.end(), p1) != np.end());
    EXPECT_TRUE(std::find(np.begin(), np.end(), p2) != np.end());
}

TEST_F(PeriodicTest, _2DParticleOnCorner3CopyTestLimitedCutoff) {
    boundary = Periodic(BoundaryLocation::LEFT, dimension, 1.0, 2);
    R3 x = {.0, 9.5, .0};
    R3 v = {.0, .0, .0};
    Particle p = Particle(x, v, 1.0, 1.0, 1.0, 0);
    const ForceSource& force_source = LennardJonesForce();
    auto new_particles = boundary.applyBoundary(p, force_source);
   
    Particle p1(p);
    Particle p2(p);
    p1.getX() = {10.0, 9.5, .0};
    p2.getX() = {10.0, -0.5, .0};
    p1.getType() = 1;
    p2.getType() = 1;
    auto np = new_particles.value_or(std::vector<Particle>());
    EXPECT_TRUE(np.size() == 2);
    EXPECT_TRUE(std::find(np.begin(), np.end(), p1) != np.end());
    EXPECT_TRUE(std::find(np.begin(), np.end(), p2) != np.end());
}
}  // namespace mol_sim
