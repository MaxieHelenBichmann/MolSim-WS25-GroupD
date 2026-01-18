#include "particles/boundaries/Periodic.h"

#include <gtest/gtest.h>
#include <physics/pairwiseforces/PairwiseForceSource.h>

#include <algorithm>

#include "particles/Particle.h"
#include "particles/boundaries/Boundary.h"
#include "physics/pairwiseforces/LennardJonesForce.h"
#include "utils/Vector.h"

namespace mol_sim {
/**
 * @brief Test Fixture for testing the Periodic boundary condition
 * Base data/setup:
 * dimension = (10, 10, 10)
 * cutoff = 1.0
 * 6 3D periodic boundaries for all sides (LEFT, RIGHT, FRONT, BACK, UPPER, LOWER)
 * 4 2D periodic boundaries for all sides (LEFT, RIGHT, FRONT, BACK)
 */
class PeriodicTest : public testing::Test {
   protected:
    R3 dimension{10.0, 10.0, 10.0};
    R3 zero{.0, .0, .0};
    double cutoff = 1.0;
    const PairwiseForceSource& force_source = *(new LennardJonesForce());
    Periodic left_boundary{BoundaryLocation::LEFT, dimension, cutoff, 3};
    Periodic right_boundary{BoundaryLocation::RIGHT, dimension, cutoff, 3};
    Periodic upper_boundary{BoundaryLocation::UPPER, dimension, cutoff, 3};
    Periodic lower_boundary{BoundaryLocation::LOWER, dimension, cutoff, 3};
    Periodic front_boundary{BoundaryLocation::FRONT, dimension, cutoff, 3};
    Periodic back_boundary{BoundaryLocation::BACK, dimension, cutoff, 3};
    Periodic left_boundary2_d{BoundaryLocation::LEFT, dimension, cutoff, 2};
    Periodic right_boundary2_d{BoundaryLocation::RIGHT, dimension, cutoff, 2};
    Periodic front_boundary2_d{BoundaryLocation::FRONT, dimension, cutoff, 2};
    Periodic back_boundary2_d{BoundaryLocation::BACK, dimension, cutoff, 2};

    PeriodicTest() = default;
    void TearDown() override { delete &force_source; }
};

//------------------------------------------General------------------------------------------------

/**
 * @brief Tests that a particle is copied correctly if its exactly on
 * the border between the normal part of the domain and the boundary region of the domain
 */
TEST_F(PeriodicTest, ParticleExactlyOnBorderInnerBoundaryCopy) {
    R3 x = {1.0, 5.0, 5.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    left_boundary.applyBoundary(p, force_source);
    Particle p1(p);
    p1.getX() = {11.0, 5.0, 5.0};
    p1.getType() = 1;
    EXPECT_EQ(p.getMirrorPositions().size(), 1);
    EXPECT_EQ(p1.getX(), p.getMirrorPositions()[0]);
}

/**
 * @brief Tests that a particle is NOT copied if it's not inside
 * the boundary region (or the halo region) of the domain.
 */
TEST_F(PeriodicTest, ParticleInInnerCellNoCopy) {
    R3 x = {5.0, 5.0, 5.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    left_boundary.applyBoundary(p, force_source);
    EXPECT_TRUE(p.getMirrorPositions().empty());
}

/**
 * @brief Tests that a particle is teleported correctly in 3D
 */
TEST_F(PeriodicTest, Particle3DTeleport) {
    R3 x = {-.5, 5.0, 5.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    left_boundary.applyBoundary(p, force_source);
    Particle p_exp(p);
    p_exp.getX() = {9.5, 5.0, 5.0};
    EXPECT_EQ(p, p_exp);
}

/**
 * @brief Tests that a particle is teleported correctly in 2D
 */
TEST_F(PeriodicTest, Particle2DTeleport) {
    R3 x = {-.5, 5.0, 0.0};
    R3 v = {.0, .0, .0};
    Particle p = Particle(x, v, 1.0, 1.0, 1.0, 0);
    const PairwiseForceSource& force_source = LennardJonesForce();
    left_boundary2_d.applyBoundary(p, force_source);
    Particle p_exp(p);
    p_exp.getX() = {9.5, 5.0, 0.0};
    EXPECT_EQ(p, p_exp);
}

//---------------------------------------------3D mirrors-----------------------------------------------
//-------------------------------------------3D LEFT Boundary-------------------------------------------
/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 0
 * of the LEFT boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DLeft0Mirror) {
    R3 x = {.0, .0, 5.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    left_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 16) | (1 << 25);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 1
 * of the LEFT boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DLeft1Mirror) {
    R3 x = {.0, 5.0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    left_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 16) | (1 << 17);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 2
 * of the LEFT boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DLeft2Mirror) {
    R3 x = {.0, 10.0, 5.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    left_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 7) | (1 << 16);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 3
 * of the LEFT boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DLeft3Mirror) {
    R3 x = {.0, 5.0, 10.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    left_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 15) | (1 << 16);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 4
 * of the LEFT boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DLeft4Mirror) {
    R3 x = {.0, .0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    left_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 16) | (1 << 17) | (1 << 25) | (1 << 26);
    EXPECT_TRUE(p.getMirrorPositions().size() == 4);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 5
 * of the LEFT boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DLeft5Mirror) {
    R3 x = {.0, 10.0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    left_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 16) | (1 << 17) | (1 << 7) | (1 << 8);
    EXPECT_TRUE(p.getMirrorPositions().size() == 4);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 6
 * of the LEFT boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DLeft6Mirror) {
    R3 x = {.0, 10.0, 10.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    left_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 6) | (1 << 7) | (1 << 15) | (1 << 16);
    EXPECT_TRUE(p.getMirrorPositions().size() == 4);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 7
 * of the LEFT boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DLeft7Mirror) {
    R3 x = {.0, .0, 10.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    left_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 24) | (1 << 25) | (1 << 15) | (1 << 16);
    EXPECT_TRUE(p.getMirrorPositions().size() == 4);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

//-------------------------------------------3D RIGHT Boundary-------------------------------------------
/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 0
 * of the RIGHT boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DRight0Mirror) {
    R3 x = {10.0, .0, 5.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    right_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 10) | (1 << 19);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 1
 * of the RIGHT boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DRight1Mirror) {
    R3 x = {10.0, 5.0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    right_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 10) | (1 << 11);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 2
 * of the RIGHT boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DRight2Mirror) {
    R3 x = {10.0, 10.0, 5.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    right_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 1) | (1 << 10);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 3
 * of the RIGHT boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DRight3Mirror) {
    R3 x = {10.0, 5.0, 10.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    right_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 9) | (1 << 10);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 4
 * of the RIGHT boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DRight4Mirror) {
    R3 x = {10.0, .0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    right_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 10) | (1 << 11) | (1 << 19) | (1 << 20);
    EXPECT_TRUE(p.getMirrorPositions().size() == 4);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 5
 * of the RIGHT boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DRight5Mirror) {
    R3 x = {10.0, 10.0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    right_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 10) | (1 << 11) | (1 << 1) | (1 << 2);
    EXPECT_TRUE(p.getMirrorPositions().size() == 4);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 6
 * of the RIGHT boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DRight6Mirror) {
    R3 x = {10.0, 10.0, 10.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    right_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 0) | (1 << 1) | (1 << 9) | (1 << 10);
    EXPECT_TRUE(p.getMirrorPositions().size() == 4);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 7
 * of the RIGHT boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DRight7Mirror) {
    R3 x = {10.0, .0, 10.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    right_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 18) | (1 << 19) | (1 << 9) | (1 << 10);
    EXPECT_TRUE(p.getMirrorPositions().size() == 4);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

//-------------------------------------------3D FRONT Boundary-------------------------------------------
/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 0
 * of the FRONT boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DFront0Mirror) {
    R3 x = {.0, .0, 5.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    front_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 22) | (1 << 25);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 1
 * of the FRONT boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DFront1Mirror) {
    R3 x = {5.0, .0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    front_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 22) | (1 << 23);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 2
 * of the FRONT boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DFront2Mirror) {
    R3 x = {10.0, .0, 5.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    front_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 22) | (1 << 19);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 3
 * of the FRONT boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DFront3Mirror) {
    R3 x = {5.0, .0, 10.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    front_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 21) | (1 << 22);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 4
 * of the FRONT boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DFront4Mirror) {
    R3 x = {.0, .0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    front_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 22) | (1 << 23) | (1 << 25) | (1 << 26);
    EXPECT_TRUE(p.getMirrorPositions().size() == 4);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 5
 * of the FRONT boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DFront5Mirror) {
    R3 x = {10.0, .0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    front_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 19) | (1 << 20) | (1 << 22) | (1 << 23);
    EXPECT_TRUE(p.getMirrorPositions().size() == 4);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 6
 * of the FRONT boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DFront6Mirror) {
    R3 x = {10.0, .0, 10.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    front_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 18) | (1 << 19) | (1 << 21) | (1 << 22);
    EXPECT_TRUE(p.getMirrorPositions().size() == 4);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 7
 * of the FRONT boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DFront7Mirror) {
    R3 x = {.0, .0, 10.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    front_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 21) | (1 << 22) | (1 << 24) | (1 << 25);
    EXPECT_TRUE(p.getMirrorPositions().size() == 4);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

//-------------------------------------------3D BACK Boundary-------------------------------------------
/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 0
 * of the BACK boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DBack0Mirror) {
    R3 x = {.0, 10.0, 5.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    back_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 4) | (1 << 7);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 1
 * of the BACK boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DBack1Mirror) {
    R3 x = {5.0, 10.0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    back_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 4) | (1 << 5);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 2
 * of the BACK boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DBack2Mirror) {
    R3 x = {10.0, 10.0, 5.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    back_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 4) | (1 << 1);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 3
 * of the BACK boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DBack3Mirror) {
    R3 x = {5.0, 10.0, 10.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    back_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 3) | (1 << 4);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 4
 * of the BACK boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DBack4Mirror) {
    R3 x = {.0, 10.0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    back_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 4) | (1 << 5) | (1 << 7) | (1 << 8);
    EXPECT_TRUE(p.getMirrorPositions().size() == 4);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 5
 * of the BACK boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DBack5Mirror) {
    R3 x = {10.0, 10.0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    back_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 1) | (1 << 2) | (1 << 4) | (1 << 5);
    EXPECT_TRUE(p.getMirrorPositions().size() == 4);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 6
 * of the BACK boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DBack6Mirror) {
    R3 x = {10.0, 10.0, 10.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    back_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 0) | (1 << 1) | (1 << 3) | (1 << 4);
    EXPECT_TRUE(p.getMirrorPositions().size() == 4);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 7
 * of the BACK boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DBack7Mirror) {
    R3 x = {.0, 10.0, 10.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    back_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 3) | (1 << 4) | (1 << 6) | (1 << 7);
    EXPECT_TRUE(p.getMirrorPositions().size() == 4);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}
//-------------------------------------------3D UPPER Boundary-------------------------------------------
/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 0
 * of the UPPER boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DUpper0Mirror) {
    R3 x = {.0, 5.0, 10.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    upper_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 12) | (1 << 15);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 1
 * of the UPPER boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DUpper1Mirror) {
    R3 x = {5.0, .0, 10.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    upper_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 12) | (1 << 21);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 2
 * of the UPPER boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DUpper2Mirror) {
    R3 x = {10.0, 5.0, 10.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    upper_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 12) | (1 << 9);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 3
 * of the UPPER boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DUpper3Mirror) {
    R3 x = {5.0, 10.0, 10.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    upper_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 12) | (1 << 3);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 4
 * of the UPPER boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DUpper4Mirror) {
    R3 x = {.0, .0, 10.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    upper_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 12) | (1 << 15) | (1 << 21) | (1 << 24);
    EXPECT_TRUE(p.getMirrorPositions().size() == 4);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 5
 * of the UPPER boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DUpper5Mirror) {
    R3 x = {10.0, .0, 10.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    upper_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 12) | (1 << 9) | (1 << 18) | (1 << 21);
    EXPECT_TRUE(p.getMirrorPositions().size() == 4);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 6
 * of the UPPER boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DUpper6Mirror) {
    R3 x = {10.0, 10.0, 10.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    upper_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 12) | (1 << 9) | (1 << 0) | (1 << 3);
    EXPECT_TRUE(p.getMirrorPositions().size() == 4);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 7
 * of the UPPER boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DUpper7Mirror) {
    R3 x = {.0, 10.0, 10.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    upper_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 12) | (1 << 15) | (1 << 3) | (1 << 6);
    EXPECT_TRUE(p.getMirrorPositions().size() == 4);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}
//-------------------------------------------3D LOWER Boundary-------------------------------------------
/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 0
 * of the LOWER boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DLower0Mirror) {
    R3 x = {.0, 5.0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    lower_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 14) | (1 << 17);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 1
 * of the LOWER boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DLower1Mirror) {
    R3 x = {5.0, .0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    lower_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 14) | (1 << 23);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 2
 * of the LOWER boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DLower2Mirror) {
    R3 x = {10.0, 5.0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    lower_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 14) | (1 << 11);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 3
 * of the LOWER boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DLower3Mirror) {
    R3 x = {5.0, 10.0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    lower_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 14) | (1 << 5);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 4
 * of the LOWER boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DLower4Mirror) {
    R3 x = {.0, .0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    lower_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 14) | (1 << 17) | (1 << 23) | (1 << 26);
    EXPECT_TRUE(p.getMirrorPositions().size() == 4);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 5
 * of the LOWER boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DLower5Mirror) {
    R3 x = {10.0, .0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    lower_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 14) | (1 << 11) | (1 << 20) | (1 << 23);
    EXPECT_TRUE(p.getMirrorPositions().size() == 4);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 6
 * of the LOWER boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DLower6Mirror) {
    R3 x = {10.0, 10.0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    lower_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 14) | (1 << 11) | (1 << 2) | (1 << 5);
    EXPECT_TRUE(p.getMirrorPositions().size() == 4);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

/**
 * @brief Tests that a 3D particle is mirrored correctly if it's in corner / edge with idx 7
 * of the LOWER boundary.
 */
TEST_F(PeriodicTest, ParticleOn3DLower7Mirror) {
    R3 x = {.0, 10.0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    lower_boundary.applyBoundary(p, force_source);
    uint32_t mirror_locations_expected = (1 << 14) | (1 << 17) | (1 << 5) | (1 << 8);
    EXPECT_TRUE(p.getMirrorPositions().size() == 4);
    EXPECT_EQ(mirror_locations_expected, p.getMirrorLocations());
}

//--------------------------------------------2D mirrors--------------------------------------------------
//---------------------------------------------2D LEFT-----------------------------------------------------
/**
 * @brief Tests that a 2D particle is mirrored correctly if on the LEFT boundary.
 */
TEST_F(PeriodicTest, ParticleOn2DLeftMirror) {
    R3 x = {.0, 5.0, 0.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    Particle p_exp(p);
    p_exp.getX() = {10.0, 5.0, 0.0};
    p_exp.getType() = 1;
    left_boundary2_d.applyBoundary(p, force_source);
    EXPECT_TRUE(p.getMirrorPositions().size() == 1);
    EXPECT_EQ(p.getMirrorPositions()[0], p_exp.getX());
}

/**
 * @brief Tests that a 2D particle is mirrored correctly if it's in corner with idx 0
 * of the LEFT boundary.
 */
TEST_F(PeriodicTest, ParticleOn2DLeft0Mirror) {
    R3 x = {.0, .0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    Particle p1(p);
    Particle p2(p);
    p1.getX() = {10.0, .0, .0};
    p2.getX() = {10.0, 10.0, .0};
    p1.getType() = 1;
    p2.getType() = 1;
    left_boundary2_d.applyBoundary(p, force_source);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_TRUE(std::ranges::find(p.getMirrorPositions(), p1.getX()) != p.getMirrorPositions().end());
    EXPECT_TRUE(std::ranges::find(p.getMirrorPositions(), p2.getX()) != p.getMirrorPositions().end());
}

/**
 * @brief Tests that a 2D particle is mirrored correctly if it's in corner with idx 3
 * of the LEFT boundary.
 */
TEST_F(PeriodicTest, ParticleOn2DLeft3Mirror) {
    R3 x = {.0, 10.0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    Particle p1(p);
    Particle p2(p);
    p1.getX() = {10.0, .0, .0};
    p2.getX() = {10.0, 10.0, .0};
    p1.getType() = 1;
    p2.getType() = 1;
    left_boundary2_d.applyBoundary(p, force_source);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_TRUE(std::ranges::find(p.getMirrorPositions(), p1.getX()) != p.getMirrorPositions().end());
    EXPECT_TRUE(std::ranges::find(p.getMirrorPositions(), p2.getX()) != p.getMirrorPositions().end());
}

//---------------------------------------------2D RIGHT-----------------------------------------------------
/**
 * @brief Tests that a 2D particle is mirrored correctly if on the RIGHT boundary.
 */
TEST_F(PeriodicTest, ParticleOn2DRightMirror) {
    R3 x = {10.0, 5.0, 0.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    Particle p_exp(p);
    p_exp.getX() = {.0, 5.0, 0.0};
    p_exp.getType() = 1;
    right_boundary2_d.applyBoundary(p, force_source);
    EXPECT_TRUE(p.getMirrorPositions().size() == 1);
    EXPECT_EQ(p.getMirrorPositions()[0], p_exp.getX());
}

/**
 * @brief Tests that a 2D particle is mirrored correctly if it's in corner with idx 1
 * of the RIGHT boundary.
 */
TEST_F(PeriodicTest, ParticleOn2DRight1Mirror) {
    R3 x = {10.0, .0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    Particle p1(p);
    Particle p2(p);
    p1.getX() = {.0, .0, .0};
    p2.getX() = {.0, 10.0, .0};
    p1.getType() = 1;
    p2.getType() = 1;
    right_boundary2_d.applyBoundary(p, force_source);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_TRUE(std::ranges::find(p.getMirrorPositions(), p1.getX()) != p.getMirrorPositions().end());
    EXPECT_TRUE(std::ranges::find(p.getMirrorPositions(), p2.getX()) != p.getMirrorPositions().end());
}

/**
 * @brief Tests that a 2D particle is mirrored correctly if it's in corner with idx 2
 * of the RIGHT boundary.
 */
TEST_F(PeriodicTest, ParticleOn2DRight2Mirror) {
    R3 x = {10.0, 10.0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    Particle p1(p);
    Particle p2(p);
    p1.getX() = {.0, .0, .0};
    p2.getX() = {.0, 10.0, .0};
    p1.getType() = 1;
    p2.getType() = 1;
    right_boundary2_d.applyBoundary(p, force_source);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_TRUE(std::ranges::find(p.getMirrorPositions(), p1.getX()) != p.getMirrorPositions().end());
    EXPECT_TRUE(std::ranges::find(p.getMirrorPositions(), p2.getX()) != p.getMirrorPositions().end());
}

//---------------------------------------------2D FRONT-----------------------------------------------------
/**
 * @brief Tests that a 2D particle is mirrored correctly if on the FRONT boundary.
 */
TEST_F(PeriodicTest, ParticleOn2DFrontMirror) {
    R3 x = {5.0, .0, 0.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    Particle p_exp(p);
    p_exp.getX() = {5.0, 10.0, 0.0};
    p_exp.getType() = 1;
    front_boundary2_d.applyBoundary(p, force_source);
    EXPECT_TRUE(p.getMirrorPositions().size() == 1);
    EXPECT_EQ(p.getMirrorPositions()[0], p_exp.getX());
}

/**
 * @brief Tests that a 2D particle is mirrored correctly if it's in corner with idx 0
 * of the FRONT boundary.
 */
TEST_F(PeriodicTest, ParticleOn2DFront0Mirror) {
    R3 x = {.0, .0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    Particle p1(p);
    Particle p2(p);
    p1.getX() = {.0, 10.0, .0};
    p2.getX() = {10.0, 10.0, .0};
    p1.getType() = 1;
    p2.getType() = 1;
    front_boundary2_d.applyBoundary(p, force_source);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_TRUE(std::ranges::find(p.getMirrorPositions(), p1.getX()) != p.getMirrorPositions().end());
    EXPECT_TRUE(std::ranges::find(p.getMirrorPositions(), p2.getX()) != p.getMirrorPositions().end());
}

/**
 * @brief Tests that a 2D particle is mirrored correctly if it's in corner with idx 1
 * of the FRONT boundary.
 */
TEST_F(PeriodicTest, ParticleOn2DFront1Mirror) {
    R3 x = {10.0, .0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    Particle p1(p);
    Particle p2(p);
    p1.getX() = {.0, 10.0, .0};
    p2.getX() = {10.0, 10.0, .0};
    p1.getType() = 1;
    p2.getType() = 1;
    front_boundary2_d.applyBoundary(p, force_source);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_TRUE(std::ranges::find(p.getMirrorPositions(), p1.getX()) != p.getMirrorPositions().end());
    EXPECT_TRUE(std::ranges::find(p.getMirrorPositions(), p2.getX()) != p.getMirrorPositions().end());
}

//---------------------------------------------2D BACK-----------------------------------------------------
/**
 * @brief Tests that a 2D particle is mirrored correctly if on the BACK boundary.
 */
TEST_F(PeriodicTest, ParticleOn2DBackMirror) {
    R3 x = {5.0, 20.0, 0.0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    Particle p_exp(p);
    p_exp.getX() = {5.0, .0, 0.0};
    p_exp.getType() = 1;
    back_boundary2_d.applyBoundary(p, force_source);
    EXPECT_TRUE(p.getMirrorPositions().size() == 1);
    EXPECT_EQ(p.getMirrorPositions()[0], p_exp.getX());
}

/**
 * @brief Tests that a 2D particle is mirrored correctly if it's in corner with idx 0
 * of the BACK boundary.
 */
TEST_F(PeriodicTest, ParticleOn2DBack0Mirror) {
    R3 x = {.0, 10.0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    Particle p1(p);
    Particle p2(p);
    p1.getX() = {.0, .0, .0};
    p2.getX() = {10.0, .0, .0};
    p1.getType() = 1;
    p2.getType() = 1;
    back_boundary2_d.applyBoundary(p, force_source);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_TRUE(std::ranges::find(p.getMirrorPositions(), p1.getX()) != p.getMirrorPositions().end());
    EXPECT_TRUE(std::ranges::find(p.getMirrorPositions(), p2.getX()) != p.getMirrorPositions().end());
}

/**
 * @brief Tests that a 2D particle is mirrored correctly if it's in corner with idx 2
 * of the BACK boundary.
 */
TEST_F(PeriodicTest, ParticleOn2DBack2Mirror) {
    R3 x = {10.0, 10.0, .0};
    Particle p = Particle(x, zero, 1.0, 1.0, 1.0, 0);
    Particle p1(p);
    Particle p2(p);
    p1.getX() = {.0, .0, .0};
    p2.getX() = {10.0, .0, .0};
    p1.getType() = 1;
    p2.getType() = 1;
    back_boundary2_d.applyBoundary(p, force_source);
    EXPECT_TRUE(p.getMirrorPositions().size() == 2);
    EXPECT_TRUE(std::ranges::find(p.getMirrorPositions(), p1.getX()) != p.getMirrorPositions().end());
    EXPECT_TRUE(std::ranges::find(p.getMirrorPositions(), p2.getX()) != p.getMirrorPositions().end());
}
}  // namespace mol_sim
