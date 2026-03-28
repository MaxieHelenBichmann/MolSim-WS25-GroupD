#include "particles/container/LinkedCellContainer.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <optional>
#include <set>
#include <vector>

#include "particles/ParticleContainer.h"
#include "utils/Vector.h"

namespace mol_sim {

// Note: Basic ParticleContainer interface tests are in ParticleContainer_test.cpp
// This file contains LinkedCellContainer-specific tests only.

static_assert(ParticleContainer<LinkedCellContainer>, "LinkedCellContainer must satisfy ParticleContainer concept");

/**
 * @brief Test Fixture for LinkedCellContainer-specific tests.
 */
class LinkedCellContainerSpecificTest : public testing::Test {
   protected:
    R3 domain_size;
    LinkedCellContainer particles_empty;
    LinkedCellContainer particles_full;
    Particle p0;
    Particle p1;
    Particle p2;
    Particle p3;
    double cutoff_radius_even = 2.5;

    LinkedCellContainerSpecificTest()
        : domain_size(R3{10.0, 10.0, 10.0}),
          particles_empty(LinkedCellContainer(domain_size, 2.5)),
          particles_full(LinkedCellContainer(domain_size, 2.5)),
          p0(Particle(0)),
          p1(Particle(1)),
          p2(Particle(2)),
          p3(Particle(3)) {}

    void SetUp() override {
        particles_empty.clear();
        particles_full.clear();
        particles_full.addParticle(p0);
        particles_full.addParticle(p1);
        particles_full.addParticle(p2);
        particles_full.addParticle(p3);
    }
};

// ============================================================================
// LinkedCellContainer-specific tests (spatial grid behavior)
// ============================================================================

/**
 * @brief Tests correct behaviour of method fitsDomain.
 */
TEST_F(LinkedCellContainerSpecificTest, FitsDomain) {
    R3 inside{5.0, 5.0, 5.0};
    EXPECT_TRUE(particles_empty.fitsDomain(inside));

    R3 far_outside{1000.0, 0.0, 0.0};
    EXPECT_FALSE(particles_empty.fitsDomain(far_outside));

    R3 close_inside{9.7, 10.0, 10.0};
    EXPECT_TRUE(particles_empty.fitsDomain(close_inside));

    R3 close_outside{10.0, 10.5, 10.0};
    EXPECT_FALSE(particles_empty.fitsDomain(close_outside));

    R3 on_min{0.0, 0.0, 0.0};
    R3 on_max{10.0, 10.0, 10.0};
    EXPECT_TRUE(particles_empty.fitsDomain(on_min));
    EXPECT_TRUE(particles_empty.fitsDomain(on_max));

    R3 outside_x{-0.1, 5.0, 5.0};
    R3 outside_y{5.0, -0.1, 5.0};
    R3 outside_z{5.0, 5.0, 10.1};
    EXPECT_FALSE(particles_empty.fitsDomain(outside_x));
    EXPECT_FALSE(particles_empty.fitsDomain(outside_y));
    EXPECT_FALSE(particles_empty.fitsDomain(outside_z));
}

/**
 * @brief Tests correct behaviour of method fitsContainer, which also considers halo cells.
 */
TEST_F(LinkedCellContainerSpecificTest, FitsContainer) {
    R3 inside{5.0, 5.0, 5.0};
    EXPECT_TRUE(particles_empty.fitsContainer(inside));

    R3 far_outside{1000.0, 0.0, 0.0};
    EXPECT_FALSE(particles_empty.fitsContainer(far_outside));

    R3 halo{-0.5, 5.0, 5.0};
    EXPECT_TRUE(particles_empty.fitsContainer(halo));

    R3 close_inside{-2.5, -2.5, -2.5};
    EXPECT_TRUE(particles_empty.fitsContainer(close_inside));

    R3 close_outside{12.5, 12.5, 12.7};
    EXPECT_FALSE(particles_empty.fitsContainer(close_outside));
}

/**
 * @brief Tests correct behaviour of the method updateParticlePosition, with different scenarios.
 */
TEST_F(LinkedCellContainerSpecificTest, UpdateParticlePositionScenarios) {
    R3 v{0.0, 0.0, 0.0};

    Particle p_inside{R3{5.0, 5.0, 5.0}, v, 1.0, 1.0, 1.0};
    Particle p_move_outside{R3{9.0, 9.0, 9.0}, v, 1.0, 1.0, 1.0};
    Particle p_move_halo{R3{1.0, 1.0, 1.0}, v, 1.0, 1.0, 1.0};

    particles_empty.addParticle(p_inside);
    particles_empty.addParticle(p_move_outside);
    particles_empty.addParticle(p_move_halo);

    ASSERT_EQ(particles_empty.size(), 3U);

    // Move inside -> inside
    particles_empty.updateParticlePosition(particles_empty.begin(), R3{6.0, 6.0, 6.0});
    EXPECT_TRUE(particles_empty.begin()->getX() == (R3{6.0, 6.0, 6.0}));

    // Move inside -> outside
    particles_empty.updateParticlePosition(particles_empty.begin() + 1, R3{15.0, 15.0, 15.0});
    EXPECT_EQ(particles_empty.size(), 2U);

    // Move inside -> halo
    particles_empty.updateParticlePosition(particles_empty.begin() + 1, R3{-1.0, -1.0, -1.0});
    EXPECT_TRUE((particles_empty.begin() + 1)->getX() == (R3{-1.0, -1.0, -1.0}));
}

/**
 * @brief Tests correct behaviour of the boundary iterator for all and a specific side.
 */
TEST_F(LinkedCellContainerSpecificTest, BoundaryIteratorDetailed) {  // NOLINT
    R3 v{0.0, 0.0, 0.0};

    Particle par_inside{R3{5.0, 5.0, 5.0}, v, 1.0, 1.0, 1.0};
    Particle par_boundary_left{R3{0.1, 5.0, 5.0}, v, 1.0, 1.0, 1.0};
    Particle par_boundary_right{R3{9.0, 5.0, 5.0}, v, 1.0, 1.0, 1.0};

    particles_empty.addParticle(par_inside);
    particles_empty.addParticle(par_boundary_left);
    particles_empty.addParticle(par_boundary_right);

    ASSERT_EQ(particles_empty.size(), 3U);

    // All boundary sides
    std::vector<Particle> boundary_all;
    auto it = particles_empty.boundaryBegin();
    auto end = particles_empty.boundaryEnd();
    for (; it != end; ++it) {
        boundary_all.push_back(*it);
    }

    EXPECT_EQ(boundary_all.size(), 2U);
    EXPECT_TRUE(std::count(boundary_all.begin(), boundary_all.end(), par_boundary_left) > 0);
    EXPECT_TRUE(std::count(boundary_all.begin(), boundary_all.end(), par_boundary_right) > 0);
    EXPECT_FALSE(std::count(boundary_all.begin(), boundary_all.end(), par_inside) > 0);

    // Only RIGHT boundary
    std::vector<Particle> boundary_right;
    for (auto it = particles_empty.boundaryBegin({BoundaryLocation::RIGHT});
         it != particles_empty.boundaryEnd({BoundaryLocation::RIGHT}); ++it) {
        boundary_right.push_back(*it);
    }

    EXPECT_EQ(boundary_right.size(), 1U);
    EXPECT_FALSE(std::count(boundary_right.begin(), boundary_right.end(), par_boundary_left) > 0);
    EXPECT_TRUE(std::count(boundary_right.begin(), boundary_right.end(), par_boundary_right) > 0);
    EXPECT_FALSE(std::count(boundary_right.begin(), boundary_right.end(), par_inside) > 0);
}

/**
 * @brief Tests correct behaviour of the halo iterator for all and a specific side.
 */
TEST_F(LinkedCellContainerSpecificTest, HaloIteratorDetailed) {  // NOLINT
    R3 v{0.0, 0.0, 0.0};

    Particle par_inside{R3{5.0, 5.0, 5.0}, v, 1.0, 1.0, 1.0};
    Particle par_boundary{R3{2.0, 2.0, 2.0}, v, 1.0, 1.0, 1.0};
    Particle par_halo_left{R3{-0.1, 5.0, 5.0}, v, 1.0, 1.0, 1.0};
    Particle par_halo_right{R3{12.0, 5.0, 5.0}, v, 1.0, 1.0, 1.0};

    particles_empty.addParticle(par_inside);
    particles_empty.addParticle(par_boundary);
    particles_empty.addParticle(par_halo_left);
    particles_empty.addParticle(par_halo_right);

    ASSERT_EQ(particles_empty.size(), 4U);

    // All boundary sides
    std::vector<Particle> halo_all;
    auto it = particles_empty.haloBegin();
    auto end = particles_empty.haloEnd();
    for (; it != end; ++it) {
        halo_all.push_back(*it);
    }

    EXPECT_EQ(halo_all.size(), 2U);
    EXPECT_TRUE(std::count(halo_all.begin(), halo_all.end(), par_halo_left) > 0);
    EXPECT_TRUE(std::count(halo_all.begin(), halo_all.end(), par_halo_right) > 0);
    EXPECT_FALSE(std::count(halo_all.begin(), halo_all.end(), par_inside) > 0);
    EXPECT_FALSE(std::count(halo_all.begin(), halo_all.end(), par_boundary) > 0);

    // Only RIGHT boundary
    std::vector<Particle> halo_right;
    for (auto it = particles_empty.haloBegin({BoundaryLocation::RIGHT});
         it != particles_empty.haloEnd({BoundaryLocation::RIGHT}); ++it) {
        halo_right.push_back(*it);
    }

    EXPECT_EQ(halo_right.size(), 1U);
    EXPECT_FALSE(std::count(halo_right.begin(), halo_right.end(), par_halo_left) > 0);
    EXPECT_TRUE(std::count(halo_right.begin(), halo_right.end(), par_halo_right) > 0);
    EXPECT_FALSE(std::count(halo_right.begin(), halo_right.end(), par_inside) > 0);
    EXPECT_FALSE(std::count(halo_right.begin(), halo_right.end(), par_boundary) > 0);
}

/**
 * @brief Tests proximity iterator with Newton's third law optimization.
 */
TEST_F(LinkedCellContainerSpecificTest, ProximityIteratorN3L) {
    LinkedCellContainer particles_one({10.0, 10.0, 10.0}, 1.0);

    R3 v{0.0, 0.0, 0.0};
    particles_one.addParticle(R3{2.5, 3.1, 3.1}, v, 1.0, 1.0, 1.0);
    particles_one.addParticle(R3{5.0, 5.0, 5.0}, v, 1.0, 1.0, 1.0);
    particles_one.addParticle(R3{6.0, 6.0, 6.0}, v, 1.0, 1.0, 1.0);
    R3 center{3.1, 3.1, 3.1};

    auto it = particles_one.proximityBegin(center, particles_one.size());
    auto end = particles_one.proximityEnd(center);

    size_t count = 0;
    while (it != end) {
        EXPECT_LE((it->getX() - center).euclidNorm(), 1.0);
        ++it;
        ++count;
    }
    // only one particle is within the radius AND considered with the Newton's third law optimization
    EXPECT_EQ(count, 1U);
}

/**
 * @brief Tests that particle neighbors are correctly invalidated after eraseParticle.
 * This is critical for membrane simulations to avoid dangling pointers.
 */
TEST_F(LinkedCellContainerSpecificTest, EraseParticleInvalidatesNeighbors) {
    LinkedCellContainer particles({10.0, 10.0, 10.0}, 5.0);

    R3 v{0.0, 0.0, 0.0};
    // Create a simple 2x2 membrane-like structure with neighbors
    particles.addParticle(R3{1.0, 1.0, 0.0}, v, 1.0, 1.0, 1.0, 2);  // idx 0
    particles.addParticle(R3{2.0, 1.0, 0.0}, v, 1.0, 1.0, 1.0, 2);  // idx 1
    particles.addParticle(R3{1.0, 2.0, 0.0}, v, 1.0, 1.0, 1.0, 2);  // idx 2
    particles.addParticle(R3{2.0, 2.0, 0.0}, v, 1.0, 1.0, 1.0, 2);  // idx 3

    ASSERT_EQ(particles.size(), 4);

    // Set up ONLY reciprocal neighbor relationships for a 2x2 grid:
    //   2  3   (top row)
    //   0  1   (bottom row)

    // Particle 0 (bottom-left) ↔ Particle 1 (bottom-right): horizontal
    particles[0].getNeighbors()[1] = 1;  // 0's right = 1
    particles[1].getNeighbors()[0] = 0;  // 1's left = 0

    // Particle 0 (bottom-left) ↔ Particle 2 (top-left): vertical
    particles[0].getNeighbors()[3] = 2;  // 0's top = 2
    particles[2].getNeighbors()[2] = 0;  // 2's bottom = 0

    // Particle 1 (bottom-right) ↔ Particle 3 (top-right): vertical
    particles[1].getNeighbors()[3] = 3;  // 1's top = 3
    particles[3].getNeighbors()[2] = 1;  // 3's bottom = 1

    // Particle 2 (top-left) ↔ Particle 3 (top-right): horizontal
    particles[2].getNeighbors()[1] = 3;  // 2's right = 3
    particles[3].getNeighbors()[0] = 2;  // 3's left = 2

    // Diagonal: Particle 0 (bottom-left) ↔ Particle 3 (top-right)
    particles[0].getNeighbors()[7] = 3;  // 0's top-right diagonal = 3
    particles[3].getNeighbors()[4] = 0;  // 3's bottom-left diagonal = 0

    // Diagonal: Particle 1 (bottom-right) ↔ Particle 2 (top-left)
    particles[1].getNeighbors()[6] = 2;  // 1's top-left diagonal = 2
    particles[2].getNeighbors()[5] = 1;  // 2's bottom-right diagonal = 1

    // Verify initial setup
    EXPECT_EQ(particles[0].getNeighbors()[1], 1);
    EXPECT_EQ(particles[1].getNeighbors()[0], 0);

    // Erase particle 1 (bottom-right)
    auto it = particles.begin();
    ++it;  // Move to particle 1
    particles.eraseParticle(it);

    EXPECT_EQ(particles.size(), 3);

    // Check that all remaining neighbor indices are valid (< 3)
    // Particle 0 should have: right=nullopt (was 1, erased), top=2, top-right diagonal=1 (was 3, swapped to 1)
    EXPECT_EQ(std::nullopt, particles[0].getNeighbors()[1]);
    EXPECT_EQ(2, particles[0].getNeighbors()[3]);
    EXPECT_EQ(1, particles[0].getNeighbors()[7]);

    // Particle 1 (was Particle 3, swapped) should have: left=2, bottom=nullopt (was 1, erased), bottom-left=0
    EXPECT_EQ(2, particles[1].getNeighbors()[0]);
    EXPECT_EQ(std::nullopt, particles[1].getNeighbors()[2]);
    EXPECT_EQ(0, particles[1].getNeighbors()[4]);

    // Particle 2 should have: bottom=0, right=1 (was 3, swapped to position 1), bottom-right=nullopt (was 1, erased)
    EXPECT_EQ(1, particles[2].getNeighbors()[1]);
    EXPECT_EQ(0, particles[2].getNeighbors()[2]);
    EXPECT_EQ(std::nullopt, particles[2].getNeighbors()[5]);
}

/**
 * @brief Tests that after erasing multiple membrane particles, remaining neighbors are valid.
 */
TEST_F(LinkedCellContainerSpecificTest, EraseMultipleMembraneParticles) {
    LinkedCellContainer particles({10.0, 10.0, 10.0}, 5.0);

    R3 v{0.0, 0.0, 0.0};
    // Create a 3x3 membrane
    for (size_t j = 0; j < 3; j++) {
        for (size_t k = 0; k < 3; k++) {
            particles.addParticle(R3{static_cast<double>(k), static_cast<double>(j), 0.0}, v, 1.0, 1.0, 1.0, 2);
        }
    }

    ASSERT_EQ(particles.size(), 9);

    // Set up some neighbor relationships for particles 3, 4, 5 (middle row)
    // Particle 4 (center): left=3, right=5
    particles[4].getNeighbors()[0] = 3;
    particles[4].getNeighbors()[1] = 5;
    particles[3].getNeighbors()[1] = 4;
    particles[5].getNeighbors()[0] = 4;

    // Erase center particle (index 4)
    auto it = particles.begin();
    std::advance(it, 4);
    particles.eraseParticle(it);

    EXPECT_EQ(particles.size(), 8);

    // After erase, particle 4's neighbors should be nullopt
    EXPECT_EQ(particles[3].getNeighbors()[1], std::nullopt);  // was pointing to index 4
    EXPECT_EQ(particles[5].getNeighbors()[0], std::nullopt);  // was pointing to index 4

    // Verify all neighbor indices are valid (< size)
    for (size_t i = 0; i < particles.size(); i++) {
        for (const auto& neighbor : particles[i].getNeighbors()) {
            if (neighbor.has_value()) {
                EXPECT_LT(neighbor.value(), particles.size()) << "Particle " << i << " has invalid neighbor index";
            }
        }
    }
}

}  // namespace mol_sim
