#include "particles/container/LinkedCellContainer.h"

#include <gtest/gtest.h>

#include <algorithm>
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

    // Set up neighbor relationships manually (as MembraneGenerator would)
    // Particle 0 (bottom-left): right=1, top=2, top-right=3
    particles[0].getNeighbors()[1] = &particles[1];  // right
    particles[0].getNeighbors()[3] = &particles[2];  // top
    particles[0].getNeighbors()[7] = &particles[3];  // top-right

    // Particle 1 (bottom-right): left=0, top=3, top-left=2
    particles[1].getNeighbors()[0] = &particles[0];  // left
    particles[1].getNeighbors()[3] = &particles[3];  // top
    particles[1].getNeighbors()[6] = &particles[2];  // top-left

    // Particle 2 (top-left): right=3, bottom=0, bottom-right=1
    particles[2].getNeighbors()[1] = &particles[3];  // right
    particles[2].getNeighbors()[2] = &particles[0];  // bottom
    particles[2].getNeighbors()[5] = &particles[1];  // bottom-right

    // Particle 3 (top-right): left=2, bottom=1, bottom-left=0
    particles[3].getNeighbors()[0] = &particles[2];  // left
    particles[3].getNeighbors()[2] = &particles[1];  // bottom
    particles[3].getNeighbors()[4] = &particles[0];  // bottom-left

    // Verify initial setup
    EXPECT_EQ(particles[0].getNeighbors()[1], &particles[1]);
    EXPECT_EQ(particles[1].getNeighbors()[0], &particles[0]);

    // Save address of particle to be erased
    Particle* addr_1 = &particles[1];

    // Erase particle 1 (bottom-right)
    auto it = particles.begin();
    ++it;  // Move to particle 1
    particles.eraseParticle(it);

    EXPECT_EQ(particles.size(), 3);

    // After erase, the last particle (old index 3) is swapped to position 1
    // So particles[1] now contains what was particles[3]
    // Verify that old neighbors pointing to erased particle 1 should be updated/invalidated
    // Note: The container doesn't automatically fix neighbor pointers - this is a manual concern

    // Check that particles still in container don't have dangling pointers to erased particle
    // In a real membrane simulation, the user would need to fix these neighbor pointers
    for (size_t i = 0; i < particles.size(); i++) {
        for (auto* neighbor : particles[i].getNeighbors()) {
            // If neighbor points to old address of particle 1, that's a dangling pointer
            if (neighbor == addr_1) {
                FAIL() << "Particle " << i << " has dangling pointer to erased particle at " << addr_1;
            }
        }
    }
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
    particles[4].getNeighbors()[0] = &particles[3];
    particles[4].getNeighbors()[1] = &particles[5];
    particles[3].getNeighbors()[1] = &particles[4];
    particles[5].getNeighbors()[0] = &particles[4];

    // Erase center particle (index 4)
    auto it = particles.begin();
    std::advance(it, 4);
    Particle* erased_addr = &(*it);
    particles.eraseParticle(it);

    EXPECT_EQ(particles.size(), 8);

    // Verify no remaining particles have dangling pointers to the erased particle
    for (size_t i = 0; i < particles.size(); i++) {
        for (auto* neighbor : particles[i].getNeighbors()) {
            if (neighbor == erased_addr) {
                FAIL() << "Particle " << i << " has dangling pointer to erased center particle";
            }
        }
    }

    // Note: In a real application, neighbor pointers would need to be manually updated
    // This test just verifies the container operation doesn't introduce new dangling pointers
}

}  // namespace mol_sim
