#include "particles/container/cells/Cell.h"

#include <gtest/gtest.h>

namespace mol_sim {

class CellTest : public ::testing::Test {
   protected:
    std::array<double, 6> unit_bounds{{0.0, 1.0, 0.0, 1.0, 0.0, 1.0}};
    Cell inner_cell{CellType::INNER, unit_bounds};
    Cell boundary_cell{CellType::BOUNDARY, unit_bounds};
    Cell halo_cell{CellType::HALO, unit_bounds};

    void SetUp() override {
        inner_cell.clear();
        boundary_cell.clear();
        halo_cell.clear();
        inner_cell.addParticle(1);
        inner_cell.addParticle(2);
        inner_cell.addParticle(3);
    }
};

TEST_F(CellTest, testAddParticle) {
    size_t old_size = inner_cell.size();
    inner_cell.addParticle(99);

    EXPECT_EQ(inner_cell.size(), old_size + 1);
    EXPECT_TRUE(inner_cell.particles().contains(99));
}

TEST_F(CellTest, testAddParticleNOP) {
    size_t old_size = inner_cell.size();
    inner_cell.addParticle(2);

    EXPECT_EQ(inner_cell.size(), old_size);
}

TEST_F(CellTest, testRemoveParticle) {
    inner_cell.removeParticle(2);

    const auto& p = inner_cell.particles();
    EXPECT_EQ(p.size(), static_cast<size_t>(2));

    EXPECT_TRUE(p.contains(1));
    EXPECT_TRUE(p.contains(3));
    EXPECT_FALSE(p.contains(2));
}

TEST_F(CellTest, testRemoveParticleNOP) {
    size_t old_size = inner_cell.size();
    inner_cell.removeParticle(999);

    EXPECT_EQ(inner_cell.size(), old_size);
    const auto& p = inner_cell.particles();
    EXPECT_TRUE(p.contains(1));
    EXPECT_TRUE(p.contains(2));
    EXPECT_TRUE(p.contains(3));
}

TEST_F(CellTest, testUpdateParticleIndex) {
    inner_cell.updateParticleIndex(2, 20);

    const auto& p = inner_cell.particles();
    EXPECT_EQ(p.size(), static_cast<size_t>(3));

    EXPECT_FALSE(p.contains(2));
    EXPECT_TRUE(p.contains(20));
}

TEST_F(CellTest, testUpdateParticleIndexNOP) {
    size_t old_size = inner_cell.size();
    auto before = inner_cell.particles();

    inner_cell.updateParticleIndex(999, 1000);

    EXPECT_EQ(inner_cell.size(), old_size);
    EXPECT_EQ(inner_cell.particles(), before);
}

TEST_F(CellTest, testClear) {
    inner_cell.clear();
    EXPECT_EQ(inner_cell.size(), static_cast<size_t>(0));
    EXPECT_TRUE(inner_cell.particles().empty());
}

TEST_F(CellTest, testParticles) {
    auto& vec = inner_cell.particles();

    vec.insert(777);

    EXPECT_EQ(inner_cell.size(), static_cast<size_t>(4));
    EXPECT_TRUE(vec.contains(777));
}

TEST_F(CellTest, testConstParticles) {
    const Cell& c_ref = inner_cell;
    const auto& vec = c_ref.particles();

    EXPECT_EQ(vec.size(), static_cast<size_t>(3));
}

TEST_F(CellTest, testFits) {
    R3 inside{0.5, 0.5, 0.5};
    EXPECT_TRUE(inner_cell.fits(inside));

    R3 on_min{0.0, 0.0, 0.0};
    R3 on_max{1.0, 1.0, 1.0};
    EXPECT_TRUE(inner_cell.fits(on_min));
    EXPECT_TRUE(inner_cell.fits(on_max));

    R3 outside_x_low{-0.1, 0.5, 0.5};
    R3 outside_x_high{1.1, 0.5, 0.5};
    R3 outside_y_low{0.5, -0.1, 0.5};
    R3 outside_z_high{0.5, 0.5, 1.1};

    EXPECT_FALSE(inner_cell.fits(outside_x_low));
    EXPECT_FALSE(inner_cell.fits(outside_x_high));
    EXPECT_FALSE(inner_cell.fits(outside_y_low));
    EXPECT_FALSE(inner_cell.fits(outside_z_high));
}

TEST_F(CellTest, testGetType) {
    EXPECT_EQ(inner_cell.getType(), CellType::INNER);
    EXPECT_EQ(boundary_cell.getType(), CellType::BOUNDARY);
    EXPECT_EQ(halo_cell.getType(), CellType::HALO);
}

TEST_F(CellTest, testSize) {
    EXPECT_EQ(halo_cell.size(), static_cast<size_t>(0));

    EXPECT_EQ(inner_cell.size(), static_cast<size_t>(3));

    inner_cell.addParticle(42);
    EXPECT_EQ(inner_cell.size(), static_cast<size_t>(4));

    inner_cell.removeParticle(2);
    EXPECT_EQ(inner_cell.size(), static_cast<size_t>(3));
}

}  // namespace mol_sim