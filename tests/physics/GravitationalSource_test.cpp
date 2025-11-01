#include <gtest/gtest.h>

#include "particles/Particle.h"
#include "physics/GravitationalForce.h"

namespace mol_sim {

class GravitationalForceTest : public testing::Test {
   protected:
    Particle p1;
    Particle p2;

    GravitationalForceTest() : p1({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0), p2({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 2.0) {}
};

TEST_F(GravitationalForceTest, TestForceCalcSimple) {
    R3 result = GravitationalForce().calculateForce(p1, p2);
    R3 expected = {2.0, 0.0, 0.0};
    EXPECT_EQ(result, expected);
}
TEST_F(GravitationalForceTest, TestForceCalcUnaligned) {
    p1 = {{1.0, 0.5, 2.0}, {0.0, 0.0, 0.0}, 1.0};
    p2 = {{1.0, 2.0, 4.0}, {0.0, 0.0, 0.0}, 2.0};
    R3 result = GravitationalForce().calculateForce(p1, p2);
    R3 expected = {0.0, 0.192, 0.256};
    EXPECT_EQ(result, expected);
}
TEST_F(GravitationalForceTest, TestForceSymmetry) {
    R3 result1 = GravitationalForce().calculateForce(p1, p2);
    R3 result2 = GravitationalForce().calculateForce(p2, p1);
    EXPECT_EQ(result1, result2 * -1);
}

}  // namespace mol_sim
