#include "utils/Vector.h"

#include <gtest/gtest.h>

namespace mol_sim {
/**
 * @brief Test Fixture for testing the GravitationalForce
 * Base data/setup:
 * p1:
 *  position: (0,0,0)
 *  velocity: (0,0,0)
 *  mass: 1
 * p2:
 *  position: (1,0,0)
 *  velocity: (0,0,0)
 *  mass: 2
 */
class VectorTest : public testing::Test {
   protected:
    Vector<double, 2> R2_1;
    Vector<double, 2> R2_2;
    Vector<double, 3> R3_1;
    Vector<double, 3> R3_2;
    Vector<double, 10> R10_1;
    Vector<double, 10> R10_2;
    Vector<int, 2> Z2_1;
    Vector<int, 2> Z2_2;
    Vector<int, 3> Z3_1;
    Vector<int, 3> Z3_2;
    Vector<int, 10> Z10_1;
    Vector<int, 10> Z10_2;

    VectorTest()
        : R2_1(Vector<double, 2>{0., 0.}),
          R2_2(Vector<double, 2>{3.0, 4.0}),
          R3_1(Vector<double, 3>{13.5, 69., .01}),
          R3_2(Vector<double, 3>{6.0, -2., 50.}),
          R10_1(Vector<double, 10>()),
          R10_2(Vector<double, 10>{1.125, 0.66, 12.3, -45., 1337., 42.8, -8.6, NAN, 66.6, 0.22}),
          Z2_1(Vector<int, 2>()),
          Z2_2(Vector<int, 2>{4, 3}),
          Z3_1(Vector<int, 3>{22, 3, 17}),
          Z3_2(Vector<int, 3>{700, 44, 222}),
          Z10_1(Vector<int, 10>()),
          Z10_2(Vector<int, 10>{20502, 1301, 5, 0, 24, 66, 69, 74, 9, 10}) {}
};

/**
 * @brief Tests correct access with the subscript operator.
 */
TEST_F(VectorTest, TestAccessSubscript) {}

/**
 * @brief Tests correct access to components using structured bindings.
 *
 */
TEST_F(VectorTest, TestAccessStructuresBindings) {}

/**
 * @brief Tests equality and ordering of different Vectors.
 */
TEST_F(VectorTest, TestComparison) {}

/**
 * @brief Tests correctness of point-wise addition of two Vectors.
 */
TEST_F(VectorTest, TestPointwiseAdditon) {}

/**
 * @brief Tests correctness of point-wise subtraction of two Vectors.
 */
TEST_F(VectorTest, TestPointwiseSubtraction) {}

/**
 * @brief Tests correctness of point-wise multiplication of two Vectors.
 */
TEST_F(VectorTest, TestPointwiseMultiplication) {}

/**
 * @brief Tests correctness of point-wise division of two Vectors.
 */
TEST_F(VectorTest, TestPointwiseDivision) {}

/**
 * @brief Tests correctness of scalar multiplication of two Vectors.
 */
TEST_F(VectorTest, TestScalarMultiplication) {}

/**
 * @brief Tests correctness of the scalar product of a Vector and a scalar in both orders.
 */
TEST_F(VectorTest, TestScalarProduct) {}

/**
 * @brief Tests correctness of the euclidean norm of a Vector.
 */
TEST_F(VectorTest, TestEuclideanNorm) {}

/**
 * @brief Tests correctness of the string representation of a Vector
 */
TEST_F(VectorTest, TestStringRepresentation) {}

}  // namespace mol_sim
