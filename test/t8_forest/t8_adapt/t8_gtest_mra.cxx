
#include <gtest/gtest.h>

#include <t8_forest/t8_adapt/multiscale.hpp>

class t8_adapt_mra_test : public ::testing::Test {
  void SetUp() {}

 public:
  size_t max_level = 5;
  size_t p = 3;
  double c_thresh = 0.1;
};

/// TODO
/// write tests for orhtogonal maskmatrices, mst -> imst -> mst, cancellation
/// property, etc.

TEST_F(t8_adapt_mra_test, cstr) {
  t8_mra::t8_multiscale<T8_ECLASS_TRIANGLE> mra(p, c_thresh, max_level);
  ASSERT_EQ(t8_mra::t8_multiscale<T8_ECLASS_TRIANGLE>::DIM, 2);

  ASSERT_EQ(mra.DIM, 2);
  ASSERT_DOUBLE_EQ(mra.c_thresh, c_thresh);
  ASSERT_EQ(max_level, mra.max_level);
  ASSERT_EQ(p, mra.polynomial_degree);
  ASSERT_EQ((p * (p + 1)) / 2, mra.dof);
  ASSERT_EQ(mra.mask_coeffs.size(), 4u);
  ASSERT_EQ(mra.inv_mask_coeffs.size(), 4u);
}
