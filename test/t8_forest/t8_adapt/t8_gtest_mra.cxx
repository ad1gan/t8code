
#include <gtest/gtest.h>

#include <t8_forest/t8_adapt/multiscale.hpp>

class t8_mra_test : public ::testing::Test {
  void SetUp() {}

 public:
  size_t max_level = 5;
  size_t p = 3;
};

TEST_F(t8_mra_test, cstr) {
  t8_mra::t8_multiscale<2> mra(p, max_level);

  ASSERT_EQ(max_level, mra.max_level);
  ASSERT_EQ(p, mra.polynomial_degree);
  ASSERT_EQ((p * (p + 1)) / 2, mra.dof);
  ASSERT_EQ(mra.mask_coeffs.size(), 4u);
  ASSERT_EQ(mra.inv_mask_coeffs.size(), 4u);
}
