#include <gtest/gtest.h>

#include <t8_forest/t8_adapt/util/mat.hpp>
#include <t8_forest/t8_adapt/util/vec.hpp>

class t8_adapt_mat_test : public ::testing::Test {
  void SetUp() {}

 public:
  size_t m = 20;
  size_t n = 10;
  double eps = 1e-14;
};

TEST_F(t8_adapt_mat_test, cstr) {
  /// cstr
  t8_mra::util::mat foo1;
  t8_mra::util::mat foo2(m, n);

  /// copy cstr/assignment
  t8_mra::util::mat bar1(foo2);
  t8_mra::util::mat bar2 = foo2;

  /// move cstr/assignment
  t8_mra::util::mat bar3(std::move(bar1));
  t8_mra::util::mat bar4 = std::move(foo2);

  for (auto j = 0u; j < bar4.cols(); ++j)
    for (auto i = 0u; i < bar4.rows(); ++i) bar4(i, j) = 0.0;
}

TEST_F(t8_adapt_mat_test, row_cols) {
  t8_mra::util::mat foo1;
  t8_mra::util::mat foo2(m, n);

  ASSERT_EQ(foo1.rows(), 0u);
  ASSERT_EQ(foo1.cols(), 0u);
  ASSERT_EQ(foo2.rows(), m);
  ASSERT_EQ(foo2.cols(), n);
}

TEST_F(t8_adapt_mat_test, resize) {
  t8_mra::util::mat foo(m, n);
  foo.resize(10u, 20u);

  ASSERT_EQ(foo.rows(), 10u);
  ASSERT_EQ(foo.cols(), 20u);

  for (auto j = 0u; j < foo.cols(); ++j)
    for (auto i = 0u; i < foo.rows(); ++i) foo(i, j) = 0.0;
}

TEST_F(t8_adapt_mat_test, access_operator) {
  t8_mra::util::mat foo(m, n);
  for (auto j = 0u; j < foo.cols(); ++j)
    for (auto i = 0u; i < foo.rows(); ++i) foo(i, j) = i + j;

  for (auto j = 0u; j < foo.cols(); ++j)
    for (auto i = 0u; i < foo.rows(); ++i) ASSERT_EQ(foo(i, j), i + j);
}

TEST_F(t8_adapt_mat_test, lr_decomposition) {
  t8_mra::util::mat foo(3u, 3u);
  std::vector<size_t> pivot;

  /// Example from https://en.wikipedia.org/wiki/LU_decomposition#Example_2
  foo(0, 0) = 0.0;
  foo(1, 0) = 4.0;
  foo(2, 0) = 2.0;

  foo(0, 1) = 5.0;
  foo(1, 1) = 2.0;
  foo(2, 1) = 7.0;

  foo(0, 2) = 22.0 / 3.0;
  foo(1, 2) = 1.0;
  foo(2, 2) = 9.0;

  t8_mra::util::lu_factors(foo, pivot);

  ASSERT_EQ(pivot[0], 1u);
  ASSERT_EQ(pivot[1], 2u);
  ASSERT_EQ(pivot[2], 0u);

  /// L Matrix
  ASSERT_NEAR(foo(1, 0), 0.5, eps * 0.5);
  ASSERT_NEAR(foo(2, 0), 0.0, 0.0);
  ASSERT_NEAR(foo(2, 1), 5.0 / 6.0, eps * 5.0 / 6.0);

  /// U Matrix
  ASSERT_NEAR(foo(0, 0), 4.0, eps * 4.0);
  ASSERT_NEAR(foo(0, 1), 2.0, eps * 2.0);
  ASSERT_NEAR(foo(0, 2), 1.0, eps * 1.0);
  ASSERT_NEAR(foo(1, 1), 6.0, eps * 6.0);
  ASSERT_NEAR(foo(1, 2), 8.5, eps * 8.5);
  ASSERT_NEAR(foo(2, 2), 0.25, eps * 0.25);

  t8_mra::util::vec x(3u);
  x(0) = 2.0;
  x(1) = 0.0;
  x(2) = 1.0;

  t8_mra::util::lu_solve(foo, pivot, x);

  ASSERT_NEAR(x(0), 37.0 / 18.0, eps * 37.0 / 18.0);
  ASSERT_NEAR(x(1), -58.0 / 9.0, eps * 58.0 / 9.0);
  ASSERT_NEAR(x(2), 14.0 / 3.0, eps * 14.0 / 3.0);
}
