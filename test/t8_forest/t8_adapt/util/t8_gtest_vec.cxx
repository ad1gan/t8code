#include <gtest/gtest.h>

#include <t8_forest/t8_adapt/util/vec.hpp>

class t8_adapt_vec_test : public ::testing::Test {
  void SetUp() {}

 public:
  size_t n = 50;
};

TEST_F(t8_adapt_vec_test, cstr) {
  /// cstr
  t8_mra::util::vec foo1;
  t8_mra::util::vec foo2(n);

  /// copy cstr/assignment
  t8_mra::util::vec bar1(foo2);
  t8_mra::util::vec bar2 = foo2;

  /// move cstr/assignment
  t8_mra::util::vec bar3(std::move(bar1));
  t8_mra::util::vec bar4 = std::move(foo2);
}

TEST_F(t8_adapt_vec_test, size) {
  t8_mra::util::vec foo1;
  t8_mra::util::vec foo2(n);

  ASSERT_EQ(foo1.size(), 0u);
  ASSERT_EQ(foo2.size(), n);
}

TEST_F(t8_adapt_vec_test, access_operator) {
  t8_mra::util::vec foo(n);
  for (auto i = 0u; i < foo.size(); ++i) foo(i) = i;
  for (auto i = 0u; i < foo.size(); ++i) ASSERT_EQ(foo(i), i);
}

TEST_F(t8_adapt_vec_test, math_operator) {
  t8_mra::util::vec foo(n);
  t8_mra::util::vec bar(n);

  for (auto i = 0u; i < foo.size(); ++i) foo(i) = i;
  bar = foo;

  foo += bar;
  for (auto i = 0u; i < foo.size(); ++i) ASSERT_EQ(foo(i), 2 * i);

  foo *= 0.5;
  for (auto i = 0u; i < foo.size(); ++i) ASSERT_EQ(foo(i), i);

  foo -= bar;
  for (auto i = 0u; i < foo.size(); ++i) ASSERT_EQ(foo(i), 0.0);
}

TEST_F(t8_adapt_vec_test, norms) {
  t8_mra::util::vec foo(3);
  foo(0) = 2.0;
  foo(1) = -3.0;
  foo(2) = 7.0;

  ASSERT_EQ(t8_mra::util::l1norm(foo), 12.0);
  ASSERT_EQ(t8_mra::util::l2norm(foo), std::sqrt(62));
  ASSERT_EQ(t8_mra::util::linftynorm(foo), 7.0);
}

TEST_F(t8_adapt_vec_test, inner_product) {
  t8_mra::util::vec foo(3);
  foo(0) = 1.0;
  foo(1) = -2.0;
  foo(2) = 3.0;

  t8_mra::util::vec bar(3);
  bar(0) = -10.0;
  bar(1) = 0.5;
  bar(2) = 4.0;

  ASSERT_EQ(t8_mra::util::inner(foo, bar), 1.0);
}
