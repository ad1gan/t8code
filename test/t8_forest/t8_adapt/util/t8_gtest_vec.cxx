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

  for (auto i = 0u; i < bar4.size(); ++i) ASSERT_EQ(bar4(i), 0.0);
}

TEST_F(t8_adapt_vec_test, size) {
  t8_mra::util::vec foo1;
  t8_mra::util::vec foo2(n);
  t8_mra::util::vec foo3({1.0, 2.0, 3.0});
  t8_mra::util::vec foo4 = {1.0, 2.0, 3.0, 0.0};

  ASSERT_EQ(foo1.size(), 0u);
  ASSERT_EQ(foo2.size(), n);
  ASSERT_EQ(foo3.size(), 3u);
  ASSERT_EQ(foo4.size(), 4u);
}

TEST_F(t8_adapt_vec_test, resize) {
  t8_mra::util::vec foo(n);
  foo.resize(20u);

  ASSERT_EQ(foo.size(), 20u);
  for (auto i = 0u; i < foo.size(); ++i) ASSERT_EQ(foo(i), 0.0);
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
  t8_mra::util::vec foo = {2.0, -3.0, 7.0};

  ASSERT_EQ(t8_mra::util::l1norm(foo), 12.0);
  ASSERT_EQ(t8_mra::util::l2norm(foo), std::sqrt(62));
  ASSERT_EQ(t8_mra::util::linftynorm(foo), 7.0);
}

TEST_F(t8_adapt_vec_test, inner_product) {
  t8_mra::util::vec foo = {1.0, -2.0, 3.0};
  t8_mra::util::vec bar = {-10.0, 0.5, 4.0};

  ASSERT_EQ(t8_mra::util::inner(foo, bar), 1.0);
}
