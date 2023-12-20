
#include <gtest/gtest.h>
#include "test_helpers.hpp"

#include "pressio/expressions.hpp"

namespace
{

template<class T>
auto create_vector(){
  Eigen::Matrix<T, -1, 1> a(6);
  a(0) = 1;
  a(1) = 5;
  a(2) = 9;
  a(3) = 13;
  a(4) = 17;
  a(5) = 21;
  return a;
}

template <typename T>
class EigenSpanTest : public testing::Test {};

using TestingTypes = ::testing::Types<double, int, float>;

} //end namespace

TYPED_TEST_SUITE(EigenSpanTest, TestingTypes);

TYPED_TEST(EigenSpanTest, baseline) {
  {
    auto v = create_vector<TypeParam>();
    auto s = pressio::span(v, 1, 0);
    ASSERT_TRUE(s.extent(0) == 0);
  }

  {
    auto v = create_vector<TypeParam>();
    auto s = pressio::span(v, 1, 3);
    ASSERT_TRUE(s.extent(0) == 3);
    ASSERT_TRUE(s(0) == TypeParam(5));
    ASSERT_TRUE(s(1) == TypeParam(9));
    ASSERT_TRUE(s(2) == TypeParam(13));
    ASSERT_TRUE(s[0] == TypeParam(5));
    ASSERT_TRUE(s[1] == TypeParam(9));
    ASSERT_TRUE(s[2] == TypeParam(13));

    ASSERT_TRUE(v(1) == TypeParam(5));
    ASSERT_TRUE(v(2) == TypeParam(9));
    ASSERT_TRUE(v(3) == TypeParam(13));
  }

  {
    auto v = create_vector<TypeParam>();
    auto s = pressio::span(v, 1, 3);
    s(0) = 66;
    s(2) = 88;
    ASSERT_TRUE(s.extent(0) == 3);
    ASSERT_TRUE(s(0) == TypeParam(66));
    ASSERT_TRUE(s(1) == TypeParam(9));
    ASSERT_TRUE(s(2) == TypeParam(88));
    ASSERT_TRUE(s[0] == TypeParam(66));
    ASSERT_TRUE(s[1] == TypeParam(9));
    ASSERT_TRUE(s[2] == TypeParam(88));

    ASSERT_TRUE(v(1) == TypeParam(66));
    ASSERT_TRUE(v(2) == TypeParam(9));
    ASSERT_TRUE(v(3) == TypeParam(88));
  }

  {
    const auto v = create_vector<TypeParam>();
    auto s = pressio::span(v, 1, 3);
    // s(0) = 66; should not compile
    ASSERT_TRUE(s.extent(0) == 3);
    ASSERT_TRUE(s(0) == TypeParam(5));
    ASSERT_TRUE(s(1) == TypeParam(9));
    ASSERT_TRUE(s(2) == TypeParam(13));

    ASSERT_TRUE(s[0] == TypeParam(5));
    ASSERT_TRUE(s[1] == TypeParam(9));
    ASSERT_TRUE(s[2] == TypeParam(13));
  }

  {
    const auto v = create_vector<TypeParam>();
    auto s = pressio::span(v, 1, 3);
    auto s2 = s;
    ASSERT_TRUE(s.data() == s2.data());
    ASSERT_TRUE(s.extent(0) == s2.extent(0));
  }

  {
    auto v = create_vector<TypeParam>();
    auto s = pressio::span(v, 1, 3);
    auto & natEx = s.native();
    EXPECT_EQ( natEx.size(), 3 );
    EXPECT_TRUE( natEx(0) == TypeParam(5) );
    EXPECT_TRUE( natEx(1) == TypeParam(9) );
    EXPECT_TRUE( natEx(2) == TypeParam(13) );
  }
}

TEST(expressions_eigen, span_traits)
{
  {
    Eigen::VectorXd o(10);
    check_span_traits<double>(o);
  }

  {
    Eigen::Matrix<int,-1,1> o(10);
    check_span_traits<int>(o);
  }

  {
    const Eigen::Matrix<int,-1,1> o(10);
    using expr_t = decltype(pressio::span(o, 0, 1));
    static_assert(pressio::Traits<expr_t>::rank == 1);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::scalar_type, int>);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::reference_type, int const &>);
  }
}
