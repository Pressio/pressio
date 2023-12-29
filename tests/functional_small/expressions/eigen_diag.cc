
#include <gtest/gtest.h>
#include "test_helpers.hpp"

#include "pressio/expressions.hpp"

namespace
{

template <typename T>
void test1(T & A)
{
  {
    const auto diagvals = pressio::diagonal(A);
    EXPECT_EQ( diagvals.extent(0), 4 );
    EXPECT_DOUBLE_EQ( diagvals(0), 1.2 );
    EXPECT_DOUBLE_EQ( diagvals(1), 6.2 );
    EXPECT_DOUBLE_EQ( diagvals(2), 11.2 );
  }
}

template <typename T>
void test2(T & A)
{
  {
    // change some entries
    auto diagvals = pressio::diagonal(A);
    EXPECT_EQ( diagvals.extent(0), 4 );
    // before changing it
    EXPECT_DOUBLE_EQ( diagvals(0), 1.2 );
    EXPECT_DOUBLE_EQ( diagvals(1), 6.2 );
    EXPECT_DOUBLE_EQ( diagvals(2), 11.2 );
    // modify
    diagvals(0) = 44.;
    diagvals(1) = 6.;
    // after
    EXPECT_DOUBLE_EQ( diagvals(0), 44. );
    EXPECT_DOUBLE_EQ( diagvals(1), 6. );
  }

  {
    // get the native expression
    const auto diagvals = pressio::diagonal(A);
    auto &natEx = diagvals.native();
    EXPECT_DOUBLE_EQ( natEx(0), 44. );
    EXPECT_DOUBLE_EQ( natEx(1), 6. );
  }
}

template <typename T>
void testConst(const T & A){
  const  auto diagvals = pressio::diagonal(A);
  EXPECT_EQ( diagvals.extent(0), 4 );
  EXPECT_DOUBLE_EQ( diagvals(0), 44. );
  EXPECT_DOUBLE_EQ( diagvals(1), 6. );
  EXPECT_DOUBLE_EQ( diagvals(2), 11.2 );

  auto & natEx = diagvals.native();
  EXPECT_DOUBLE_EQ( natEx(0), 44. );
  EXPECT_DOUBLE_EQ( natEx(1), 6. );
}

TEST(expressions_eigen, diag)
{
  // col-major matrix (which is default in Eigen)
  using eigmat_t = Eigen::MatrixXd;

  eigmat_t A(4,4);
  fill_matrix(A);
  test1(A);
  test2(A);
  testConst(A);
}

TEST(expressions_eigen, diagRowMajor)
{
  using eigmat_t = Eigen::Matrix<double, -1, -1, Eigen::RowMajor>;
  eigmat_t A(4,4);
  fill_matrix(A);
  test1(A);
  test2(A);
  testConst(A);
}

TEST(expressions_eigen, diag_traits)
{
  {
    using T = Eigen::MatrixXd;
    T o(10, 10);
    using expr_t = decltype(pressio::diagonal(o));
    static_assert(pressio::Traits<expr_t>::rank == 1);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::scalar_type, double>);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::reference_type, double &>);
  }

  {
    using T = Eigen::Matrix<int,-1,-1>;
    T o(10, 10);
    using expr_t = decltype(pressio::diagonal(o));
    static_assert(pressio::Traits<expr_t>::rank == 1);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::scalar_type, int>);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::reference_type, int &>);
  }

  {
    using T = Eigen::Matrix<int,-1,-1>;
    const T o(10,10);
    using expr_t = decltype(pressio::diagonal(o));
    static_assert(pressio::Traits<expr_t>::rank == 1);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::scalar_type, int>);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::reference_type, int const &>);
  }
}

} // namespace
