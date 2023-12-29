
#include <gtest/gtest.h>
#include "pressio/expressions.hpp"
#include "test_helpers.hpp"

namespace
{
template<class T>
auto create_matrix(){
  Eigen::MatrixXd A(6,4);
  helpers::fill_matrix_seq(A);
  return A;
}

template<class T>
void execute()
{
  {
    auto A = create_matrix<T>();
    auto c = pressio::column(A, 0);
    ASSERT_TRUE(c.extent(0) == 6);
    ASSERT_TRUE(c.extent(1) == 1);
  }

  {
    auto A = create_matrix<T>();
    for (int j=0; j<4; ++j){
      auto c = pressio::column(A, j);
      ASSERT_TRUE(c.extent(0) == 6);
      ASSERT_TRUE(c.extent(1) == 1);

      ASSERT_TRUE(c(0) == A(0,j));
      ASSERT_TRUE(c(1) == A(1,j));
      ASSERT_TRUE(c(2) == A(2,j));
      ASSERT_TRUE(c(3) == A(3,j));
      ASSERT_TRUE(c(4) == A(4,j));
      ASSERT_TRUE(c(5) == A(5,j));
      ASSERT_TRUE(c[0] == A(0,j));
      ASSERT_TRUE(c[1] == A(1,j));
      ASSERT_TRUE(c[2] == A(2,j));
      ASSERT_TRUE(c[3] == A(3,j));
      ASSERT_TRUE(c[4] == A(4,j));
      ASSERT_TRUE(c[5] == A(5,j));
    }
  }

  {
    auto A = create_matrix<T>();
    auto c = pressio::column(A, 1);
    ASSERT_TRUE(c.extent(0) == 6);
    ASSERT_TRUE(c.extent(1) == 1);

    ASSERT_TRUE(c(2) == T(10));
    ASSERT_TRUE(c(3) == T(14));
    ASSERT_TRUE(A(2,1) == T(10));
    ASSERT_TRUE(A(3,1) == T(14));

    c(2) = 181;
    c(3) = 201;
    ASSERT_TRUE(A(2,1) == T(181));
    ASSERT_TRUE(A(3,1) == T(201));
  }

  {
    const auto A = create_matrix<T>();
    auto c = pressio::column(A, 1);
    ASSERT_TRUE(c.extent(0) == 6);
    ASSERT_TRUE(c.extent(1) == 1);
    // the following should not compile because of const
    // c(2) = 181;
    // c(3) = 201;
  }
}
} //end namespace

TEST(expressions_eigen, column_double){
  execute<double>();
}

TEST(expressions_eigen, column_int){
  execute<int>();
}

TEST(expressions_eigen, column_traits)
{
  {
    using T = Eigen::MatrixXd;
    T o(10, 10);
    using expr_t = decltype(pressio::column(o, 0));
    static_assert(pressio::Traits<expr_t>::rank == 1);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::scalar_type, double>);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::reference_type, double &>);
  }

  {
    using T = Eigen::Matrix<int,-1,-1>;
    T o(10, 10);
    using expr_t = decltype(pressio::column(o, 0));
    static_assert(pressio::Traits<expr_t>::rank == 1);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::scalar_type, int>);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::reference_type, int &>);
  }


  {
    using T = Eigen::Matrix<int,-1,-1>;
    const T o(10,10);
    using expr_t = decltype(pressio::column(o, 0));
    static_assert(pressio::Traits<expr_t>::rank == 1);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::scalar_type, int>);
    static_assert(std::is_same_v<pressio::Traits<expr_t>::reference_type, int const &>);
  }
}
