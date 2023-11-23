
#include <gtest/gtest.h>
#include "pressio/expressions.hpp"

namespace
{
template<class T>
auto create_matrix(){
  using eigmat_t = Eigen::MatrixXd;
  eigmat_t A(6,4);
  A(0,0) = 1.;  A(0,1) = 2.;  A(0,2) = 3.;  A(0,3) = 4.;
  A(1,0) = 5.;  A(1,1) = 6.;  A(1,2) = 7.;  A(1,3) = 8.;
  A(2,0) = 9.;  A(2,1) = 10.; A(2,2) = 11.; A(2,3) = 12.;
  A(3,0) = 13.; A(3,1) = 14.; A(3,2) = 15.; A(3,3) = 16.;
  A(4,0) = 17.; A(4,1) = 18.; A(4,2) = 19.; A(4,3) = 20.;
  A(5,0) = 21.; A(5,1) = 22.; A(5,2) = 23.; A(5,3) = 24.;
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
