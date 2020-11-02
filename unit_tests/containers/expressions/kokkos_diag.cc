
#include <gtest/gtest.h>
#include "pressio_containers.hpp"

namespace
{
template <typename T>
void fillMatrix(T & A)
{
  A(0,0) = 1.2; A(0,1) = 2.;  A(0,2) = 3.;   A(0,3) = 4.;
  A(1,0) = 5.;  A(1,1) = 6.2; A(1,2) = 7.;   A(1,3) = 8.;
  A(2,0) = 9.;  A(2,1) = 10.; A(2,2) = 11.2; A(2,3) = 12.;
  A(3,0) = 13.; A(3,1) = 14.; A(3,2) = 15.;  A(3,3) = 16.;
}
};

TEST(containers_expressions_kokkos, diag1)
{
  using n_t = Kokkos::View<double **, Kokkos::HostSpace>;
  using my_t = pressio::containers::DenseMatrix<n_t>;
  my_t A(4,4);
  fillMatrix(A);

  const auto d = pressio::containers::diag(A);
  static_assert
    (std::is_const<
     typename std::remove_reference<decltype(d(0))>::type
     >::value, "");

  EXPECT_EQ( d.extent(0), 4 );
  EXPECT_DOUBLE_EQ( d(0), 1.2 );
  EXPECT_DOUBLE_EQ( d(1), 6.2 );
  EXPECT_DOUBLE_EQ( d(2), 11.2 );
  EXPECT_DOUBLE_EQ( d(3), 16. );

  EXPECT_DOUBLE_EQ( A(0,0), 1.2 );
  EXPECT_DOUBLE_EQ( A(1,1), 6.2 );
  EXPECT_DOUBLE_EQ( A(2,2), 11.2 );
  EXPECT_DOUBLE_EQ( A(3,3), 16. );

  //change a, should change span too
  A(1,1)=43.;
  EXPECT_DOUBLE_EQ( d(1), 43. );
}


TEST(containers_expressions_kokkos, diag2)
{
  using n_t = Kokkos::View<double **, Kokkos::HostSpace>;
  using my_t = pressio::containers::DenseMatrix<n_t>;
  my_t A(4,4);
  fillMatrix(A);
  EXPECT_DOUBLE_EQ( A(0,0), 1.2 );
  EXPECT_DOUBLE_EQ( A(1,1), 6.2 );
  EXPECT_DOUBLE_EQ( A(2,2), 11.2 );
  EXPECT_DOUBLE_EQ( A(3,3), 16. );

  auto d = pressio::containers::diag(A);
  static_assert
    (!std::is_const<
     typename std::remove_reference<decltype(d(0))>::type
     >::value, "");

  d(0) = 44.;
  d(2) = 22.;

  // both A and d should be changed
  EXPECT_DOUBLE_EQ( d(0), 44. );
  EXPECT_DOUBLE_EQ( d(1), 6.2 );
  EXPECT_DOUBLE_EQ( d(2), 22. );
  EXPECT_DOUBLE_EQ( d(3), 16. );
  EXPECT_DOUBLE_EQ( A(0,0), 44. );
  EXPECT_DOUBLE_EQ( A(1,1), 6.2 );
  EXPECT_DOUBLE_EQ( A(2,2), 22. );
  EXPECT_DOUBLE_EQ( A(3,3), 16. );
}

TEST(containers_expressions_kokkos, diag3)
{
  using n_t = Kokkos::View<double **, Kokkos::HostSpace>;
  using my_t = pressio::containers::DenseMatrix<n_t>;
  n_t A0("A0", 4,4);
  fillMatrix(A0);

  const my_t A(A0);
  EXPECT_DOUBLE_EQ( A(0,0), 1.2 );
  EXPECT_DOUBLE_EQ( A(1,1), 6.2 );
  EXPECT_DOUBLE_EQ( A(2,2), 11.2 );
  EXPECT_DOUBLE_EQ( A(3,3), 16. );

  auto d = pressio::containers::diag(A);
  // since the matrix A is const, diag should be read-only
  static_assert
    (std::is_const<
     typename std::remove_reference<decltype(d(0))>::type
     >::value, "");

  EXPECT_DOUBLE_EQ( d(0), 1.2 );
  EXPECT_DOUBLE_EQ( d(1), 6.2 );
  EXPECT_DOUBLE_EQ( d(2), 11.2 );
  EXPECT_DOUBLE_EQ( d(3), 16. );
}
