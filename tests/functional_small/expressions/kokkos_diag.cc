
#include <gtest/gtest.h>
#include "pressio_expressions.hpp"

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

TEST(expressions_kokkos, diag1)
{
  using n_t = Kokkos::View<double **, Kokkos::HostSpace>;
  n_t A("A", 4,4);
  fillMatrix(A);

  auto d = pressio::diag(A);
  d(0) = 5.;
  EXPECT_EQ( d.extent(0), 4 );
  EXPECT_DOUBLE_EQ( d(0), 5. );
  EXPECT_DOUBLE_EQ( d(1), 6.2 );
  EXPECT_DOUBLE_EQ( d(2), 11.2 );
  EXPECT_DOUBLE_EQ( d(3), 16. );
  EXPECT_DOUBLE_EQ( A(0,0), 5. );
  EXPECT_DOUBLE_EQ( A(1,1), 6.2 );
  EXPECT_DOUBLE_EQ( A(2,2), 11.2 );
  EXPECT_DOUBLE_EQ( A(3,3), 16. );

  const auto d2 = pressio::diag(A);
  d2(2) = 55.;
  EXPECT_DOUBLE_EQ( d2(0), 5. );
  EXPECT_DOUBLE_EQ( d2(1), 6.2 );
  EXPECT_DOUBLE_EQ( d2(2), 55. );
  EXPECT_DOUBLE_EQ( d2(3), 16. );
  EXPECT_DOUBLE_EQ( A(0,0), 5. );
  EXPECT_DOUBLE_EQ( A(1,1), 6.2 );
  EXPECT_DOUBLE_EQ( A(2,2), 55. );
  EXPECT_DOUBLE_EQ( A(3,3), 16. );
}

TEST(expressions_kokkos, diag2)
{
  using n_t = Kokkos::View<double**, Kokkos::HostSpace>;
  n_t A("A", 4,4);
  fillMatrix(A);

  using T = Kokkos::View<const double**, Kokkos::HostSpace>;
  T B = A;

  auto d = pressio::diag(B);
  static_assert
    (std::is_const<
     typename std::remove_reference<decltype(d(0))>::type
     >::value, "");
}
