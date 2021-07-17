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
  A(4,0) = 12.; A(4,1) = 14.; A(4,2) = 5.;   A(4,3) = 6.;
}
};

TEST(expressions_kokkos, subspan0)
{
  using n_t  = Kokkos::View<double**, Kokkos::HostSpace>;
  const n_t A("A",5,4);
  auto ss = pressio::subspan(A,
					 std::make_pair(0,2),
					 std::make_pair(1,2) );
  static_assert
    (!std::is_const<
     typename std::remove_reference<decltype(ss(0,0))>::type
     >::value, "");
}

TEST(expressions_kokkos, subspan1)
{
  using T  = Kokkos::View<double **, Kokkos::HostSpace>;
  T A("A",5,4);
  fillMatrix(A);

  const auto ss = pressio::subspan(A,
					       std::make_pair(0,2),
					       std::make_pair(2,4) );
  EXPECT_EQ( ss.extent(0), 2 );
  EXPECT_EQ( ss.extent(1), 2 );

  EXPECT_DOUBLE_EQ( ss(0,0), 3.); EXPECT_DOUBLE_EQ( ss(0,1), 4.);
  EXPECT_DOUBLE_EQ( ss(1,0), 7);  EXPECT_DOUBLE_EQ( ss(1,1), 8.);

  //change A, should change expr too
  A(0,3)=43.;
  EXPECT_DOUBLE_EQ( ss(0,1), 43. );
}


TEST(expressions_kokkos, subspan2)
{
  using T  = Kokkos::View<double**, Kokkos::HostSpace>;
  T A("A",5,4);
  fillMatrix(A);

  auto ss = pressio::subspan(A,
                 std::make_pair(2,4),
                 std::make_pair(2,4) );
  EXPECT_EQ( ss.extent(0), 2 );
  EXPECT_EQ( ss.extent(1), 2 );
  ss(0,0) = -1.;
  ss(1,1) = -2.;
  EXPECT_DOUBLE_EQ( A(2,2), -1.); 
  EXPECT_DOUBLE_EQ( A(3,3), -2.);  

  const auto ss2 = pressio::subspan(A,
                 std::make_pair(2,3),
                 std::make_pair(2,4) );
  EXPECT_EQ( ss2.extent(0), 1);
  EXPECT_EQ( ss2.extent(1), 2);
  ss2(0,0) = -3.;
  ss2(0,1) = -4.;
  EXPECT_DOUBLE_EQ( A(2,2), -3.); 
  EXPECT_DOUBLE_EQ( A(2,3), -4.); 
  EXPECT_DOUBLE_EQ( A(3,3), -2.); 
}
