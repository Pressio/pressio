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
  A(4,0) = 12.; A(4,1) = 14.; A(4,2) = 5.;   A(4,3) = 6.;
}
};

TEST(containers_expressions_kokkos, subspan0)
{
  using n_t  = Kokkos::View<double **, Kokkos::HostSpace>;
  using my_t = pressio::containers::DenseMatrix<n_t>;
  n_t A("A",5,4);
  const my_t Aw(A);
  auto ss = pressio::containers::subspan(Aw,
					 std::make_pair(0,2),
					 std::make_pair(1,2) );
  static_assert
    (std::is_const<
     typename std::remove_reference<decltype(ss(0,0))>::type
     >::value, "");
}

TEST(containers_expressions_kokkos, span1)
{
  using n_t  = Kokkos::View<double **, Kokkos::HostSpace>;
  using my_t = pressio::containers::DenseMatrix<n_t>;
  my_t A(5,4);
  fillMatrix(A);

  const auto ss = pressio::containers::subspan(A,
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

TEST(containers_expressions_kokkos, span2)
{
  using n_t  = Kokkos::View<double **, Kokkos::HostSpace>;
  using my_t = pressio::containers::DenseMatrix<n_t>;
  my_t A(5,4);
  fillMatrix(A);

  auto ss = pressio::containers::subspan(A,
					 std::make_pair(0,2),
					 std::make_pair(2,4) );
  static_assert
    (!std::is_const<
     typename std::remove_reference<decltype(ss(0,0))>::type
     >::value, "");

  EXPECT_DOUBLE_EQ( ss(0,0), 3.); EXPECT_DOUBLE_EQ( ss(0,1), 4.);
  EXPECT_DOUBLE_EQ( ss(1,0), 7);  EXPECT_DOUBLE_EQ( ss(1,1), 8.);

  ss(1,0) = 43.;
  ss(1,1) = 23.;

  EXPECT_DOUBLE_EQ( ss(1,0), 43.);
  EXPECT_DOUBLE_EQ( ss(1,1), 23.);
  EXPECT_DOUBLE_EQ( A(1,2), 43.);
  EXPECT_DOUBLE_EQ( A(1,3), 23.);
}

TEST(containers_expressions_kokkos, span3)
{
  using n_t  = Kokkos::View<double **, Kokkos::HostSpace>;
  using my_t = pressio::containers::DenseMatrix<n_t>;
  n_t A0("A",5,4);
  fillMatrix(A0);

  const my_t A(A0);
  auto ss = pressio::containers::subspan(A,
					 std::make_pair(0,2),
					 std::make_pair(2,4) );
  static_assert
    (std::is_const<
     typename std::remove_reference<decltype(ss(0,0))>::type
     >::value, "");

  EXPECT_DOUBLE_EQ( ss(0,0), 3.); EXPECT_DOUBLE_EQ( ss(0,1), 4.);
  EXPECT_DOUBLE_EQ( ss(1,0), 7);  EXPECT_DOUBLE_EQ( ss(1,1), 8.);
}
