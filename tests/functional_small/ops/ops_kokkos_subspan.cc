
#include <gtest/gtest.h>
#include "pressio/ops.hpp"

using vec_t = Kokkos::View<double*>;
using mat_t = Kokkos::View<double**>;

TEST(ops_kokkos, subspan_extent)
{
  mat_t A("A",5,5);
  std::pair<int,int> r(1,3);
  std::pair<int,int> c(2,5);

  auto ex = pressio::subspan(A,r,c);
  ASSERT_TRUE(pressio::ops::extent(ex,0)==2);
  ASSERT_TRUE(pressio::ops::extent(ex,1)==3);
}

TEST(ops_kokkos, subspan_scale)
{
  mat_t A("A",4,5);
  KokkosBlas::fill(A,1.);

  std::pair<int,int> r(1,3);
  std::pair<int,int> c(2,4);
  auto sp = pressio::subspan(A,r,c);
  pressio::ops::scale(sp, 3.);

  auto A_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A);
  ASSERT_DOUBLE_EQ(A_h(0,0),1.);
  ASSERT_DOUBLE_EQ(A_h(0,1),1.);
  ASSERT_DOUBLE_EQ(A_h(0,2),1.);
  ASSERT_DOUBLE_EQ(A_h(0,3),1.);
  ASSERT_DOUBLE_EQ(A_h(0,4),1.);

  ASSERT_DOUBLE_EQ(A_h(1,0),1.);
  ASSERT_DOUBLE_EQ(A_h(1,1),1.);
  ASSERT_DOUBLE_EQ(A_h(1,2),3.);
  ASSERT_DOUBLE_EQ(A_h(1,3),3.);
  ASSERT_DOUBLE_EQ(A_h(1,4),1.);

  ASSERT_DOUBLE_EQ(A_h(2,0),1.);
  ASSERT_DOUBLE_EQ(A_h(2,1),1.);
  ASSERT_DOUBLE_EQ(A_h(2,2),3.);
  ASSERT_DOUBLE_EQ(A_h(2,3),3.);
  ASSERT_DOUBLE_EQ(A_h(2,4),1.);

  ASSERT_DOUBLE_EQ(A_h(3,0),1.);
  ASSERT_DOUBLE_EQ(A_h(3,1),1.);
  ASSERT_DOUBLE_EQ(A_h(3,2),1.);
  ASSERT_DOUBLE_EQ(A_h(3,3),1.);
  ASSERT_DOUBLE_EQ(A_h(3,4),1.);
}

TEST(ops_kokkos, subspan_set_zero)
{
  mat_t A("A",4,5);
  KokkosBlas::fill(A,1.);

  std::pair<int,int> r(1,3);
  std::pair<int,int> c(2,4);
  auto sp = pressio::subspan(A,r,c);
  pressio::ops::set_zero(sp);

  auto A_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A);
  ASSERT_DOUBLE_EQ(A_h(0,0),1.);
  ASSERT_DOUBLE_EQ(A_h(0,1),1.);
  ASSERT_DOUBLE_EQ(A_h(0,2),1.);
  ASSERT_DOUBLE_EQ(A_h(0,3),1.);
  ASSERT_DOUBLE_EQ(A_h(0,4),1.);

  ASSERT_DOUBLE_EQ(A_h(1,0),1.);
  ASSERT_DOUBLE_EQ(A_h(1,1),1.);
  ASSERT_DOUBLE_EQ(A_h(1,2),0.);
  ASSERT_DOUBLE_EQ(A_h(1,3),0.);
  ASSERT_DOUBLE_EQ(A_h(1,4),1.);

  ASSERT_DOUBLE_EQ(A_h(2,0),1.);
  ASSERT_DOUBLE_EQ(A_h(2,1),1.);
  ASSERT_DOUBLE_EQ(A_h(2,2),0.);
  ASSERT_DOUBLE_EQ(A_h(2,3),0.);
  ASSERT_DOUBLE_EQ(A_h(2,4),1.);

  ASSERT_DOUBLE_EQ(A_h(3,0),1.);
  ASSERT_DOUBLE_EQ(A_h(3,1),1.);
  ASSERT_DOUBLE_EQ(A_h(3,2),1.);
  ASSERT_DOUBLE_EQ(A_h(3,3),1.);
  ASSERT_DOUBLE_EQ(A_h(3,4),1.);
}

TEST(ops_kokkos, subspan_fill)
{
  mat_t A("A",4,5);
  KokkosBlas::fill(A,1.);

  std::pair<int,int> r(1,3);
  std::pair<int,int> c(2,4);
  auto sp = pressio::subspan(A,r,c);
  pressio::ops::fill(sp, 44.);

  auto A_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A);
  ASSERT_DOUBLE_EQ(A_h(0,0),1.);
  ASSERT_DOUBLE_EQ(A_h(0,1),1.);
  ASSERT_DOUBLE_EQ(A_h(0,2),1.);
  ASSERT_DOUBLE_EQ(A_h(0,3),1.);
  ASSERT_DOUBLE_EQ(A_h(0,4),1.);

  ASSERT_DOUBLE_EQ(A_h(1,0),1.);
  ASSERT_DOUBLE_EQ(A_h(1,1),1.);
  ASSERT_DOUBLE_EQ(A_h(1,2),44.);
  ASSERT_DOUBLE_EQ(A_h(1,3),44.);
  ASSERT_DOUBLE_EQ(A_h(1,4),1.);

  ASSERT_DOUBLE_EQ(A_h(2,0),1.);
  ASSERT_DOUBLE_EQ(A_h(2,1),1.);
  ASSERT_DOUBLE_EQ(A_h(2,2),44.);
  ASSERT_DOUBLE_EQ(A_h(2,3),44.);
  ASSERT_DOUBLE_EQ(A_h(2,4),1.);

  ASSERT_DOUBLE_EQ(A_h(3,0),1.);
  ASSERT_DOUBLE_EQ(A_h(3,1),1.);
  ASSERT_DOUBLE_EQ(A_h(3,2),1.);
  ASSERT_DOUBLE_EQ(A_h(3,3),1.);
  ASSERT_DOUBLE_EQ(A_h(3,4),1.);
}
